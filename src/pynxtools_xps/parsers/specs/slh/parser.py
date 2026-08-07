# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
Metadata parser for reading for metadata exported from SpecsLab Prodigy.
"""

import csv
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, ClassVar

import pandas as pd

from pynxtools_xps.logging import _logger
from pynxtools_xps.parsers.base import ParsedSpectrum, _XPSMetadataParser
from pynxtools_xps.parsers.specs.sle.parser import SPECSSLEParser
from pynxtools_xps.parsers.versioning import (
    VersionRange,
    VersionTuple,
    normalize_version,
)

_BASE_TIME_FORMAT = "%Y-%b-%d %H:%M:%S.%f"


class SPECSMetadataCSVParser(_XPSMetadataParser):
    """
    Metadata-only parser for SpecsLab Prodigy metadata CSV exports.

    Parses a ``.csv`` result file exported from SpecsLab .slh files and
    injects metadata into the ``metadata`` of matching
    ``ParsedSpectrum`` objects produced by ``SPECSSLEParser``.

    Compatible primary parser: ``SPECSSLEParser``.
    """

    compatible_primary_parser: ClassVar[type[SPECSSLEParser]] = SPECSSLEParser
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".csv",)

    requires_version: ClassVar[bool] = True
    supported_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid
        ((0, 6), (0, 7)),  # 0.6, not 0.7
    )
    supported_primary_parser_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid??
        ((1, 2), (1, 3)),  # 1.2, not 1.3
        ((1, 8), (1, 14)),
    )

    def __init__(self):
        super().__init__()
        self._dataframe = pd.DataFrame()

    def matches_file(self, file: Path) -> bool:
        """Return True for SpecsLab Prodigy metadata CSV exports."""
        # TODO: not yet implemented -- CSV history exports are a separate,
        # not-yet-scoped follow-up to the .slh work (see SPECSMetadataSLHParser).
        return False

    def _read_csv(self):
        pass

    def _build_entries(self):
        # TODO: work with self._dataframe
        pass

    def _parse(self, file: Path, **kwargs) -> None:
        """Parse quantification CSV and populate ``self._data``."""
        with open(file) as f:
            all_rows = list(csv.reader(f, delimiter="\t"))

        data = self._read_csv()
        self._build_entries()

    def update_main_file_data(self, main_file_data: dict[str, ParsedSpectrum]) -> None:
        """
        Merge ``self._data`` metadata into matching entries of *main_file_data*.

        Args:
            main_file_data: Mapping from NeXus entry name to ``ParsedSpectrum``,
                as produced by the compatible primary parser.
        """
        if not self._data:
            return


class SPECSMetadataSLHParser(_XPSMetadataParser):
    """
    Metadata-only parser for SpecsLab Prodigy parameter-history (.slh) files.

    An .slh file is a SQLite database logging one group of instrument
    parameters (e.g. "xray", "nap_parameters", "phoibos_voltages") over a
    whole session, not any single spectrum's own acquisition window.
    Schema: ``ParameterHistory(ID, Name, UniqueDev, Device, Parameter,
    BaseTime)`` joined 1:1 on ``ID`` with ``ParameterInfo(ID, DeviceType,
    Command, ReadablePar, Prototype, Unit, Scaling, Representation)``, plus
    ``NumericalHistoryData(ID, Offset_s, Observation)``. Schema-driven, not
    hardcoded to specific parameter names -- ``ParameterHistory`` is the
    source of truth for what a given file logs.

    Because an .slh file isn't entry-keyed, ``_parse`` leaves ``self._data``
    empty (just to satisfy the base class) and instead stores every
    parameter's full time series in ``self._parameter_logs``.
    ``update_main_file_data`` is overridden to slice that per entry, using
    each entry's own ``time_stamp``/``dwell_time``/``n_values``/
    ``total_scans`` metadata to compute its acquisition window.

    Handles two real-world quirks: ``ParameterInfo.Unit`` can contain
    non-UTF-8 bytes, and ``NumericalHistoryData`` isn't guaranteed to be
    stored in chronological order.

    Compatible primary parser: ``SPECSSLEParser``.
    """

    compatible_primary_parser: ClassVar[type[SPECSSLEParser]] = SPECSSLEParser
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".slh",)

    requires_version: ClassVar[bool] = True
    supported_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid
        ((0, 6), (0, 7)),  # 0.6, not 0.7
    )
    supported_primary_parser_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid??
        ((1, 2), (1, 3)),  # 1.2, not 1.3
        ((1, 8), (1, 14)),
    )

    _SQLITE_MAGIC = b"SQLite format 3\x00"
    _REQUIRED_TABLES = frozenset(
        {"ParameterHistory", "NumericalHistoryData", "ParameterInfo"}
    )

    def __init__(self):
        super().__init__()
        self.con: sqlite3.Connection
        self.cur: sqlite3.Cursor
        self._parameter_logs: pd.DataFrame = pd.DataFrame()

    def matches_file(self, file: Path) -> bool:
        """SQLite magic bytes, then the specific table set this format
        always has -- .sle files are SQLite too, so the magic bytes alone
        don't distinguish them."""
        try:
            with open(file, "rb") as f:
                if f.read(16) != self._SQLITE_MAGIC:
                    return False
            conn = sqlite3.connect(str(file))
            try:
                cur = conn.cursor()
                cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
                found = {row[0] for row in cur.fetchall()}
                return self._REQUIRED_TABLES <= found
            finally:
                conn.close()
        except Exception:
            return False

    def detect_version(self, file: Path) -> VersionTuple | None:
        """Read `Configuration.Value WHERE Key='Version'`, same convention
        as `SPECSSLEParser._get_version`."""
        try:
            self._initiate_file_connection(file)
            try:
                self.cur.execute("SELECT Value FROM Configuration WHERE Key='Version'")
                row = self.cur.fetchone()
            finally:
                self._close_con()
            return normalize_version(row[0]) if row else None
        except Exception:
            return None

    def _initiate_file_connection(self, file: str | Path) -> None:
        """Open the SQLite connection to *file*."""
        self.con = sqlite3.connect(str(file))
        # ParameterInfo.Unit can contain non-UTF-8 bytes (e.g. mis-encoded "°C").
        self.con.text_factory = lambda raw: raw.decode("utf-8", errors="replace")
        self.cur = self.con.cursor()

    def _close_con(self) -> None:
        self.con.close()

    def _parse(self, file: Path, **kwargs: Any) -> None:
        """Read every logged parameter's full time series from *file*.

        kwargs
        ------
        slh_resample_interval : str | None
            Pandas offset alias (default ``"1s"``) to resample each
            parameter's raw, irregularly-sampled observations to,
            forward-filling gaps. Pass ``None`` to keep raw timestamps
            instead. Passed through from ``params.yaml``, same mechanism
            as ``SPECSSLEParser``'s ``remove_align``/``sum_channels``.
        """
        self._initiate_file_connection(file)
        try:
            param_df = self._get_parameter_metadata()
            info_df = self._get_parameter_info()
            param_df = param_df.join(info_df, how="left")

            resample_interval = kwargs.get("slh_resample_interval", "1s")
            self._parameter_logs = self._get_data_for_param_df(
                param_df, resample_interval=resample_interval
            )
        finally:
            self._close_con()
        # self._data deliberately stays empty -- see class docstring.

    def _get_parameter_metadata(self) -> pd.DataFrame:
        """Read ``ParameterHistory``: which parameters this file logs."""
        query = "SELECT ID, Name, UniqueDev, Device, Parameter, BaseTime FROM ParameterHistory"
        results = self.cur.execute(query).fetchall()

        columns = [
            "id",
            "name",
            "unique_device_name",
            "device",
            "measured_param",
            "base_time",
        ]
        rows = []
        for result in results:
            row = dict(zip(columns, result))
            row["base_time"] = datetime.strptime(row["base_time"], _BASE_TIME_FORMAT)
            rows.append(row)

        df = pd.DataFrame(rows, columns=columns)
        df.set_index("id", drop=True, inplace=True)
        return df

    def _get_parameter_info(self) -> pd.DataFrame:
        """Read ``ParameterInfo``: unit + scaling for each parameter."""
        query = (
            "SELECT ID, DeviceType, Command, ReadablePar, Unit, Scaling, "
            "Representation FROM ParameterInfo"
        )
        results = self.cur.execute(query).fetchall()

        columns = [
            "id",
            "device_type",
            "command",
            "readable_parameter",
            "unit",
            "scaling",
            "representation",
        ]
        df = pd.DataFrame(results, columns=columns)
        df.set_index("id", drop=True, inplace=True)
        return df

    def _get_data_for_param_df(
        self, param_df: pd.DataFrame, resample_interval: str | None
    ) -> pd.DataFrame:
        """Read every parameter's raw observations, scaled and timestamped.

        `Observation` is stored unscaled, so it's multiplied by `Scaling`.
        """
        rows = []
        for param_id, row in param_df.iterrows():
            self.cur.execute(
                "SELECT Offset_s, Observation FROM NumericalHistoryData "
                "WHERE ID = ? ORDER BY Offset_s ASC",
                (int(param_id),),
            )
            raw = self.cur.fetchall()
            if not raw:
                _logger.debug(
                    "SLH parameter %r (id=%s) has no numerical history data.",
                    row["name"],
                    param_id,
                )
                continue

            data = pd.DataFrame(raw, columns=["Offset_s", "Observation"])
            scaling = row.get("scaling")
            scaling = 1.0 if scaling in (None, "") else float(scaling)
            data["Observation"] = data["Observation"] * scaling
            data["date_time"] = row["base_time"] + pd.to_timedelta(
                data["Offset_s"], unit="s"
            )

            if resample_interval:
                data = self._resample_df(data, interval=resample_interval)
            else:
                data = data.set_index("date_time")

            rows.append(
                {
                    "name": row["name"],
                    "device": row["device"] or "",
                    "measured_param": row["measured_param"] or "",
                    "unit": row.get("unit") or "",
                    "timestamps": data.index.strftime("%Y-%m-%dT%H:%M:%S.%fZ").tolist(),
                    "values": data["Observation"].tolist(),
                }
            )
        return pd.DataFrame(
            rows,
            columns=[
                "name",
                "device",
                "measured_param",
                "unit",
                "timestamps",
                "values",
            ],
        )

    def _resample_df(self, df: pd.DataFrame, interval: str) -> pd.DataFrame:
        """Resample to *interval*, forward-filling gaps between observations."""
        resampled = df.set_index("date_time").resample(interval, label="right").mean()
        resampled["Observation"] = resampled["Observation"].ffill()
        resampled["Offset_s"] = resampled["Offset_s"].ffill()
        return resampled.dropna(subset=["Observation"])

    def update_main_file_data(self, main_file_data: dict[str, ParsedSpectrum]) -> None:
        """Slice this file's parameter logs into each matching entry's
        acquisition window, computed from the entry's own metadata.

        Overrides the default entry-name-matching merge -- an .slh file
        isn't entry-keyed.
        """
        if self._parameter_logs.empty:
            return

        for spectrum in main_file_data.values():
            window = self._entry_time_window(spectrum)
            if window is None:
                continue
            start, end = window

            for _, log in self._parameter_logs.iterrows():
                if not log["timestamps"]:
                    continue
                key_base = f"{log['measured_param']}_{log['device']}".strip(
                    "_"
                ).replace(" ", "_")
                timestamps = pd.to_datetime(log["timestamps"])
                mask = (timestamps >= start) & (timestamps <= end)
                if not mask.any():
                    continue

                spectrum.metadata[key_base] = [
                    v for v, m in zip(log["values"], mask) if m
                ]
                spectrum.metadata[f"{key_base}/@units"] = log["unit"]
                spectrum.metadata[f"{key_base}/time"] = [
                    t for t, m in zip(log["timestamps"], mask) if m
                ]

    @staticmethod
    def _entry_time_window(
        spectrum: ParsedSpectrum,
    ) -> tuple[pd.Timestamp, pd.Timestamp] | None:
        """Compute one entry's acquisition [start, end] from its own metadata.

        Duration is estimated as `dwell_time * n_values * total_scans` --
        a generous first-order estimate, not exact.
        """
        start_raw = spectrum.metadata.get("time_stamp")
        if not start_raw:
            return None
        try:
            start = pd.Timestamp(start_raw)
        except (ValueError, TypeError):
            return None

        dwell_time = spectrum.metadata.get("dwell_time") or 0.0
        n_values = spectrum.metadata.get("n_values") or 0
        total_scans = spectrum.metadata.get("total_scans") or 1
        try:
            duration_s = float(dwell_time) * float(n_values) * float(total_scans)
        except (TypeError, ValueError):
            duration_s = 0.0

        end = start + pd.Timedelta(seconds=duration_s) if duration_s > 0 else start
        return start, end
