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

import re
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

# .csv export timestamp columns -- see SPECSMetadataCSVParser docstring.
_UNIVERSAL_TIME_COLUMN = "Universal Time"
_LOCAL_DATE_COLUMN = "Local Date"
_LOCAL_TIME_COLUMN = "Local Time of Day"
_RELATIVE_TIME_COLUMN = "Relative Time [s]"
_TIME_COLUMNS: tuple[str, ...] = (
    _UNIVERSAL_TIME_COLUMN,
    _LOCAL_DATE_COLUMN,
    _LOCAL_TIME_COLUMN,
)
_NON_PARAMETER_COLUMNS = frozenset({*_TIME_COLUMNS, _RELATIVE_TIME_COLUMN})
_LOCAL_DATE_TIME_FORMAT = "%Y-%m-%d %H:%M:%S.%f"

# e.g. "CoilCurrent Current [µA] (Phoibos)": name, unit in [], device in ().
_HEADER_RE = re.compile(
    r"^(?P<param>.*?)\s*\[(?P<unit>[^\]]*)\]\s*(?:\((?P<device>[^)]*)\))?\s*$"
)


class _ParameterLogWindowMixin:
    """Shared ``update_main_file_data`` for parsers storing a session's
    logged parameters in ``self._parameter_logs`` (columns: ``name``,
    ``device``, ``measured_param``, ``unit``, ``timestamps``, ``values``).
    Used by both ``SPECSMetadataSLHParser`` and ``SPECSMetadataCSVParser``,
    which read the same underlying data from different file formats.
    """

    _parameter_logs: pd.DataFrame

    def update_main_file_data(self, main_file_data: dict[str, ParsedSpectrum]) -> None:
        """Attach each entry's parameter readings within its own time window."""
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
        """Acquisition [start, end] from an entry's own time_stamp/dwell_time/
        n_values/total_scans metadata."""
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


class SPECSMetadataCSVParser(_ParameterLogWindowMixin, _XPSMetadataParser):
    """
    Metadata-only parser for SpecsLab Prodigy parameter-history CSV exports.

    Same underlying data as an ``.slh`` file (see ``SPECSMetadataSLHParser``),
    exported as a flat, semicolon-delimited table. Header shape::

        "Universal Time";"Local Date";"Local Time of Day";"Relative Time [s]";"CoilCurrent Current [µA] (Phoibos)";...
        2025-Apr-03 06:57:15.625973;2025-04-03;08:57:15.625973;0.000135;0;...

    - Delimiter is ``;``, not ``,``, despite the ``.csv`` extension.
    - Timestamp: prefers ``Universal Time``; falls back to ``Local Date``
      + ``Local Time of Day``; skips gracefully if neither is usable.
    - Sparse table: blank cells mean "unchanged", not zero. Each
      parameter's series only uses its own non-blank rows, then resamples
      with forward-fill (``csv_resample_interval``), same as the .slh path.
    - Unit/device come from the header (``"<param> [<unit>] (<device>)"``).
    - Read with ``pandas.read_csv`` (files can be very large).

    ``supported_versions`` is empty: nothing in a CSV export carries a
    version marker to gate on.

    Compatible primary parser: ``SPECSSLEParser``.
    """

    compatible_primary_parser: ClassVar[type[SPECSSLEParser]] = SPECSSLEParser
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".csv",)

    requires_version: ClassVar[bool] = False
    supported_versions: ClassVar[tuple[VersionRange, ...]] = ()
    supported_primary_parser_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid??
        ((1, 2), (1, 3)),  # 1.2, not 1.3
        ((1, 8), (1, 14)),
    )

    def __init__(self):
        super().__init__()
        self._parameter_logs: pd.DataFrame = pd.DataFrame()

    def matches_file(self, file: Path) -> bool:
        """Return True for SpecsLab Prodigy metadata CSV exports.

        Peeks at the header row: confirms the ``;`` delimiter and that at
        least one of the three time columns is present. Rejects
        comma-delimited CSVs and binary formats such as .sle/.slh.
        """
        try:
            with open(file, encoding="utf-8") as f:
                header_line = f.readline()
        except Exception:
            return False

        if ";" not in header_line:
            return False

        fields = [
            field.strip().strip('"') for field in header_line.strip("\r\n").split(";")
        ]
        return any(column in fields for column in _TIME_COLUMNS)

    def _parse(self, file: Path, **kwargs: Any) -> None:
        """Read every logged parameter's full time series from *file*.

        kwargs
        ------
        csv_resample_interval : str | None
            Pandas offset alias (default ``"1s"``) to resample each
            parameter's non-blank observations to, forward-filling gaps.
            A separate kwarg from ``slh_resample_interval`` so the two
            file types can be tuned independently.
        """
        resample_interval = kwargs.get("csv_resample_interval", "1s")

        df = pd.read_csv(file, sep=";")
        df.columns = [str(column).strip() for column in df.columns]

        date_time = self._extract_date_time(df, file)
        if date_time is None:
            self._parameter_logs = pd.DataFrame()
            return

        df["_date_time"] = date_time
        df = df[df["_date_time"].notna()]

        self._parameter_logs = self._build_parameter_logs(df, resample_interval)
        # self._data deliberately stays empty -- see class docstring.

    @staticmethod
    def _extract_date_time(df: pd.DataFrame, file: Path) -> pd.Series | None:
        """Timestamp Series from Universal Time, else Local Date + Local
        Time of Day, else None (logged) if neither is present."""
        if _UNIVERSAL_TIME_COLUMN in df.columns:
            return pd.to_datetime(
                df[_UNIVERSAL_TIME_COLUMN], format=_BASE_TIME_FORMAT, errors="coerce"
            )

        if _LOCAL_DATE_COLUMN in df.columns and _LOCAL_TIME_COLUMN in df.columns:
            combined = (
                df[_LOCAL_DATE_COLUMN].astype(str)
                + " "
                + df[_LOCAL_TIME_COLUMN].astype(str)
            )
            return pd.to_datetime(
                combined, format=_LOCAL_DATE_TIME_FORMAT, errors="coerce"
            )

        _logger.warning(
            "SPECS CSV metadata file %s has none of the expected time columns "
            "(%s); skipping.",
            file,
            ", ".join(_TIME_COLUMNS),
        )
        return None

    def _build_parameter_logs(
        self, df: pd.DataFrame, resample_interval: str | None
    ) -> pd.DataFrame:
        """Build one row per parameter column, same shape as the .slh path."""
        parameter_columns = [
            c
            for c in df.columns
            if c not in _NON_PARAMETER_COLUMNS and c != "_date_time"
        ]

        rows = []
        for column in parameter_columns:
            measured_param, unit, device = self._parse_header(column)

            values = pd.to_numeric(df[column], errors="coerce")
            valid = values.notna()
            if not valid.any():
                _logger.debug(
                    "SPECS CSV parameter column %r has no observations.", column
                )
                continue

            data = pd.DataFrame(
                {
                    "date_time": df.loc[valid, "_date_time"],
                    "Observation": values.loc[valid],
                }
            ).sort_values("date_time")

            if resample_interval:
                data = self._resample_df(data, interval=resample_interval)
            else:
                data = data.set_index("date_time")

            rows.append(
                {
                    "name": f"{measured_param} ({device})"
                    if device
                    else measured_param,
                    "device": device,
                    "measured_param": measured_param,
                    "unit": unit,
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

    @staticmethod
    def _parse_header(column: str) -> tuple[str, str, str]:
        """Split ``"<param> [<unit>] (<device>)"`` into its three parts.
        Falls back to the whole header as the name if it doesn't match."""
        match = _HEADER_RE.match(column.strip())
        if not match:
            return column.strip(), "", ""
        measured_param = match.group("param").strip()
        unit = match.group("unit").strip()
        device = (match.group("device") or "").strip()
        return measured_param, unit, device

    @staticmethod
    def _resample_df(df: pd.DataFrame, interval: str) -> pd.DataFrame:
        """Resample to *interval*, forward-filling gaps between observations."""
        resampled = df.set_index("date_time").resample(interval, label="right").mean()
        resampled["Observation"] = resampled["Observation"].ffill()
        return resampled.dropna(subset=["Observation"])


class SPECSMetadataSLHParser(_ParameterLogWindowMixin, _XPSMetadataParser):
    """
    Metadata-only parser for SpecsLab Prodigy parameter-history (.slh) files.

    An .slh file is a SQLite database logging one group of instrument
    parameters (e.g. "xray", "nap_parameters", "phoibos_voltages") over a
    whole session, not any single spectrum's own acquisition window.
    Schema: ``ParameterHistory(ID, Name, UniqueDev, Device, Parameter,
    BaseTime)`` + ``ParameterInfo(ID, DeviceType, Command, ReadablePar,
    Prototype, Unit, Scaling, Representation)`` (joined 1:1 on ``ID``) +
    ``NumericalHistoryData(ID, Offset_s, Observation)``. Schema-driven, not
    hardcoded to specific parameter names.

    ``_parse`` leaves ``self._data`` empty (an .slh file isn't
    entry-keyed) and stores every parameter's time series in
    ``self._parameter_logs`` instead; ``update_main_file_data`` slices
    that per entry using each entry's own ``time_stamp``/``dwell_time``/
    ``n_values``/``total_scans`` metadata.

    Notes:
    - ``ParameterInfo.Unit`` can contain non-UTF-8 bytes -- decoded with
      ``errors="replace"``.
    - ``NumericalHistoryData`` rows aren't guaranteed chronological order
      -- always queried with ``ORDER BY``.
    - ``Observation`` is stored unscaled; multiplied by ``Scaling`` here.

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
        """Return True for SpecsLab Prodigy history (.slh) files.

        Positive structural identification, same convention as
        ``SPECSSLEParser.matches_file``: SQLite magic bytes first (cheap),
        then confirm the specific tables this format always has -- .sle
        files are SQLite too, so the magic bytes alone don't distinguish
        them; the table set does (.sle has no ``ParameterHistory`` table
        at all).
        """
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
        """Read `Configuration.Value WHERE Key='Version'` (same convention
        as `SPECSSLEParser._get_version`)."""
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
        # Unit can contain non-UTF-8 bytes; a strict decode would crash the parse.
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
            forward-filling gaps. Pass ``None`` to keep the raw,
            irregular timestamps instead.
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

        ``Observation`` is stored unscaled; multiplied by ``Scaling`` here.
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

    # update_main_file_data / _entry_time_window: see _ParameterLogWindowMixin,
    # shared verbatim with SPECSMetadataCSVParser.
