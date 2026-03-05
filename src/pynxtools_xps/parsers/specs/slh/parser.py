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
from pathlib import Path
from typing import ClassVar

import pandas as pd

from pynxtools_xps.parsers.base import ParsedSpectrum, _XPSMetadataParser, _XPSParser
from pynxtools_xps.parsers.specs.sle.parser import SpecsSLEParser
from pynxtools_xps.parsers.versioning import VersionRange


class SPECSMetadataCSVParser(_XPSMetadataParser):
    """
    Metadata-only parser for SpecsLab Prodigy metadata CSV exports.

    Parses a ``.csv`` result file exported from SpecsLab .slh files and
    injects metadata into the ``metadata`` of matching
    ``ParsedSpectrum`` objects produced by ``SpecsSLEParser``.

    Compatible primary parser: ``SpecsSLEParser``.
    """

    compatible_primary_parser: ClassVar[type[SpecsSLEParser]] = SpecsSLEParser
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".csv",)

    requires_version: ClassVar[bool] = True
    supported_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid
        ((0, 6), (0, 7)),  # 0.6, not 0.7
    )
    supported_parser_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid??
        # TODO: generalize this
        ((1, 2), (1, 3)),  # 1.2, not 1.3
        ((1, 8), (1, 14)),  # 4.1 – 4.100
    )

    def __init__(self):
        super().__init__()
        self._dataframe = pd.DataFrame()

    def matches_file(self, file: Path) -> bool:
        """Return True for SpecsLab Prodigy metadata CSV exports."""
        return False

    def _supports_parser(cls, parser: _XPSParser) -> bool:
        if not super()._supports_parser(parser):
            return False
        # TODO: check that the supported prodigy versions match
        return True

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
    Metadata-only parser for SpecsLab Prodigy history (.slh) files.

    Parses a ``.slh`` history file exported from SpecsLab Prodigy and
    injects metadata into the ``metadata`` of matching
    ``ParsedSpectrum`` objects produced by ``SpecsSLEParser``.

    Compatible primary parser: ``SpecsSLEParser``.
    """

    compatible_primary_parser: ClassVar[type[SpecsSLEParser]] = SpecsSLEParser
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".slh",)

    requires_version: ClassVar[bool] = True
    supported_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid
        ((0, 6), (0, 7)),  # 0.6, not 0.7
    )
    supported_primary_parser_versions: ClassVar[tuple[VersionRange, ...]] = (
        # TODO: are these valid??
        # TODO: generalize this
        ((1, 2), (1, 3)),  # 1.2, not 1.3
        ((1, 8), (1, 14)),  # 4.1 – 4.100
    )

    def __init__(self):
        super().__init__()
        self._dataframe = pd.DataFrame()
        self.con: sqlite3.Connection
        self.cur: sqlite3.Cursor

        self.encoding = ["f", 4]

    def matches_file(self, file: Path) -> bool:
        """Return True for SpecsLab Prodigy metadata CSV exports."""
        # TODO: missing implementation
        return False

    def _supports_parser(cls, parser: _XPSParser) -> bool:
        if not super()._supports_parser(parser):
            return False
        # TODO: check that the supported prodigy versions match
        return True

    def _initiate_file_connection(self, file: str | Path):
        """Set the SQLlite connection of the file to be opened."""
        sql_connection = file
        self.con = sqlite3.connect(sql_connection)
        self.cur = self.con.cursor()

    def _close_con(self):
        """
        Close the database connection.

        Returns
        -------
        None.

        """
        self.con.close()

    def _parse(self, file, **kwargs):
        # initiate connection to sql file
        self._initiate_file_connection(file)

        # read and parse sle file
        param_df = self._get_parameter_metadata()
        # parameter_info = self._get_parameter_info()
        data = self._get_data_for_param_df(param_df)

        df = param_df.update(data)

        return data  # df.to_dict(orient="records")

    def _get_parameter_metadata(self):
        """

        Returns
        -------
        None.

        """
        query = "SELECT * FROM 'ParameterHistory'"
        results = self.cur.execute(query).fetchall()

        columns = [
            "id",
            "name",
            "unique_device_name",
            "device",
            "measured_param",
            "base_time",
        ]
        df = pd.DataFrame(columns=columns)

        for result in results:
            new_row = {}
            for key, value in zip(columns, result):
                if key == "base_time":
                    # convert base time to time stamp
                    value = datetime.strptime(value, "%Y-%b-%d %H:%M:%S.%f")
                new_row[key] = value

            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

        df.set_index("id", drop=True, inplace=True)

        return df

    def _get_parameter_info(self):
        """

        Returns
        -------
        None.

        """
        parameter_info = []
        query = "SELECT ID, DeviceType, Command, ReadablePar, Unit, Scaling, Representation FROM ParameterInfo"
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

        for result in results:
            param_dict = {}
            for key, value in zip(columns, result):
                param_dict[key] = value

            parameter_info += [param_dict]
        return parameter_info

    def _get_data_for_param_df(self, param_df):
        data_df = pd.DataFrame()

        for ID, row in param_df.iterrows():
            query = f"SELECT Observation, Offset_s FROM 'NumericalHistoryData' WHERE ID={ID}"
            data = self.cur.execute(query).fetchall()

            data = pd.DataFrame(data)
            data.columns = ["Observation", "Offset_s"]
            data["date_time"] = (
                pd.TimedeltaIndex(data.Offset_s, unit="s") + row["base_time"]
            )
            data = self._resample_df(data, interval="1s").reindex()
            data = data.to_dict(orient="list")

            data_df = pd.concat([data_df, pd.DataFrame([data])], ignore_index=True)

        return data_df

    def _resample_df(self, df, interval="1s"):
        df_resampled = df.copy()
        # Resample to 1s intervals

        df_resampled = df_resampled.resample(
            interval, on="date_time", label="right"
        ).mean()
        # If a value is not given, take the previous value
        df_resampled["Observation"] = df_resampled["Observation"].ffill()
        df_resampled["Offset_s"] = df_resampled["Offset_s"].ffill()

        return df_resampled

    def _check_encoding(self):
        """
        Check whether the binary data should be decoded float or double.

        Returns
        -------
        None.

        """
        # TODO: check if the same version in the SLE PARSER can be reused
        query = "SELECT LENGTH(Data),ChunkSize FROM CountRateData LIMIT 1"
        self.cur.execute(query)
        data, chunksize = self.cur.fetchall()[0]

        encodings_map = {
            "double": ["d", 8],
            "float": ["f", 4],
        }

        if data / chunksize == 4:
            self.encoding = encodings_map["float"]
        elif data / chunksize == 8:
            self.encoding = encodings_map["double"]
        else:
            print("This binary encoding is not supported.")

    def _build_entries(self, data):
        # TODO: work with self._dataframe
        pass

    # def _convert_to_common_format(self):
    #     """
    #     Reformat spectra into the format needed for the Converter object
    #     """
    #     for spec in self.spectra:
    #         spec["data"] = {}
    #         spec["data"]["x"] = self._get_energy_data(spec)

    #         channels = [
    #             key
    #             for key in spec
    #             if any(name in key for name in ["cps_ch_", "cps_calib"])
    #         ]

    #         for channel_key in channels:
    #             spec["data"][channel_key] = np.array(spec[channel_key])
    #         for channel_key in channels:
    #             spec.pop(channel_key)

    #         spec["y_units"] = "Counts per Second"
