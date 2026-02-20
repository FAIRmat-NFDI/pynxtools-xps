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
Parser for reading XPS (X-ray Photoelectron Spectroscopy) metadata from
Kratos instruments (currently only after exporting to .vms format), to be
passed to MPES nxdl (NeXus Definition Language) template.
"""

from pathlib import Path

from pynxtools_xps.mapping import _Value
from pynxtools_xps.parsers.base import _XPSParser
from pynxtools_xps.parsers.kratos.data_model import KratosMetadata
from pynxtools_xps.parsers.kratos.metadata import _context


class KratosParser(_XPSParser):
    """
    A parser for reading in data from Kratos spectrometers.
    """

    def __init__(self):
        """
        Construct the parser.

        """
        self.metadata = KratosMetadata()

    def parse_file(self, file: str | Path, **kwargs):
        # TODO: parse actual data, not just metadata!
        """
        Parse the data file into a list of dictionaries.

        Parsed data is stored in the attribute 'self.data'.
        Each dictionary in the data list is a grouping of related
        attributes. The dictionaries are later re-structured into a
        nested dictionary that more closely resembles the domain logic.

        Parameters
        ----------
        file : Union[str, Path]
            XPS data filepath.

        Returns
        -------
        list[dict[str, _Value]]
            Flat list of dictionaries containing one spectrum each.

        """
        # TODO: write parser for actual Kratos data file
        # header, data = self._separate_header_and_data()
        # self.parse_header_into_metadata(header)

        return self.data

    def parse_header_into_metadata(self, header: list[str]):
        """
        Parse header into KratosMetadata dataclass.

        Parameters
        ----------
        header : List[str]
            Header data for one spectrum as a String.

        """
        datacls_fields = list(self.metadata.__dataclass_fields__.keys())
        datacls_fields = [field for field in datacls_fields if "_units" not in field]

        for line in header:
            try:
                key, raw_value = line.split(" : ")

            except ValueError:
                continue

            key = _context.normalize_key(key)

            if key in datacls_fields:
                field_type = type(getattr(self.metadata, key))

                key, value, unit = _context.format(key, raw_value)

                setattr(self.metadata, key, field_type(value))

                if unit:
                    setattr(self.metadata, f"{key}_units", unit)

        self.metadata.validate_types()

    def flatten_metadata(self) -> dict[str, _Value]:
        """
        Flatten metadata so that nested dictionaries (to depth of 1) are
        lifted to the top level using underscore-separated keys.

        Keys ending in '_units' are converted to the NOMAD-style
        '/@units' suffix.

        Parameters
        ----------
        metadata_dict : dict
            Metadata dict with KratosMetadata fields as keys.

        Returns
        -------
        flattened_dict : dict
            Flatted metadata_dict without any nested substructure.

        """

        flattened: dict[str, _Value] = {}

        def change_unit_path(flattened: dict[str, _Value], key: str):
            """Sets up unit for in flattened_dict a given key ."""
            if "_units" in key:
                new_key = key.replace("_units", "/@units")
                flattened[new_key] = flattened.pop(key)

        for key, value in self.metadata.dict().items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    flattened[f"{key}_{sub_key}"] = sub_value
                    change_unit_path(flattened, f"{key}_{sub_key}")
            else:
                flattened[key] = value
            change_unit_path(flattened, key)

        return flattened
