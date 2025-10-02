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

import re
from pathlib import Path
from typing import Any, Union

from pynxtools_xps.kratos.kratos_data_model import KratosMetadata
from pynxtools_xps.reader_utils import convert_pascal_to_snake, extract_unit
from pynxtools_xps.value_mappers import convert_bool, parse_datetime

SETTINGS_MAP: dict[str, str] = {
    "location_i_d": "location_id",
    "tilt": "sample_tilt",
    "lens": "lens_mode",
    "anode_library": "anode_material",
    "start": "energy_start",
    "end": "energy_end",
    "centre": "energy_centre",
    "width": "energy_width",
    "x-ray_power": "x_ray_power",
}

KEYS_WITH_UNITS: list[str] = [
    "filament_current",
    "filament_bias",
    "charge_balance",
    "emission_current",
    "resolution",
    "energy_start",
    "energy_end",
    "energy_centre",
    "energy_width",
    "step_size",
    "dwell_time",
    "sweep_time",
    "x_ray_power",
]

UNIT_MISSING: dict[str, str] = {
    "sample_tilt": "degree",
}

POSSIBLE_DATE_FORMATS: list[str] = ["%d.%m.%Y %H:%M", "%d/%m/%Y %H:%M"]


class KratosParser:
    """
    A parser for reading in data from Kratos spectrometers.
    """

    def __init__(self):
        """
        Construct the parser.

        """
        self.raw_data: str = ""
        self.spectra: list[dict[str, Any]] = []

        self.metadata = KratosMetadata()

        self.value_function_map: dict[str, Any] = {
            "date_created": parse_datetime,
            "description": _convert_description,
            "charge_neutraliser": convert_bool,
            "deflection": _convert_xray_deflection,
        }

    def parse_file(self, file: str | Path, **kwargs):
        """
        TODO: parse actual data, not just metadata!

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
        list[dict[str, Any]]
            Flat list of dictionaries containing one spectrum each.

        """
        # TODO: write parser for actual Kratos data file
        # header, data = self._separate_header_and_data()
        # self.parse_header_into_metadata(header)

        return self.spectra

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
                key, value = line.split(" : ")

            except ValueError:
                continue

            key = convert_pascal_to_snake(key)
            key = SETTINGS_MAP.get(key, key)

            if key in datacls_fields:
                field_type = type(getattr(self.metadata, key))

                if key in KEYS_WITH_UNITS:
                    value, unit = extract_unit(key, value, UNIT_MISSING)  # type: ignore[assignment]
                    setattr(self.metadata, f"{key}_units", unit)

                value = self.map_values(key, value, field_type)

                setattr(self.metadata, key, value)

        self.metadata.validate_types()

    def map_values(self, key: str, value, field_type):
        """
        Map values to corresponding structure and field type.

        Parameters
        ----------
        key : str
            KratosMetadata dataclass key.
        value : str
            Original value.
        field_type : class
            Class to be mapped onto.

        Returns
        -------
        value
            Value of correct type and internal structure.

        """

        if key in self.value_function_map:
            map_fn = self.value_function_map[key]
            if key == "date_created":
                value = map_fn(value, POSSIBLE_DATE_FORMATS)
            else:
                value = map_fn(value)
        return field_type(value)

    def flatten_metadata(self) -> dict[str, Any]:
        """
        Flatten metadata dict so that key-value pairs of nested
        dictionaries are at the top level.

        Parameters
        ----------
        metadata_dict : dict
            Metadata dict with KratosMetadata fields as keys.

        Returns
        -------
        flattened_dict : dict
            Flatted metadata_dict without any nested substructure.

        """

        flattened_dict: dict[str, Any] = {}

        def setup_unit(flattened_dict: dict[str, Any], unit_key: str):
            """Sets up unit for in flattened_dict a given key ."""
            if "_units" in unit_key:
                new_key = unit_key.replace("_units", "/@units")
                flattened_dict[new_key] = flattened_dict.pop(unit_key)

            elif unit_key in KEYS_WITH_UNITS and not flattened_dict.get(
                f"{unit_key}/@units", None
            ):
                try:
                    flattened_dict[f"{unit_key}/@units"] = UNIT_MISSING[unit_key]
                except KeyError:
                    pass

        for key, value in self.metadata.dict().items():
            if isinstance(value, dict):
                for subkey, subvalue in value.items():
                    flattened_dict[f"{key}_{subkey}"] = subvalue
            else:
                flattened_dict[key] = value
                setup_unit(flattened_dict, key)

        return flattened_dict


def _convert_description(value: str) -> dict[str, Any]:
    """Map all items in description to a dictionary."""
    pattern = re.compile(
        r"\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)\s*([a-zA-Z]+)"
    )
    match = pattern.match(value)

    if match:
        x, y, z, unit = match.groups()
        return {
            "x": x,
            "x_units": unit,
            "y": y,
            "y_units": unit,
            "z": z,
            "z_units": unit,
        }
    else:
        raise ValueError(f"Invalid input string '{value}' for description.")


def _separate_val_and_unit(value: str) -> tuple[Any, str]:
    """Map all items in energy_referencing to a dictionary."""

    pattern = re.compile(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)([a-zA-Z]+)")
    match = pattern.match(value)

    if match:
        value, unit = match.groups()
        return value, unit
    else:
        raise ValueError(f"Input string '{value}' does not contain a value and unit.")


def _convert_xray_deflection(value: str) -> dict[str, Any]:
    """Convert deflection like (0.000, 0.000)mm to dict."""

    pattern = re.compile(
        r"\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)\s*([a-zA-Z]+)"
    )
    match = pattern.match(value)

    if match:
        x, y, unit = match.groups()
        return {
            "x": float(x),
            "y": float(y),
            "x_units": unit,
            "y_units": unit,
        }
    else:
        raise ValueError(f"Invalid input string: '{value}' for X-ray deflection.")
