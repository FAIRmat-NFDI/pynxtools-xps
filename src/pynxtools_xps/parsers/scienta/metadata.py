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
Metadata mapping for the Scienta parser.
"""

import datetime
import re
from functools import partial
from typing import Any

import numpy as np

# from pynxtools_xps.reader_utils import _re_map_single_value, convert_pascal_to_snake
from pynxtools_xps.mapping import (
    _convert_detector_acquisition_mode,
    _convert_energy_scan_mode,
    _convert_energy_type,
    _convert_slit_type,
    _MetadataContext,
    _ValueMap,
    parse_datetime,
)

# TODO: define these!
_POSSIBLE_DATE_FORMATS: list[str] = []


def _get_key_value_unit(line: str) -> tuple[str, Any, str | None]:
    """
    Split the line at the '=' sign and return a normalized key-value pair.
    """
    try:
        key_raw, value_str = line.split("=", 1)
        key, value, unit = _context.format(key_raw.strip(), value_str.strip())
    except ValueError:
        key, value, unit = "", "", ""
    return key, value, unit


def _extract_energy_units(energy_units: str) -> str:
    """
    Extract energy units from the strings for energy_units.
    Binding Energy [eV] -> eV

    """
    units = re.search(r"\[(.*?)\]", energy_units)
    return units.group(1) if units is not None else "eV"


def _separate_dimension_scale(scale: str) -> np.typing.ArrayLike:
    """
    Separate the str of the dimension scale into a numpy array

    Parameters
    ----------
    scale : str
        Str of the form "600 599.5 599 598.5 5".

    Returns
    -------
    np.ndarray
        Dimension scale as a numpy array.

    """
    return np.array([float(s) for s in scale.split(" ")])


def _construct_date_time(date_string: str, time_string: str) -> str | None:
    """
    Convert the native time format to the datetime string
    in the ISO 8601 format: '%Y-%b-%dT%H:%M:%S.%fZ'.

    """

    def _parse_date(date_string: str) -> datetime.datetime:
        possible_date_formats = ["%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y"]
        for date_fmt in possible_date_formats:
            try:
                return datetime.datetime.strptime(date_string, date_fmt)
            except ValueError as err:
                continue
        raise ValueError("Date format not recognized")

    def _parse_time(time_string: str) -> datetime.time:
        possible_time_formats = ["%H:%M:%S", "%I:%M:%S %p"]
        for time_fmt in possible_time_formats:
            try:
                return datetime.datetime.strptime(time_string, time_fmt).time()
            except ValueError as err:
                continue
        raise ValueError("Time format not recognized")

    try:
        date = _parse_date(date_string)
        time = _parse_time(time_string)
        datetime_obj = datetime.datetime.combine(date, time)

        return datetime_obj.astimezone().isoformat()

    except ValueError as err:
        raise ValueError(
            "Date and time could not be converted to ISO 8601 format."
        ) from err


# TODO: do we need this? If so, can it be more general?
def _check_valid_value(value: str | int | float | bool | np.ndarray) -> bool:
    """
    Check if a value is valid.

    Strings and arrays are considered valid if they are non-empty.
    Numbers and booleans are always considered valid.

    Args:
        value (str | int | float | bool | np.ndarray):
            The value to check. Can be a scalar, string, or NumPy array.

    Returns:
        bool: True if the value is valid, False otherwise.
    """
    if isinstance(value, str | int | float) and value is not None:
        return True
    if isinstance(value, bool):
        return True
    if isinstance(value, np.ndarray) and value.size != 0:
        return True
    return False


_KEY_MAP: dict[str, str] = {
    "number_of_regions": "no_of_regions",
    "version": "software_version",
    "dimension_name": "energy_units",
    "dimension_size": "energy_size",
    "dimension_scale": "energy_axis",
    "number_of_sweeps": "no_of_scans",
    "energy_unit": "energy_scale_2",
    "low_energy": "start_energy",
    "high_energy": "stop_energy",
    "energy_step": "step_size",
    "step_time": "dwell_time",
    "detector_first_x-_channel": "detector_first_x_channel",
    "detector_last_x-_channel": "detector_last_x_channel",
    "detector_first_y-_channel": "detector_first_y_channel",
    "detector_last_y-_channel": "detector_last_y_channel",
    "file": "data_file",
    "sequence": "sequence_file",
    "spectrum_name": "spectrum_type",
    "instrument": "instrument_name",
    "location": "vendor",
    "user": "user_name",
    "sample": "sample_name",
    "comments": "spectrum_comment",
    "date": "start_date",
}

_VALUE_MAP: _ValueMap = {
    "no_of_regions": int,
    "energy_size": int,
    "pass_energy": float,
    "no_of_scans": int,
    "excitation_energy": float,
    "center_energy": float,
    "start_energy": float,
    "stop_energy": float,
    "step_size": float,
    "dwell_time": float,
    "detector_first_x_channel": int,
    "detector_last_x_channel": int,
    "detector_first_y_channel": int,
    "detector_last_y_channel": int,
    "time_per_spectrum_channel": float,
    "energy_units": _extract_energy_units,
    "energy_axis": _separate_dimension_scale,
    "energy_scale": _convert_energy_type,
    "energy_scale_2": _convert_energy_type,
    "acquisition_mode": _convert_energy_scan_mode,
    "time_per_spectrum_channel": float,
    "manipulator_r1": float,
    "manipulator_r2": float,
    "start_time": partial(
        parse_datetime,
        possible_date_formats=_POSSIBLE_DATE_FORMATS,
    ),
    "stop_time": partial(
        parse_datetime,
        possible_date_formats=_POSSIBLE_DATE_FORMATS,
    ),
    "preset_type": lambda x: x.lower(),
    "source_type": lambda x: x.lower(),
    "energy_mode": _convert_energy_type,
    "duration": int,
    "acquisition/spectrum_definition/acquisition_mode": _convert_detector_acquisition_mode,
    "instrument/analyser/slit/type": _convert_slit_type,
}

_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {
    "pass_energy": "eV",
    "excitation_energy": "eV",
    "energy_axis": "eV",
    "center_energy": "eV",
    "start_energy": "eV",
    "stop_energy": "eV",
    "step_size": "eV",
    "dwell_time": "s",
    "time_per_spectrum_channel": "s",
}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
