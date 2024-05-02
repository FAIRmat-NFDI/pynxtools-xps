"""
Helper functions for populating NXmpes template
"""

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
import re
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Union
from pathlib import Path
from scipy.interpolate import interp1d
import numpy as np

from dataclasses import dataclass


@dataclass
class XpsDataclass:
    """Generic class to hold a data model and a type validation method."""

    def validate_types(self):
        ret = True
        for field_name, field_def in self.__dataclass_fields__.items():
            actual_type = type(getattr(self, field_name))
            if actual_type != field_def.type:
                print(f"\t{field_name}: '{actual_type}' instead of '{field_def.type}'")
                ret = False
        return ret

    def __post_init__(self):
        if not self.validate_types():
            raise ValueError("Wrong types")

    def dict(self):
        return self.__dict__.copy()


class XPSMapper(ABC):
    """Abstract base class from mapping from a parser to NXmpes template"""

    def __init__(self):
        self.file: Union[str, Path] = ""
        self.raw_data: List[str] = []
        self._xps_dict: Dict[str, Any] = {}
        self._root_path = "/ENTRY[entry]"

        self.parser = None

    @abstractmethod
    def _select_parser(self):
        """
        Select the correct parser for the file extension and format.

        Should be implemented by the inheriting mapper.

        Returns
        -------
        Parser

        """

    @property
    def data_dict(self) -> dict:
        """Getter property."""
        return self._xps_dict

    def parse_file(self, file, **kwargs):
        """
        Parse the file using the Scienta TXT parser.

        """
        self.file = file
        self.parser = self._select_parser()
        self.raw_data = self.parser.parse_file(file, **kwargs)

        file_key = f"{self._root_path}/File"
        self._xps_dict[file_key] = file

        self.construct_data()

        return self.data_dict

    @abstractmethod
    def construct_data(self):
        """
        Map from individual parser format to NXmpes-ready dict.

        Should be implemented by the inheriting mapper.

        """


def convert_snake_to_pascal(str_value: str):
    """Convert snakecase text to pascal case."""
    return str_value.replace("_", " ").title().replace(" ", "")


def convert_pascal_to_snake(str_value: str):
    """Convert pascal case text to snake case."""
    # Convert CamelCase to snake_case
    snake_case = re.sub(r"(?<!^)(?=[A-Z])", "_", str_value)

    # Convert whitespace to underscores and remove extra underscores
    snake_case_cleaned = re.sub(r"\s+", "_", snake_case).replace("__", "_")

    return snake_case_cleaned.lower()


def safe_arange_with_edges(start: float, stop: float, step: float):
    """
    In order to avoid float point errors in the division by step.

    Parameters
    ----------
    start : float
        Smallest value.
    stop : float
        Biggest value.
    step : float
        Step size between points.

    Returns
    -------
    ndarray
        1D array with values in the interval (start, stop),
        incremented by step.

    """
    return step * np.arange(start / step, (stop + step) / step)


def check_uniform_step_width(lst: List[float]):
    """
    Check to see if a non-uniform step width is used in an list.

    Parameters
    ----------
    lst : list
        List of data points.

    Returns
    -------
    bool
        False if list is non-uniformally spaced.

    """
    start = lst[0]
    stop = lst[-1]
    step = get_minimal_step(lst)

    if step != 0.0 and np.abs((stop - start) / step) > len(lst):
        return False
    return True


def get_minimal_step(lst):
    """
    Return the minimal difference between two consecutive values
    in a list. Used for extracting minimal difference in a
    list with non-uniform spacing.

    Parameters
    ----------
    lst : list or np.ndarray
        List of data points.

    Returns
    -------
    step : float
        Non-zero, minimal distance between consecutive data
        points in lst.

    """
    lst1 = np.roll(lst, -1)
    diff = np.abs(np.subtract(lst, lst1))
    step = round(np.min(diff[diff != 0]), 2)

    return step


def _resample_array(y, x0, x1):
    """
    Resample an array (y) which has the same initial spacing
    of another array(x0), based on the spacing of a new
    array(x1).

    Parameters
    ----------
    y : array
        Lineshape array or list.
    x0 : array
        x array with old spacing.
    x1 : array
        x array with new spacing.

    Returns
    -------
    list
        Interpolated y array.

    """
    # pylint: disable=invalid-name
    interp_fn = interp1d(x0, y, axis=0, fill_value="extrapolate")
    return interp_fn(x1)


def interpolate_arrays(x: List[float], array_list: List[np.ndarray]):
    """
    Interpolate data points in case a non-uniform step width was used.

    Parameters
    ----------
    x : list
        List of non-uniformally spaced data points.
    array_list : list
        List of arrays to be interpolated according to new x axis.

    Returns
    -------
    x, array_list
        Interpolated x axis and list of arrays

    """
    # pylint: disable=invalid-name
    if not isinstance(array_list, list):
        array_list = [array_list]
    start = x[0]
    stop = x[-1]
    step = get_minimal_step(x)
    if start > stop:
        # pylint: disable=arguments-out-of-order
        new_x = np.flip(safe_arange_with_edges(stop, start, step))
    else:
        new_x = safe_arange_with_edges(start, stop, step)

    output_list = [_resample_array(arr, x, new_x) for arr in array_list]

    return new_x, output_list


def check_for_allowed_in_list(value, allowed_values: List[Any]):
    """
    Check if a value is a list of values.
    If not, raise Exception.
    """

    if value not in allowed_values:
        raise Exception(f"{value} not in allowed values: {allowed_values}.")
    return value


def re_map_keys(dictionary: Dict[str, Any], key_map: Dict[str, str]):
    """
    Map the keys in a dictionary such that they are replaced by the values
    in key_map and return the dictionary with replaced keys.

    This is often used to map some metadata keys in a vendor file format
    to the common metadata names. used in NXmpes/NXxps.

    Parameters
    ----------
    dictionary : Dict[str, Any]
        Dictionary with XPS metadata.
    key_map : Dict[str, str]
        Mapping of keys from vendor file format to common metadata names.
        Example from the VMS parser:
            key_map = {
               "block_id": "region",
               "sample_id": "sample_name",
               "technique": "analysis_method",
               "source_energy": "excitation_energy",
            }

    Returns
    -------
    dictionary : Dict[str, Any]
        Dictionary with changed keys.

    """
    for k in key_map.keys():
        if k in dictionary:
            dictionary[key_map[k]] = dictionary.pop(k)
    return dictionary


def re_map_values(dictionary: Dict[str, Any], map_functions: Dict[str, Any]):
    """
    Map the values in a dicitionary using functions defined in a
    mapping functions dictionary.

    This is often used to map some metadata in a vendor file format
    to the enumerations used in NXmpes/NXxps.

    Parameters
    ----------
    dictionary : Dict[str, Any]
        Dictionary with XPS metadata.
    map_functions : Dict[str, Any]
       Mapping functions for keys in the dictionary.
       Example from the SLE parser:
           map_functions = {
               "energy/@type": self._change_energy_type,
               "excitation_energy": self._convert_excitation_energy,
               "time_stamp": self._convert_date_time,
               "energy_scan_mode": self._convert_energy_scan_mode,
           }

    Returns
    -------
    dictionary : Dict[str, Any]
        Dictionary with changed values.

    """
    for key, map_fn in map_functions.items():
        if key in dictionary:
            dictionary[key] = map_fn(dictionary[key])
    return dictionary


def drop_unused_keys(dictionary: Dict[str, Any], keys_to_drop: List[str]):
    """
    Remove any keys parsed from sle that are not needed

    Parameters
    ----------
    dictionary : dict
        Dictionary with data and metadata for a spectrum.
    keys_to_drop : list
        List of metadata keys that are not needed.

    Returns
    -------
    None.

    """
    for key in keys_to_drop:
        if key in dictionary:
            dictionary.pop(key)


def update_dict_without_overwrite(d1: Dict[str, Any], d2: Dict[str, Any]):
    """Update d1 with d2, but don't overwrite existing keys."""
    d1.update({k: v for k, v in d2.items() if k not in d1})


def construct_data_key(spectrum: Dict[str, Any]):
    """
    Construct a key for the 'data' field of the xps_dict.
    Output example: cycle0_scan0.

    """
    if "loop_no" in spectrum:
        cycle_key = f'cycle{spectrum["loop_no"]}'
    else:
        cycle_key = "cycle0"

    if "scan_no" in spectrum:
        scan_key = f'scan{spectrum["scan_no"]}'
    else:
        scan_key = "scan0"

    return f"{cycle_key}_{scan_key}"


def construct_detector_data_key(spectrum: Dict[str, Any]):
    """
    Construct a key for the detector data fields of the xps_dict.
    Output example: 'cycles/Cycle_0/scans/Scan_0'

    """
    if "loop_no" in spectrum:
        cycle_key = f'cycles/Cycle_{spectrum["loop_no"]}'
    else:
        cycle_key = "cycles/Cycle_0"

    if "scan_no" in spectrum:
        scan_key = f'scans/Scan_{spectrum["scan_no"]}'
    else:
        scan_key = "scans/Scan_0"

    key = f"{cycle_key}/{scan_key}"

    if "channel_no" in spectrum:
        key += f'/channels/Channel_{spectrum["channel_no"]}'

    return key


KEY_PATTERNS = [
    re.compile(rf"{key_part}(.*?)(?=\/|$)")
    for key_part in ["Group_", "Region_", "RegionData_"]
]


def align_name_part(name_part: str):
    """Make one part of the entry name compliant with NeXus standards."""
    replacements = {
        " ": "_",
        ",": "",
        ".": "_",
    }

    for key, val in replacements.items():
        name_part = name_part.replace(key, val)

    return name_part


def construct_entry_name(key: str):
    """Construct entry name."""
    name_parts = []

    for key_pattern in KEY_PATTERNS:
        match = re.search(key_pattern, key)
        if match:
            name_part = align_name_part(match.group(1))
            name_parts.append(name_part)
    return "__".join(name_parts)
