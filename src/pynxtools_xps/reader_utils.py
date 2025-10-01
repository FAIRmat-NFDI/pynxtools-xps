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
Helper functions for populating NXmpes template
"""

import logging
import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
import pint
from pynxtools.units import ureg
from scipy.interpolate import interp1d

logger = logging.getLogger(__name__)


@dataclass
class XpsDataclass:
    """Generic class to hold a data model and a type validation method."""

    def validate_types(self):
        ret = True
        for field_name, field_def in self.__dataclass_fields__.items():
            actual_type = type(getattr(self, field_name))
            if actual_type != field_def.type:
                logger.warning(
                    f"Type mismatch in dataclass {type(self).__name__}. {field_name}: '{actual_type}' instead of '{field_def.type}'"
                )
                ret = False
        return ret

    def __post_init__(self):
        if not self.validate_types():
            raise ValueError(f"Type mismatch in dataclass {type(self).__name__}")

    def dict(self):
        return self.__dict__.copy()


class XPSMapper(ABC):
    """Abstract base class from mapping from a parser to NXmpes template"""

    def __init__(self):
        self.file: str | Path = ""
        self.raw_data: list[str] = []
        self._xps_dict: dict[str, Any] = {}

        self.parser = None

    @abstractmethod
    def _select_parser(self):
        """
        Select the correct parser for the file extension and format.

        Should be implemented by the inheriting mapper.
        """

    @property
    def data_dict(self) -> dict:
        """Getter property."""
        return self._xps_dict

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the file using the Scienta TXT parser.

        """
        self.file = file
        self.parser = self._select_parser()
        self.raw_data = self.parser.parse_file(file, **kwargs)

        self._xps_dict["File"] = file
        self._xps_dict["file_ext"] = os.path.splitext(file)[1]

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


def safe_arange_with_edges(start: float, stop: float, step: float) -> np.ndarray:
    """
    In order to avoid float point errors in the division by step.

    Args:
        start (float): Lower limit.
        stop (float): Upper limit.
        step (float): Step size between points.
    Returns:
        ndarray
            1D array with values in the interval (start, stop),
            incremented by step.

    """
    return step * np.arange(start / step, (stop + step) / step)


def check_uniform_step_width(lst: list[float]) -> bool:
    """
    Check to see if a non-uniform step width is used in an list.

    Args:
        lst (list): List of data points.

    Returns:
        bool: False if list is non-uniformally spaced.

    """
    start = lst[0]
    stop = lst[-1]
    step = get_minimal_step(lst)

    if step != 0.0 and np.abs((stop - start) / step) > len(lst):
        return False
    return True


def get_minimal_step(lst: list[float] | np.ndarray) -> float:
    """
    Return the minimal difference between two consecutive values
    in a list. Used for extracting minimal difference in a
    list with non-uniform spacing.

    Args:
        lst (list): List of data points.

    Returns:
        step (float): Non-zero, minimal distance between consecutive data
        points in lst.

    """
    lst1 = np.roll(lst, -1)
    diff = np.abs(np.subtract(lst, lst1))
    step = round(np.min(diff[diff != 0]), 2)

    return step


def _resample_array(y: np.ndarray, x0: np.ndarray, x1: np.ndarray) -> np.ndarray:
    """
    Resample an array (y) which has the same initial spacing
    of another array(x0), based on the spacing of a new
    array(x1).

    Args:
        y (np.ndarray): Lineshape array or list.
        x0 (np.ndarray): x array with old spacing.
        x1 (np.ndarray): x array with new spacing.

    Returns:
        np.ndarray: Interpolated y array.

    """
    # pylint: disable=invalid-name
    interp_fn = interp1d(x0, y, axis=0, fill_value="extrapolate")
    return interp_fn(x1)


def interpolate_arrays(x: list[float], array_list: list[np.ndarray]):
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

    output_list = [_resample_array(arr, np.array(x), new_x) for arr in array_list]

    return new_x, output_list


def check_for_allowed_in_list(value: Any, allowed_values: list[Any]):
    """
    Check if a value is a list of values.
    If not, raise Exception.
    """

    if value not in allowed_values:
        raise Exception(f"{value} not in allowed values: {allowed_values}.")
    return value


def re_map_keys(dictionary: dict[str, Any], key_map: dict[str, str]) -> dict[str, Any]:
    """
    Map the keys in a dictionary such that they are replaced by the values
    in key_map and return the dictionary with replaced keys.

    This is often used to map some metadata keys in a vendor file format
    to the common metadata names. used in NXmpes/NXxps.

    Args:
        dictionary (dict[str, Any]): Dictionary with XPS metadata.
        key_map (dict[str, str]): Mapping of keys from vendor file format
        to common metadata names.

        Example from the VMS parser:
            key_map = {
               "block_id": "region",
               "sample_id": "sample_name",
               "technique": "analysis_method",
               "source_energy": "excitation_energy",
            }

    Returns:
        dictionary (dict[str, Any]): Dictionary with changed keys.

    """
    for key in key_map:
        if key in dictionary:
            map_key = key_map.get(key, key)
            dictionary[map_key] = dictionary.pop(key)
    return dictionary


def re_map_values(
    dictionary: dict[str, Any], map_functions: dict[str, Any]
) -> dict[str, Any]:
    """
    Map the values in a dicitionary using functions defined in a
    mapping functions dictionary.

    This is often used to map some metadata in a vendor file format
    to the enumerations used in NXmpes/NXxps.

    Args:
        dictionary (dict[str, Any]): Dictionary with XPS metadata.
        map_functions (dict[str, Any]): Mapping functions for keys in the dictionary.

        Example from the SLE parser:
        map_functions = {
            "energy/@type": self._change_energy_type,
            "excitation_energy": self._convert_excitation_energy,
            "time_stamp": self._convert_date_time,
            "energy_scan_mode": self._convert_energy_scan_mode,
        }

    Returns:
        dictionary (dict[str, Any]): Dictionary with changed values.

    """
    for key, map_fn in map_functions.items():
        if key in dictionary:
            dictionary[key] = map_fn(dictionary[key])
    return dictionary


def _re_map_single_value(
    input_key: str,
    value: str | int | float | bool | np.ndarray | None,
    map_functions: dict[str, Any],
    **kwargs,
) -> str | int | float | bool | np.ndarray | None:
    """
    Map the values returned from the file to the preferred format for
    the parser output.

    """
    if isinstance(value, str) and value is not None:
        value = value.rstrip("\n")

    for key, map_method in map_functions.items():
        if key == input_key:
            map_method = map_functions[key]
            return map_method(value, **kwargs) if kwargs else map_method(value)  # type: ignore[operator]

    return value


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


def drop_unused_keys(dictionary: dict[str, Any], keys_to_drop: list[str]) -> None:
    """
    Remove unwanted keys from a dictionary.

    Args:
        dictionary (dict[str, Any]):
            Dictionary containing data and metadata for a spectrum.
        keys_to_drop (list[str]):
            List of keys that should be removed from the dictionary.

    Returns:
        None
    """
    for key in keys_to_drop:
        dictionary.pop(key, None)


def update_dict_without_overwrite(d1: dict[str, Any], d2: dict[str, Any]):
    """Update d1 with d2, but don't overwrite existing keys."""
    d1.update({k: v for k, v in d2.items() if k not in d1})


def construct_data_key(spectrum: dict[str, Any]) -> str:
    """
    Construct a key for the 'data' field of the xps_dict.
    Output example: cycle0_scan0.

    """
    if "loop_no" in spectrum:
        cycle_key = f"cycle{spectrum['loop_no']}"
    else:
        cycle_key = "cycle0"

    if "scan_no" in spectrum:
        scan_key = f"scan{spectrum['scan_no']}"
    else:
        scan_key = "scan0"

    return f"{cycle_key}_{scan_key}"


def align_name_part(name_part: str):
    """Make one part of the entry name compliant with NeXus standards."""
    translation_table = str.maketrans(
        {
            " ": "_",
            ",": "",
            ".": "_",
            "-": "_",
            ":": "_",
            "+": "_",
            "/": "_",
            "=": "",
        }
    )

    return name_part.translate(translation_table)


def construct_entry_name(parts: list[str]) -> str:
    """Construct name for the NXentry instances."""
    if len(parts) == 1:
        return align_name_part(parts[0])
    return "__".join([align_name_part(part) for part in parts])


def _format_value(value: int | float | str) -> int | float | str:
    """
    Formats the input value as an int or float if it's a numeric string.

    If a string represents a float (e.g., "5.0"), it remains a float even if it
    has no decimals.

    Args:
        value (Union[str, float, int]): The input value to format.

    Returns:
        Union[int, float, str]: The formatted value as int or float if numeric;
                                otherwise, returns the original value.

    """
    if isinstance(value, str) and re.match(r"^-?\d*\.?\d+(?:[eE][-+]?\d+)?$", value):
        # Check for decimal to ensure float, even if the decimal part is zero
        return float(value) if "." in value or "e" in value.lower() else int(value)
    return value


UNIT_PATTERN = re.compile(r"^([-+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([a-zA-Z/\s]+)$")


def split_value_and_unit(
    value_str: str,
) -> tuple[int | float | str, str]:
    """
    Splits a string into a numerical value and its associated unit.

    If the string contains a pattern where a number (which can be negative and
    in scientific notation) is followed by a unit, it returns a tuple containing
    the value (as a float) and the unit (as a string). If the string does not match
    this pattern, it returns the original string.

    Args:
        value_str (str): The input string to split.

    Returns:
        Tuple[Union[int, float, str], str]:
            - (value, unit) if a numeric value is detected.
            - (original string, "") if no numeric value is detected.

    """
    match = UNIT_PATTERN.match(value_str)
    if match:
        value = _format_value(match.group(1))
        unit = match.group(2).replace(" ", "") if match.group(2) else ""
        return value, unit
    return _format_value(value_str), ""


def extract_unit(
    key: str, value_str: str, unit_missing: dict[str, str] | None = None
) -> tuple[int | float | str, str]:
    """
    Extract a numeric value and its associated unit from a metadata string.

    The function identifies and separates numerical and unit components from
    the `value_str` string. If no unit is found in `value_str`, it checks the
    `unit_missing` dictionary for a default unit.

    Example:
        analyzer_work_function = "4.506eV"
        -> (4.506, "eV")

    Args:
        key (str): Key associated with the value.
        value_str (str): Combined numeric value and unit as a string.
        unit_missing (Optional[dict[str, str]]): Optional dictionary with default units
            for keys missing units. Defaults to None.

    Returns:
        tuple[Union[int, float, str], str]:
            - A tuple with the numeric value (int, float, or str) and the unit.
            - If no unit is found in `value_str`, it uses `unit_missing` if available.
            - If no unit is found in either, returns an empty string for the unit.

    """
    value, unit = split_value_and_unit(value_str)

    if not unit:
        unit = unit_missing.get(key, "")

    return value, unit


def check_units(template_path: str, unit: str) -> None:
    """
    Check that the unit is a valid pint unit.

    Args:
        template_path (str): Path of a Template object.
        unit (str): String representation of a unit.

    """
    error_txt = f"Invalid unit '{unit}' at path: {template_path}"

    if unit is not None:
        error_txt = f"Invalid unit '{unit}' at path: {template_path}"
        try:
            ureg.Unit(unit)
        except pint.errors.UndefinedUnitError as pint_err:
            raise pint.errors.UndefinedUnitError(error_txt) from pint_err
        except TypeError as type_err:
            logger.warning(f"WARNING: {error_txt}")
