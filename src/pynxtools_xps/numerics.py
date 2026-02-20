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
Numeric helper functions
"""

from typing import Any

import numpy as np
import pint
from pynxtools.units import ureg
from scipy.interpolate import interp1d

from pynxtools_xps.logging import _logger


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
            _logger.warning(f"WARNING: {error_txt}")


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
        bool: False if list is non-uniformly spaced.

    """
    start = lst[0]
    stop = lst[-1]
    step = _get_minimal_step(lst)

    if step != 0.0 and np.abs((stop - start) / step) > len(lst):
        return False
    return True


def _get_minimal_step(lst: list[float] | np.ndarray) -> float:
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
    interp_fn = interp1d(x0, y, axis=0, fill_value="extrapolate")
    return interp_fn(x1)


def interpolate_arrays(x: list[float], array_list: list[np.ndarray]):
    """
    Interpolate data points in case a non-uniform step width was used.

    Parameters
    ----------
    x : list
        List of non-uniformly spaced data points.
    array_list : list
        List of arrays to be interpolated according to new x axis.

    Returns
    -------
    x, array_list
        Interpolated x axis and list of arrays

    """
    if not isinstance(array_list, list):
        array_list = [array_list]
    start = x[0]
    stop = x[-1]
    step = _get_minimal_step(x)
    if start > stop:
        # pylint: disable=arguments-out-of-order
        new_x = np.flip(safe_arange_with_edges(stop, start, step))
    else:
        new_x = safe_arange_with_edges(start, stop, step)

    output_list = [_resample_array(arr, np.array(x), new_x) for arr in array_list]

    return new_x, output_list


def _check_for_allowed_in_list(value: Any, allowed_values: list[Any]):
    """
    Check if a value is a list of values.
    If not, raise Exception.
    """

    if value not in allowed_values:
        raise Exception(f"{value} not in allowed values: {allowed_values}.")
    return value
