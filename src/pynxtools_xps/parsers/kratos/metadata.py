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
Metadata mapping for the Kratos parser.
"""

import re
from functools import partial
from typing import Any

from pynxtools_xps.mapping import (
    _convert_bool,
    _MetadataContext,
    _ValueMap,
    parse_datetime,
)

_POSSIBLE_DATE_FORMATS: list[str] = ["%d.%m.%Y %H:%M", "%d/%m/%Y %H:%M"]

_DESCRIPTION_PATTERN: re.Pattern = re.compile(
    r"\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)\s*([a-zA-Z]+)"
)


def _convert_description(value: str) -> dict[str, Any]:
    """Map all items in description to a dictionary."""

    match = _DESCRIPTION_PATTERN.match(value)

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


_KEY_MAP: dict[str, str] = {
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

_VALUE_MAP: _ValueMap = {
    "date_created": partial(
        parse_datetime,
        possible_date_formats=_POSSIBLE_DATE_FORMATS,
    ),
    "description": _convert_description,
    "charge_neutraliser": _convert_bool,
    "deflection": _convert_xray_deflection,
}

# TODO: what to do with these, should be handled by generic format_key_value_and_unit
# KEYS_WITH_UNITS: list[str] = [
#     "filament_current",
#     "filament_bias",
#     "charge_balance",
#     "emission_current",
#     "energy_start",
#     "energy_end",
#     "energy_centre",
#     "energy_width",
#     "step_size",
#     "dwell_time",
#     "sweep_time",
#     "x_ray_power",
# ]

_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {
    "sample_tilt": "degree",
    "resolution": "eV",
}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)

# TODO: is this covered already?
# _UNIT_PATTERN  = re.compile(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)([a-zA-Z]+)")
