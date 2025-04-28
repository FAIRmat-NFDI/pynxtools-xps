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
# pylint: disable=too-many-lines,too-few-public-methods
"""
Utility function for mapping keys and values in the pynxtools template.
"""

import re
import datetime
from typing import Dict, List, Any

ENERGY_TYPE_MAP = {
    "BE": "binding",
    "KE": "kinetic",
    "Binding": "binding",
    "Binding Energy": "binding",
    "binding energy": "binding",
    "Kinetic": "kinetic",
    "Kinetic Energy": "kinetic",
    "kinetic energy": "kinetic",
    "Analyser Energy": "binding",
}

ENERGY_SCAN_MODE_MAP = {
    "Fixed": "fixed_energy",
    "fixed": "fixed_energy",
    "FixedEnergies": "fixed_energy",
    "Swept": "fixed_analyzer_transmission",
    "FixedAnalyzerTransmission": "fixed_analyzer_transmission",
    "FAT": "fixed_analyzer_transmission",
    "FixedRetardationRatio": "fixed_retardation_ratio",
    "FRR": "fixed_retardation_ratio",
    "Snapshot": "snapshot",
    "SnapshotFAT": "snapshot",
}

MEASUREMENT_METHOD_MAP = {
    "XPS": "X-ray photoelectron spectroscopy (XPS)",
    "UPS": "ultraviolet photoelectron spectroscopy (UPS)",
    "ElectronSpectroscopy": "electron spectroscopy for chemical analysis (ESCA)",
    "NAPXPS": "near ambient pressure X-ray photoelectron spectroscopy (NAPXPS)",
    "ARXPS": "angle-resolved X-ray photoelectron spectroscopy (ARXPS)",
}

BOOL_MAP = {
    "yes": True,
    "Yes": True,
    "no": False,
    "No": False,
    "On": True,
    "Off": False,
}

UNIT_MAP = {
    "Counts": "counts_per_second",
    "counts/s": "counts_per_second",
    "CPS": "counts_per_second",
    "u": "um",
    "KV": "kV",
    "seconds": "s",
    "(min)": "min",
    "Percent": "",  # should be changed back to percent once pint is updated
    "atom": "dimensionless",
    "eV/atom": "eV",
    "microm": "micrometer",
    "d": "degree",
    "nU": "V",  # workaround for SPECS SLE reader
    "s-1": "1/s",  # workaround for SPECS XY reader
    "norm": None,  # workaround for SPECS XY reader
}


def _replace_from_map(value: Any, value_map: Dict[str, Any]):
    """
    For a given value, return a new value if the value is
    part of the value_map.
    """
    if value in value_map:
        return value_map[value]
    return value


def convert_energy_type(energy_type: str):
    """
    Change the strings for energy type to the allowed
    values in NXmpes.

    """
    return _replace_from_map(energy_type, ENERGY_TYPE_MAP)


def convert_energy_scan_mode(energy_scan_mode: str):
    """
    Change the strings for energy scan mode to
    the allowed values in NXmpes.

    """
    return _replace_from_map(energy_scan_mode, ENERGY_SCAN_MODE_MAP)


def convert_measurement_method(measurement_method: str):
    """
    Change the strings for measurement method to
    the allowed values in NXmpes.

    """
    return _replace_from_map(measurement_method, MEASUREMENT_METHOD_MAP)


def convert_bool(bool_like: str):
    """Convert "yes", "no" to actual boooleans."""
    return _replace_from_map(bool_like, BOOL_MAP)


def convert_units(units: str):
    """Map y_units to shortened values."""
    return _replace_from_map(units, UNIT_MAP)


def get_units_for_key(unit_key: str, unit_map: Dict[str, str]) -> str:
    """
    Get correct units for a given key from a dictionary with unit map.
    Parameters
    ----------
    unit_key : str
       Key of type <mapping>:<spectrum_key>, e.g.
       detector/detector_voltage
    Returns
    -------
    str
        Unit for that unit_key.
    """
    regex_match = re.search(r"\[([A-Za-z0-9_]+)\]", unit_key)
    if regex_match is None:
        return unit_map.get(unit_key, None)
    return regex_match.group(1)


def parse_datetime(
    datetime_string: str,
    possible_date_formats: List[str],
    tzinfo: datetime.tzinfo = datetime.timezone.utc,
) -> str:
    """
    Convert a date string to ISO 8601 format with optional timezone handling.

    Convert the native time format to the datetime string
    in the ISO 8601 format: '%Y-%b-%dT%H:%M:%S.%fZ'.
    For different vendors, there are different possible date formats,
    all of which can be checked with this method.
    Optionally, a timezone (tzinfo) can be applied to the datetime object if provided.

    Parameters
    ----------
    datetime_string : str
        String representation of the date
    possible_date_formats : List[str]
        List of possible date time formats to attempt for parsing.
    tzinfo: datetime.tzinfo
        A tzinfo object specifying the desired timezone to apply to the datetime object.
        Defaults to UTC (datetime.timezone.utc).

    Raises
    ------
    ValueError
        If the time format cannot be converted, a ValueError is raised.

    Returns
    -------
    str
        Datetime in ISO 8601 format.
    """
    for date_fmt in possible_date_formats:
        try:
            datetime_obj = datetime.datetime.strptime(datetime_string, date_fmt)

            if tzinfo is not None:
                # Apply the specified timezone to the datetime object
                datetime_obj = datetime_obj.replace(tzinfo=tzinfo)

            # Convert to ISO 8601 format
            return datetime_obj.isoformat()

        except ValueError:
            continue

    raise ValueError(
        f"Datetime {datetime_string} could not be converted to ISO 8601 format, "
        f"as it does not match any of these formats: {possible_date_formats}."
    )
