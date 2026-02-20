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

import datetime
import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, TypeAlias

import numpy as np

BOOL_MAP: dict[str, bool] = {
    "yes": True,
    "Yes": True,
    "no": False,
    "No": False,
    "On": True,
    "Off": False,
}

MEASUREMENT_METHOD_MAP: dict[str, tuple[str, str]] = {
    "XPS": ("XPS", "X-ray photoelectron spectroscopy"),
    "UPS": ("UPS", "ultraviolet photoelectron spectroscopy"),
    "ESCA": ("XPS", "electron spectroscopy for chemical analysis"),
    "NAPXPS": ("NAPXPS", "near ambient pressure X-ray photoelectron spectroscopy"),
    "ARXPS": ("ARXPS", "angle-resolved X-ray photoelectron spectroscopy"),
}

ENERGY_TYPE_MAP: dict[str, str] = {
    "BE": "binding",
    "KE": "kinetic",
    "Binding": "binding",
    "Binding Energy": "binding",
    "binding energy": "binding",
    "binding_energy": "binding",
    "Kinetic": "kinetic",
    "Kinetic Energy": "kinetic",
    "kinetic energy": "kinetic",
    "kinetic_energy": "kinetic",
    "Analyser Energy": "binding",
    "analyser_energy": "binding",
}

ENERGY_SCAN_MODE_MAP: dict[str, str] = {
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

ACQUSITION_MODE_MAP: dict[str, str] = {
    "Image": "pulse counting",
    "Events": "pulse counting",
}

SLIT_TYPE_MAP: dict[str, str] = {"Straight": "straight slit", "Curved": "curved slit"}


def _replace_from_map(value: Any, value_map: dict[str, Any]):
    """
    For a given value, return a new value if the value is
    part of the value_map.
    """
    return value_map.get(value, value)


def _make_converter(value_map: dict[str, Any]) -> Callable[[Any], Any]:
    return lambda value: _replace_from_map(value, value_map)


_convert_bool = _make_converter(BOOL_MAP)
_convert_detector_acquisition_mode = _make_converter(ACQUSITION_MODE_MAP)
_convert_energy_scan_mode = _make_converter(ENERGY_SCAN_MODE_MAP)
_convert_energy_type = _make_converter(ENERGY_TYPE_MAP)
_convert_measurement_method = _make_converter(MEASUREMENT_METHOD_MAP)
_convert_slit_type = _make_converter(SLIT_TYPE_MAP)


def parse_datetime(
    datetime_string: str,
    possible_date_formats: list[str],
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
        if date_fmt == "%Y-%m-%dT%H:%M:%S.%f%z":
            # strptime only supports six digits for microseconds
            datetime_string = datetime_string[:-7] + datetime_string[-6:]

        try:
            datetime_obj = datetime.datetime.strptime(datetime_string, date_fmt)

            if tzinfo is not None:
                datetime_obj = datetime_obj.replace(tzinfo=tzinfo)

            return datetime_obj.isoformat()

        except ValueError:
            continue

    raise ValueError(
        f"Datetime {datetime_string} could not be converted to ISO 8601 format, "
        f"as it does not match any of these formats: {possible_date_formats}."
    )


def convert_snake_to_pascal(value: str):
    """Convert snakecase text to pascal case."""
    return value.replace("_", " ").title().replace(" ", "")


def convert_pascal_to_snake(value: str) -> str:
    """Convert PascalCase text to snake_case, preserving bracketed content."""

    def replace_non_bracketed(match):
        content = match.group(0)

        # If already ALL CAPS (with optional underscores), just lowercase
        if content.isupper():
            return content.lower()

        # split acronym -> word (TFCParameters -> TFC_Parameters)
        snake_case = re.sub(
            r"(?<=[A-Z])(?=[A-Z][a-z])",
            "_",
            content,
        )

        # split camel/pascal case (myVariable -> my_Variable)
        snake_case = re.sub(
            r"(?<=[a-z0-9])(?=[A-Z])",
            "_",
            snake_case,
        )

        # normalize separators
        snake_case = re.sub(r"[\s\-]+", "_", snake_case)

        # collapse duplicate underscores
        return re.sub(r"_+", "_", snake_case)

    pattern = r"(\[.*?\]|[^[]+)"
    parts = re.sub(
        pattern,
        lambda m: (
            m.group(0) if m.group(0).startswith("[") else replace_non_bracketed(m)
        ),
        value,
    )

    return parts.lower()


# TODO: is datetime proper here?
_Value: TypeAlias = str | int | float | np.ndarray | datetime.datetime | None
_ValueMap: TypeAlias = dict[str, Callable[[Any], Any]]

_UNIT_IN_KEY_RE: re.Pattern[str] = re.compile(r"^(?P<key>.+?)\s*\[(?P<unit>[^\]]+)\]$")
_UNIT_IN_VALUE_PATTERN: re.Pattern[str] = re.compile(
    r"^([-+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([a-zA-Z/\s]+)$"
)


def _split_key_and_unit(
    key: str,
) -> tuple[str, str | None]:
    """
    Extract a unit embedded in a key name.

    The function detects units of the form ``"<key> [<unit>]"`` where the unit
    is enclosed in square brackets at the end of the key. If found, the unit
    is separated from the key.

    Special handling is applied for ambiguous units such as ``"nU"``, where
    the resolved unit depends on the semantic meaning of the key.

    Args:
        key (str): Formatted key name, potentially containing an embedded unit.

    Returns:
        tuple[str, str | None]:
            - Key with any embedded unit removed.
            - Normalized unit, or None if no unit was found or it was disabled.
    """
    match = _UNIT_IN_KEY_RE.match(key)
    if not match:
        return key, None

    base_key = match.group("key").strip()
    unit = match.group("unit")

    # Exceptional case: ambiguous "nU" unit
    if unit == "nU":
        if base_key.endswith("voltage"):
            return base_key, "V"
        if base_key.endswith("current"):
            return base_key, "A"

    return base_key, unit


def _split_value_and_unit(
    value: str,
) -> tuple[_Value, str | None]:
    """
    Splits a Value string into a numerical value and its associated unit.

    If the string contains a pattern where a number (which can be negative and
    in scientific notation) is followed by a unit, it returns a tuple containing
    the value and the unit (as a string). If the string does not match
    this pattern, it returns the original string.

    Args:
        value_str (str): The input string to split.

    Returns:
        Tuple[Union[int, float, str], str]:
            - (value, unit) if a numeric value is detected.
            - (original string, "") if no numeric value is detected.

    """
    match = _UNIT_IN_VALUE_PATTERN.match(value)
    if match:
        value = match.group(1)
        unit = match.group(2).replace(" ", "") if match.group(2) else ""
        return value, unit
    return value, None


@dataclass(slots=True)
class _MetadataContext:
    """
    Context object encapsulating key, value, and unit normalization logic.

    The class defines a deterministic normalization pipeline applied to
    metadata extracted from vendor files. The processing order is:

        1. Normalize the key using ``key_map`` or snake_case conversion.
        2. Extract value and unit from the value string.
        3. If no unit is present, extract a unit embedded in the key.
        4. If still undefined, assign a unit from ``default_units``.
        5. Normalize the unit using ``unit_map``.
        6. Apply value mapping functions defined in ``value_map``.
        7. Normalize numeric values (int/float conversion).

    The class contains no runtime state beyond configuration and can be
    reused across multiple parsing operations.
    """

    key_map: dict[str, str]
    value_map: _ValueMap
    unit_map: dict[str, str | None]
    default_units: dict[str, str]

    def normalize_key(self, key: str) -> str:
        """Map key and convert to snake_case if unmapped."""
        key = convert_pascal_to_snake(key)
        return self.key_map.get(key, key)

    def parse_value_and_unit(
        self,
        key: str,
        value: _Value,
    ) -> tuple[_Value, str | None]:
        """
        Extract value and unit from input.

        Unit extraction is only attempted if the value is a string.
        Non-string values pass through unchanged.
        """

        if value is None:
            return value, None

        if isinstance(value, str):
            if not value:
                return "", None

            value, unit = _split_value_and_unit(value)
            return value, unit

        # already typed value → no unit extraction possible
        return value, None

    def resolve_unit_from_key(
        self,
        key: str,
        unit: str | None,
    ) -> tuple[str, str | None]:
        """
        Extract unit from key if not present in value (fallback that is only used
        if unit is not in value)
        """

        if unit is not None:
            return key, unit

        return _split_key_and_unit(key)

    def _format_value(self, value: _Value) -> _Value:
        """
        Formats the input value as an int or float if it's a numeric string.

        If a string represents a float (e.g., "5.0"), it remains a float even if it
        has no decimals.

        Args:
            value (Value): The input value to format.

        Returns:
            Value: The formatted value as int or float if numeric;
            otherwise, returns the original value.

        """
        if isinstance(value, str) and re.match(
            r"^-?\d*\.?\d+(?:[eE][-+]?\d+)?$", value
        ):
            # Check for decimal to ensure float, even if the decimal part is zero
            return float(value) if "." in value or "e" in value.lower() else int(value)
        return value

    def map_value(self, key: str, value: _Value) -> _Value:
        """
        Apply value mapping functions if defined.

        This is often used to map some metadata in a vendor file format
        to the enumerations used in NXmpes/NXxps.

        self.value_map (dict[str, Any]): Mapping functions for keys in the dictionary.

        Example from the SLE parser:
        self.value_map = {
            "energy/@type": self._change_energy_type,
            "excitation_energy": self._convert_excitation_energy,
            "time_stamp": self._convert_date_time,
            "energy_scan_mode": self.__convert_energy_scan_mode,
        }
        """
        # TODO: what is this needed, was in reader_utils._re_map_single_value
        # if isinstance(value, str) and value is not None:
        #     value = value.rstrip("\n")

        map_fn: Callable[[_Value], _Value] | None = self.value_map.get(key)
        if map_fn is None:
            return value
        return map_fn(value)

    def format(
        self,
        key: str,
        value: _Value,
    ) -> tuple[str, _Value, str | None]:

        key = self.normalize_key(key)

        value, unit = self.parse_value_and_unit(key, value)

        key, unit = self.resolve_unit_from_key(key, unit)

        if unit is None:
            unit = self.default_units.get(key)

        unit = self.unit_map.get(unit, unit) if unit else None

        value = self.map_value(key, value)
        value = self._format_value(value)

        return key, value, unit


# TODO: this is broken, needs to be properly fixed because it is used
# def _setattr_with_datacls_type(self, datacls: dataclass, attr: str, value: _Value):
#     """
#     Set attribute of dataclass instance with the field_type
#     defined in the dataclass.
#     """
#     field_type = type(getattr(datacls, attr))
#     setattr(datacls, attr, field_type(value))


# ToDO: is the regex still needed? If not, remove entirely
# def get_units_for_key(unit_key: str, unit_map: dict[str, str]) -> str:
#     """
#     Get correct units for a given key from a dictionary with unit map.
#     Parameters
#     ----------
#     unit_key : str
#        Key of type <mapping>:<spectrum_key>, e.g.
#        detector/detector_voltage
#     Returns
#     -------
#     str
#         Unit for that unit_key.
#     """
#     regex_match = re.search(r"\[([A-Za-z0-9_]+)\]", unit_key)
#     if regex_match is None:
#         return unit_map.get(unit_key, None)
#     return regex_match.group(1)


# ToDO: these should be not needed eventually
def _re_map_keys(
    dictionary: dict[str, Any], context: _MetadataContext
) -> dict[str, Any]:
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
    for key in list(dictionary.keys()):
        new_key = context.normalize_key(key)
        if new_key != key:
            dictionary[new_key] = dictionary.pop(key)

    return dictionary


def _re_map_values(
    dictionary: dict[str, Any], context: _MetadataContext
) -> dict[str, Any]:
    """
    Map the values in a dictionary using _MetadataContext mapper. defined in a
    mapping functions dictionary.

    Returns:
        dictionary (dict[str, Any]): Dictionary with changed values.

    """
    for key, value in dictionary.items():
        dictionary[key] = context.map_value(key, value)
    return dictionary
