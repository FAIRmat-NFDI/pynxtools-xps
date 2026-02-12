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
Mappings for Specs Lab Prodigy SLE format reader.
"""

import re
from typing import Any

import numpy as np
from lxml import etree as ET

from pynxtools_xps.reader_utils import (
    _format_value,
    _re_map_single_value,
    convert_pascal_to_snake,
    extract_unit,
)
from pynxtools_xps.value_mappers import (
    UNIT_MAP,
    convert_energy_scan_mode,
    convert_energy_type,
    parse_datetime,
)

KEY_MAP: dict[str, str | dict[str, str]] = {
    # SQL spectrum metadata
    "SpectrumID": "spectrum_id",
    "EnergyChns": "energy_channels",
    "NonEnergyChns": "non_energy_channels",
    "Samples": "n_values",
    "ElectronEnergy": "electron_energy",
    "EnergyType": "energy/@type",
    "Step": "step_size",
    "EpassOrRR": "pass_energy_or_retardation_ratio",
    "Ubias": "bias_voltage",
    "Udet": "detector_voltage",
    "Wf": "work_function",
    "Timestamp": "time_stamp",
    # Spectrum group settings
    "ScanMode": {"Name": "energy_scan_mode"},
    "SlitInfo": {"Entrance": "entrance_slit", "Exit": "exit_slit"},
    "Lens": {},
    "EnergyChannelCalibration": {
        "Dir": "calibration_file/dir",
        "File": "calibration_file/path",
    },
    "Transmission": {"File": "transmission_function/file"},
    "Iris": {"Diameter": "iris_diameter"},
    # Spectrum settings
    "Ebin": "binding_energy",
    "Ekin": "kinetic_energy",
    "End": "end_energy",
    "DwellTime": "dwell_time",
    "NumScans": "total_scans",
    "LensMode": "lens_mode",
    "Timestamp": "time_stamp",
    "Entrance": "entrance_slit",
    "Exit": "exit_slit",
    "Epass": "pass_energy",
    "VoltageRange": "voltage_energy_range",
    # spectrometer settings
    "Coil Current [mA]": "coil_current [mA]",
    "Pre Defl Y [nU]": "pre_deflector_y_current [nU]",
    "Pre Defl X [nU]": "pre_deflector_x_current [nU]",
    "L1 [nU]": "lens1_voltage [nU]",
    "L2 [nU]": "lens2_voltage [nU]",
    "Focus Displacement 1 [nu]": "focus_displacement_current [nU]",
    "Detector Voltage [V]": "detector_voltage [V]",
    "Bias Voltage Electrons [V]": "bias_voltage_electrons [V]",
    "Bias Voltage Ions [V]": "bias_voltage_ions [V]",
    # source settings
    "anode": "source_label",
    "uanode": "source_voltage",
    "iemission": "emission_current",
    "ihv": "source_high_voltage",
    "ufilament": "filament_voltage",
    "ifilament": "filament_current",
    "DeviceExcitationEnergy": "excitation_energy",
    "panode": "anode_power",
    "temperature": "source_temperature",
    # scan metadata
    "Node": "node_id",
    "ScanDate": "timestamp",
    "Eexc": "excitation_energy",
    "DeviceExcitationEnergy": "excitation_energy",
    "Trace": "trace",
    "RawID": "raw_id",
    "Channel": "channel",
    # transmission_data_map
    "Data": "data",
}

VALUE_MAP: dict[str, Any] = {
    "energy/@type": convert_energy_type,
    "excitation_energy": float,
    "time_stamp": (parse_datetime, {"possible_date_formats": ["%Y-%b-%d %H:%M:%S.%f"]}),
    "energy_scan_mode": convert_energy_scan_mode,
}

KEYS_TO_DROP: list[str] = [
    "Work Function",
]

UNITS: dict[str, str] = {
    "work_function": "eV",
    "excitation_energy": "eV",
    "iris_diameter": "mm",
    "step_size": "eV",
    "detector_voltage": "V",
    "dwell_time": "s",
    "raw_data/raw": "counts_per_second ",
    "polar_angle": "degree ",
    "azimuth_angle": "degree",
    "pass_energy": "eV",
    "start_energy": "eV",
    "emission_current": "A",
    "source_voltage": "V",
    "transmission_function/kinetic_energy": "eV",
}

POSSIBLE_DATE_FORMATS: list[str] = ["%Y-%b-%d %H:%M:%S.%f"]

_UNIT_IN_KEY_RE = re.compile(r"^(?P<key>.+?)\s*\[(?P<unit>[^\]]+)\]$")


def extract_unit_from_key(
    key: str,
    unit_map: dict[str, str | None],
) -> tuple[str, str | None]:
    """
    Extract and normalize a unit embedded in a key name.

    The function detects units of the form ``"<key> [<unit>]"`` where the unit
    is enclosed in square brackets at the end of the key. If found, the unit
    is removed from the key and normalized using ``unit_map``.

    Special handling is applied for ambiguous units such as ``"nU"``, where
    the resolved unit depends on the semantic meaning of the key.

    Args:
        key (str): Formatted key name, potentially containing an embedded unit.
        unit_map (dict[str, str | None]): Mapping from raw units to normalized
            units. A mapped value of ``None`` explicitly disables the unit.

    Returns:
        tuple[str, str | None]:
            - Key with any embedded unit removed.
            - Normalized unit, or None if no unit was found or it was disabled.
    """
    match = _UNIT_IN_KEY_RE.match(key)
    if not match:
        return key, None

    base_key = match.group("key")
    raw_unit = match.group("unit")

    # Exceptional case: ambiguous "nU" unit
    if raw_unit == "nU":
        if base_key.endswith("voltage"):
            return base_key, "V"
        if base_key.endswith("current"):
            return base_key, "A"

    # Default normalization via unit map
    unit = unit_map.get(raw_unit, raw_unit)

    return base_key, unit


def format_key_value_and_unit(
    key: str, value_str: str
) -> tuple[str, str | int | float | bool | np.ndarray | None, str | None]:
    """
    Formats a key and value string and extracts an associated unit.

    This reader-level helper:
    1. Maps the key using `KEY_MAP` and converts it to snake_case.
    2. Extracts the numeric value and unit from the value string.
    3. Falls back to extracting the unit from the key if not present in the value.
    4. Formats the numeric part of the value (int / float when possible).
    5. Remaps the value if specified in `VALUE_MAP`.

    Args:
        key (str): Raw key name as found in the input.
        value_str (str): Raw value string, potentially containing a unit.

    Returns:
        tuple[str, int | float | str, str | None]:
            - Formatted key (snake_case, mapped if applicable).
            - Formatted value (numeric when possible, otherwise unchanged).
            - Extracted unit, or None if no unit could be determined.
    """
    key = KEY_MAP.get(key, convert_pascal_to_snake(key))  # type: ignore[assignment]

    value, unit = extract_unit(key, value_str, UNITS)

    if not unit:
        key, unit = extract_unit_from_key(key, UNIT_MAP)

    value = _format_value(value)
    value = _re_map_single_value(key, value, VALUE_MAP)  # type: ignore[assignment]

    return key, value, unit


def iterate_xml_at_tag(xml_elem: ET.Element, tag: str) -> dict[str, str | float | int]:
    """
    Iterates through XML elements at the specified tag and formats their attributes.

    Parameters
    ----------
    xml_elem : ET.Element
        The XML element to search within.

    tag : str
        The tag name to find in the XML structure.

    Returns
    -------
    dict[str, Union[str, float, int]]
        A dictionary containing formatted attribute values keyed by their corresponding names.
    """

    sub_elem = xml_elem.find(tag)

    settings: dict[str, Any] = {}

    special_key_map: str | dict[str, str] = KEY_MAP.get(tag, {})

    if sub_elem is not None and isinstance(special_key_map, dict):
        for param in sub_elem.iter():
            for key, value in param.attrib.items():
                key = special_key_map.get(key, key)
                key, value, unit = format_key_value_and_unit(key, value)
                settings[key] = value
                if unit:
                    settings[f"{key}/@units"] = unit

    return settings
