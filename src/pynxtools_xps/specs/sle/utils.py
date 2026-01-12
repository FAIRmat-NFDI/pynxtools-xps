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

from typing import Any

from lxml import etree as ET

from pynxtools_xps.reader_utils import (
    _format_value,
    _re_map_single_value,
    convert_pascal_to_snake,
    extract_unit,
)
from pynxtools_xps.value_mappers import (
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
    "VoltageRange": "voltage_range",
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


def format_key_and_value(key: str, value_str: str) -> tuple[str, Any]:
    """
    Formats a key and a corresponding value string according to a series of transformations.

    This function:
    1. Maps the key based on a predefined dictionary (`KEY_MAP`).
    2. Converts the key from PascalCase to snake_case.
    3. Extracts the numeric value and unit from the value string.
    4. Formats the numeric part of the value according to its expected type.
    5. Remaps the value to a new format if specified in `VALUE_MAP`.

    Args:
        key (str): The key associated with the value, which may need mapping and formatting.
        value_str (str): The value string to format and separate into numeric value and unit.

    Returns:
        tuple[Any, str]:
            - The formatted key (converted to snake_case and remapped if needed).
            - The formatted value, with numeric value processed and remapped according to `VALUE_MAP`.
    """
    key = KEY_MAP.get(key, convert_pascal_to_snake(key))  # type: ignore[assignment]

    value, unit = extract_unit(key, value_str)
    value = _format_value(value)
    value = _re_map_single_value(key, value, VALUE_MAP)  # type: ignore[assignment]

    return key, value


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
                key, value = format_key_and_value(key, value)
                settings[key] = value

    return settings
