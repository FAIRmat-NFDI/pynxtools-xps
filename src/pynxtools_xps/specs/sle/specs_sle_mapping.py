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

from typing import List, Dict, Any


# from pynxtools_xps.reader_utils import (
#     convert_pascal_to_snake,
#     _re_map_single_value,
# )

from pynxtools_xps.value_mappers import (
    # convert_measurement_method,
    convert_energy_type,
    convert_energy_scan_mode,
    # MEASUREMENT_METHOD_MAP,
    convert_units,
    parse_datetime,
)


KEY_MAP: Dict[str, str] = {
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
    "Trace": "trace",
    "RawID": "raw_id",
    "Channel": "channel",
    # transmission_data_map
    "Data": "data",
}

VALUE_MAP: Dict[str, Any] = {
    "energy/@type": convert_energy_type,
    "excitation_energy": float,
    "time_stamp": parse_datetime,
    "energy_scan_mode": convert_energy_scan_mode,
}

KEYS_TO_DROP: List[str] = [
    "Work Function",
]

POSSIBLE_DATE_FORMATS: List[str] = ["%Y-%b-%d %H:%M:%S.%f"]
