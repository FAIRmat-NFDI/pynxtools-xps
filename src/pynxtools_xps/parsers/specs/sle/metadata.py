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
Metadata mapping for the SPECS SLE reader.
"""

from functools import partial

from pynxtools_xps.mapping import (
    _convert_energy_scan_mode,
    _convert_energy_type,
    _convert_measurement_method,
    _MetadataContext,
    _ValueMap,
    parse_datetime,
)

_POSSIBLE_DATE_FORMATS: list[str] = ["%Y-%b-%d %H:%M:%S.%f"]

# TODO: should this be generalized?
# _KEYS_TO_DROP: list[str] = [
#     "Work Function",
# ]

# TODO: wrong type, does this need to be a dict of dicts?
_KEY_MAP: dict[str, str | dict[str, str]] = {
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

_VALUE_MAP: _ValueMap = {
    "energy/@type": _convert_energy_type,
    "excitation_energy": float,
    "time_stamp": partial(
        parse_datetime,
        possible_date_formats=_POSSIBLE_DATE_FORMATS,
    ),
    "energy_scan_mode": _convert_energy_scan_mode,
    "measurement_type": _convert_measurement_method,
}

_UNIT_MAP: dict[str, str | None] = {
    "a.u.": "counts",
    "Counts": "counts",
    "counts/s": "counts_per_second",
    "CPS": "counts_per_second",
    "u": "um",
    "KV": "kV",
    "seconds": "s",
    "(min)": "min",
    # "Percent": "",  # should be changed back to percent once pint is updated
    "atom": "dimensionless",
    "eV/atom": "eV",
    "microm": "micrometer",
    "d": "degree",
    "nU": "",  # workaround for SPECS SLE reader
}

_DEFAULT_UNITS: dict[str, str] = {
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

# TODO: remove type-ignore
_context = _MetadataContext(
    key_map=_KEY_MAP,  # type: ignore[arg-type]
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
