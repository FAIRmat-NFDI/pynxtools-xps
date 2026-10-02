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
Metadata mapping for the SPECS metadata reader.
"""

from pynxtools_xps.mapping import _MetadataContext, _ValueMap

# TODO: these are almost certainly incorrect
_POSSIBLE_DATE_FORMATS: list[str] = ["%Y-%b-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S.%f"]

# TODO: is this needed?
_KEY_MAP: dict[str, str] = {
    # spectrometer_setting_map
    "Coil Current [mA]": "coil_current [mA]",
    "Pre Defl Y [nU]": "pre_deflector_y_current [nU]",
    "Pre Defl X [nU]": "pre_deflector_x_current [nU]",
    "L1 [nU]": "lens1_voltage [nU]",
    "L2 [nU]": "lens2_voltage [nU]",
    "Focus Displacement 1 [nu]": "focus_displacement_current [nU]",
    "Detector Voltage [V]": "detector_voltage [V]",
    "Bias Voltage Electrons [V]": "bias_voltage_electrons [V]",
    "Bias Voltage Ions [V]": "bias_voltage_ions [V]",
    # source_setting_map
    "anode": "source_label",
    "uanode": "source_voltage",
    "iemission": "emission_current",
    "ihv": "source_high_voltage",
    "ufilament": "filament_voltage",
    "ifilament": "filament_current",
    "DeviceExcitationEnergy": "excitation_energy",
    "panode": "anode_power",
    "temperature": "source_temperature",
    # sql_metadata_map
    "EnergyType": "x_units",
    "EpassOrRR": "pass_energy",
    "Wf": "workfunction",
    "Timestamp": "time_stamp",
    "Samples": "n_values",
    "ElectronEnergy": "start_energy",
    "Step": "step_size",
}

_VALUE_MAP: _ValueMap = {
    # TODO: are these valid
    # "x_units": self._change_energy_type,
    # "time_stamp": self._convert_date_time,
}
_UNIT_MAP: dict[str, str | None] = {}
_DEFAULT_UNITS: dict[str, str] = {}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
