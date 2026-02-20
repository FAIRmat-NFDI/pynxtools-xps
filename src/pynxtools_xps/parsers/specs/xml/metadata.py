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
Metadata mapping for the SPECS XY parser.
"""

from pynxtools_xps.mapping import (
    _convert_energy_scan_mode,
    _convert_energy_type,
    _convert_measurement_method,
    _MetadataContext,
    _ValueMap,
)

_KEY_MAP: dict[str, str] = {
    "group": "group_id",
    "analyzer_slit": "entrance_slit",
    "energy_axis": "x_units",
    "source": "source_label",
    "polar_angle": "source_polar_angle",
    "azimuth_angle": "source_azimuth_angle",
    "ex_energy": "excitation_energy",
    # "acquisition_date": "time_stamp",
    "analyzer": "analyzer_name",
    "analyzer_lens": "lens_mode",
    "analyzer_lens_mode": "lens_mode",
    "curves/scan": "curves_per_scan",
    "values/curve": "n_values",
    "bias_voltage": "bias_voltage_electrons",
    "binding_energy": "start_energy",
    "eff._workfunction": "work_function",
    "comment": "comments",
}

_VALUE_MAP: _ValueMap = {
    "analysis_method": _convert_measurement_method,
    "scan_mode": _convert_energy_scan_mode,
    "bias_voltage_electrons": float,
    "n_values": int,
    "excitation_energy": float,
    "kinetic_energy": float,
    "work_function": float,
    "dwell_time": float,
    "detector_voltage": float,
    "curves_per_scan": int,
    "pass_energy": float,
    "spectrum_id": int,
    "x_units": _convert_energy_type,
}

_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {
    "work_function": "eV",
    "excitation_energy": "eV",
    "pass_energy": "eV",
    "bias_voltage_electrons": "V",
    "dwell_time": "s",
    "step_size": "eV",
}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
