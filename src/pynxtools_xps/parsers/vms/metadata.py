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
Metadata mapping for the VAMAS parser.
"""

from pynxtools_xps.mapping import (
    _convert_energy_scan_mode,
    _convert_energy_type,
    _convert_measurement_method,
    _MetadataContext,
    _ValueMap,
)

EXP_MODES = [
    "MAP",
    "MAPDP",
    "MAPSV",
    "MAPSVDP",
    "NORM",
    "SDP",
    "SDPSC",
    "SEM",
]

ALLOWED_TECHNIQUES = [
    "AES",
    "AES diff",
    "AES dir",
    "EDX",
    "ELS",
    "FABMS",
    "FABMS energy spec",
    "ISS",
    "SIMS",
    "SIMS energy spec",
    "SNMS",
    "SNMS energy spec",
    "UPS",
    "XPS",
    "XRF",
]

_KEY_MAP: dict[str, str] = {
    "block_id": "region",
    "sample_id": "sample_name",
    "technique": "analysis_method",
    "source_energy": "excitation_energy",
    "analyzer_mode": "scan_mode",
    "resolution": "pass_energy",
    "transition_label": "transition",
    "abscissa_label": "energy_label",
    "abscissa_units": "energy_units",
    "abscissa_start": "start_energy",
    "abscissa_step": "step_size",
    "variable_label_1": "y_labels_1",
    "variable_units_1": "y_units_1",
    "variable_label_2": "y_labels_2",
    "variable_units_2": "y_units_2",
    "species_label": "element",
}

_VALUE_MAP: _ValueMap = {
    "analysis_method": _convert_measurement_method,
    "energy_label": _convert_energy_type,
    "scan_mode": _convert_energy_scan_mode,
}

_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {
    "excitation_energy": "eV",
    "extent": "um",
    "particle_charge": "C",
    "source_beam_width_x": "um",
    "source_beam_width_y": "um",
    "source_polar_angle": "degree",
    "source_azimuth_angle": "degree",
    "source_power": "W",
    "sputter_ion_charge": "C",
    "sputter_source_energy": "eV",
    "sputter_source_beam_current": "A",
    "sputter_source_width_x": "um",
    "sputter_source_width_y": "um",
    "sputter_source_incidence_polar_angle": "degree",
    "sputter_source_azimuth_angle": "degree",
    "analyzer_take_off_polar_angle": "degree",
    "analyzer_take_off_azimuth_angle": "degree",
    "analysis_width_x": "m",
    "analysis_width_y": "m",
    "target_bias": "V",
    "time_correction": "s",
    "work_function": "eV",
    "spatial_acceptance": "um",
    "pass_energy": "eV",
    "differential_width_aes": "eV",
    "dwell_time": "s",
    "sample_normal_polar_angle_of_tilt": "degree ",
    "sample_normal_tilt_azimuth_angle": "degree",
    "sample_rotation_angle": "degree",
    "start_energy": "eV",
    "step_size": "eV",
}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
