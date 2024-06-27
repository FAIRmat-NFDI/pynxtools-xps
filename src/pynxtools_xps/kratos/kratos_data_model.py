"""
Data model for Kratos spectrometers for export to CasaXPS.
"""
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
# pylint: disable=too-many-instance-attributes

from dataclasses import dataclass
from dataclasses import field

from pynxtools_xps.reader_utils import XpsDataclass


@dataclass
class KratosMetadata(XpsDataclass):
    """An object to store the Kratos metadata."""

    created_by: str = ""
    date_created: str = ""
    computer: str = ""
    location_id: str = ""
    sample: str = ""
    description: dict = field(default_factory=dict)
    sample_tilt: str = ""
    charge_neutraliser: str = ""
    filament_current: float = 0.0
    filament_current_units: str = ""
    filament_bias: float = 0.0
    filament_bias_units: str = ""
    charge_balance: float = 0.0
    charge_balance_units: str = ""
    mode: str = ""
    tuning: str = ""
    emission_current: float = 0.0
    emission_current_units: str = ""
    anode_material: str = ""
    collimation: str = ""
    lens_mode: str = ""
    resolution: float = 0.0
    resolution_units: str = "eV"
    deflection: dict = field(default_factory=dict)

    energy_start: float = 0.0
    energy_start_units: str = "eV"
    energy_end: float = 0.0
    energy_end_units: str = "eV"
    energy_centre: float = 0.0
    energy_centre_units: str = "eV"
    energy_width: float = 0.0
    energy_width_units: str = "eV"
    step_size: float = 0.0
    step_size_units: str = "eV"
    number_steps: int = 0
    dwell_time: float = 0.0
    dwell_time_units: str = "s"
    sweep_time: float = 0.0
    sweep_time_units: str = "s"
    sweeps: int = 0

    x_ray_power: float = 0.0
    x_ray_power_units: str = "W"
