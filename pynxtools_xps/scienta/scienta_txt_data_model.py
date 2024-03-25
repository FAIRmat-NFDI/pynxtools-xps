"""
Data model for Vamas ISO standard.
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

from dataclasses import dataclass, field
import numpy as np

from pynxtools_xps.reader_utils import XpsDataclass


@dataclass
class ScientaHeader(XpsDataclass):
    no_of_regions: int = 0
    software_version: str = ""


@dataclass
class ScientaRegion(XpsDataclass):
    region_id: int = 0
    region_name: str = ""
    energy_type: str = ""
    energy_size: int = 0
    energy_axis: np.ndarray = np.array([])
    lens_mode: str = ""
    pass_energy: float = 0.0
    pass_energy_units: str = "eV"
    no_of_scans: int = 0
    excitation_energy: float = 0.0
    excitation_energy_units: str = "eV"
    energy_type: str = "binding"  # energy_scale
    acquisition_mode: str = ""
    energy_units: str = "kinetic"  # ???
    center_energy: float = 0.0
    center_energy_units: str = "eV"
    start_energy: float = 0.0
    start_energy_units: str = "eV"
    stop_energy: float = 0.0
    stop_energy_units: str = "eV"
    step_size: float = 0.0
    stop_energy_units: str = "eV"
    dwell_time: float = 0.0
    dwell_time_units: str = "eV"
    detector_first_x_channel: int = 0
    detector_last_x_channel: int = 0
    detector_first_y_channel: int = 0
    detector_last_y_channel: int = 0
    number_of_slices: str = ""
    data_file: str = ""
    sequence_file: str = ""
    spectrum_type: str = ""
    instrument_name: str = ""
    vendor: str = ""
    user_name: str = ""
    sample_name: str = ""
    spectrum_comment: str = ""
    start_date: str = ""
    start_time: str = ""
    time_per_spectrum_channel: float = 0.0
    time_per_spectrum_channel_units: str = "s"
    detector_mode: str = ""
    data: dict = field(default_factory=dict)
