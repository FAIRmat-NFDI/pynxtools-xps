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
from pynxtools_xps.reader_utils import XpsDataclass


@dataclass(order=True)
class VamasHeader(XpsDataclass):
    """An object to store the Vamas header information."""

    format_id: str = (
        "VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4"
    )
    institute_id: str = "Not Specified"
    instrument_model_id: str = "Not Specified"
    operator_id: str = "Not Specified"
    experiment_id: str = "Not Specified"
    num_comment_lines: int = 0
    comment_lines: list = field(
        default_factory=list
    )  # ["Casa Info Follows CasaXPS Version 2.3.22PR1.0\n0"]
    exp_mode: str = "NORM"
    scan_mode: str = "REGULAR"
    num_spectral_regions: int = 0
    num_analysis_positions: int = 0
    num_x_coords: int = 0
    num_y_coords: int = 0
    num_exp_var: int = 1
    num_entries_in_inclusion_list: int = 0
    inclusion_list: list = field(default_factory=list)
    num_manually_entered_items_in_block: int = 0
    manually_entered_items_in_block: list = field(default_factory=list)
    num_future_upgrade_exp_entries: int = 0
    num_future_upgrade_block_entries: int = 0
    future_upgrade_exp_entries: list = field(default_factory=list)
    num_blocks: int = 1


@dataclass
class VamasBlock(XpsDataclass):
    """An object to store a block of spectrum data and meta-data."""

    block_id: str = ""
    sample_id: str = ""
    year: int = 0
    month: int = 0
    day: int = 0
    hour: int = 0
    minute: int = 0
    second: int = 0
    no_hrs_in_advance_of_gmt: int = 0
    num_comment_lines: int = 0
    # This list should contain one element per for each
    # line in the comment block
    comment_lines: list = field(default_factory=list)
    technique: str = ""
    x_coord: float = 0.0
    y_coord: float = 0.0
    exp_var_value: str = ""
    source_label: str = ""
    sputter_ion_atomic_number: int = 0
    sputter_ion_num_atoms: int = 0
    sputter_ion_charge: int = 1
    source_energy: float = 0.0
    source_power: str = "0"
    source_beam_width_x: str = "0"
    source_beam_width_y: str = "0"
    field_of_view_x: float = 0.0
    field_of_view_y: float = 0.0
    first_linescan_x_start: float = 0.0
    first_linescan_y_start: float = 0.0
    first_linescan_x_end: float = 0.0
    first_linescan_y_end: float = 0.0
    last_linescan_x_end: float = 0.0
    last_linescan_y_end: float = 0.0
    source_polar_angle: float = 0.0
    source_azimuth_angle: float = 180.0
    analyser_mode: str = ""
    resolution: float = 0.0
    differential_width_aes: float = 0.0
    magnification: float = 1.0
    work_function: float = 0.0
    target_bias: float = 0.0
    # analyser slit length divided by the magnification
    # of the analyser transfer lens
    analysis_width_x: float = 0.0
    analysis_width_y: float = 0.0
    # degrees from upward z-direction,
    # defined by the sample stage
    analyser_take_off_polar_angle: float = 0.0
    analyser_take_off_azimuth_angle: float = 0.0
    species_label: str = ""
    transition_label: str = ""
    particle_charge: int = -1
    abscissa_label: str = "kinetic energy"
    abscissa_units: str = "eV"
    abscissa_start: float = 0.0
    abscissa_step: float = 0.0
    no_variables: int = 2
    variable_label_1: str = "counts"
    variable_units_1: str = "d"
    variable_label_2: str = "Transmission"
    variable_units_2: str = "d"
    signal_mode: str = "pulse counting"
    dwell_time: float = 0.0
    no_scans: int = 0
    time_correction: float = 0.0

    sputter_source_energy: float = 0.0
    sputter_source_beam_current: float = 0.0
    sputter_source_width_x: float = 0.0
    sputter_source_width_y: float = 0.0
    sputter_source_incidence_polar_angle: float = 0.0
    sputter_source_azimuth_angle: float = 0.0

    # degrees from upward z-direction,
    # defined by the sample stage
    sample_normal_polar_angle_of_tilt: float = 0.0
    # degrees clockwise from the y-direction towards the
    # operator, defined by the sample stage
    sample_normal_tilt_azimuth_angle: float = 0.0
    sample_rotation_angle: float = 0.0
    no_additional_params: int = 2
    num_ord_values: int = 0
    future_upgrade_block_entries: list = field(default_factory=list)


@dataclass
class ExpVariable(XpsDataclass):
    label: str = ""
    unit: str = ""


@dataclass
class VamasAdditionalParam(XpsDataclass):
    label: str = ""
    unit: str = ""
    value: float = 0.0


@dataclass
class OrdinateValue(XpsDataclass):
    min_ord_value: float = 0.0
    max_ord_value: float = 0.0
