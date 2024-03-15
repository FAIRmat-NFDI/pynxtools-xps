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


@dataclass
class VamasHeader(XpsDataclass):
    """An object to store the Vamas header information."""

    format_id: str = (
        "VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4"
    )
    institute_id: str = "Not Specified"
    instrument_model_id: str = "Not Specified"
    operator_id: str = "Not Specified"
    experiment_id: str = "Not Specified"
    no_comment_lines: str = "2"
    comment_lines: list = field(
        default_factory=list
    )  # ["Casa Info Follows CasaXPS Version 2.3.22PR1.0\n0"]
    exp_mode: str = "NORM"
    scan_mode: str = "REGULAR"
    nr_regions: str = "0"
    nr_exp_var: str = "1"
    exp_var_label: str = "Exp Variable"
    exp_var_unit: str = "d"
    unknown_3: str = "0"
    unknown_4: str = "0"
    unknown_5: str = "0"
    unknown_6: str = "0"
    no_blocks: str = "1"


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
    no_comment_lines: int = 0
    # This list should contain one element per for each
    # line in the comment block
    comment_lines: list = field(default_factory=list)
    technique: str = ""
    exp_var_value: str = ""
    source_label: str = ""
    source_energy: float = 0.0
    unknown_1: str = "0"
    unknown_2: str = "0"
    unknown_3: str = "0"
    source_analyzer_angle: float = 0.0
    unknown_4: str = "180"
    analyzer_mode: str = ""
    resolution: float = 0.0
    magnification: str = "1"
    work_function: float = 0.0
    target_bias: float = 0.0
    # analyser slit length divided by the magnification
    # of the analyser transfer lens
    analyzer_width_x: float = 0.0
    analyzer_width_y: float = 0.0
    # degrees from upward z-direction,
    # defined by the sample stage
    analyzer_take_off_polar_angle: float = 0.0
    analyzer_azimuth: float = 0.0
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
    time_correction: str = "0"
    # degrees from upward z-direction,
    # defined by the sample stage
    sample_angle_tilt: float = 0.0
    # degrees clockwise from the y-direction towards the
    # operator, defined by the sample stage
    sample_tilt_azimuth: float = 0.0
    sample_rotation: float = 0.0
    no_additional_params: int = 2
    param_label_1: str = "ESCAPE DEPTH TYPE"
    param_unit_1: str = "d"
    param_value_1: str = "0"
    param_label_2: str = "MFP Exponent"
    param_unit_2: str = "d"
    param_value_2: str = "0"
    num_ord_values: int = 0
    min_ord_value_1: float = 0.0
    max_ord_value_1: float = 0.0
    min_ord_value_2: float = 0.0
    max_ord_value_2: float = 0.0
    data_string: str = ""
