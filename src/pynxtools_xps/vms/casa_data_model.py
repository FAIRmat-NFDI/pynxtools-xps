"""
Data model for CasaXPS.
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

import re
from typing import List, Dict, Any
from dataclasses import dataclass, field

from pynxtools_xps.reader_utils import XpsDataclass
from pynxtools_xps.value_mappers import convert_energy_type


class CasaProcess:
    """Processing and fitting information from CasaXPS."""

    def __init__(self):
        self.casa_data: Dict[str, List[Any]] = {
            "energy_calibrations": [],
            "intensity_calibrations": [],
            "smoothings": [],
            "regions": [],
            "components": [],
        }

    def process_comments(self, comment_list: List[str]):
        line_map = {
            "Calib": self.process_energy_calibration,
            "KECalib": self.process_intensity_calibration,
            "TransmissionSet": self.process_intensity_calibration,
            "InCalib": self.process_intensity_calibration,
            "Sm": self.process_smoothing,
            "CASA region": self.process_region,
            "CASA comp": self.process_component,
        }

        n_processes, n_annotations, n_regions, n_components = [
            int(element)
            for element in comment_list
            if isinstance(element, str) and element.isdigit()
        ][:4]
        n_annotations *= 2  # account for empty line

        self.no_of_casa_lines = (
            5 + n_processes + n_annotations + n_regions + n_components
        )

        for line in comment_list[: self.no_of_casa_lines]:
            for token, func in line_map.items():
                if line.startswith(token):
                    func(line)

        return self.casa_data

    def flatten_metadata(self):
        # Write process data
        process_key_map: Dict[str, List[str]] = {
            "process": ["energy_calibrations", "intensity_calibrations", "smoothings"],
            "peak_fitting": ["regions", "components"],
        }

        flattened_dict: Dict[str, Any] = {}

        for grouping, process_key_list in process_key_map.items():
            for spectrum_key in process_key_list:
                processes = self.casa_data[spectrum_key]
                for i, process in enumerate(processes):
                    process_key = f"{spectrum_key}/{spectrum_key.rstrip('s')}{i}"
                    for key, value in process.dict().items():
                        key = key.replace("_units", "/@units")
                        flattened_dict[f"{process_key}/{key}"] = value

        return flattened_dict

    def process_energy_calibration(self, calib_str: str):
        """
        Process one alignment.

        Parameters
        ----------
        align_str : str
            String with alignment definition in a VMS file.
            Typical example:
                'Calib M = 6.8 A = 0 BE ADD'


        Returns
        -------
        None.

        """
        split_line = self._split_line(calib_str)

        energy_calib = CasaEnergyCalibration()

        measured_index = split_line.index("M")
        aligned_index = split_line.index("A")
        try:
            energy_calib.measured_energies = [float(split_line[measured_index + 2])]
        except ValueError:
            energies = [
                float(energy) for energy in split_line[measured_index + 2].split(",")
            ]
            energy_calib.measured_energies = energies
            energy_calib.range_calibration = True
        energy_calib.aligned_energy = float(split_line[aligned_index + 2])

        energy_calib.energy_type = convert_energy_type(split_line[7])
        energy_calib.operation = split_line[8]

        if energy_calib.operation == "ADD":
            energy_calib.energy_offset = (
                min(energy_calib.measured_energies) - energy_calib.aligned_energy
            )
            energy_calib.operation = "addition"

        energy_calib.validate_types()

        self.casa_data["energy_calibrations"] += [energy_calib]

    def process_intensity_calibration(self, calib_str: str):
        """Process different intensity calibrations."""

        intensity_calib = CasaIntensityCalibration()

        def get_power_law_coefficient(line: str):
            intensity_calib.power_law_coefficient = float(line.split("KECalib ")[1])

        def get_transmission_function_ordinate(line: str):
            intensity_calib.tf_function_ordinate_no = int(
                line.split("TransmissionSet 1.0 ")[1]
            )

        def get_both_calibrations(line: str):
            split_line = line.split(" ")
            intensity_calib.power_law_coefficient = float(split_line[1])
            intensity_calib.tf_function_ordinate_no = int(split_line[2])

        calib_map = {
            "KECalib": get_power_law_coefficient,
            "TransmissionSet": get_transmission_function_ordinate,
            "InCalib": get_both_calibrations,
        }

        for token, func in calib_map.items():
            if calib_str.startswith(token):
                func(calib_str)

        intensity_calib.validate_types()

        self.casa_data["intensity_calibrations"] += [intensity_calib]

    def process_smoothing(self, region_str: str):
        # TODO: parse smoothing of data
        pass

    def process_region(self, region_str: str):
        """
        Process one fit region.

        Parameters
        ----------
        region_str : str
            String with region definition in a VMS file.
            Typical example:
                'CASA region (*O1s*) (*Shirley*) 718.88221 726.6 9.74 2 0 9 1450.5943 -450 0 0 (*Mo 3d*) 95.9219 0 9.74'

        Returns
        -------
        None.

        """
        split_line = self._split_line(region_str)

        region = CasaRegion()
        region.name = split_line[2]
        region.bg_type = split_line[3]
        region.start = float(split_line[4])
        region.end = float(split_line[5])
        region.rsf = float(split_line[6])
        region.av_width = float(split_line[7])
        region.start_offset = float(split_line[8])
        region.end_offset = float(split_line[9])
        region.cross_section = split_line[10:14]
        region.tag = split_line[14]
        region.unknown_0 = float(split_line[15])
        region.unknown_1 = float(split_line[16])
        region.rsf_effective = float(split_line[17])

        region.validate_types()

        self.casa_data["regions"] += [region]

    def process_component(self, component_str: str):
        """
        Process one fit component if it is defined in the

        Parameters
        ----------
        component_str : str
            String with component definition in a VMS file.
            Typical example:
                'CASA comp (*O1s, lattice*) (*GL(50)*) Area 308996.39 1e-020 15401311 -1 1 MFWHM 1.5315216 0.32 2 -1 1 Position 723.30224 719.3797 726.3 -1 1 RSF 0.063 MASS 15.9994 INDEX -1 (*O 1s*) CONST (**) UNCORRECTEDRSF 2.85'


        Returns
        -------
        None.

        """
        split_line = self._split_line(component_str)

        component = CasaComponent()
        component.name = split_line[2]
        component.lineshape = split_line[3]

        area_index = split_line.index("Area")
        component.area = float(split_line[area_index + 1])
        component.area_min = float(split_line[area_index + 2])
        component.area_max = float(split_line[area_index + 3])
        component.area_ref_comp_id = int(split_line[area_index + 4])
        component.area_ref_comp_factor = float(split_line[area_index + 4])

        width_index = split_line.index("MFWHM")
        component.width = float(split_line[width_index + 1])
        component.width_min = float(split_line[width_index + 2])
        component.width_max = float(split_line[width_index + 3])
        component.width_ref_comp_id = int(split_line[width_index + 4])
        component.width_ref_comp_factor = float(split_line[width_index + 5])

        position_index = split_line.index("Position")
        component.position = float(split_line[position_index + 1])
        component.position_min = float(split_line[position_index + 2])
        component.position_max = float(split_line[position_index + 3])
        component.position_ref_comp_id = int(split_line[position_index + 4])
        component.position_ref_comp_factor = float(split_line[position_index + 5])

        rsf_index = split_line.index("RSF")
        component.rsf = float(split_line[rsf_index + 1])
        mass_index = split_line.index("MASS")
        component.mass = float(split_line[mass_index + 1])
        comp_index = split_line.index("INDEX")
        component.index = int(split_line[comp_index + 1])
        component.tag = str(split_line[comp_index + 2])
        const_index = split_line.index("CONST")
        component.const = str(split_line[const_index + 1])
        uncorrected_rsf_index = split_line.index("UNCORRECTEDRSF")
        component.uncorrected_rsf = float(split_line[uncorrected_rsf_index + 1])

        component.validate_types()

        self.casa_data["components"] += [component]

    def _split_line(self, line):
        """Split line in Casa processing."""
        regex_pattern = r"\s+|\(\*(?:(?!\*\)).)*?\*\)|\S+"
        result = [match for match in re.findall(regex_pattern, line) if match.strip()]
        return [re.sub(r"[\(*?\*)]", "", text) for text in result]


@dataclass
class CasaEnergyCalibration(XpsDataclass):
    """An object to store one CasaXPS energy calibration."""

    energy_type: str = "binding"
    energy_offset: float = 0.0
    energy_offset_units: str = "eV"
    measured_energies: list = field(default_factory=list)
    measured_energies_units: str = "eV"
    aligned_energy: float = 0.0
    aligned_energy_units: str = "eV"
    operation: str = "eV"
    range_calibration: bool = False


@dataclass
class CasaIntensityCalibration(XpsDataclass):
    """An object to store one CasaXPS energy calibration."""

    # TODO: enable more types of intensity calibrations
    power_law_coefficient: float = 0.0
    tf_function_ordinate_no: int = 0


@dataclass
class CasaRegion(XpsDataclass):
    """An object to store one CasaXPS peak fit region."""

    name: str = ""
    rsf: float = 0.0
    start: float = 0.0
    start_units: str = "eV"
    end: float = 0.0
    end_units: str = "eV"
    bg_type: str = ""
    av_width: float = 0.0
    av_width_units: str = "eV"
    start_offset: float = 0.0
    start_offset_units: str = "counts_per_second"
    end_offset: float = 0.0
    end_offset_units: str = "counts_per_second"
    cross_section: list = field(default_factory=list)
    tag: str = ""
    unknown_0: float = 0.0
    unknown_1: float = 0.0
    rsf_effective: float = 0.0


@dataclass
class CasaComponent(XpsDataclass):
    """An object to store one CasaXPS peak fit component."""

    name: str = ""
    index: int = -1
    lineshape: str = ""

    area: float = 0.0
    area_min: float = 0.0
    area_max: float = 0.0
    area_ref_comp_id: int = -1
    area_ref_comp_factor: float = 0.0

    width: float = 0.0
    width_units: str = "eV"
    width_min: float = 0.0
    width_min_units: str = "eV"
    width_max: float = 0.0
    width_max_units: str = "eV"
    width_ref_comp_id: int = -1
    width_ref_comp_factor: float = 0.0

    position: float = 0.0
    position_units: str = "eV"
    position_min: float = 0.0
    position_min_units: str = "eV"
    position_max: float = 0.0
    position_max_units: str = "eV"
    position_ref_comp_id: int = -1
    position_ref_comp_factor: float = 0.0

    rsf: float = 0.0
    uncorrected_rsf: float = 0.0
    mass: float = 0.0
    tag: str = ""
    const: str = ""  # CONST
