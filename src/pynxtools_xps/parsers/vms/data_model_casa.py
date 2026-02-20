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
Data model for CasaXPS.
"""

import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pynxtools_xps.logging import _logger
from pynxtools_xps.mapping import _convert_energy_type
from pynxtools_xps.models.backgrounds import (
    LinearBackground,
    Shirley,
    StepDown,
    StepUp,
    TougaardU3,
    TougaardU4,
)
from pynxtools_xps.models.lineshapes import (
    DoniachSunjic,
    GaussianLorentzianProduct,
    GaussianLorentzianSum,
    LorentzianAsymmetric,
    LorentzianFinite,
)
from pynxtools_xps.parsers.base import _XPSDataclass

_LINESHAPES: dict[str, tuple[Any, ...]] = {
    "GL": ("Gaussian-Lorentzian Product", GaussianLorentzianProduct),
    "SGL": ("Gaussian-Lorentzian Sum", GaussianLorentzianSum),
    "LA": ("Asymmetric Lorentzian", LorentzianAsymmetric),
    "DS": ("Doniach-Sunjic", DoniachSunjic),
    "LF": ("Asymmetric Finite", LorentzianFinite),
}

_BACKGROUNDS: dict[str, tuple[Any, ...]] = {
    "Linear": ("Linear", LinearBackground),
    "Shirley": ("Shirley Sum", Shirley),
    "Step Up": ("Step Up", StepUp),
    "Step Down": ("Step Down", StepDown),
    "U 2 Tougaard": ("Tougaard", TougaardU3),
    "U 3 Tougaard": ("Tougaard", TougaardU3),
    "U 4 Tougaard": ("Tougaard", TougaardU4),
}


def _split_after_letters(string: str):
    """
    Splits the part of the string after the leading letters by commas.

    Parameters:
        string (str): The input string.

    Returns:
        tuple: A tuple with the leading letters and the list of values after splitting by commas.
    """
    # Extract the leading letters
    match = re.match(r"([A-Za-z\s]+)", string)
    leading_letters = match.group(0) if match else ""

    # Remove the leading letters and split the remaining part by commas
    remaining_part = string[len(leading_letters) :]
    split_values = remaining_part.split(",")

    if split_values == [""]:
        return leading_letters, []

    return leading_letters, [float(val) for val in split_values]


@dataclass
class CasaRegion(_XPSDataclass):
    """An object to store one CasaXPS peak fit region."""

    name: str = ""
    rsf: float = 0.0
    start: float = 0.0
    start_units: str = "eV"
    end: float = 0.0
    end_units: str = "eV"
    bg_type: str = ""
    function_type: str = ""
    av_width: float = 0.0
    av_width_units: str = "eV"
    start_offset: float = 0.0
    start_offset_units: str = ""
    end_offset: float = 0.0
    end_offset_units: str = ""
    cross_section: list = field(default_factory=list)
    tag: str = ""
    unknown_0: float = 0.0
    unknown_1: float = 0.0
    rsf_effective: float = 0.0

    def calculate_background(self, x: np.ndarray, y: np.ndarray):
        background_params = [
            float(param) for param in self.cross_section if float(param)
        ]

        try:
            background_class = _BACKGROUNDS[self.bg_type][-1]

            try:
                background = background_class(*background_params)
            except TypeError:
                background = background_class()

            if self.end > self.start:
                min_x = self.start
                max_x = self.end
            else:
                min_x = self.end
                max_x = self.start

            region = np.argwhere((x >= min_x) & (x <= max_x))
            fit_region = slice(region[0, 0], region[-1, 0], 1)

            self.start_offset = 100

            y_start_offset = y[0] * (self.start_offset / 100.0)
            y_end_offset = y[-1] * (self.end_offset / 100.0)
            y[0] -= y_start_offset
            y[-1] -= y_end_offset

            x, y = x[fit_region], y[fit_region]

            self.data = background.calc_background(x, y)

            self.formula = background.formula()
            self.description = str(background)

        except KeyError:
            _logger.warning(
                f"Background {self.name} (type {self.bg_type}) could not be parsed because no model exists for it."
            )


@dataclass
class CasaComponent(_XPSDataclass):
    """An object to store one CasaXPS peak fit component."""

    name: str = ""
    index: int = -1
    lineshape: str = ""
    function_type: str = ""

    area: float = 0.0
    area_min: float = 0.0
    area_max: float = 0.0
    area_ref_comp_id: int = 1
    area_ref_comp_factor: int = 1

    width: float = 0.0
    width_units: str = "eV"
    width_min: float = 0.0
    width_min_units: str = "eV"
    width_max: float = 0.0
    width_max_units: str = "eV"
    width_ref_comp_id: int = -1
    width_ref_comp_factor: int = 1

    position: float = 0.0
    position_units: str = "eV"
    position_min: float = 0.0
    position_min_units: str = "eV"
    position_max: float = 0.0
    position_max_units: str = "eV"
    position_ref_comp_id: int = -1
    position_ref_comp_factor: int = 1

    rsf: float = 0.0
    uncorrected_rsf: float = 0.0
    mass: float = 0.0
    tag: str = ""
    const: str = ""  # CONST

    atomic_concentration: float = 0.0

    def calculate_lineshape(self, x: np.ndarray):
        leading_letters, params = _split_after_letters(self.lineshape)

        try:
            peak_class = _LINESHAPES[leading_letters][-1]

            peak_parameters = [self.position, self.width, self.area] + params

            peak = peak_class(*peak_parameters)

            self.data = peak.calc_lineshape(x)
            self.formula = peak.formula()
            self.description = str(peak)

        except KeyError:
            _logger.warning(
                f"Component {self.name} (index {self.index}) could not be parsed because no model exists for it."
            )


class CasaProcess:
    """Processing and fitting information from CasaXPS."""

    def __init__(self):
        self.casa_data: dict[str, list[Any]] = {
            "energy_calibrations": [],
            "intensity_calibrations": [],
            "smoothings": [],
            "regions": [],
            "components": [],
        }

    def process_comments(self, comment_list: list[str]):
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
        """Flatten all process metadata."""
        process_keys: list[str] = [
            "energy_calibrations",
            "intensity_calibrations",
            "smoothings",
            "regions",
            "components",
        ]

        flattened: dict[str, Any] = {}

        for process_key in process_keys:
            processes = self.casa_data[process_key]
            for i, process in enumerate(processes):
                spectrum_key = f"{process_key.rstrip('s')}{i}"
                for key, value in process.dict().items():
                    key = key.replace("_units", "/@units")
                    flattened[f"{spectrum_key}/{key}"] = value

        return flattened

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

        measured_index = split_line.index("M")
        aligned_index = split_line.index("A")

        # --- measured energies -------------------------------------------------
        raw_measured = split_line[measured_index + 2]

        try:
            measured_energies = [float(raw_measured)]
            range_calibration = False
        except ValueError:
            measured_energies = [float(e) for e in raw_measured.split(",")]
            range_calibration = True

        aligned_energy = float(split_line[aligned_index + 2])

        # --- operation handling ------------------------------------------------
        operation = split_line[8]
        energy_offset = None

        if operation == "ADD":
            energy_offset = min(measured_energies) - aligned_energy
            operation = "addition"

        energy_calib = CasaEnergyCalibration(
            measured_energies=measured_energies,
            aligned_energy=aligned_energy,
            energy_type=_convert_energy_type(split_line[7]),
            operation=operation,
            range_calibration=range_calibration,
            energy_offset=energy_offset,
        )

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

        bg_type = split_line[3]

        region = CasaRegion(
            name=split_line[2],
            bg_type=bg_type,
            function_type=_BACKGROUNDS[bg_type][0],
            start=float(split_line[4]),
            end=float(split_line[5]),
            rsf=float(split_line[6]),
            av_width=float(split_line[7]),
            start_offset=float(split_line[8]),
            end_offset=float(split_line[9]),
            cross_section=split_line[10:14],
            tag=split_line[14],
            unknown_0=float(split_line[15]),
            unknown_1=float(split_line[16]),
            rsf_effective=float(split_line[17]),
        )
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

        def idx(token: str) -> int:
            return split_line.index(token)

        split_line = self._split_line(component_str)

        lineshape = split_line[3]
        function_key = _split_after_letters(lineshape)[0]

        area_index = idx("Area")
        width_index = idx("MFWHM")
        position_index = idx("Position")
        rsf_index = idx("RSF")
        mass_index = idx("MASS")
        comp_index = idx("INDEX")
        const_index = idx("CONST")
        uncorrected_rsf_index = idx("UNCORRECTEDRSF")

        component = CasaComponent(
            name=split_line[2],
            lineshape=lineshape,
            function_type=_LINESHAPES[function_key][0],
            area=float(split_line[area_index + 1]),
            area_min=float(split_line[area_index + 2]),
            area_max=float(split_line[area_index + 3]),
            area_ref_comp_id=int(split_line[area_index + 4]),
            area_ref_comp_factor=int(split_line[area_index + 5]),
            width=float(split_line[width_index + 1]),
            width_min=float(split_line[width_index + 2]),
            width_max=float(split_line[width_index + 3]),
            width_ref_comp_id=int(split_line[width_index + 4]),
            width_ref_comp_factor=int(split_line[width_index + 5]),
            position=float(split_line[position_index + 1]),
            position_min=float(split_line[position_index + 2]),
            position_max=float(split_line[position_index + 3]),
            position_ref_comp_id=int(split_line[position_index + 4]),
            position_ref_comp_factor=int(split_line[position_index + 5]),
            rsf=float(split_line[rsf_index + 1]),
            mass=float(split_line[mass_index + 1]),
            index=int(split_line[comp_index + 1]),
            tag=str(split_line[comp_index + 2]),
            const=str(split_line[const_index + 1]),
            uncorrected_rsf=float(split_line[uncorrected_rsf_index + 1]),
        )

        component.validate_types()

        self.casa_data["components"] += [component]

    def _split_line(self, line):
        """Split line in Casa processing."""
        regex_pattern = r"\s+|\(\*(?:(?!\*\)).)*?\*\)|\S+"
        result = [match for match in re.findall(regex_pattern, line) if match.strip()]
        return [re.sub(r"[\(*?\*)]", "", text) for text in result]


@dataclass
class CasaEnergyCalibration(_XPSDataclass):
    """An object to store one CasaXPS energy calibration."""

    energy_type: str = "binding"
    energy_offset: float = 0.0
    energy_offset_units: str = "eV"
    measured_energies: list = field(default_factory=list)
    measured_energies_units: str = "eV"
    aligned_energy: float = 0.0
    aligned_energy_units: str = "eV"
    operation: str = "ADD"
    range_calibration: bool = False

    def apply_energy_shift(self, x: float):
        if self.operation == "addition":
            if self.energy_type == "binding":
                return x - self.energy_offset

            elif self.energy_type == "kinetic":
                return x + self.energy_offset

        if self.range_calibration:
            pass  # ToDo: apply range calibration


@dataclass
class CasaIntensityCalibration(_XPSDataclass):
    """An object to store one CasaXPS energy calibration."""

    # TODO: enable more types of intensity calibrations
    power_law_coefficient: float = 0.0
    tf_function_ordinate_no: int = 0
