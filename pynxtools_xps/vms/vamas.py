"""
Class for reading XPS files from raw VMS data.
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
# pylint: disable=too-many-lines

import re
from copy import deepcopy
import datetime
from pathlib import Path
from typing import Any, Dict, List, Union

from itertools import groupby
import xarray as xr
import numpy as np

from pynxtools_xps.vms.vamas_data_model import VamasHeader, VamasBlock
from pynxtools_xps.vms.casa_data_model import CasaProcess
# from pynxtools_xps.phi.spe_pro_phi import PhiParser

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    construct_data_key,
    construct_detector_data_key,
    to_snake_case,
    get_minimal_step,
)


class VamasMapper(XPSMapper):
    """
    Class for restructuring .txt data file from
    Vamas format into python dictionary.
    """

    config_file = "config_vms.json"

    def __init__(self):
        self.file: Union[str, Path] = ""
        self.parsers: List[Any] = []

        self.units: dict = {
            "instrument/sample_normal_polarangle_tilt": "degree ",
            "instrument/sample_tilt_azimuth": "degree",
            "instrument/sample_rotation_angle": "degree",
            "source/source_analyzer_angle": "degree",
            "source/excitation_energy": "eV",
            "source/particle_charge": "C",
            "analyser/analyzer_take_off_azimuth": "degree",
            "analyser/analyzer_take_off_polar": "degree",
            "analyser/analysis_width_x": "m",
            "analyser/analysis_width_y": "m",
            "analyser/target_bias": "V",
            "analyser/time_correction": "s",
            "analyser/work_function": "eV",
            "energydispersion/pass_energy": "eV",
            "detector/dwell_time": "s",
            "data/start_energy": "eV",
            "data/step_size": "eV",
        }

        super().__init__()

    def _select_parser(self):
        """
        Select parser based on the structure of the Vamas file,
        i.e., whether it is regular or irregular.

        Returns
        -------
        VamasParserVMS
            Vamas parser for reading this file structure.

        """
        return VamasParser()

    def construct_data(self):
        """Map VMS format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code
        spectra = deepcopy(self.raw_data)

        self._xps_dict["data"]: Dict[str, Any] = {}

        key_map = {
            "user": [],
            "instrument": [
                "sample_normal_polarangle_tilt",
                "sample_tilt_azimuth",
                "sample_rotation_angle",
            ],
            "source": [
                "source_label",
                "source_analyzer_angle",
            ],
            "beam": ["excitation_energy", "particle_charge"],
            "analyser": [
                "analyzer_take_off_azimuth",
                "analyzer_take_off_polar",
                "analysis_width_x",
                "analysis_width_y",
                "target_bias",
                "time_correction",
                "work_function",
            ],
            "collectioncolumn": [
                "magnification",
            ],
            "energydispersion": [
                "scan_mode",
                "pass_energy",
            ],
            "detector": [
                "signal_mode",
                "dwell_time",
            ],
            "manipulator": [],
            "sample": ["sample_name"],
            "data": [
                "x_label",
                "x_units",
                "y_labels_1",
                "y_units_1",
                "y_labels_2",
                "y_units_2",
                "n_values",
                "start_energy",
                "step_size",
            ],
            "region": [
                "analysis_method",
                "element",
                "group_id",
                "group_name",
                "region",
                "scan_no",
                "scans",
                "spectrum_id",
                "spectrum_type",
                "transition",
                "time_stamp",
            ],
        }

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum, key_map)

    def _update_xps_dict_with_spectrum(
        self, spectrum: Dict[str, Any], key_map: Dict[str, str]
    ):
        """Map one spectrum from raw data to NXmpes-ready dict."""
        # pylint: disable=too-many-locals,duplicate-code
        group_parent = f'{self._root_path}/Group_{spectrum["group_name"]}'
        region_parent = f'{group_parent}/regions/Region_{spectrum["spectrum_type"]}'
        instrument_parent = f"{region_parent}/instrument"
        analyser_parent = f"{instrument_parent}/analyser"

        path_map: Dict[str, str] = {
            "user": f"{region_parent}/user",
            "instrument": f"{instrument_parent}",
            "source": f"{instrument_parent}/source",
            "beam": f"{instrument_parent}/beam",
            "analyser": f"{analyser_parent}",
            "collectioncolumn": f"{analyser_parent}/collectioncolumn",
            "energydispersion": f"{analyser_parent}/energydispersion",
            "detector": f"{analyser_parent}/detector",
            "manipulator": f"{instrument_parent}/manipulator",
            "sample": f"{region_parent}/sample",
            "data": f"{region_parent}/data",
            "energy_referencing": f"{region_parent}/calibrations/energy_referencing",
            "peak_fitting": f"{region_parent}/peak_fitting",
            "region": f"{region_parent}",
        }

        used_keys = []

        for grouping, spectrum_keys in key_map.items():
            root = path_map[str(grouping)]
            for spectrum_key in spectrum_keys:
                mpes_key = spectrum_key.rsplit(" ", 1)[0]
                self._xps_dict[f"{root}/{mpes_key}"] = spectrum[spectrum_key]

                unit_key = f"{grouping}/{spectrum_key}"
                units = self._get_units_for_key(unit_key)
                if units:
                    self._xps_dict[f"{root}/{mpes_key}/@units"] = units

        # Write process data
        process_key_map: Dict[str, List[str]] = {
            "energy_referencing": ["alignments"],
            "peak_fitting": ["regions", "components"],
        }

        for grouping, process_key_list in process_key_map.items():
            root = path_map[str(grouping)]
            for spectrum_key in process_key_list:
                try:
                    processes = spectrum[spectrum_key]
                    for i, process in enumerate(processes):
                        process_key = (
                            f"{root}/{spectrum_key}/{spectrum_key.rstrip('s')}{i}"
                        )
                        for key, value in process.dict().items():
                            key = key.replace("_units", "/@units")
                            self._xps_dict[f"{process_key}/{key}"] = value
                    used_keys += [spectrum_key]
                except KeyError:
                    pass

        # Create keys for writing to data and detector
        entry = construct_entry_name(region_parent)
        scan_key = construct_data_key(spectrum)
        detector_data_key_child = construct_detector_data_key(spectrum)
        detector_data_key = f'{path_map["detector"]}/{detector_data_key_child}/counts'

        energy = np.array(spectrum["data"]["x"])
        intensity_raw = np.array(spectrum["data"]["y"])
        intensity_cps = np.array(spectrum["data"]["y_cps"])

        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        # Write averaged cycle data to 'data'.
        all_scan_data = [
            np.array(value)
            for key, value in self._xps_dict["data"][entry].items()
            if scan_key.split("_")[0] in key
        ]

        # Write averaged cycle data to 'data'.
        averaged_scans = np.mean(all_scan_data, axis=0)
        if averaged_scans.size == 1:
            # on first scan in cycle
            averaged_scans = intensity_cps

        try:
            self._xps_dict["data"][entry][scan_key.split("_")[0]] = xr.DataArray(
                data=averaged_scans,
                coords={"energy": energy},
            )
        except ValueError:
            pass

        self._xps_dict["data"][entry][scan_key] = xr.DataArray(
            data=intensity_cps, coords={"energy": energy}
        )
        used_keys += ["data"]

        # Write raw intensities to 'detector'.
        self._xps_dict[detector_data_key] = intensity_raw

        # Write additional keys to region parent.
        for spectrum_key, value in spectrum.items():
            if spectrum_key not in used_keys:
                self._xps_dict[f"{region_parent}/{spectrum_key}"] = value

    def _get_units_for_key(self, unit_key: str):
        """
        Get correct units for a given key.

        Parameters
        ----------
        unit_key : str
           Key of type <mapping>:<spectrum_key>, e.g.
           detector/detector_voltage

        Returns
        -------
        str
            Unit for that unit_key.

        """
        try:
            return re.search(r"\[([A-Za-z0-9_]+)\]", unit_key).group(1)
        except AttributeError:
            try:
                return self.units[unit_key]
            except KeyError:
                return ""


class VamasParser:
    """A parser for reading vamas files."""

    def __init__(self):
        """Construct the vamas parser.

        Class attributes are a VamasHeader, which stores the vamas header
        attributes, blocks, which store the individual Block objects. Each
        block represents one spectrum, then there are several kinds of
        vamas attribute keys, which are used, depending on how the
        vamas file is formatted.
        """
        self.data: List[str] = []

        self.header = VamasHeader()
        self.blocks: List[VamasBlock] = []

        self.attrs = {
            "common_header": [
                "format_id",
                "institute_id",
                "instrumentModel_id",
                "operator_id",
                "experiment_id",
                "no_comment_lines",
            ],
            "exp_var": ["exp_var_label", "exp_var_unit"],
            "norm_header": [
                "scan_mode",
                "nr_regions",
                "nr_exp_var",
                "unknown_3",
                "unknown_4",
                "unknown_5",
                "unknown_6",
                "no_blocks",
            ],
            "map_header": [
                "scan_mode",
                "nr_regions",
                "nr_positions",
                "nr_x_coords",
                "nr_y_coords",
                "nr_exp_var",
                "unknown_3",
                "unknown_4",
                "unknown_5",
                "unknown_6",
                "no_blocks",
            ],
            "norm_block": [
                "block_id",
                "sample_id",
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "second",
                "no_hrs_in_advance_of_gmt",
                "no_comment_lines",
                "comment_lines",
                "technique",
                "exp_var_value",
                "source_label",
                "source_energy",
                "unknown_1",
                "unknown_2",
                "unknown_3",
                "source_analyzer_angle",
                "unknown_4",
                "analyzer_mode",
                "resolution",
                "magnification",
                "work_function",
                "target_bias",
                "analyzer_width_x",
                "analyzer_width_y",
                "analyzer_take_off_polar_angle",
                "analyzer_azimuth",
                "species_label",
                "transition_label",
                "particle_charge",
                "abscissa_label",
                "abscissa_units",
                "abscissa_start",
                "abscissa_step",
                "no_variables",
                "variable_label1",
                "variable_units1",
                "variable_label2",
                "variable_units2",
                "signal_mode",
                "dwell_time",
                "no_scans",
                "time_correction",
                "sample_angle_tilt",
                "sample_tilt_azimuth",
                "sample_rotation",
                "no_additional_params",
                "param_label_1",
                "param_unit_1",
                "param_value_1",
                "param_label_2",
                "param_unit_2",
                "param_value_2",
                "num_ord_values",
                "min_ord_value_1",
                "max_ord_value_1",
                "min_ord_value_2",
                "max_ord_value_2",
                "data_string",
            ],
            "map_block": [
                "block_id",
                "sample_id",
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "second",
                "no_hrs_in_advance_of_gmt",
                "no_comment_lines",
                "comment_lines",
                "technique",
                "x_coord",
                "y_coord",
                "exp_var_value",
                "source_label",
                "source_energy",
                "unknown_1",
                "unknown_2",
                "unknown_3",
                "fov_x",
                "fov_y",
                "source_analyzer_angle",
                "unknown_4",
                "analyzer_mode",
                "resolution",
                "magnification",
                "work_function",
                "target_bias",
                "analyzer_width_x",
                "analyzer_width_y",
                "analyzer_take_off_polar_angle",
                "analyzer_azimuth",
                "species_label",
                "transition_label",
                "particle_charge",
                "abscissa_label",
                "abscissa_units",
                "abscissa_start",
                "abscissa_step",
                "no_variables",
                "variable_label_1",
                "variable_units_1",
                "variable_label_2",
                "variable_units_2",
                "signal_mode",
                "dwell_time",
                "no_scans",
                "time_correction",
                "sample_angle_tilt",
                "sample_tilt_azimuth",
                "sample_rotation",
                "no_additional_params",
                "param_label_1",
                "param_unit_1",
                "param_value_1",
                "param_label_2",
                "param_unit_2",
                "param_value_2",
                "num_ord_values",
                "min_ord_value_1",
                "max_ord_value_1",
                "min_ord_value_2",
                "max_ord_value_2",
                "data_string",
            ],
        }

    def parse_file(self, file: Union[str, Path]):
        """Parse the vamas file into a list of dictionaries.

        Parameters
        ----------
        file: str
           The location and name of the vamas file to be parsed.
        """
        self._read_lines(file)
        self._parse_header()
        self._parse_blocks()
        return self.build_list()

    def _read_lines(self, file: Union[str, Path]):
        """Read in vamas text file."""
        with open(file, "rb") as vms_file:
            for line in vms_file:
                if line.endswith(b"\r\n"):
                    self.data += [line.decode("utf-8", errors="ignore").strip()]

    def _parse_header(self):
        """Parse the vama header into a VamasHeader object.

        The common_header_attr are the header attributes that are common
        to both types of Vamas format (NORM and MAP).
        Returns
        -------
        None.
        """
        for attr in self.attrs["common_header"]:
            setattr(self.header, attr, self.data.pop(0).strip())
        no_comment_lines = int(self.header.no_comment_lines)
        comments = []
        for _ in range(no_comment_lines):
            comments += [self.data.pop(0)]
        self.header.comment_lines = comments
        self.header.exp_mode = self.data.pop(0).strip()
        if self.header.exp_mode == "NORM":
            for attr in self.attrs["norm_header"]:
                setattr(self.header, attr, self.data.pop(0).strip())
                if attr == "nr_exp_var":
                    self._add_exp_var()

        elif self.header.exp_mode == "MAP":
            for attr in self.attrs["map_header"]:
                setattr(self.header, attr, self.data.pop(0).strip())
                if attr == "nr_exp_var":
                    self._add_exp_var()

        self.header.validate_types()

    def _add_exp_var(self):
        """Add experimental variable to header."""
        for _ in range(int(self.header.nr_exp_var)):
            for attr in self.attrs["exp_var"]:
                setattr(self.header, attr, self.data.pop(0).strip())

    def _parse_blocks(self):
        """Parse all (metadata) of Vamas blocks."""
        for _ in range(int(self.header.no_blocks)):
            self._parse_one_block()

    def _parse_one_block(self):
        """Parse one Vamas Block."""
        if self.header.exp_mode == "NORM":
            self.blocks += [self._parse_norm_block()]
        elif self.header.exp_mode == "MAP":
            self.blocks += [self._parse_map_block()]

    def _parse_norm_block(self):
        """
        Use this method when the NORM keyword is present.

        Returns
        -------
        block : vamas.Block object.
            A block represents one spectrum with its metadata.

        """
        # pylint: disable=too-many-statements
        block = VamasBlock()
        block.block_id = self.data.pop(0).strip()
        block.sample_id = self.data.pop(0).strip()
        block.year = int(self.data.pop(0).strip())
        block.month = int(self.data.pop(0).strip())
        block.day = int(self.data.pop(0).strip())
        block.hour = int(self.data.pop(0).strip())
        block.minute = int(self.data.pop(0).strip())
        block.second = int(self.data.pop(0).strip().split(".")[0])
        block.no_hrs_in_advance_of_gmt = int(self.data.pop(0).strip())
        block.no_comment_lines = int(self.data.pop(0).strip())
        for _ in range(block.no_comment_lines):
            block.comment_lines += [self.data.pop(0)]
        block.technique = self.data.pop(0).strip()
        for _ in range(int(self.header.nr_exp_var)):
            block.exp_var_value = self.data.pop(0).strip()
        block.source_label = self.data.pop(0).strip()
        block.source_energy = float(self.data.pop(0).strip())
        block.unknown_1 = self.data.pop(0).strip()
        block.unknown_2 = self.data.pop(0).strip()
        block.unknown_3 = self.data.pop(0).strip()
        block.source_analyzer_angle = float(self.data.pop(0).strip())
        block.unknown_4 = self.data.pop(0).strip()
        block.analyzer_mode = self.data.pop(0).strip()
        block.resolution = float(self.data.pop(0).strip())
        block.magnification = self.data.pop(0).strip()
        block.work_function = float(self.data.pop(0).strip())
        block.target_bias = float(self.data.pop(0).strip())
        block.analyzer_width_x = float(self.data.pop(0).strip())
        block.analyzer_width_y = float(self.data.pop(0).strip())
        block.analyzer_take_off_polar_angle = float(self.data.pop(0).strip())
        block.analyzer_azimuth = float(self.data.pop(0).strip())
        block.species_label = self.data.pop(0).strip()
        block.transition_label = self.data.pop(0).strip()
        block.particle_charge = int(self.data.pop(0).strip())

        if self.header.scan_mode == "REGULAR":
            block.abscissa_label = self.data.pop(0).strip()
            block.abscissa_units = self.data.pop(0).strip()
            block.abscissa_start = float(self.data.pop(0).strip())
            block.abscissa_step = float(self.data.pop(0).strip())

            block.no_variables = int(self.data.pop(0).strip())
            for var in range(block.no_variables):
                name = "variable_label_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())
                name = "variable_units_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())

        else:
            block.no_variables = int(self.data.pop(0).strip()) - 1
            block.abscissa_label = self.data.pop(0).strip()
            block.abscissa_units = self.data.pop(0).strip()
            for var in range(block.no_variables):
                name = "variable_label_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())
                name = "variable_units_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())

        block.signal_mode = self.data.pop(0).strip()
        block.dwell_time = float(self.data.pop(0).strip())
        block.no_scans = int(self.data.pop(0).strip())
        block.time_correction = self.data.pop(0).strip()
        block.sample_angle_tilt = float(self.data.pop(0).strip())
        block.sample_tilt_azimuth = float(self.data.pop(0).strip())
        block.sample_rotation = float(self.data.pop(0).strip())
        block.no_additional_params = int(self.data.pop(0).strip())
        for param in range(block.no_additional_params):
            name = "param_label_" + str(param + 1)
            setattr(block, name, self.data.pop(0))
            name = "param_unit_" + str(param + 1)
            setattr(block, name, self.data.pop(0))
            name = "param_value_" + str(param + 1)
            setattr(block, name, self.data.pop(0))
        block.num_ord_values = int(self.data.pop(0).strip())
        if self.header.scan_mode == "IRREGULAR":
            del self.data[:2]
        for var in range(block.no_variables):
            name = "min_ord_value_" + str(var + 1)
            setattr(block, name, float(self.data.pop(0).strip()))
            name = "max_ord_value_" + str(var + 1)
            setattr(block, name, float(self.data.pop(0).strip()))

        self._add_data_values(block)

        block.validate_types()
        return block

    def _parse_map_block(self):
        """
        Use this method when the MAP keyword is present.

        Returns
        -------
        block : vamas.Block object.
            A block represents one spectrum with its metadata.

        """
        # pylint: disable=too-many-statements
        block = VamasBlock()
        block.block_id = self.data.pop(0).strip()
        block.sample_id = self.data.pop(0).strip()
        block.year = int(self.data.pop(0).strip())
        block.month = int(self.data.pop(0).strip())
        block.day = int(self.data.pop(0).strip())
        block.hour = int(self.data.pop(0).strip())
        block.minute = int(self.data.pop(0).strip())
        block.second = int(self.data.pop(0).strip())
        block.no_hrs_in_advance_of_gmt = int(self.data.pop(0).strip())
        block.no_comment_lines = int(self.data.pop(0).strip())
        for _ in range(block.no_comment_lines):
            self.data.pop(0)
            block.comment_lines += [self.data.pop(0)]
        block.technique = self.data.pop(0).strip()
        block.x_coord = self.data.pop(0).strip()
        block.y_coord = self.data.pop(0).strip()
        block.exp_var_value = self.data.pop(0).strip()
        block.source_label = self.data.pop(0).strip()
        block.source_energy = float(self.data.pop(0).strip())
        block.unknown_1 = self.data.pop(0).strip()
        block.unknown_2 = self.data.pop(0).strip()
        block.unknown_3 = self.data.pop(0).strip()
        block.fov_x = self.data.pop(0).strip()
        block.fov_y = self.data.pop(0).strip()
        block.source_analyzer_angle = float(self.data.pop(0).strip())
        block.unknown_4 = self.data.pop(0).strip()
        block.analyzer_mode = self.data.pop(0).strip()
        block.resolution = float(self.data.pop(0).strip())
        block.magnification = self.data.pop(0).strip()
        block.work_function = float(self.data.pop(0).strip())
        block.target_bias = float(self.data.pop(0).strip())
        block.analyzer_width_x = float(self.data.pop(0).strip())
        block.analyzer_width_y = float(self.data.pop(0).strip())
        block.analyzer_take_off_polar_angle = float(self.data.pop(0).strip())
        block.analyzer_azimuth = float(self.data.pop(0).strip())
        block.species_label = self.data.pop(0).strip()
        block.transition_label = self.data.pop(0).strip()
        block.particle_charge = int(self.data.pop(0).strip())

        if self.header.scan_mode == "REGULAR":
            block.abscissa_label = self.data.pop(0).strip()
            block.abscissa_units = self.data.pop(0).strip()
            block.abscissa_start = float(self.data.pop(0).strip())
            block.abscissa_step = float(self.data.pop(0).strip())

            block.no_variables = int(self.data.pop(0).strip())
            for var in range(block.no_variables):
                name = "variable_label_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())
                name = "variable_units_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())

        else:
            block.no_variables = int(self.data.pop(0).strip()) - 1
            block.abscissa_label = self.data.pop(0).strip()
            block.abscissa_units = self.data.pop(0).strip()
            for var in range(block.no_variables):
                name = "variable_label_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())
                name = "variable_units_" + str(var + 1)
                setattr(block, name, self.data.pop(0).strip())

        block.signal_mode = self.data.pop(0).strip()
        block.dwell_time = float(self.data.pop(0).strip())
        block.no_scans = int(self.data.pop(0).strip())
        block.time_correction = self.data.pop(0).strip()
        block.sample_angle_tilt = float(self.data.pop(0).strip())
        block.sample_tilt_azimuth = float(self.data.pop(0).strip())
        block.sample_rotation = float(self.data.pop(0).strip())
        block.no_additional_params = int(self.data.pop(0).strip())
        for param in range(block.no_additional_params):
            name = "param_label_" + str(param + 1)
            setattr(block, name, self.data.pop(0))
            name = "param_unit_" + str(param + 1)
            setattr(block, name, self.data.pop(0))
            name = "param_value_" + str(param + 1)
            setattr(block, name, self.data.pop(0))
        block.num_ord_values = int(self.data.pop(0).strip())
        if self.header.scan_mode == "IRREGULAR":
            del self.data[:2]
        for var in range(block.no_variables):
            name = "min_ord_value_" + str(var + 1)
            setattr(block, name, float(self.data.pop(0).strip()))
            name = "max_ord_value_" + str(var + 1)
            setattr(block, name, float(self.data.pop(0).strip()))

        self._add_data_values(block)

        return block

    def _add_data_values(self, block: VamasBlock):
        """Add data values to a Vamas data block."""
        if self.header.scan_mode == "REGULAR":
            self._add_regular_data(block)
        elif self.header.scan_mode == "IRREGULAR":
            self._add_irregular_data(block)

    def _add_regular_data(self, block: VamasBlock):
        """Parse data with regularly spaced energy axis."""
        data_dict: Dict[str, List] = {}

        start = float(block.abscissa_start)
        step = float(block.abscissa_step)
        num = int(block.num_ord_values / block.no_variables)
        energy = [round(start + i * step, 2) for i in range(num)]

        if block.abscissa_label == "binding energy":
            energy.reverse()

        setattr(block, "x", energy)

        for var in range(block.no_variables):
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)
            data_dict[name] = []

        data_array = list(np.array(self.data[: block.num_ord_values], dtype=float))

        self.data = self.data[block.num_ord_values :]

        for var in range(block.no_variables):
            max_var = block.no_variables
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)
            data_array_slice = data_array[var::max_var]
            data_dict[name] = data_array_slice
            setattr(block, name, data_dict[name])

    def _add_irregular_data(self, block: VamasBlock):
        """Parse data with regularly spaced energy axis."""
        data_dict: Dict[str, List] = {}

        block_data = list(np.array(self.data[: block.num_ord_values], dtype=float))

        energy = block_data[:: block.no_variables + 1]
        if block.abscissa_label == "binding energy":
            energy.reverse()

        setattr(block, "x", energy)
        block.abscissa_start = float(min(energy))
        block.abscissa_step = float(get_minimal_step(energy))

        for var in range(block.no_variables):
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)
            data_dict[name] = []

        for var in range(block.no_variables):
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)
            data_array_slice = block_data[var + 1 :: block.no_variables + 1]
            data_dict[name] = data_array_slice
            setattr(block, name, data_dict[name])

        self.data = self.data[block.num_ord_values :]  # + block.no_variables :]

    def _get_scan_numbers_for_spectra(self, spectra: List[Dict]):
        """
        For a flat list of spectra dictionaries, group the spectra
        by group name and spectrum type and iteratively give them
        scan numbers.

        Parameters
        ----------
        spectra : list
            List of dicts with each dict containing data and metadata
            for one spectrum.

        Returns
        -------
        flattened_spectra : list
            Same list of dicts, but each spectrum gets a scan number.

        """
        grouped_spectra = [
            list(y)
            for x, y in groupby(
                sorted(spectra, key=lambda x: (x["group_name"], x["spectrum_type"])),
                lambda x: (x["group_name"], x["spectrum_type"]),
            )
        ]

        for group in grouped_spectra:
            for i, spectrum in enumerate(group):
                spectrum["scan_no"] = i

        flattened_spectra = [
            spectrum for group in grouped_spectra for spectrum in group
        ]

        return flattened_spectra

    def handle_header_comments(self, comment_list: List[str]):
        """Handle comments (incl. Casa info) for the header."""
        comments = {}

        special_keys = {
            "Casa Info Follows": self._handle_casa_header,
            "SpecsLab Prodigy": self._handle_prodigy_header,
            # "SOFH": self._handle_phi_header,
        }

        for keyword, handle_func in special_keys.items():
            if any(keyword in line for line in comment_list):
                index = [i for i, line in enumerate(comment_list) if keyword in line][0]

                if keyword == "Casa Info Follows":
                    special_comments = comment_list[index]
                    comment_list = comment_list[index + 1 :]

                if keyword == "SpecsLab Prodigy":
                    special_comments = comment_list[index]
                    comment_list = comment_list[index + 1 :]

                if keyword == "SOFH":
                    end_index = [
                        i for i, line in enumerate(comment_list) if "EOFH" in line
                    ][0]
                    special_comments = comment_list[index : end_index + 1]  # type: ignore[assignment]
                    del comment_list[index : end_index + 1]

                comments.update(handle_func(special_comments))

        # Handle non-special comments.
        for line in comment_list:
            for sep in ("=", ":"):
                try:
                    key, value = [part.strip(" ") for part in line.split(sep, 1)]
                    comments[to_snake_case(key)] = value
                except ValueError:
                    continue

        return comments

    def _handle_casa_header(self, comment_line: str):
        """Get information about CasaXPS version."""
        return {
            "casa_version": comment_line.split("Casa Info Follows CasaXPS Version")[
                1
            ].strip()
        }

    def _handle_prodigy_header(self, comment_line: str):
        """Get information about SpecsLab Prodigy version."""
        return {"prodigy_version": comment_line.split("Version")[1].strip()}

    # =============================================================================
    #     def _handle_phi_header(self, comment_list: List[str]):
    #         """Get metadta from Phi system."""
    #         phi_parser = PhiParser()
    #         phi_parser.parse_header_into_metadata(comment_list)
    #
    #         phi_comments = phi_parser.metadata.dict()
    #
    #         regions = phi_parser.parse_spectral_regions(comment_list)
    #         areas = phi_parser.parse_spatial_areas(comment_list)
    #
    #         for region in regions:
    #             for area in areas:
    #                 concatenated = {**region.dict(), **area.dict()}
    #
    #             phi_comments.update(concatenated)
    #
    #         return phi_comments
    # =============================================================================

    def handle_block_comments(self, comment_list: List[str]):
        """Handle comments (incl. Casa fitting) for each block."""
        comments = {}

        if "Casa Info Follows" in comment_list[0]:
            # Get all processing and fitting data from Casa comments.
            casa = CasaProcess()
            casa_data = casa.process_comments(comment_list)

            comments.update(casa_data)

            no_of_casa_lines = 1

            for number in (
                "n_alignments",
                "n_unknown_processes",
                "n_regions",
                "n_components",
            ):
                occurence = getattr(casa, number)
                no_of_casa_lines += 1
                if occurence >= 1:
                    no_of_casa_lines += occurence

            non_casa_comments = comment_list[no_of_casa_lines:]

        else:
            non_casa_comments = comment_list

        for line in non_casa_comments:
            for sep in ("=", ":"):
                try:
                    key, value = [part.strip(" ") for part in line.split("=", 1)]
                    comments[to_snake_case(key)] = value
                except ValueError:
                    continue
        return comments

    def build_list(self):
        """
        Construct a list of dictionaries from the Vamas objects

        Returns
        -------
        List
            Each list element is a dictionary with the data and
            metadata of one spectrum.

        """
        group_id = -1
        temp_group_name = ""
        spectra = []

        header_dict = {to_snake_case(k): v for (k, v) in self.header.dict().items()}
        del header_dict["comment_lines"]
        header_dict.update(self.handle_header_comments(self.header.comment_lines))

        for spectrum_id, block in enumerate(self.blocks):
            group_name = block.sample_id
            # This set of conditions detects if the group name has changed.
            # If it has, then it increments the group_idx.
            if group_name != temp_group_name:
                temp_group_name = group_name
                group_id += 1

            spectrum_type = str(block.species_label + block.transition_label)

            settings = {
                "region": block.block_id,
                "sample_name": block.sample_id,
                "analysis_method": block.technique,
                "source_label": block.source_label,
                "excitation_energy": block.source_energy,
                "source_analyzer_angle": block.source_analyzer_angle,
                "scan_mode": block.analyzer_mode,
                "pass_energy": block.resolution,
                "magnification": block.magnification,
                "work_function": block.work_function,
                "target_bias": block.target_bias,
                "analysis_width_x": block.analyzer_width_x,
                "analysis_width_y": block.analyzer_width_y,
                "analyzer_take_off_polar": block.analyzer_take_off_polar_angle,
                "analyzer_take_off_azimuth": block.analyzer_azimuth,
                "element": block.species_label,
                "transition": block.transition_label,
                "particle_charge": block.particle_charge,
                "x_label": block.abscissa_label,
                "x_units": block.abscissa_units,
                "start_energy": block.abscissa_start,
                "step_size": block.abscissa_step,
                "y_labels_1": block.variable_label_1,
                "y_units_1": block.variable_units_1,
                "y_labels_2": block.variable_label_2,
                "y_units_2": block.variable_units_2,
                "signal_mode": block.signal_mode,
                "dwell_time": block.dwell_time,
                "time_correction": block.time_correction,
                "sample_normal_polarangle_tilt": block.sample_angle_tilt,
                "sample_tilt_azimuth": block.sample_tilt_azimuth,
                "sample_rotation_angle": block.sample_rotation,
                "n_values": int(block.num_ord_values / block.no_variables),
            }
            settings.update(header_dict)

            comment_dict = self.handle_block_comments(block.comment_lines)
            settings.update(comment_dict)

            # Convert the native time format to the datetime string
            # in the ISO 8601 format
            tzinfo = datetime.timezone(
                datetime.timedelta(hours=block.no_hrs_in_advance_of_gmt)
            )
            try:
                date_time = datetime.datetime(
                    block.year,
                    block.month,
                    block.day,
                    block.hour,
                    block.minute,
                    block.second,
                    tzinfo=tzinfo,
                )
            except ValueError:
                date_time = datetime.datetime.min

            data = {"x": block.x}
            for var in range(int(block.no_variables)):
                if var == 0:
                    key = "y"

                    data["y"] = getattr(block, "y")

                    if block.variable_label_1 in ["Intensity", "counts"]:
                        y_cps = [np.round(y / block.dwell_time, 2) for y in block.y]
                        data["y_cps"] = y_cps

                else:
                    key = "y" + str(var)
                    data[key] = getattr(block, key)

            spec_dict = {
                "time_stamp": date_time,
                "group_name": group_name,
                "group_id": group_id,
                "spectrum_type": spectrum_type,
                "spectrum_id": spectrum_id,
                "scans": block.no_scans,
                "data": data,
            }
            spec_dict.update(settings)
            spectra += [spec_dict]

        spectra = self._get_scan_numbers_for_spectra(spectra)

        return spectra
