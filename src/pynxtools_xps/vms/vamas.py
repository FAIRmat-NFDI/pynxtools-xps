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

from copy import deepcopy
import datetime
from pathlib import Path
from typing import Any, Dict, List, Union
import warnings
from itertools import groupby
import xarray as xr
import numpy as np

from pynxtools_xps.vms.vamas_data_model import (
    VamasHeader,
    VamasBlock,
    ExpVariable,
    VamasAdditionalParam,
    OrdinateValue,
)
from pynxtools_xps.vms.vamas_comment_handler import handle_comments

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    construct_data_key,
    convert_pascal_to_snake,
    get_minimal_step,
    check_for_allowed_in_list,
    re_map_keys,
    re_map_values,
    drop_unused_keys,
    update_dict_without_overwrite,
)
from pynxtools_xps.value_mappers import (
    convert_measurement_method,
    convert_energy_type,
    convert_energy_scan_mode,
    get_units_for_key,
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

UNITS: dict = {
    "beam_xray/excitation_energy": "eV",
    "beam_xray/extent": "um",
    "beam_xray/particle_charge": "C",
    "beam_xray/source_beam_width_x": "um",
    "beam_xray/source_beam_width_y": "um",
    "beam_xray/source_polar_angle": "degree",
    "beam_xray/source_azimuth_angle": "degree",
    "source_xray/source_power": "W",
    "source_sputter/sputter_ion_charge": "C",
    "source_sputter/sputter_source_energy": "eV",
    "source_sputter/sputter_source_beam_current": "A",
    "source_sputter/sputter_source_width_x": "um",
    "source_sputter/sputter_source_width_y": "um",
    "source_sputter/sputter_source_incidence_polar_angle": "degree",
    "source_sputter/sputter_source_azimuth_angle": "degree",
    "electronanalyser/analyser_take_off_polar_angle": "degree",
    "electronanalyser/analyser_take_off_azimuth_angle": "degree",
    "electronanalyser/analysis_width_x": "m",
    "electronanalyser/analysis_width_y": "m",
    "electronanalyser/target_bias": "V",
    "electronanalyser/time_correction": "s",
    "electronanalyser/work_function": "eV",
    "collectioncolumn/spatial_acceptance": "um",
    "energydispersion/pass_energy": "eV",
    "energydispersion/differential_width_aes": "eV",
    "detector/dwell_time": "s",
    "sample/sample_normal_polar_angle_of_tilt": "degree ",
    "sample/sample_normal_tilt_azimuth_angle": "degree",
    "sample/sample_rotation_angle": "degree",
    "data/start_energy": "eV",
    "data/step_size": "eV",
}


class VamasMapper(XPSMapper):
    """
    Class for restructuring .txt data file from
    Vamas format into python dictionary.
    """

    config_file = "config_vms.json"

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

        template_key_map = {
            "user": [],
            "instrument": [],
            "beam_xray": [
                "excitation_energy",
                "extent",
                "particle_charge",
                "source_polar_angle",
                "source_azimuth_angle",
            ],
            "source_xray": [
                "source_label",
                "source_power",
                "field_of_view_x",
                "field_of_view_y",
            ],
            "source_sputter": [
                "sputter_ion_atomic_number",
                "sputter_ion_num_atoms",
                "sputter_ion_charge",
                "sputter_source_energy",
                "sputter_source_beam_current",
                "sputter_source_width_x",
                "sputter_source_width_y",
                "sputter_source_incidence_polar_angle",
                "sputter_source_azimuth_angle",
            ],
            "electronanalyser": [
                "analyser_take_off_polar_angle",
                "analyser_take_off_azimuth_angle",
                "analysis_width_x",
                "analysis_width_y",
                "target_bias",
                "time_correction",
                "work_function",
            ],
            "collectioncolumn": [
                "magnification",
                "spatial_acceptance",
            ],
            "energydispersion": [
                "scan_mode",
                "pass_energy",
                "differential_width_aes",
            ],
            "detector": [
                "signal_mode",
                "dwell_time",
            ],
            "manipulator": [],
            "sample": [
                "sample_name",
                "sample_normal_polar_angle_of_tilt",
                "sample_normal_tilt_azimuth_angle",
                "sample_rotation_angle",
            ],
            "data": [
                "energy_label",
                "energy_units",
                "y_labels_1",
                "y_units_1",
                "y_labels_2",
                "y_units_2",
                "n_values",
                "start_energy",
                "step_size",
            ],
            "profiling": [
                "first_linescan_x_start",
                "first_linescan_y_start",
                "first_linescan_x_end",
                "first_linescan_y_end",
                "last_linescan_x_end",
                "last_linescan_y_end",
            ],
            "region": [
                "format_id",
                "institute_id",
                "instrument_model_id",
                "operator_id",
                "experiment_id",
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
                "x_coord",
                "y_coord",
            ],
        }

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum, template_key_map)

    def _update_xps_dict_with_spectrum(
        self, spectrum: Dict[str, Any], template_key_map: Dict[str, str]
    ):
        """Map one spectrum from raw data to NXmpes-ready dict."""
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []
        for part in ["group_name", "spectrum_type"]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        if len(entry_parts) == 1:
            entry_parts += [spectrum["time_stamp"]]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        instrument_parent = f"{entry_parent}/instrument"
        analyser_parent = f"{instrument_parent}/electronanalyser"

        path_map: Dict[str, str] = {
            "user": f"{entry_parent}/user",
            "instrument": f"{instrument_parent}",
            "beam_xray": f"{instrument_parent}/beam_xray",
            "source_xray": f"{instrument_parent}/source_xray",
            "source_sputter": f"{instrument_parent}/source_sputter",
            "electronanalyser": f"{analyser_parent}",
            "collectioncolumn": f"{analyser_parent}/collectioncolumn",
            "energydispersion": f"{analyser_parent}/energydispersion",
            "detector": f"{analyser_parent}/detector",
            "manipulator": f"{instrument_parent}/manipulator",
            "sample": f"{entry_parent}/sample",
            "data": f"{entry_parent}/data",
            "process": f"{entry_parent}/process",
            "peak_fitting": f"{entry_parent}/peak_fitting",
            "profiling": f"{entry_parent}/profiling",
            "region": f"{entry_parent}/region",
        }

        used_keys = []

        for grouping, spectrum_keys in template_key_map.items():
            root = path_map[grouping]
            for spectrum_key in spectrum_keys:
                mpes_key = spectrum_key.rsplit(" ", 1)[0]
                try:
                    self._xps_dict[f"{root}/{mpes_key}"] = spectrum[spectrum_key]
                except KeyError:
                    pass

                unit_key = f"{grouping}/{spectrum_key}"
                units = get_units_for_key(unit_key, UNITS)
                if units is not None:
                    self._xps_dict[f"{root}/{mpes_key}/@units"] = units
                used_keys += [spectrum_key]

        # Write process data
        process_key_map: Dict[str, List[str]] = {
            "process": ["energy_calibrations", "intensity_calibrations", "smoothings"],
            "peak_fitting": ["regions", "components"],
        }

        for grouping, process_key_list in process_key_map.items():
            root = path_map[grouping]

            for spectrum_key in process_key_list:
                process_keys = [key for key in spectrum if key.startswith(spectrum_key)]
                for process_key in process_keys:
                    try:
                        process_key = process_key.replace("_units", "/@units")
                        self._xps_dict[f"{root}/{process_key}"] = spectrum[process_key]
                        used_keys += [process_key]
                    except KeyError:
                        pass

        # Create key for writing to data.
        scan_key = construct_data_key(spectrum)

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
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

        # Write raw intensities to 'cycle0'.
        self._xps_dict["data"][entry][f"{scan_key}_chan0"] = xr.DataArray(
            data=intensity_raw, coords={"energy": energy}
        )

        # Write additional keys to region parent.
        for spectrum_key, value in spectrum.items():
            if spectrum_key not in used_keys:
                self._xps_dict[f"{entry_parent}/misc/{spectrum_key}"] = value


KEY_MAP = {
    "block_id": "region",
    "sample_id": "sample_name",
    "technique": "analysis_method",
    "source_energy": "excitation_energy",
    "analyser_mode": "scan_mode",
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

VALUE_MAP = {
    "analysis_method": convert_measurement_method,
    "energy_label": convert_energy_type,
    "scan_mode": convert_energy_scan_mode,
}


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

    def parse_file(self, file: Union[str, Path], **kwargs):
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
                if line.endswith((b"\r\n", b"\n")):
                    self.data += [line.decode("utf-8", errors="ignore").strip()]

    def _extract_n_lines_to_list(self, number_of_lines):
        """Ectract n number of lines to a list."""
        extracted = []
        for _ in range(number_of_lines):
            extracted += [self.data.pop(0)]

        return extracted

    def _setattr_with_datacls_type(self, datacls, attr, value):
        """
        Set attribute of dataclass instance with the field_type
        defined in the dataclass.
        """
        field_type = type(getattr(datacls, attr))
        setattr(datacls, attr, field_type(value))

    def _parse_header(self):
        """Parse the Vamas header into a VamasHeader object.

        The common_header_attr are the header attributes that are common
        to all types of Vamas experiment modes.
        Returns
        -------
        None.
        """
        common_attrs = [
            "format_id",
            "institute_id",
            "instrument_model_id",
            "operator_id",
            "experiment_id",
            "num_comment_lines",
        ]
        map_attrs = ["num_analysis_positions", "num_x_coords", "num_y_coords"]

        for attr in common_attrs:
            self._setattr_with_datacls_type(self.header, attr, self.data.pop(0).strip())

        num_comment_lines = int(self.header.num_comment_lines)
        self.header.comment_lines = self._extract_n_lines_to_list(num_comment_lines)

        self.header.exp_mode = check_for_allowed_in_list(
            self.data.pop(0).strip(), EXP_MODES
        )
        self.header.scan_mode = self.data.pop(0).strip()

        if self.header.exp_mode in ["MAP", "MAPDP", "NORM", "SDP"]:
            self.header.num_spectral_regions = int(self.data.pop(0).strip())
        else:
            delattr(self.header, "num_spectral_regions")

        for attr in map_attrs:
            if self.header.exp_mode in ["MAP", "MAPDP"]:
                self._setattr_with_datacls_type(
                    self.header, attr, self.data.pop(0).strip()
                )
            else:
                delattr(self.header, attr)

        self.header.num_exp_var = int(self.data.pop(0).strip())

        for exp_var_no in range(int(self.header.num_exp_var)):
            exp_var = ExpVariable()
            for attr in ["label", "unit"]:
                self._setattr_with_datacls_type(exp_var, attr, self.data.pop(0).strip())
                setattr(
                    self.header, f"exp_var_{exp_var_no}_{attr}", getattr(exp_var, attr)
                )

        self.header.num_entries_in_inclusion_list = int(self.data.pop(0).strip())
        self.header.inclusion_list = self._extract_n_lines_to_list(
            self.header.num_entries_in_inclusion_list
        )

        self.header.num_manually_entered_items_in_block = int(self.data.pop(0).strip())
        self.header.manually_entered_items_in_block = self._extract_n_lines_to_list(
            self.header.num_manually_entered_items_in_block
        )

        self.header.num_future_upgrade_exp_entries = int(self.data.pop(0).strip())
        self.header.num_future_upgrade_block_entries = int(self.data.pop(0).strip())
        self.header.future_upgrade_exp_entries = self._extract_n_lines_to_list(
            self.header.num_future_upgrade_exp_entries
        )

        self.header.num_blocks = int(self.data.pop(0).strip())
        self.header.validate_types()

    def _parse_blocks(self):
        """Parse all (metadata) of Vamas blocks."""
        for _ in range(int(self.header.num_blocks)):
            self.blocks += [self._parse_one_block()]

    def _parse_one_block(self):
        """
        Parse one Vamas Block.

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
        block.num_comment_lines = int(self.data.pop(0).strip())
        for _ in range(block.num_comment_lines):
            block.comment_lines += [self.data.pop(0)]
        block.technique = check_for_allowed_in_list(
            self.data.pop(0).strip(), ALLOWED_TECHNIQUES
        )

        for attr in ["x_coord", "y_coord"]:
            if self.header.exp_mode in ["MAP", "MAPDP"]:
                self._setattr_with_datacls_type(block, attr, self.data.pop(0).strip())
            else:
                delattr(block, attr)

        for _ in range(int(self.header.num_exp_var)):
            block.exp_var_value = self.data.pop(0).strip()
        block.source_label = self.data.pop(0).strip()

        for attr in [
            "sputter_ion_atomic_number",
            "sputter_ion_num_atoms",
            "sputter_ion_charge",
        ]:
            if self.header.exp_mode in [
                "MAPDP",
                "MAPSVDP",
                "SDP",
                "SDPSV",
            ] or block.technique in [
                "FABMS",
                "FABMS energy spec",
                "ISS",
                "SIMS",
                "SIMS energy spec",
                "SNMS",
                "SNMS energy spec",
            ]:
                self._setattr_with_datacls_type(block, attr, self.data.pop(0).strip())
            else:
                delattr(block, attr)
        block.source_energy = float(self.data.pop(0).strip())
        block.source_power = self.data.pop(0).strip()
        block.source_beam_width_x = self.data.pop(0).strip()
        block.source_beam_width_x = self.data.pop(0).strip()

        for attr in ["field_of_view_x", "field_of_view_y"]:
            if self.header.exp_mode in ["MAP", "MAPDP", "MAPSV", "MAPSVDP", "SEM"]:
                self._setattr_with_datacls_type(block, attr, self.data.pop(0).strip())
            else:
                delattr(block, attr)

        for attr in [
            "first_linescan_x_start",
            "first_linescan_y_start",
            "first_linescan_x_end",
            "first_linescan_y_end",
            "last_linescan_x_end",
            "last_linescan_y_end",
        ]:
            if self.header.exp_mode in ["MAP", "MAPSVDP", "SEM"]:
                self._setattr_with_datacls_type(block, attr, self.data.pop(0).strip())
            else:
                delattr(block, attr)

        block.source_polar_angle = float(self.data.pop(0).strip())
        block.source_azimuth_angle = float(self.data.pop(0).strip())
        block.analyser_mode = self.data.pop(0).strip()
        block.resolution = float(self.data.pop(0).strip())

        if block.technique == "AES diff":
            self.header.differential_width_aes = float(self.data.pop(0).strip())
        else:
            delattr(block, "differential_width_aes")

        block.magnification = float(self.data.pop(0).strip())
        block.work_function = float(self.data.pop(0).strip())
        block.target_bias = float(self.data.pop(0).strip())
        block.analysis_width_x = float(self.data.pop(0).strip())
        block.analysis_width_y = float(self.data.pop(0).strip())
        block.analyser_take_off_polar_angle = float(self.data.pop(0).strip())
        block.analyser_take_off_azimuth_angle = float(self.data.pop(0).strip())
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
        block.time_correction = float(self.data.pop(0).strip())

        for attr in [
            "sputter_source_energy",
            "sputter_source_beam_current",
            "sputter_source_width_x",
            "sputter_source_width_y",
            "sputter_source_incidence_polar_angle",
            "sputter_source_azimuth_angle",
        ]:
            if self.header.exp_mode in [
                "MAPDP",
                "MAPSVDP",
                "SDP",
                "SDPSV",
            ] and block.technique in [
                "AES",
                "AES diff",
                "AES dir",
                "EDX",
                "ELS",
                "UPS",
                "XPS",
                "XRF",
            ]:
                self._setattr_with_datacls_type(block, attr, self.data.pop(0).strip())
            else:
                delattr(block, attr)

        block.sample_normal_polar_angle_of_tilt = float(self.data.pop(0).strip())
        block.sample_normal_tilt_azimuth_angle = float(self.data.pop(0).strip())
        block.sample_rotation_angle = float(self.data.pop(0).strip())
        block.no_additional_params = int(self.data.pop(0).strip())

        for param_no in range(block.no_additional_params):
            param = VamasAdditionalParam()
            for attr in ["label", "unit", "value"]:
                self._setattr_with_datacls_type(param, attr, self.data.pop(0).strip())
                setattr(block, f"param_{param_no}_{attr}", getattr(param, attr))

        block.future_upgrade_block_entries = self._extract_n_lines_to_list(
            self.header.num_future_upgrade_block_entries
        )

        block.num_ord_values = int(self.data.pop(0).strip())
        if self.header.scan_mode == "IRREGULAR":
            del self.data[:2]

        for var_no in range(block.no_variables):
            var = OrdinateValue()
            for attr in ["min_ord_value", "max_ord_value"]:
                self._setattr_with_datacls_type(var, attr, self.data.pop(0).strip())
                setattr(block, f"{attr}_{var_no + 1}", getattr(var, attr))

        self._add_data_values(block)

        block.validate_types()

        return block

    def _add_data_values(self, block: VamasBlock):
        """Add data values to a Vamas data block."""
        if self.header.scan_mode == "REGULAR":
            self._add_regular_data(block)
        elif self.header.scan_mode == "IRREGULAR":
            self._add_irregular_data(block)
        elif self.header.scan_mode == "MAPPING":
            self._add_mapping_data(block)

    def _add_regular_data(self, block: VamasBlock):
        """Parse data with regularly spaced energy axis."""
        data_dict: Dict[str, np.ndarray] = {}

        start = float(block.abscissa_start)
        step = float(block.abscissa_step)
        num = int(block.num_ord_values / block.no_variables)
        energy = np.array([round(start + i * step, 2) for i in range(num)])

        if block.abscissa_label == "binding energy":
            energy = np.flip(energy)

        setattr(block, "x", energy)

        for var in range(block.no_variables):
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)

        data_array = np.array(self.data[: block.num_ord_values], dtype=float)

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
        data_dict: Dict[str, np.ndarray] = {}

        block_data = np.array(self.data[: block.num_ord_values], dtype=float)

        energy = block_data[:: block.no_variables + 1]
        if block.abscissa_label == "binding energy":
            energy = np.flip(energy)

        setattr(block, "x", energy)
        block.abscissa_start = float(min(energy))
        block.abscissa_step = float(get_minimal_step(energy))

        for var in range(block.no_variables):
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)

        for var in range(block.no_variables):
            if var == 0:
                name = "y"
            else:
                name = "y" + str(var)
            data_array_slice = block_data[var + 1 :: block.no_variables + 1]
            data_dict[name] = data_array_slice
            setattr(block, name, data_dict[name])

        self.data = self.data[block.num_ord_values :]

    def _add_mapping_data(self, block: VamasBlock):
        """
        Parse data in mapping format.

        TBD!
        """
        pass

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

        header_dict = {
            convert_pascal_to_snake(k): v for (k, v) in self.header.dict().items()
        }
        del header_dict["comment_lines"]

        header_comments = handle_comments(
            self.header.comment_lines, comment_type="header"
        )

        update_dict_without_overwrite(header_dict, header_comments)

        for spectrum_id, block in enumerate(self.blocks):
            group_name = block.sample_id
            # This set of conditions detects if the group name has changed.
            # If it has, then it increments the group_idx.
            if group_name != temp_group_name:
                temp_group_name = group_name
                group_id += 1

            spectrum_type = str(block.species_label + block.transition_label)

            settings = {
                convert_pascal_to_snake(k): v for (k, v) in block.dict().items()
            }

            update_dict_without_overwrite(settings, header_dict)

            settings["n_values"] = int(block.num_ord_values / block.no_variables)

            # Remap to the MPES-prefered keys and values .
            re_map_keys(settings, KEY_MAP)
            re_map_values(settings, VALUE_MAP)

            comment_dict = handle_comments(block.comment_lines, comment_type="block")

            re_map_keys(comment_dict, KEY_MAP)
            re_map_values(comment_dict, VALUE_MAP)

            update_dict_without_overwrite(settings, comment_dict)

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
                date_time = datetime.datetime(1, 1, 1, 0, 0, 0, tzinfo=tzinfo)

            if isinstance(date_time, datetime.datetime):
                date_time = date_time.isoformat()

            # Map x-y values to 2D lists.
            settings["extent"] = np.array(
                [
                    settings["source_beam_width_x"],
                    settings["source_beam_width_y"],
                ],
                dtype=float,
            )
            settings["spatial_acceptance"] = np.array(
                [
                    settings["analysis_width_x"],
                    settings["analysis_width_y"],
                ],
                dtype=float,
            )

            data = {"x": block.x}

            for var in range(int(block.no_variables)):
                if var == 0:
                    key = "y"

                    data["y"] = getattr(block, "y")
                    del settings["y"]

                    if block.variable_label_1 in ["Intensity", "counts"]:
                        y_cps = [np.round(y / block.dwell_time, 2) for y in block.y]
                        data["y_cps"] = np.array(y_cps)

                else:
                    key = "y" + str(var)
                    data[key] = getattr(block, key)
                    del settings[key]

            spec_dict = {
                "time_stamp": date_time,
                "group_name": group_name,
                "group_id": group_id,
                "spectrum_type": spectrum_type,
                "spectrum_id": spectrum_id,
                "scans": block.no_scans,
                "data": data,
            }

            remove_keys = [
                "comment_lines",
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "second",
                "no_hrs_in_advance_of_gmt",
                "source_beam_width_x",
                "source_beam_width_y",
                "x",
            ]

            drop_unused_keys(settings, remove_keys)

            spec_dict.update(settings)
            spectra += [spec_dict]

        spectra = self._get_scan_numbers_for_spectra(spectra)

        return spectra
