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
Parser for reading XPS (X-ray Photoelectron Spectroscopy) data from
Phi PHI VersaProbe 4 instruments (.spe or .pro format), to be passed to
MPES nxdl (NeXus Definition Language) template.
"""

import re
import warnings
import copy
import datetime
import logging
import struct
from typing import Any, Dict, List, Union
from pathlib import Path
import pytz
import xarray as xr
import numpy as np

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    safe_arange_with_edges,
    convert_pascal_to_snake,
    split_value_and_unit,
)

from pynxtools_xps.value_mappers import (
    convert_energy_scan_mode,
    convert_measurement_method,
    convert_bool,
    convert_units,
)

from pynxtools_xps.phi.phi_data_model import (
    PhiMetadata,
    PhiSpectralRegion,
    PhiSpatialArea,
)

logger = logging.getLogger(__name__)

SETTINGS_MAP: Dict[str, str] = {
    "FileDesc": "file_description",
    "acq_filename": "acquisition_filename",
    "acq_file_date": "acquisition_file_date",
    "institution": "vendor",
    "operator": "user_name",
    "experiment_i_d": "experiment_id",
    "analyzer_work_fcn": "analyzer_work_function",
    "analyzer_retard_gain": "analyzer_retardation_gain",
    "reg_image_interval": "register_image_interval",
    "reg_image_mode": "register_image_mode",
    "reg_image_last": "register_image_last",
    "platen_i_d": "platen_id",
    "s_x_i_filename": "sxi_filename",
    "intensity_recal": "intensity_recalibration",
    "intensity_cal_coeff": "intensity_calibration_coefficients",
    "energy_recal": "energy_recalibration",
    "s_c_a_multiplier_voltage": "sca_multiplier_voltage",
    "C60IonGun": "c60_ion_gun",
    "t_f_c_parameters": "tfc_parameters",
    "image_size_x_y": "image_size",
    "float_volt": "float_voltage",
    "float_enable": "float_enabled",
    "grid_volt": "grid_voltage",
    "condensor_volt": "condenser_lens_voltage",
    "objective_volt": "objective_lens_voltage",
    "bend_volt": "bend_voltage",
    "neutral_float_volt": "neutral_float_voltage",
    "neutral_float_enable": "neutral_float_enabled",
    "neutral_grid_volt": "neutral_grid_voltage",
    "neutral_condensor_volt": "neutral_condenser_lens_voltage",
    "neutral_objective_volt": "neutral_objective_lens_voltage",
    "neutral_bend_volt": "neutral_bend_voltage",
    "no_d_p_data_cyc": "no_depth_profile_cycles",
    "no_pre_sputter_cyc": "no_pre_sputter_cycles",
    "sample_rotation": "profiling_sample_rotation",
    "depth_recal": "profiling_depth_recalibration",
    "sputter_mode": "profiling_sputter_mode",
    "no_depth_reg": "profiling_no_depth_regions",
    "depth_cal_def": "depth_calibration_definition",
    "analyser_mode": "energy_scan_mode",
    "surv_num_cycles": "survey_num_of_cycles",
    "surv_time_per_step": "survey_dwell_time",
    "no_spectral_reg_full": "no_spectral_regions_full",
    "no_spectral_region": "no_spectral_regions",
    "spectral_reg_def_full": "spectral_region_definition_full",
    "spectral_reg_def2_full": "spectral_region_definition2_full",
    "spectral_reg_background_full": "spectral_region_background_full",
    "spectral_reg_hero_full": "spectral_region_hero_full",
    "spectral_reg_i_r_full": "spectral_region_ir_full",
    "no_spectral_reg": "no_spectral_regions",
    "spectral_reg_def": "spectral_region_definition",
    "spectral_reg_def2": "spectral_region_definition2",
    "spectral_reg_background": "spectral_region_background",
    "spectral_reg_hero": "spectral_region_hero",
    "spectral_reg_i_r": "spectral_region_ir",
    "no_spatial_area": "no_spatial_areas",
    "spatial_area_def": "spatial_area_definition",
    "spatial_h_r_photo_cor": "spatial_hr_photo_correction",
    "xray_offset_in_um": "xray_offset",
    "xray_mag_factor": "xray_magnification_factor",
    "xray_rotation_in_deg": "xray_rotation",
    "xray_setting": "xray_settings",
    "neutralizer_current": "flood_gun_current",
    "neutralizer_energy": "flood_gun_energy",
    "flood_gun_filament": "flood_gun_filament_current",
    "flood_gun_extractor": "flood_gun_extractor_voltage",
    "detector _acq _time": "detector_acquisition_time",
    "number _of _channels": "number_of_channels",
    "stage_position": "stage_positions",
    "defect_positioner_i_d": "defect_positioner_id",
    "defect_positioner_aligment": "defect_positioner_alignment",
    "gcib_beam": "gcib_high_voltage",
    "gcib_wien": "gcib_wien_filter_voltage",
    "gcib_bend": "gcib_bend_voltage",
    "gcib_magnet": "gcib_magnet_current",
    "auto_e_gun_neut": "auto_flood_gun",
    "auto_ion_neut": "auto_neutral_ion_source",
}

KEYS_WITH_UNITS: List[str] = [
    "analyzer_work_function",
    "source_analyzer_angle",
    "analyzer_solid_angle",
    "scan_deflection_span_x",
    "scan_deflection_span_y",
    "scan_deflection_offset_x",
    "scan_deflection_offset_y",
    "sca_multiplier_voltage",
    "delay_before_acquire",
    "sputter_current",
    "sputter_rate",
    "sputter_energy",
    "float_voltage",
    "target_sputter_time",
    "sputter_emission",
    "grid_voltage",
    "condenser_lens_voltage",
    "objective_lens_voltage",
    "bend_voltage",
    "deflection_bias",
    "ion_gun_gas_pressure",
    "sputter_emission",
    "deflection_bias",
    "neutral_current",
    "neutral_rate",
    "neutral_energy",
    "neutral_float_voltage",
    "neutral_grid_voltage",
    "neutral_condenser_lens_voltage",
    "neutral_objective_lens_voltage",
    "neutral_bend_voltage",
    "neutral_target_timed_on_time",
    "neutral_emission",
    "neutral_deflection_bias",
    "neutral_ion_gun_gas_pressure",
    "profiling_sputter_delay",
    "survey_dwell_time",
    "xray_anode_power",
    "xray_power",
    "xray_beam_voltage",
    "xray_beam_diameter",
    "xray_condenser_lens_voltage",
    "xray_objective_coil_current",
    "xray_blanking_voltage",
    "xray_filament_current",
    "xray_rotation",
    "xray_emission_current",
    "xray_max_filament_current",
    "xray_stigmator_x",
    "xray_stigmator_y",
    "xray_offset_x",
    "xray_offset_y",
    "flood_gun_current",
    "flood_gun_energy",
    "flood_gun_extractor_voltage",
    "flood_gun_filament_current",
    "flood_gun_pulse_length",
    "flood_gun_pulse_frequency",
    "flood_gun_time_per_step",
    "flood_gun_ramp_rate",
    "flood_gun_x_steering",
    "flood_gun_y_steering",
    "sxi_binding_energy",
    "sxi_pass_energy",
    "sxi_lens2_voltage",
    "sxi_lens3_voltage",
    "sxi_lens4_voltage",
    "sxi_lens5_voltage",
    "sxi_rotator",
    "sxi_lens_bias_voltage",
    "sxi_shutter_bias_voltage",
    "detector_acquisition_time",
    "stage_current_rotation_speed",
    "defect_positioner_u",
    "defect_positioner_v",
    "defect_positioner_x",
    "defect_positioner_y",
    "defect_positioner_z",
    "defect_positioner_tilt",
    "defect_positioner_rotation",
    "gcib_sputter_rate",
    "gcib_high_voltage",
    "gcib_ionization",
    "gcib_extractor",
    "gcib_wien_filter_voltage",
    "gcib_bend_voltage",
    "gcib_emission",
    "gcib_magnet_current",
    "gcib_focus",
    "gcib_objective",
    "gcib_focus",
    "gcib_gas_pressure",
    "gcib_cluster_size",
    "gcib_energy_per_atom",
    "deconvolution_pass_energy",
]

UNIT_MISSING: Dict[str, str] = {
    "intensity": "counts_per_second",
    "grid_voltage": "V",
    "condenser_lens_voltage": "V",
    "objective_lens_voltage": "V",
    "bend_voltage": "V",
    "neutral_grid_voltage": "V",
    "neutral_condenser_lens_voltage": "V",
    "neutral_objective_lens_voltage": "V",
    "stage_current_rotation_speed": "degree/s",
    "neutral_bend_voltage": "V",
    "defect_positioner_u": "mm",
    "defect_positioner_v": "mm",
    "defect_positioner_x": "mm",
    "defect_positioner_y": "mm",
    "defect_positioner_z": "mm",
    "defect_positioner_tilt": "degree",
    "defect_positioner_rotation": "degree",
    "flood_gun_pulse_frequency": "1/s",
    "flood_gun_x_steering": "mm",
    "flood_gun_y_steering": "mm",
    "profiling_sputter_delay": "s",
    "scan_deflection_span_x": "mm",
    "scan_deflection_span_y": "mm",
    "scan_deflection_offset_x": "mm",
    "scan_deflection_offset_y": "mm",
    "survey_dwell_time": "s",
    "xray_offset_x": "um",
    "xray_offset_y": "um",
    "xray_rotation": "degree",
    "xray_stigmator_x": "mm",
    "xray_stigmator_y": "mm",
}


class MapperPhi(XPSMapper):
    """
    Class for restructuring .xy data file from
    Phi vendor into python dictionary.
    """

    config_file = "config_phi.json"

    def __init__(self):
        super().__init__()
        self.write_channels_to_data = True

    def _select_parser(self):
        """Select the proper Phi data parser."""
        return PhiParser()

    def construct_data(self):
        """Map Phi format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: Dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []

        for part in ["group_name", "spectrum_type"]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            if key.startswith("entry"):
                entry_parent = "/ENTRY[entry]"
                key = key.replace("entry/", "", 1)
            mpes_key = f"{entry_parent}/{key}"
            self._xps_dict[mpes_key] = value

            # units = get_units_for_key(key, UNITS)
            # if units is not None:
            #     self._xps_dict[f"{mpes_key}/@units"] = units

        # Create key for writing to data
        cycle_key = "cycle0"

        energy = np.array(spectrum["energy"])

        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        # Write averaged data to 'data'.
        all_scan_data = np.array(
            [np.array(value) for key, value in spectrum["data"].items()]
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            averaged_scans = np.mean(all_scan_data, axis=0)

        self._xps_dict["data"][entry][cycle_key] = xr.DataArray(
            data=averaged_scans,
            coords={"energy": energy},
        )

        # Write scan data to 'data'.
        for scan_no, intensity in spectrum["data"].items():
            xarr_key = f"{cycle_key}_{scan_no}"
            self._xps_dict["data"][entry][xarr_key] = xr.DataArray(
                data=intensity,
                coords={"energy": energy},
            )


class PhiParser:  # pylint: disable=too-few-public-methods
    """
    A parser for reading in PHI VersaProbe 4 data in the .spe or
    .pro format.
    Tested with Software version SS 3.3.3.2.
    """

    def __init__(self):
        """
        Construct the parser.

        """
        self.raw_data: str = ""
        self.spectra: List[Dict[str, Any]] = []

        self.metadata = PhiMetadata()

        self.binary_header_length = 4
        self.spectra_header_length = 24
        self.encoding: List[str, int] = ["<f", 4]

        self.binary_header: np.ndarray = None
        self.spectra_header: np.ndarray = None

        self.value_function_map = {
            "technique": convert_measurement_method,
            "technique_ex": convert_measurement_method,
            "file_type": _map_file_type,
            "file_date": _parse_datetime,
            "acquisition_file_date": _parse_datetime,
            "energy_reference": _convert_energy_referencing,
            "intensity_calibration_coefficients": _map_to_list,
            "energy_recalibration": convert_bool,
            "scan_deflection_span": _map_to_xy,
            "scan_deflection_offset": _map_to_xy,
            "tfc_parameters": _map_to_list,
            "image_size": _map_to_xy,
            "float_enabled": convert_bool,
            "sputter_raster": _map_to_xy_with_units,
            "sputter_raster_offset": _map_to_xy_with_units,
            "neutral_raster": _map_to_xy_with_units,
            "neutral_raster_offset": _map_to_xy_with_units,
            "profiling_xray_off_during_sputter": convert_bool,
            "profiling_source_blank_during_sputter": convert_bool,
            "profiling_depth_recalibration": convert_bool,
            "energy_scan_mode": convert_energy_scan_mode,
            "xray_source": _convert_xray_source_params,
            "xray_stigmator": _map_to_xy,
            "xray_offset": _map_to_xy,
            "xray_magnification_factor": _map_to_xy,
            "xray_delay_factor": _map_to_xy,
            "xray_high_power": convert_bool,
            "xray_emission_control": convert_bool,
            "xray_settings": _convert_xray_source_settings,
            "sxi_auto_contrast": convert_bool,
            "sxi_shutter_bias": convert_bool,
            "stage_positions": _convert_stage_positions,
            "gcib_raster_size": _map_to_xy_with_units,
            "gcib_raster_offset": _map_to_xy_with_units,
            "auto_flood_gun": convert_bool,
            "auto_neutral_ion_source": convert_bool,
            "presputter": convert_bool,
        }

    def parse_file(self, file: Union[str, Path], **kwargs):
        """
        Parse the .spe, .pro file into a list of dictionaries.

        Parsed data is stored in the attribute 'self.data'.
        Each dictionary in the data list is a grouping of related
        attributes. The dictionaries are later re-structured into a
        nested dictionary that more closely resembles the domain logic.

        Parameters
        ----------
        file : str
            XPS data filepath.

        Returns
        -------
        list
            Flat list of dictionaries containing one spectrum each.

        """
        self.raw_data = self._read_lines(file)
        header, data = self._separate_header_and_data()
        # l = b'\x01'.join(data)

        self.parse_header_into_metadata(header)
        regions = self.parse_spectral_regions(header)
        areas = self.parse_spatial_areas(header)

        self.add_regions_and_areas_to_spectra(regions, areas)

        self.parse_binary_header(data)
        self.parse_data_into_spectra(data)

        self.add_metadata_to_each_spectrum()

        return self.spectra

    def _read_lines(self, file: Union[str, Path]):
        """
        Read in all lines from the file as a list of strings.

        Parameters
        ----------
        file : str
            XPS data filepath.

        Returns
        -------
        None.

        """
        with open(file, mode="rb") as txt_file:  # "utf-8"
            data = txt_file.read()

        return data

    def _separate_header_and_data(self):
        """
        Separate header (with metadata) from data.

        Returns
        -------
        None.

        """
        try:  # Windows
            header, data = self.raw_data.split(b"EOFH\r\n")
            header = header.decode("ascii").split("\r")
            header = [line.strip() for line in header if line.strip()]

        except ValueError:  # Linux, MacOS
            header, data = self.raw_data.split(b"EOFH\n")
            header = header.decode("ascii").split("\n")
            header = [line.strip() for line in header if line.strip()]

        return header, data

    def parse_header_into_metadata(self, header: List[str]):
        """
        Parse header into PHiMetadata dataclass.

        Parameters
        ----------
        header : str
            Header data for one spectrum as a String.

        """

        def map_keys(key: str, channel_count: int):
            """Maps key names for better identification of fields."""
            prefix_map = {
                "neut": ("neut_", "neutral_"),
                "prof": ("prof_", "profiling_"),
                "x_ray": ("x_ray", "xray"),
                "egun_neut": ("egun_neut", "flood_gun"),
            }

            for prefix, (initial, replacement) in prefix_map.items():
                if key.startswith(prefix):
                    key = key.replace(initial, replacement)

            if key.startswith("sxi_lens"):
                key += "_voltage"

            if key.startswith("channel _info"):
                key = f"channel_{channel_count}_info"
                channel_count += 1

            key_map: Dict[str, str] = {
                "desc": "description",
                "defect_pos": "defect_positioner",
                "g_c_i_b": "gcib",
            }

            for old_key, new_key in key_map.items():
                if old_key in key:
                    key = key.replace(old_key, new_key)

            key = SETTINGS_MAP.get(key, key)

            return key, channel_count

        channel_count = 1
        datacls_fields = list(self.metadata.__dataclass_fields__.keys())
        datacls_fields = [field for field in datacls_fields if "_units" not in field]

        for line in header:
            try:
                key, value = line.split(": ")

            except ValueError:
                key = line.strip(": ")
                value = ""

            key = convert_pascal_to_snake(key)
            key, channel_count = map_keys(key, channel_count)

            if key in datacls_fields:
                field_type = type(getattr(self.metadata, key))

                if key in KEYS_WITH_UNITS:
                    value, unit = self.extract_unit(key, value)
                    setattr(self.metadata, f"{key}_units", unit)

                if not key.startswith("channel_"):
                    value = self.map_values(key, value, field_type)

                else:
                    value = _convert_channel_info(value)

                setattr(self.metadata, key, value)

        self.metadata.validate_types()

    def parse_spectral_regions(self, header: List[str]):
        """
        Parse spectral regions definitions.

        Parameters
        ----------
        header : list
            List of header strings.

        Returns
        -------
        regions_full : list
            List of regions with `full` keyword.
        regions : list
            List of regions without `full` keyword.

        """
        spectral_defs = [line for line in header if line.startswith("SpectralReg")]

        regions = []
        regions_full = []
        for n in range(self.metadata.no_spectral_regions_full):
            regions_full += [PhiSpectralRegion(region_id=n + 1)]
        for n in range(self.metadata.no_spectral_regions):
            regions += [PhiSpectralRegion(region_id=n + 1)]

        concept_map = {
            "SpectralRegDef": "region_definition",
            "SpectralRegDef2": "region_definition2",
            "SpectralRegBackground": "region_background",
            "SpectralRegHero": "region_hero",
            "SpectralRegIR": "region_ir",
        }

        for region in regions_full:
            region.full_region = True
            for file_key, metadata_key in concept_map.items():
                if spectral_defs and spectral_defs[0].startswith(file_key):
                    setattr(region, metadata_key, spectral_defs[0].split(" ", 1)[1])
                    spectral_defs.pop(0)

        for region in regions:
            region.full_region = False

            for file_key, metadata_key in concept_map.items():
                if spectral_defs and spectral_defs[0].startswith(file_key):
                    setattr(region, metadata_key, spectral_defs[0].split(" ", 1)[1])
                    spectral_defs.pop(0)

        for region in regions:
            def_split = region.region_definition.split(" ")
            region.spectrum_type = def_split[2]
            region.n_values = int(def_split[4])
            step = -float(def_split[5])
            start = float(def_split[6])
            stop = float(def_split[7])
            region.dwell_time = float(def_split[-3])
            region.dwell_time_units = "s"
            region.pass_energy = float(def_split[-2])
            region.energy_type = "binding"
            region.energy_units = "eV"

            region.energy = np.flip(safe_arange_with_edges(stop, start, step))

            region.validate_types()

        return regions

    def parse_spatial_areas(self, header: List[str]):
        """
        Parse spatial areas definitions.

        Parameters
        ----------
        header : list
            List of header strings.

        Returns
        -------
        areas: list
            List of spatial areas and their definitions.

        """
        spatial_defs = [line for line in header if line.startswith("Spatial")]

        areas = []
        for n in range(self.metadata.no_spatial_areas):
            areas += [PhiSpatialArea(area_id=n)]

        concept_map = {
            "SpatialAreaDef": "area_definition",
            "SpatialAreaDesc": "area_description",
            "SpatialHRPhotoCor": "area_hr_photo_correction",
        }

        for area in areas:
            for file_key, metadata_key in concept_map.items():
                if spatial_defs and spatial_defs[0].startswith(file_key):
                    setattr(area, metadata_key, spatial_defs[0].split(" ", 1)[1])
                    spatial_defs.pop(0)

        return areas

    def add_regions_and_areas_to_spectra(self, regions: list, areas: list):
        """
        Define each spectra by its region and area defintions.

        Parameters
        ----------
        regions : list
            List of PhiSpectralRegion objects.
        areas : list
            List of PhiSpatialArea objects

        Returns
        -------
        None.

        """
        for region in regions:
            for area in areas:
                concatenated = {**region.dict(), **area.dict()}

                region_and_areas: Dict[str, Any] = {}

                for key, value in concatenated.items():
                    replacement_map = {
                        "_units": "/@units",
                        "energy_type": "energy/@type",
                    }

                    new_key = key

                    for suffix, replacement in replacement_map.items():
                        if suffix in key:
                            new_key = key.replace(suffix, replacement)

                    region_and_areas[new_key] = value

                self.spectra += [region_and_areas]

    def parse_binary_header(self, binary_data):
        """
        Read the binary headers
        Assuming the headers are 4 bytes unsigned integers
        Each spectrum gets 24 unsigned 4 byte integers

        Parameters
        ----------
        binary_data : bytes
            Binary XPS data, format is 64 bit float.

        """
        binary_header = struct.unpack("I", binary_data[: self.binary_header_length])[0]

        for i, spectrum in enumerate(self.spectra):
            start = (
                self.binary_header_length * self.spectra_header_length * i
                + self.binary_header_length
            )
            stop = start + self.binary_header_length * self.spectra_header_length
            spectrum_header = struct.unpack(
                "I" * self.spectra_header_length, binary_data[start:stop]
            )

            spectrum.update(
                {
                    "binary_header": binary_header,
                    "spectrum_header": np.array(spectrum_header),
                    "n_values": spectrum_header[8],
                    "n_scans": spectrum_header[9],
                    "binary_start": spectrum_header[23],
                    "binary_len": spectrum_header[22],
                }
            )

        self._check_encoding()

    def parse_data_into_spectra(self, binary_data):
        """
        Parse the data of all spectra.

        Parameters
        ----------
        binary_data : bytes
            Binary XPS data, format is 64 bit float.

        """
        for i, spectrum in enumerate(self.spectra):
            float_buffer = self.encoding[1]

            # spectrum["n_scans"] = 2
            spectrum["data"] = {}
            start = spectrum["binary_start"]
            stop = start + spectrum["binary_len"]

            binary_spectrum_data = binary_data[start:stop]
            n_values_binary = spectrum["n_values"] * float_buffer

            for scan_no in range(spectrum["n_scans"]):
                spec_start = scan_no * n_values_binary
                spec_stop = spec_start + n_values_binary
                scan_data = binary_spectrum_data[spec_start:spec_stop]

                spectrum["data"][f"scan_{scan_no}"] = self._parse_binary_scan(scan_data)

    def _parse_binary_scan(self, binary_scan_data):
        """
        For each spectrum scan, parse the XPS data by unpacking the
        64 bit floats.

        Parameters
        ----------
        binary_spectrum_data : bytes
            Binary data containing the intensity data for one
            spectrum.

        Returns
        -------
        parsed_data : np.ndarray
            One-dimensional array with spectrum intensities.

        """
        float_buffer = self.encoding[1]
        encoding = self.encoding[0]
        parsed_data = []

        n_values = int(len(binary_scan_data) / float_buffer)

        for i in range(n_values):
            start = i * float_buffer
            stop = start + float_buffer
            parsed_data.append(
                struct.unpack_from(encoding, binary_scan_data[start:stop])[0]
            )

        return parsed_data

    def _check_encoding(self):
        """Check if the binary data is single or double encoded."""
        datasize = sum(
            [s["spectrum_header"][8] * s["spectrum_header"][9] for s in self.spectra]
        )
        binary_size = sum([s["spectrum_header"][-2] for s in self.spectra])

        encodings_map = {
            "double": ["<d", 8],
            "float": ["<f", 4],
        }
        if binary_size / datasize == 4:
            self.encoding = encodings_map["float"]
        elif binary_size / datasize == 8:
            self.encoding = encodings_map["double"]
        else:
            logger.error("This binary encoding is not supported.")

    def extract_unit(self, key: str, value):
        """
        Extract units for the metadata containing unit information.

        Example:
            analyzer_work_function: 4.506 eV
            -> analyzer_work_function: 4.506,
               analyzer_work_function_units: eV,

        Parameters
        ----------
        key : str
            Key of the associated value.
        value : str
            Combined unit and value information.

        Returns
        -------
        value :
            value with units.
        unit : str
            Associated unit.

        """
        try:
            value, unit = split_value_and_unit(value)
        except ValueError:
            unit = ""

        unit = convert_units(unit)

        if key in UNIT_MISSING:
            unit = UNIT_MISSING[key]

        return value, unit

    def map_values(self, key: str, value, field_type):
        """
        Map values to corresponding structure and field type.

        Parameters
        ----------
        key : str
            PhiMetadata dataclass key.
        value : str
            Original value.
        field_type : class
            Class to be mapped onto.

        Returns
        -------
        value
            Value of correct type and internal structure.

        """
        if key in self.value_function_map:
            map_fn = self.value_function_map[key]
            value = map_fn(value)
        return field_type(value)

    def add_metadata_to_each_spectrum(self):
        """
        Add flattened dict with PhiMetadata fields to each spectrum.

        """
        flattened_metadata = self.flatten_metadata()

        for spectrum in self.spectra:
            spectrum.update(flattened_metadata)
            spectrum["intensity/@units"] = UNIT_MISSING["intensity"]

    def flatten_metadata(self):
        """
        Flatten metadata dict so that key-value pairs of nested
        dictionaries are at the top level.

        Parameters
        ----------
        metadata_dict : dict
            Metadata dict with PhiMetadata fields as keys.

        Returns
        -------
        flattened_dict : dict
            Flatted metadata_dict without any nested substructure.

        """

        def shorten_supkey(supkey: str):
            """Shorted the key for some nested dicts."""
            shortened_key_map = {
                "xray_source": "xray",
                "xray_settings": "xray",
                "stage_positions": "stage",
            }
            if supkey in shortened_key_map:
                return shortened_key_map[supkey]
            return supkey

        flattened_dict = {}

        for key, value in self.metadata.dict().items():
            if isinstance(value, dict):
                for subkey, subvalue in value.items():
                    supkey = shorten_supkey(key)
                    flattened_dict[f"{supkey}_{subkey}"] = subvalue
            else:
                flattened_dict[key] = value

        for key in flattened_dict.copy().keys():
            if "_units" in key:
                new_key = key.replace("_units", "/@units")
                flattened_dict[new_key] = flattened_dict.pop(key)

        for key in flattened_dict.copy().keys():
            if key in KEYS_WITH_UNITS and f"{key}/@units" not in flattened_dict:
                flattened_dict[f"{key}/@units"] = UNIT_MISSING[key]

        return flattened_dict


def _map_file_type(value: str):
    """Map file_type to easily understandable values."""
    value = value.strip()
    file_type_map = {
        "SPECTRUM": "single_spectrum",
        "DEPTHPRO": "depth_profile",
    }

    if value in file_type_map:
        return file_type_map[value]
    return value


def _parse_datetime(value: str):
    """
    Parse datetime into a datetime.datetime object.

    Parameters
    ----------
    value : str
        String representation of the date in the format
        "%m/%d/%y.

    Returns
    -------
    date_object : str
        Datetime in ISO8601 format.

    """
    year, month, day = value.strip().split(" ")
    date_object = datetime.datetime(
        year=int(year), month=int(month), day=int(day), tzinfo=pytz.timezone("UTC")
    )

    return date_object.isoformat()


def _map_to_list(value: str):
    """Map all items in value to a list."""
    try:
        sep = ", " if ", " in value else " "
        values = [float(val) for val in value.split(sep)]
    except ValueError:
        sep = ","
        values = [float(val) for val in value.split(sep)]

    return values


def _map_to_xy(value: str):
    """Map items in value to a dictionary with keys x and y."""
    x, y = value.split(" ")

    return {"x": float(x), "y": float(y)}


def _map_to_xy_with_units(value: str):
    """Map items in value to a dictionary with keys x ."""
    x, y, unit = value.split(" ")
    return {"x": float(x), "x_units": unit, "y": float(y), "y_units": unit}


def _convert_energy_referencing(value: str):
    """Map all items in energy_referencing to a dictionary."""
    peak, energy = value.split(" ")
    return {"peak": peak, "energy": float(energy), "energy_units": "eV"}


def _convert_channel_info(value: str):
    """Split channel information into list."""
    channel_number, setting_a, setting_b = value.split(" ")
    return [int(setting_a), float(setting_b)]


def _convert_xray_source_params(value: str):
    """Map all items in xray_source_params to a dictionary."""
    label, energy, mono = value.split(" ")

    return {
        "anode_material": label,
        "energy": float(energy),
        "energy_units": "eV",
        "monochromatized": mono,
    }


def _convert_xray_source_settings(value: str):
    """Map all items in xray_source_settings to a dictionary."""
    (xray_settings) = re.split(r"(\d+)", value)

    for i, setting in enumerate(xray_settings):
        xray_settings[i] = convert_units(setting)

    return {
        "spot_size": float(xray_settings[1]),
        "spot_size_units": xray_settings[2],
        "power": float(xray_settings[3]),
        "power_units": xray_settings[4],
        "high_voltage": float(xray_settings[5]),
        "high_voltage_units": xray_settings[6],
    }


def _convert_stage_positions(value: str):
    """Map all items in stage_positions to a dictionary."""
    x, y, z, azimuth, polar = value.split(" ")

    return {
        "x": float(x),
        "x_units": "mm",
        "y": float(y),
        "y_units": "mm",
        "z": float(z),
        "z_units": "mm",
        "azimuth": float(azimuth),
        "azimuth_units": "degree",
        "polar": float(polar),
        "polar_units": "degree",
    }
