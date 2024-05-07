"""
Class for reading XPS files from TXT export of Scienta.
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
import copy
import datetime
from pathlib import Path
from typing import Any, Dict, List, Union, Optional
import xarray as xr
import numpy as np

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    construct_data_key,
    construct_detector_data_key,
    convert_pascal_to_snake,
    _re_map_single_key,
    _re_map_single_value,
    _check_valid_value,
)
from pynxtools_xps.value_mappers import get_units_for_key

from pynxtools_xps.scienta.scienta_txt_data_model import ScientaHeader, ScientaRegion
from pynxtools_xps.scienta.scienta_mappings import UNITS, KEY_MAP, VALUE_MAP


def _extract_energy_units(energy_units: str):
    """
    Extract energy units from the strings for energy_units.
    Binding Energy [eV] -> eV

    """
    return re.search(r"\[(.*?)\]", energy_units).group(1)


def _separate_dimension_scale(scale: str):
    """
    Seperate the str of the dimension scale into a numpy array

    Parameters
    ----------
    scale : str
        Str of the form "600 599.5 599 598.5 5".

    Returns
    -------
    np.ndarray
        Dimension scale as a numpy array.

    """
    return np.array([float(s) for s in scale.split(" ")])


def _construct_date_time(date_string: str, time_string: str) -> Optional[str]:
    """
    Convert the native time format to the datetime string
    in the ISO 8601 format: '%Y-%b-%dT%H:%M:%S.%fZ'.

    """

    def _parse_date(date_string: str) -> datetime.datetime:
        possible_date_formats = ["%Y-%m-%d", "%m/%d/%Y"]
        for date_fmt in possible_date_formats:
            try:
                return datetime.datetime.strptime(date_string, date_fmt)
            except ValueError:
                pass
        raise ValueError("Date format not recognized")

    def _parse_time(time_string: str) -> datetime.time:
        possible_time_formats = ["%H:%M:%S", "%I:%M:%S %p"]
        for time_fmt in possible_time_formats:
            try:
                return datetime.datetime.strptime(time_string, time_fmt).time()
            except ValueError:
                pass
        raise ValueError("Time format not recognized")

    try:
        date = _parse_date(date_string)
        time = _parse_time(time_string)
        datetime_obj = datetime.datetime.combine(date, time)

        return datetime_obj.astimezone().isoformat()

    except ValueError as e:
        raise ValueError(
            "Date and time could not be converted to ISO 8601 format."
        ) from e


class TxtMapperScienta(XPSMapper):
    """
    Class for restructuring .txt data file from
    Scienta TXT export into python dictionary.
    """

    config_file = "config_scienta_txt.json"

    def _select_parser(self):
        """
        Select Scienta TXT parser.
        Currently, there is only one parser.

        Returns
        -------
        ScientaTxtParser
            Parser for reading .txt file exported by Scienta.

        """
        return ScientaTxtParser()

    def construct_data(self):
        """Map TXT format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        template_key_map = {
            "file_info": ["data_file", "sequence_file"],
            "user": [
                "user_name",
            ],
            "instrument": [
                "instrument_name",
                "vendor",
            ],
            "source_xray": [],
            "beam_xray": [
                "excitation_energy",
            ],
            "analyser": [],
            "collectioncolumn": [
                "lens_mode",
            ],
            "energydispersion": [
                "acquisition_mode",
                "pass_energy",
            ],
            "detector": [
                "detector_first_x_channel",
                "detector_first_y_channel",
                "detector_last_x_channel",
                "detector_last_y_channel",
                "detector_mode",
                "dwell_time",
                "time_per_spectrum_channel",
            ],
            "manipulator": [
                "manipulator_r1",
                "manipulator_r2",
            ],
            "calibration": [],
            "sample": ["sample_name"],
            "region": [
                "center_energy",
                "energy_axis",
                "energy_scale",
                "energy_scale_2",
                "energy_size",
                "no_of_scans",
                "region_id",
                "spectrum_comment",
                "start_energy",
                "step_size",
                "stop_energy",
                "time_stamp",
            ],
            # 'unused': [
            #     'energy_unit',
            #     'number_of_slices',
            #     'software_version',
            #     'spectrum_comment',
            #     'start_date',
            #     'start_time',
            #     'time_per_spectrum_channel'
            # ]
        }

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum, template_key_map)

    def _update_xps_dict_with_spectrum(
        self, spectrum: Dict[str, Any], template_key_map: Dict[str, List[str]]
    ):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        group_parent = f'{self._root_path}/Group_{spectrum["spectrum_type"]}'
        region_parent = f'{group_parent}/Region_{spectrum["region_name"]}'
        file_parent = f"{region_parent}/file_info"
        instrument_parent = f"{region_parent}/instrument"
        analyser_parent = f"{instrument_parent}/analyser"

        path_map = {
            "file_info": f"{file_parent}",
            "user": f"{region_parent}/user",
            "instrument": f"{instrument_parent}",
            "source_xray": f"{instrument_parent}/source_xray",
            "beam_xray": f"{instrument_parent}/beam_xray",
            "analyser": f"{analyser_parent}",
            "collectioncolumn": f"{analyser_parent}/collectioncolumn",
            "energydispersion": f"{analyser_parent}/energydispersion",
            "detector": f"{analyser_parent}/detector",
            "manipulator": f"{instrument_parent}/manipulator",
            "calibration": f"{instrument_parent}/calibration",
            "sample": f"{region_parent}/sample",
            "data": f"{region_parent}/data",
            "region": f"{region_parent}/region",
        }

        for grouping, spectrum_keys in template_key_map.items():
            root = path_map[str(grouping)]

            for spectrum_key in spectrum_keys:
                mpes_key = spectrum_key.rsplit(" ", 1)[0]
                try:
                    self._xps_dict[f"{root}/{mpes_key}"] = spectrum[spectrum_key]
                except KeyError:
                    pass

                unit_key = f"{grouping}/{spectrum_key}"
                units = get_units_for_key(unit_key, UNITS)
                if units:
                    self._xps_dict[f"{root}/{mpes_key}/@units"] = units

        # Create keys for writing to data and detector
        entry = construct_entry_name(region_parent)
        scan_key = construct_data_key(spectrum)
        detector_data_key_child = construct_detector_data_key(spectrum)
        detector_data_key = f'{path_map["detector"]}/{detector_data_key_child}/counts'

        # Write raw data to detector.
        self._xps_dict[detector_data_key] = spectrum["data"]["intensity"]

        # If multiple spectra exist to entry, only create a new
        # xr.Dataset if the entry occurs for the first time.
        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        energy = np.array(spectrum["data"]["energy"])
        intensity = spectrum["data"]["intensity"]

        # Write to data in order: scan, cycle, channel

        # Write averaged cycle data to 'data'.
        all_scan_data = [
            value
            for key, value in self._xps_dict["data"][entry].items()
            if scan_key.split("_")[0] in key
        ]
        averaged_scans = np.mean(all_scan_data, axis=0)
        if averaged_scans.size == 1:
            # on first scan in cycle
            averaged_scans = intensity

        self._xps_dict["data"][entry][scan_key.split("_")[0]] = xr.DataArray(
            data=averaged_scans,
            coords={"energy": energy},
        )

        # Write scan data to 'data'.
        self._xps_dict["data"][entry][scan_key] = xr.DataArray(
            data=intensity, coords={"energy": energy}
        )

        # Write channel data to 'data'.
        channel_key = f"{scan_key}_chan0"
        self._xps_dict["data"][entry][channel_key] = xr.DataArray(
            data=intensity, coords={"energy": energy}
        )


class ScientaTxtParser:
    """Parser for Scienta TXT exports."""

    # pylint: disable=too-few-public-methods

    def __init__(self):
        self.lines: List[str] = []
        self.header = ScientaHeader()
        self.spectra: Dict[str, Any] = []

    def parse_file(self, file: Union[str, Path]):
        """
        Parse the file's data and metadata into a flat
        list of dictionaries.


        Parameters
        ----------
        file : str
            Filepath of the TXT file to be read.

        Returns
        -------
        self.spectra
            Flat list of dictionaries containing one spectrum each.

        """
        self._read_lines(file)
        self._parse_header()

        for region_id in range(1, self.header.no_of_regions + 1):
            self._parse_region(region_id)

        return self.spectra

    def _read_lines(self, file: Union[str, Path]):
        """
        Read all lines from the input txt files.


        Parameters
        ----------
        file : str
            Filepath of the TXT file to be read.

        Returns
        -------
        None.

        """
        with open(file, encoding="utf-8") as txt_file:
            for line in txt_file:
                self.lines += [line]

    def _parse_header(self):
        """
        Parse header with information about the software version
        and the number of spectra in the file.

        Returns
        -------
        None.

        """
        n_headerlines = 4
        headerlines = self.lines[:n_headerlines]
        self.lines = self.lines[n_headerlines:]

        for line in headerlines:
            key, value = self._get_key_value_pair(line)
            if key:
                value = _re_map_single_value(key, value, VALUE_MAP)
                setattr(self.header, key, value)
        self.header.validate_types()

    def _parse_region(self, region_id):
        """
        Parse data from one region (i.e., one measured spectrum)
        into a dictioanry and append to all spectra.

        Parameters
        ----------
        region_id : int
            Number of the region in the file.

        Returns
        -------
        None.

        """
        region = ScientaRegion()
        region.region_id = region_id

        bool_variables = {
            "in_region": False,
            "in_info": False,
            "in_run_mode_info": False,
            "in_ui_info": False,
            "in_manipulator": False,
            "in_data": False,
        }

        energies: List[float] = []
        intensities: List[float] = []

        line_start_patterns = {
            "in_region": f"[Region {region_id}",
            "in_info": f"[Info {region_id}",
            "in_run_mode_info": f"[Run Mode Information {region_id}",
            "in_ui_info": f"[User Interface Information {region_id}",
            "in_manipulator": f"[Manipulator {region_id}",
            "in_data": f"[Data {region_id}",
        }

        for line in self.lines:
            for bool_key, line_start in line_start_patterns.items():
                if line.startswith(line_start):
                    bool_variables[bool_key] = True
                if line.startswith("\n"):
                    bool_variables[bool_key] = False

            if bool_variables["in_region"]:
                # Read region meta data.
                key, value = self._get_key_value_pair(line)
                if "dimension" in key:
                    key_part = f"dimension_{key.rsplit('_')[-1]}"
                    key = _re_map_single_key(key_part, KEY_MAP)

                    value = _re_map_single_value(key, value, VALUE_MAP)
                if _check_valid_value(value):
                    setattr(region, key, value)

            if bool_variables["in_info"] or bool_variables["in_run_mode_info"]:
                # Read instrument meta data for this region.
                key, value = self._get_key_value_pair(line)
                if _check_valid_value(value):
                    setattr(region, key, value)

            if bool_variables["in_ui_info"]:
                key, value = self._get_key_value_pair(line)
                if _check_valid_value(value):
                    if bool_variables["in_manipulator"]:
                        key = f"manipulator_{key}"
                        value = _re_map_single_value(key, value, VALUE_MAP)
                    setattr(region, key, value)

            if bool_variables["in_data"]:
                # Read XY data for this region.
                try:
                    [energy, intensity] = [float(s) for s in line.split(" ") if s != ""]
                    energies.append(energy)
                    intensities.append(intensity)
                except ValueError:
                    # First line
                    pass

        region.data = {"energy": np.array(energies), "intensity": np.array(intensities)}

        # Convert date and time to ISO8601 date time.
        region.time_stamp = _construct_date_time(region.start_date, region.start_time)

        region.validate_types()

        region_dict = {**self.header.dict(), **region.dict()}

        for key in region_dict.copy().keys():
            if "_units" in key:
                new_key = key.replace("_units", "/@units")
                region_dict[new_key] = region_dict.pop(key)

        self.spectra.append(region_dict)

    def _get_key_value_pair(self, line: str):
        """
        Split the line at the '=' sign and return a
        key-value pair. The values are mapped according
        to the desired format.

        Parameters
        ----------
        line : str
            One line from the input file.

        Returns
        -------
        key : str
            Anything before the '=' sign, mapped to the desired
            key format.
        value : obj
            Anything after the '=' sign, mapped to the desired
            value format and type.

        """
        try:
            [key, value] = line.split("=")
            key = convert_pascal_to_snake(key)
            key = _re_map_single_key(key, KEY_MAP)
            value = _re_map_single_value(key, value, VALUE_MAP)

        except ValueError:
            key, value = "", ""

        return key, value
