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
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Union
import pytz
import xarray as xr
import numpy as np

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    construct_data_key,
    construct_detector_data_key,
    convert_pascal_to_snake,
)

from pynxtools_xps.scienta.scienta_txt_data_model import ScientaHeader, ScientaRegion

KEY_MAP: Dict[str, str] = {
    "number_of_regions": "no_of_regions",
    "version": "software_version",
    "dimension_name": "energy_type",
    "dimension_size": "energy_size",
    "dimension_scale": "energy_axis",
    "number_of_sweeps": "no_of_scans",
    "low_energy": "start_energy",
    "high_energy": "stop_energy",
    "energy_step": "step_size",
    "step_time": "dwell_time",
    "detector_first_x-_channel": "detector_first_x_channel",
    "detector_last_x-_channel": "detector_last_x_channel",
    "detector_first_y-_channel": "detector_first_y_channel",
    "detector_last_y-_channel": "detector_last_y_channel",
    "file": "data_file",
    "sequence": "sequence_file",
    "spectrum_name": "spectrum_type",
    "instrument": "instrument_name",
    "location": "vendor",
    "user": "user_name",
    "sample": "sample_name",
    "comments": "spectrum_comment",
    "date": "start_date",
    "time": "start_time",
    "name": "energy_type",
    "size": "energy_size",
    "scale": "energy_axis",
    "Transmission": "fixed analyzer transmission",
}


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
        ScientaTxtHelper
            Parser for reading .txt file exported by Scienta.

        """
        return ScientaTxtHelper()

    def construct_data(self):
        """Map TXT format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        key_map = {
            "file_info": ["data_file", "sequence_file"],
            "user": [
                "user_name",
            ],
            "instrument": [
                "instrument_name",
                "vendor",
            ],
            "source": [],
            "beam": [
                "excitation_energy",
                "excitation_energy/@units",
            ],
            "analyser": [],
            "collectioncolumn": [
                "lens_mode",
            ],
            "energydispersion": [
                "acquisition_mode",
                "pass_energy",
                "pass_energy/@units",
            ],
            "detector": [
                "detector_first_x_channel",
                "detector_first_y_channel",
                "detector_last_x_channel",
                "detector_last_y_channel",
                "detector_mode",
                "dwell_time",
                "dwell_time/@units",
                "time_per_spectrum_channel",
                "time_per_spectrum_channel/@units",
            ],
            "manipulator": [],
            "calibration": [],
            "sample": ["sample_name"],
            "data": [
                "energy_unit",
                "energy/@units",
                "energy_axis",
                "energy_type",
                "center_energy",
                "center_energy/@units",
                "start_energy",
                "start_energy/@units",
                "stop_energy",
                "stop_energy/@units",
                "step_size",
                "step_size/@units",
            ],
            "region": [
                "energy_scale",
                "energy_size",
                "no_of_scans",
                "region_id",
                "spectrum_comment",
                "start_energy",
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
            self._update_xps_dict_with_spectrum(spectrum, key_map)

    def _update_xps_dict_with_spectrum(self, spectrum, key_map):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        group_parent = f'{self._root_path}/Region_{spectrum["spectrum_type"]}'
        region_parent = f'{group_parent}/regions/Region_{spectrum["region_name"]}'
        file_parent = f"{region_parent}/file_info"
        instrument_parent = f"{region_parent}/instrument"
        analyser_parent = f"{instrument_parent}/analyser"

        path_map = {
            "file_info": f"{file_parent}",
            "user": f"{region_parent}/user",
            "instrument": f"{instrument_parent}",
            "source": f"{instrument_parent}/source",
            "beam": f"{instrument_parent}/beam",
            "analyser": f"{analyser_parent}",
            "collectioncolumn": f"{analyser_parent}/collectioncolumn",
            "energydispersion": f"{analyser_parent}/energydispersion",
            "detector": f"{analyser_parent}/detector",
            "manipulator": f"{instrument_parent}/manipulator",
            "calibration": f"{instrument_parent}/calibration",
            "sample": f"{region_parent}/sample",
            "data": f"{region_parent}/data",
            "region": f"{region_parent}",
        }

        for grouping, spectrum_keys in key_map.items():
            root = path_map[str(grouping)]
            for spectrum_key in spectrum_keys:
                try:
                    units = re.search(r"\[([A-Za-z0-9_]+)\]", spectrum_key).group(1)
                    mpes_key = spectrum_key.rsplit(" ", 1)[0]
                    self._xps_dict[f"{root}/{mpes_key}/@units"] = units
                    self._xps_dict[f"{root}/{mpes_key}"] = spectrum[spectrum_key]
                except AttributeError:
                    mpes_key = spectrum_key
                    self._xps_dict[f"{root}/{mpes_key}"] = spectrum[spectrum_key]

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


class ScientaTxtHelper:
    """Parser for Scienta TXT exports."""

    # pylint: disable=too-few-public-methods

    def __init__(self):
        self.lines: List[str] = []
        self.header = ScientaHeader()
        self.spectra: Dict[str, Any] = []

    def parse_file(self, file):
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

    def _read_lines(self, file):
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
                value = self._re_map_values(key, value)
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

        begin_region = False
        begin_info = False
        begin_data = False

        energies: List[float] = []
        intensities: List[float] = []

        for line in self.lines:
            if line.startswith(f"[Region {region_id}"):
                begin_region = True
            if line.startswith(f"[Info {region_id}"):
                begin_info = True
            if line.startswith(f"[Data {region_id}"):
                begin_data = True
                continue

            if line.startswith("\n"):
                begin_region = False
                begin_info = False
                begin_data = False

            if begin_region:
                # Read region meta data.
                key, value = self._get_key_value_pair(line)
                if "dimension" in key:
                    key_part = f"dimension_{key.rsplit('_')[-1]}"
                    key = self._re_map_keys(key_part)

                    value = self._re_map_values(key, value)
                if self._check_valid_value(value):
                    setattr(region, key, value)

            if begin_info:
                # Read instrument meta data for this region.
                key, value = self._get_key_value_pair(line)
                if self._check_valid_value(value):
                    setattr(region, key, value)

            if begin_data:
                # Read XY data for this region.
                [energy, intensity] = [float(s) for s in line.split(" ") if s != ""]

                energies.append(energy)
                intensities.append(intensity)

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

    def _check_valid_value(self, value):
        """
        Check if a string or an array is empty.

        Parameters
        ----------
        value : obj
            For testing, this can be a str or a np.ndarray.

        Returns
        -------
        bool
            True if the string or np.ndarray is not empty.

        """
        for datatype in [str, int, float]:
            if isinstance(value, datatype) and value:
                return True
        if isinstance(value, bool):
            return True
        if isinstance(value, np.ndarray) and value.size != 0:
            return True
        return False

    def _get_key_value_pair(self, line):
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
            key = self._re_map_keys(key)
            value = self._re_map_values(key, value)

        except ValueError:
            key, value = "", ""

        return key, value

    def _re_map_keys(self, key):
        """
        Map the keys returned from the file to the preferred keys for
        the parser output.

        """
        if key in KEY_MAP:
            return KEY_MAP[key]
        return key

    def _re_map_values(self, input_key: str, value: str):
        """
        Map the values returned from the file to the preferred format for
        the parser output.

        """
        try:
            value = value.rstrip("\n")
        except AttributeError:
            pass

        value_map = {
            "no_of_regions": int,
            "energy_size": int,
            # "energy_axis",
            "pass_energy": float,
            "no_of_scans": int,
            "excitation_energy": float,
            "center_energy": float,
            "start_energy": float,
            "stop_energy": float,
            "step_size": float,
            "dwell_time": float,
            "detector_first_x_channel": int,
            "detector_last_x_channel": int,
            "detector_first_y_channel": int,
            "detector_last_y_channel": int,
            "time_per_spectrum_channel": float,
            "energy_unit": _change_energy_type,
            "energy_axis": _separate_dimension_scale,
            "energy_scale": _change_energy_type,
            "acquisition_mode": _change_scan_mode,
            "time_per_spectrum_channel": float,
        }

        keys = list(value_map.keys())

        for k in keys:
            if k in input_key:
                map_method = value_map[k]
                value = map_method(value)  # type: ignore[operator]
        return value


def str_to_int(str_value: str):
    """
    Make str_value an integer.
    Maps e.g. from '0001' string to 1.
    """
    return int(str_value)


def str_to_float(str_value: str):
    """
    Make str_value an integer.
    Maps e.g. from '0001' string to 1.
    """
    return int(str_value)


def _change_energy_type(energy_type: str):
    """
    Change the strings for energy type to the preferred format.

    """
    if "Binding" in energy_type:
        return "binding"
    if "Kinetic" in energy_type:
        return "kinetic"
    return None


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


def _change_scan_mode(acquisition_mode: str):
    """
    Change the strings for acquisition mode type to
    the allowed values in NXmpes.

    """
    mode_map = {"Fixed": "fixed", "Swept": "sweep"}
    if acquisition_mode in mode_map:
        return mode_map[acquisition_mode]
    return None


def _construct_date_time(date: str, time: str):
    """
    Convert the native time format to the datetime string
    in the ISO 8601 format: '%Y-%b-%dT%H:%M:%S.%fZ'.

    """
    date_time = datetime.combine(
        datetime.strptime(date, "%Y-%m-%d"),
        datetime.strptime(time, "%H:%M:%S").time(),
    )

    localtz = pytz.timezone("Europe/Berlin")
    date_time = localtz.localize(date_time)

    return date_time.isoformat()
