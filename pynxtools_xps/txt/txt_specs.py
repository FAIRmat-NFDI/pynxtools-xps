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
import pytz
import xarray as xr
import numpy as np

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    construct_data_key,
    construct_detector_data_key,
)


class TxtMapperSpecs(XPSMapper):
    """
    Class for restructuring .txt data file from
    Specs TXT export into python dictionary.
    """

    config_file = "config_txt_specs.json"

    def _select_parser(self):
        """
        Select Specs TXT parser.
        Currently, there is only one parser.

        Returns
        -------
        ScientaTxtHelper
            Parser for reading .txt file exported by Specs.

        """
        return SpecsTxtHelper()

    def construct_data(self):
        """Map TXT format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        key_map = {
            "file_info": [],
            "user": [],
            "instrument": [],
            "source": [],
            "beam": [],
            "analyser": [],
            "collectioncolumn": [
                "lens_mode",
            ],
            "energydispersion": [
                "acquisition_mode",
                "pass_energy",
            ],
            "detector": [],
            "manipulator": [],
            "calibration": [],
            "sample": [],
            "data": [],
            "region": [],
            # 'unused': []
        }

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum, key_map)

    def _update_xps_dict_with_spectrum(self, spectrum, key_map):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        group_parent = f'{self._root_path}/RegionGroup_{spectrum["spectrum_type"]}'
        region_parent = f'{group_parent}/regions/RegionData_{spectrum["region_name"]}'
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
        self._xps_dict[detector_data_key] = spectrum["data"]["y"]

        # If multiple spectra exist to entry, only create a new
        # xr.Dataset if the entry occurs for the first time.
        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        energy = np.array(spectrum["data"]["x"])
        intensity = spectrum["data"]["y"]

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


class SpecsTxtHelper:
    """Parser for Specs TXT exports."""

    # pylint: disable=too-few-public-methods

    def __init__(self):
        self.lines: list = []
        self.spectra: list = []
        self.no_of_regions = 0

        keys_map = {}

        lens_mode_map = {"Transmission": "fixed analyzer transmission"}

        self.key_maps = [
            keys_map,
            lens_mode_map,
        ]

        self.value_map = {"key": "func"}

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
        self._parse_data()

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
        header = self.lines[:n_headerlines]
        self.lines = self.lines[n_headerlines:]

        for line in header:
            key, value = self._get_key_value_pair(line)
            if key:
                value = self._re_map_values(key, value)
                setattr(self, key, value)

    def _parse_data(self, region_id):
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
        # to be implemented
        pass

    def _re_map_keys(self, key):
        """
        Map the keys returned from the file to the preferred keys for
        the parser output.

        """
        maps = {}
        for key_map in self.key_maps:
            maps.update(key_map)

        keys = list(maps.keys())

        if key in keys:
            key = maps[key]

        return key

    def _change_energy_type(self, energy):
        """
        Change the strings for energy type to the preferred format.

        """
        if "Binding" in energy:
            return "binding energy"
        if "Kinetic" in energy:
            return "kinetic energy"
        return None

    def _construct_date_time(self, date, time):
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

    def _re_map_values(self, input_key, value):
        """
        Map the values returned from the file to the preferred format for
        the parser output.

        """
        try:
            value = value.rstrip("\n")
        except AttributeError:
            pass

        keys = list(self.value_map.keys())

        for k in keys:
            if k in input_key:
                if not isinstance(self.value_map[k], list):
                    map_methods = [self.value_map[k]]
                else:
                    map_methods = self.value_map[k]
                for method in map_methods:
                    value = method(value)
        return value
