"""
Parser for reading XPS (X-ray Photoelectron Spectroscopy) metadata
from CSV exported from Specs Lab Prodigy, to be passed to
mpes nxdl (NeXus Definition Language) template.
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
# pylint: disable=too-many-instance-attributes

import re
import pandas as pd
from copy import copy
import sqlite3
from datetime import datetime
import numpy as np


class CsvMapperSpecs:
    """
    Class for restructuring .sle data file from
    specs vendor into python dictionary.
    """

    def __init__(self):
        self.parsers = [
            CsvProdigyParser,
        ]

        self.file = None

        self.raw_data: list = []
        self._xps_dict: dict = {}

        self._root_path = "/ENTRY[entry]"

    def _select_parser(self):
        """
        Select the correct parser for the file extension and format.

        Should be implemented by the inheriting mapper.

        Returns
        -------
        Parser

        """
        return self.parsers[0]

    def construct_data(self):
        """Map CSV format to NXmpes-ready dict."""
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        key_map = {
            "user": [],
            "instrument": [],
            "source": [],
            "beam": [],
            "analyser": [],
            "collectioncolumn": [],
            "energydispersion": [],
            "detector": [],
            "manipulator": [],
            "calibration": [],
            "data": [],
            "region": [],
        }

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum, key_map)

    def _update_xps_dict_with_spectrum(self, spectrum, key_map):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals
        group_parent = f'{self._root_path}/RegionGroup_{spectrum["group_name"]}'
        region_parent = f'{group_parent}/regions/RegionData_{spectrum["spectrum_type"]}'
        instrument_parent = f"{region_parent}/instrument"
        analyser_parent = f"{instrument_parent}/analyser"

        path_map = {
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

        self._xps_dict[f'{path_map["analyser"]}/name'] = spectrum["devices"][0]
        self._xps_dict[f'{path_map["source"]}/name'] = spectrum["devices"][1]


class CsvProdigyParser:
    """
    Generic parser without reading capabilities,
    to be used as template for implementing parsers for different versions.
    """

    supported_versions = ["0.6"]
    supported_prodigy_versions = ["1.2", "1.8", "1.9", "1.10", "1.11", "1.12", "1.13"]

    def __init__(self):
        self.con = ""
        self.data = pd.DataFrame()
        self.con = None

        keys_map = {}

        spectrometer_setting_map = {
            "Coil Current [mA]": "coil_current [mA]",
            "Pre Defl Y [nU]": "pre_deflector_y_current [nU]",
            "Pre Defl X [nU]": "pre_deflector_x_current [nU]",
            "L1 [nU]": "lens1_voltage [nU]",
            "L2 [nU]": "lens2_voltage [nU]",
            "Focus Displacement 1 [nu]": "focus_displacement_current [nU]",
            "Detector Voltage [V]": "detector_voltage [V]",
            "Bias Voltage Electrons [V]": "bias_voltage_electrons [V]",
            "Bias Voltage Ions [V]": "bias_voltage_ions [V]",
        }

        self.sql_metadata_map = {
            "EnergyType": "x_units",
            "EpassOrRR": "pass_energy",
            "Wf": "workfunction",
            "Timestamp": "time_stamp",
            "Samples": "n_values",
            "ElectronEnergy": "start_energy",
            "Step": "step_size",
        }

        self.key_maps = [
            keys_map,
            spectrometer_setting_map,
            source_setting_map,
            self.sql_metadata_map,
        ]

    def parse_file(self, file, **kwargs):
        """
        Parse the file's metadata into a flat list of dictionaries.

        Parameters
        ----------
        file : str
            Filepath of the SLE file to be read.

        Returns
        -------
        self.spectra
            Flat list of dictionaries containing measured
            metadata of one parameter each.

        """
        # read and parse sle file
        data = self._read_csv()
        data = restructure_data(data)

        return data

    def _read_csv(self):
        pass

    def restructure_data(self, data):
        return data

    def _convert_date_time(self, timestamp):
        """
        Convert the native time format to the one we decide to use.
        Returns datetime string in the format '%Y-%b-%d %H:%M:%S.%f'.

        """
        date_time = datetime.strptime(timestamp, "%Y-%b-%d %H:%M:%S.%f")
        date_time = datetime.strftime(date_time, "%Y-%m-%d %H:%M:%S.%f")
        return date_time

    def _re_map_keys(self, dictionary, key_map):
        """
        Map the keys returned from the SQL table to the preferred keys for
        the parser output.

        """
        keys = list(key_map.keys())
        for k in keys:
            if k in dictionary.keys():
                dictionary[key_map[k]] = dictionary.pop(k)
        return dictionary

    def _drop_unused_keys(self, dictionary, keys_to_drop):
        """
        Remove any keys parsed from sle that are not needed

        Parameters
        ----------
        dictionary : dict
            Dictionary with data and metadata for a spectrum.
        keys_to_drop : list
            List of metadata keys that are not needed.

        Returns
        -------
        None.

        """
        for key in keys_to_drop:
            if key in dictionary.keys():
                dictionary.pop(key)

    def _re_map_values(self, dictionary):
        """
        Map the values returned from the SQL table to the preferred format.

        Parameters
        ----------
        dictionary : dict
            Dictionary with data and metadata for a spectrum.

        Returns
        -------
        dictionary : dict
            Dictionary with data and metadata for a spectrum with
            preferred keys for values.

        """
        for key, values in self.value_map.items():
            dictionary[key] = values(dictionary[key])
        return dictionary

    def _convert_to_common_format(self):
        """
        Reformat spectra into the format needed for the Converter object
        """
        maps = {}
        for key_map in self.key_maps:
            maps.update(key_map)
        for spec in self.spectra:
            self._re_map_keys(spec, maps)
            self._re_map_values(spec)
            self._drop_unused_keys(spec, self.keys_to_drop)
            spec["data"] = {}
            spec["data"]["x"] = self._get_energy_data(spec)

            channels = [
                key
                for key in spec
                if any(name in key for name in ["cps_ch_", "cps_calib"])
            ]

            for channel_key in channels:
                spec["data"][channel_key] = np.array(spec[channel_key])
            for channel_key in channels:
                spec.pop(channel_key)

            spec["y_units"] = "Counts per Second"
