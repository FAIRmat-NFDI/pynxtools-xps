"""
Class for reading XPS files from Scienta spectrometers.
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

import re
import copy
import warnings
from pathlib import Path
from typing import Any, Dict, List, Union, Tuple
import xarray as xr
import numpy as np
from igor2 import binarywave

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_entry_name,
    construct_data_key,
    _re_map_single_value,
    _check_valid_value,
)
from pynxtools_xps.value_mappers import get_units_for_key, convert_units

from pynxtools_xps.scienta.scienta_data_model import ScientaHeader, ScientaRegion
from pynxtools_xps.scienta.scienta_mappings import (
    UNITS,
    VALUE_MAP,
    _get_key_value_pair,
    _construct_date_time,
)


class MapperScienta(XPSMapper):
    """
    Class for restructuring data from
    Scienta spectrometers into a structured python
    dictionaries.
    """

    config_file = "config_scienta.json"

    __prmt_file_ext__ = [
        "ibw",
        "txt",
    ]

    __file_err_msg__ = (
        "The Scienta reader currently only allows files with "
        "the following extensions: "
        f"{__prmt_file_ext__}."
    )

    def _select_parser(self):
        """
        Select Scienta parser based on the file extension.

        Returns
        -------
        ScientaParser
            Parser for reading .txt or .ibw files exported by Scienta.

        """
        if str(self.file).endswith(".txt"):
            return ScientaTxtParser()
        elif str(self.file).endswith(".ibw"):
            return ScientaIgorParser()
        raise ValueError(MapperScienta.__file_err_msg__)

    def construct_data(self):
        """Map Parser data to NXmpes-ready dict."""
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
            "electronanalyser": [],
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
                "intensity/@units",
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
        entry_parts = []
        for part in ["spectrum_type", "region_name"]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        file_parent = f"{entry_parent}/file_info"
        instrument_parent = f"{entry_parent}/instrument"
        analyser_parent = f"{instrument_parent}/electronanalyser"

        path_map = {
            "file_info": f"{file_parent}",
            "user": f"{entry_parent}/user",
            "instrument": f"{instrument_parent}",
            "source_xray": f"{instrument_parent}/source_xray",
            "beam_xray": f"{instrument_parent}/beam_xray",
            "electronanalyser": f"{analyser_parent}",
            "collectioncolumn": f"{analyser_parent}/collectioncolumn",
            "energydispersion": f"{analyser_parent}/energydispersion",
            "detector": f"{analyser_parent}/detector",
            "manipulator": f"{instrument_parent}/manipulator",
            "calibration": f"{instrument_parent}/calibration",
            "sample": f"{entry_parent}/sample",
            "data": f"{entry_parent}/data",
            "region": f"{entry_parent}/region",
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
                if units is not None:
                    self._xps_dict[f"{root}/{mpes_key}/@units"] = units

        # Create key for writing to data
        scan_key = construct_data_key(spectrum)

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
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
        self.spectra: List[Dict[str, Any]] = []

    def parse_file(self, file: Union[str, Path], **kwargs):
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
            key, value = _get_key_value_pair(line)
            if key:
                setattr(self.header, key, value)
        self.header.validate_types()

    def _parse_region(self, region_id: int):
        """
        Parse data from one region (i.e., one measured spectrum)
        into a dictionary and append to all spectra.

        Parameters
        ----------
        region_id : int
            Number of the region in the file.

        Returns
        -------
        None.

        """
        region = ScientaRegion(region_id=region_id)

        bool_variables = {
            "in_region": False,
            "in_region_info": False,
            "in_run_mode_info": False,
            "in_ui_info": False,
            "in_manipulator": False,
            "in_data": False,
        }

        energies: List[float] = []
        intensities: List[float] = []

        line_start_patterns = {
            "in_region": f"[Region {region_id}",
            "in_region_info": f"[Info {region_id}",
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

            if any(
                [
                    bool_variables["in_region"],
                    bool_variables["in_region_info"],
                    bool_variables["in_run_mode_info"],
                ]
            ):
                # Read instrument meta data for this region.
                key, value = _get_key_value_pair(line)
                if _check_valid_value(value):
                    setattr(region, key, value)

            if bool_variables["in_ui_info"]:
                key, value = _get_key_value_pair(line)
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
        region_dict["intensity/@units"] = "counts"

        self.spectra.append(region_dict)


class ScientaIgorParser:
    """Parser for Scienta TXT exports."""

    # pylint: disable=too-few-public-methods

    def __init__(self):
        self.lines: List[str] = []
        self.spectra: List[Dict[str, Any]] = []

    def parse_file(self, file: Union[str, Path], **kwargs):
        """
        Reads the igor binarywave files and returns a list of
        dictionary containing the wave data.

        Parameters
        ----------
        file : str
            Filepath of the TXT file to be read.

        Returns
        -------
        self.spectra
            Flat list of dictionaries containing one spectrum each.

        """
        ibw = binarywave.load(file)
        ibw_version, wave = ibw["version"], ibw["wave"]

        notes = self._parse_note(wave["note"])

        data_unit_label, data_unit = self._parse_unit(wave["data_units"])
        dimension_unit_label, dimension_unit = self._parse_unit(wave["dimension_units"])

        wave_header = wave["wave_header"]
        data = wave["wData"]

        # Not needed at the moment.
        # TODO: Add support for formulas if they are written by the
        # measurement software.
        # formula = wave["formula"]
        # labels = wave["labels"]
        # spectrum_indices = wave["sIndices"]
        # bin_header = wave["bin_header"]

        if len(data.shape) == 1:
            self.no_of_regions = 1
        else:
            self.no_of_regions = data.shape[0]

        for region_id in range(0, self.no_of_regions):
            region = ScientaRegion(region_id=region_id)
            region_fields = list(region.__dataclass_fields__.keys())
            overwritten_fields = ["region_id", "time_stamp", "data"]
            unused_notes_keys = []

            for key, note in notes.items():
                if _check_valid_value(note):
                    if key in region_fields:
                        setattr(region, key, note)
                        overwritten_fields += [key]
                    else:
                        unused_notes_keys += [key]

            energies = self.axis_for_dim(wave_header, dim=region_id)
            axs_unit = self.axis_units_for_dim(wave_header, dim=region_id)

            if data.ndim == 1:
                intensities = data
            else:
                intensities = data[region_id]

            # Convert date and time to ISO8601 date time.
            region.time_stamp = _construct_date_time(
                region.start_date, region.start_time
            )

            region.energy_size = len(energies)
            region.energy_axis = energies

            region.data = {
                "energy": np.array(energies),
                "intensity": np.array(intensities),
            }

            region.validate_types()

            spectrum_dict = region.dict()

            for key in unused_notes_keys:
                spectrum_dict[key] = notes[key]

            spectrum_dict["igor_binary_wave_format_version"] = ibw_version
            spectrum_dict["intensity/@units"] = convert_units(data_unit_label)

            self.spectra.append(spectrum_dict)

        return self.spectra

    def _parse_note(self, bnote: bytes) -> Dict[str, Any]:
        """
        Parses the note field of the igor binarywave file.
        It assumes that the note field contains key-value pairs
        of the form 'key=value' separated by newlines.

        Parameters
        ----------
        bnote : bytes
            The bytes of the binarywave note field.

        Returns
        -------
        Dict[str, Any]
            The dictionary of the parsed note field.

        """
        note = bnote.decode("utf-8").replace("\r", "\n")

        notes = {}

        for line in note.split("\n"):
            key, value = _get_key_value_pair(line)
            if key:
                notes[key] = value

        return notes

    def _parse_unit(self, bunit: bytes) -> Tuple[str, object]:
        """
        Extracts the label and unit from a string containing a label
        followed by a unit enclosed in square brackets.

        Parameters
        ----------
        bunit: bytes
            The input string containing the label and unit.

        Returns
        -------
        tuple
            A tuple containing:
            - unit_label : str
                The extracted label.
            - unit : str
                The extracted unit.
        """
        unit = bunit.decode("utf-8").replace("\r", "\n")

        pattern = r"([\w\s]+)\s*\[([\w\s.]+)\]"
        matches = re.match(pattern, unit)
        if matches is not None:
            label = matches.group(1).strip()
            unit = matches.group(2).strip()
            return label, unit
        return "", ""

    def axis_for_dim(self, wave_header: Dict[str, Any], dim: int) -> np.ndarray:
        """
        Returns the axis values for a given dimension from the wave header.

        Parameters
        ----------
        wave_header : Dict[str, Any]
            The wave_header of the ibw file.
        dim : int
            The dimension to return the axis for..

        Returns
        -------
        np.ndarray
            Axis values for a given dimension.

        """
        return (
            wave_header["sfA"][dim] * np.arange(wave_header["nDim"][dim])
            + wave_header["sfB"][dim]
        )

    def axis_units_for_dim(self, wave_header: Dict[str, Any], dim: int) -> str:
        """
        Returns the unit for a given dimension from the wave header.

        Parameters
        ----------
        wave_header : Dict[str, Any]
            The wave_header of the ibw file.
        dim : int
            The dimension to return the axis for..

        Returns
        -------
        str:
            The axis units.

        """
        unit_arr = wave_header["dimUnits"][dim]

        unit = ""
        for elem in unit_arr:
            unit += elem.decode("utf-8")

        return unit
