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
# pylint: disable=too-many-lines,too-few-public-methods
"""
Classes for reading XPS files from TXT export of CasaXPS.
"""

import matplotlib.pyplot as plt


import re
import itertools
import operator
import copy
import warnings
from typing import Any, Dict, List
from collections import Counter
from abc import ABC, abstractmethod
import csv
import xarray as xr
import numpy as np

from pynxtools_xps.reader_utils import (
    XPSMapper,
    check_uniform_step_width,
    get_minimal_step,
    interpolate_arrays,
    construct_entry_name,
    construct_data_key,
    convert_pascal_to_snake,
    _format_value,
)
from pynxtools_xps.value_mappers import get_units_for_key, convert_units

UNITS: Dict[str, str] = {
    "step_size": "eV",
}

KEY_MAP: Dict[str, str] = {
    "K.E.": "kinetic_energy",
    "B.E.": "binding_energy",
    "Counts": "counts",
    "CPS": "counts_per_second",
    "Background": "background_intensity",
    "Background CPS": "background_intensity_cps",
    "Envelope": "fit_sum",
    "Envelope CPS": "fit_sum_cps",
    "%At Conc": "atomic_concentration",
}


def handle_repetitions(input_list: List[str]) -> List[str]:
    """
    Process a list of strings to handle repeated items by appending a suffix
    to each duplicate item. The suffix is in the format '_n', where 'n' is the
    occurrence number of that item in the list.

    Parameters:
    - input_list (List[str]): A list of strings where repeated items are
      identified and renamed with a suffix.

    Returns:
    - List[str]: A new list where repeated items are modified by appending
      a suffix to make them unique.
    """
    counts = Counter(input_list)
    result = []
    occurrences = {}

    for item in input_list:
        if counts[item] > 1:
            # If the item has been seen before, add a suffix
            if item not in occurrences:
                occurrences[item] = 0
            occurrences[item] += 1
            result.append(f"{item}_{occurrences[item]}")
        else:
            result.append(item)

    return result


def select_from_list(input_list: List[str], skip: int, keep_middle: int) -> List[str]:
    """
    Select items from a list according to the specified pattern:
    - Extract the first (2 + count + keep_middle) items,
    - Skip the next 'count' number of items,
    - Extract everything after the skipped items.

    Parameters:
    - input_list (List[str]): The list of strings to process.
    - skip (int): The number of items to skip after the initial selection.
    - keep_middle (int): The number of items to keep in the middle.

    Returns:
    - List[str]: The processed list after extracting and skipping items.
    """
    first_part = input_list[: (2 + skip + keep_middle)]
    skip_part = input_list[(2 + skip + keep_middle) : (2 + skip + keep_middle + skip)]
    remaining_part = input_list[(2 + skip + keep_middle + skip) :]

    return first_part + remaining_part


def get_dict_keys(header_lines: List[str]) -> List[str]:
    """
    Maps a list of header strings to their corresponding keys based on a predefined mapping.

    Args:
        header_lines (List[str]): A list of header strings to be mapped.

    Returns:
        List[str]: A list of keys, where each header is replaced by its mapped value
                   or left unchanged if no mapping is found.
    """
    return [KEY_MAP.get(header, header) for header in header_lines if header]


class TxtMapperVamasExport(XPSMapper):
    """
    Class for restructuring .txt data file from
    Casa TXT export (from Vamas) into python dictionary.
    """

    config_file = "config_vms.json"

    def __init__(self):
        self.parser_map = {
            "rows_of_tables": TextParserRows,
            "columns_of_tables": TextParserColumns,
        }
        super().__init__()

    def _get_file_type(self, file):
        """
        Check which export option was used in CasaXPS.

        Parameters
        ----------
        file : str
            XPS data filepath.

        Returns
        -------
        str
            Either columns_of_tables or rows_of_tables.

        """
        with open(file, encoding="utf-8") as txt_file:
            first_line = txt_file.readline()
            if first_line.startswith("Cycle"):
                return "columns_of_tables"
            return "rows_of_tables"

    def _select_parser(self):
        """
        Select parser based on the structure of the text file

        Returns
        -------
        TextParser
            TextParser for CasaXPS export from Vamas files.

        """
        return self.parser_map[self._get_file_type(self.file)]()

    def construct_data(self):
        """Map TXT format to NXmpes-ready dict."""
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: Dict[str, Any] = {}

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: Dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []
        for part in ["group_name", "spectrum_type"]:
            val = spectrum.get(part)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            if key.startswith("entry"):
                entry_parent = "/ENTRY[entry]"
                key = key.replace("entry/", "", 1)
            mpes_key = f"{entry_parent}/{key}"
            self._xps_dict[mpes_key] = value

            units = get_units_for_key(key, UNITS)
            if units is not None:
                self._xps_dict[f"{mpes_key}/@units"] = units

        # Create key for writing to data.
        scan_key = construct_data_key(spectrum)

        energy = np.array(spectrum["kinetic_energy/data"])
        # energy = np.array(spectrum["binding_energy/data"])
        intensity = np.array(spectrum["counts_per_second/data"])
        # intensity = np.array(spectrum["counts/data"])

        # If multiple spectra exist to entry, only create a new
        # xr.Dataset if the entry occurs for the first time.
        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

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

        self._xps_dict["data"][entry][scan_key] = xr.DataArray(
            data=intensity, coords={"energy": energy}
        )


class TextParser(ABC):  # pylint: disable=too-few-public-methods
    """
    Parser for ASCI files exported from CasaXPS.
    """

    def __init__(self):
        self.lines: List[str] = []
        self.n_headerlines: int = 7
        self.uniform_energy_steps: bool = True

    def parse_file(self, file, uniform_energy_steps=True, **kwargs):
        """
        Parse the file into a list of dictionaries.

        Parsed data stored in the attribute 'self.data'.

        Parameters
        ----------
        file : str
            XPS data filepath.
        uniform_energy_steps : bool, optional
            If true, the spectra are interpolate to have uniform
            energy steps. The default is True.
        **kwargs : dict
            n_headerlines: number of header_lines in each data block.

        Returns
        -------
        dict
            DESCRIPTION.

        """
        if "n_headerlines" in kwargs:
            self.n_headerlines = kwargs["n_headerlines"]

        self.uniform_energy_steps = uniform_energy_steps

        self._read_lines(file)
        blocks = self._parse_blocks()
        return self._build_list_of_dicts(blocks)

    def _read_lines(self, file):
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
        with open(file, encoding="utf-8") as txt_file:
            for line in txt_file:
                self.lines += [line]

    @abstractmethod
    def _parse_blocks(self):
        """
        Extract spectrum blocks from full data string.

        This method has to be implemented in the inherited parsers.

        Returns
        -------
        blocks : list
            List of strings, with each string containing one spectrum's
            data and metadata.

        """
        return []

    @abstractmethod
    def _build_list_of_dicts(self, blocks):
        """
        Build list of dictionaries, with each dict containing data
        and metadata of one spectrum (block).

        This method has to be implemented in the inherited parsers.

        Parameters
        ----------
        blocks : list
            List of data blocks containing one spectrum each.

        Returns
        -------
        spectra : list
            List of dicts with spectrum data and metadata.

        """
        return []

    def _separate_header_and_data(self, block):
        """
        Separate header (with metadata) from data for one measurement
        block

        Returns
        -------
        None.

        """
        header = block[: self.n_headerlines]
        data = block[self.n_headerlines :]

        return header, data


class TextParserRows(TextParser):
    """
    Parser for ASCI files exported from CasaXPS using the
    'Rows of Tables' option.
    """

    def __init__(self):
        super().__init__()
        self.n_headerlines: int = 7

    def _parse_blocks(self):
        """
        With the 'Rows of Tables' option, there is only one block
        with common metadata.

        """
        return self.lines

    def _build_list_of_dicts(self, blocks):
        """
        Build list of dictionaries, with each dict containing data
        and metadata of one spectrum.

        Parameters
        ----------
        blocks : list
            List of data blocks containing one spectrum each.

        Returns
        -------
        spectra : list
            List of dicts with spectrum data and metadata.

        """
        header, data_lines = self._separate_header_and_data(blocks)
        settings = self._parse_header(header)
        data = self._parse_data(data_lines)

        return [
            {**spec_settings, **spec_data}
            for spec_settings, spec_data in zip(settings, data)
        ]

    def _parse_header(self, header):
        """
        Parse header into metadata dictionary.

        Parameters
        ----------
        header : str
            Header data for one spectrum as a String.

        Returns
        -------
        settings : list
            List of dicts with measurement settings for
            one spectrum each.

        """
        settings = []
        for spec_header in header[-1].split("\t")[1::3]:
            try:
                group_name = spec_header.split(":")[1]
                region = spec_header.split(":")[2]
                y_units = convert_units(spec_header.split(":")[-1])
            except IndexError:
                group_name = spec_header if spec_header.strip() else "group"
                region = spec_header if spec_header.strip() else "region"
                y_units = "counts_per_second"

            spectrum_settings = {
                "group_name": group_name,
                "spectrum_type": region,
                "energy_type": "binding",
                "y_units": y_units,
            }
            settings += [spectrum_settings]

        return settings

    def _parse_data(self, data_lines):
        """
        Extract energy and intensity data.

        Parameters
        ----------
        data_lines : list
            List of lines with measurement data.

        Returns
        -------
        list
            List of dicts containing the binding energy
            and the intensity axes of one spectrum each.

        """
        data_lines = [x.split("\t") for x in data_lines]
        for line in data_lines:
            del line[2::3]
            del line[-1]

        lines = [[] for _ in range(max(len(line) for line in data_lines))]

        for line in data_lines:
            for i, data_point in enumerate(line):
                try:
                    lines[i].append(float(data_point.strip("\n")))
                except ValueError:
                    pass

        data = []

        for x_bin, intensity in zip(lines[::2], lines[1::2]):
            x_bin, intensity = np.array(x_bin), np.array(intensity)

            if self.uniform_energy_steps and not check_uniform_step_width(x_bin):
                x_bin, intensity = interpolate_arrays(x_bin, intensity)

            spectrum = {
                "data": {
                    "binding_energy": np.array(x_bin),
                    "intensity": np.array(intensity).squeeze(),
                },
                "start_energy": x_bin[0],
                "stop_energy": x_bin[-1],
                "energy_type": "binding",
                "y_units": "counts_per_second",
            }

            if check_uniform_step_width(x_bin):
                spectrum["step_size"] = get_minimal_step(x_bin)

            data += [spectrum]

        return data


class TextParserColumns(TextParser):
    """
    Parser for ASCI files exported from CasaXPS using the
    'Columns of Tables' option.
    """

    def __init__(self):
        super().__init__()
        self.n_headerlines = 8

    def _parse_blocks(self):
        """
        Extract spectrum blocks from full data string.

        Returns
        -------
        blocks : list
            List of strings, with each string containing one spectrum's
            data and metadata.

        """
        blocks = [
            list(g) for _, g in itertools.groupby(self.lines, lambda i: "Cycle " in i)
        ]
        blocks = [operator.add(*blocks[i : i + 2]) for i in range(0, len(blocks), 2)]

        return blocks

    def _parse_block_data(self, block_lines):
        """
        Parses a block of spectral data into metadata and a DataFrame of measurements.

        Args:
            block (list of str): The raw lines of the spectral data block.

        Returns:
            dict: A dictionary with metadata and a DataFrame of measurements.
        """

        metadata = {}
        data = {}
        fit_data = {}

        in_data_section = False

        for line in block_lines:
            line = line.strip()
            if line.startswith("Cycle"):
                # Extract cycle and scan type
                metadata["cycle"], metadata["source"], metadata["spectrum_type"] = map(
                    str.strip, line.split(":")
                )

            elif line.startswith("Characteristic Energy"):
                # Parse characteristic energy and acquisition time
                metadata_match = re.findall(r"([\w\s]+)\t([\deE\+\-.]+)", line)
                matches = {key.strip(): float(value) for key, value in metadata_match}
                for key, value in matches.items():
                    key, unit = key.rsplit(" ", 1)
                    key = convert_pascal_to_snake(key)
                    metadata.update({key: value, f"{key}/@units": unit})

            elif line.startswith("Name"):
                comp_names = handle_repetitions(get_dict_keys(line.split("\t")[1:]))
                for comp_index, comp_name in enumerate(comp_names):
                    data[f"component{comp_index}"] = {
                        "name": comp_name,
                        "data": [],
                        "data_cps": [],
                    }
                n_components = len(comp_names)

            elif line.startswith(("Area", "Width", "Position", "data")):
                line_split = line.split("\t")
                fit_data[line_split[0]] = [_format_value(val) for val in line_split[1:]]

            elif "K.E." in line and "Counts" in line:
                # Parse column headers
                all_names = get_dict_keys(line.split("\t"))

                keep_middle = 2

                for name in ["background_intensity", "fit_sum"]:
                    if name in all_names:
                        keep_middle += 1

                names = select_from_list(
                    all_names, skip=n_components, keep_middle=keep_middle
                )

                names = handle_repetitions(names)

                non_comp_names = []

                for name in names:
                    if (
                        not any(
                            subdict.get("name") == name for subdict in data.values()
                        )
                        and "CPS" not in name
                    ):
                        data[name] = {"name": name, "data": [], "data_cps": []}

                        non_comp_names += [name]

                # Circumvents the problem that there are two columns for
                # each component, but two components can also have the
                # same name.
                new_names = (
                    non_comp_names[:2]
                    + comp_names
                    + non_comp_names[2 : 2 + keep_middle]
                    + comp_names
                    + non_comp_names[2 + keep_middle :]
                )

                in_data_section = True

            elif in_data_section:
                values = [val for val in line.split("\t") if val]

                assert len(values) == len(new_names), f"{new_names}"

                lineshape_in = []
                for name, value in zip(new_names, values):
                    matching_key = name
                    for key, subdict in data.items():
                        if subdict.get("name") == name:
                            matching_key = key
                            break

                    if name not in lineshape_in:
                        data[matching_key]["data"].append(_format_value(value))
                        lineshape_in += [name]
                    else:
                        data[matching_key]["data_cps"].append(_format_value(value))
        flattened = {}
        for i, (supkey, subdict) in enumerate(data.items()):
            for subkey, value in subdict.items():
                if supkey == value:
                    continue
                if value and any(str(val).strip() for val in value):
                    if "data" in subkey:
                        value = np.array(value)
                    flattened[f"{supkey}/{subkey}"] = value

            for param in ("Area", "FWHM", "Position", "data"):
                if param in fit_data and i < len(fit_data[param]):
                    param_value = _format_value(fit_data[param][i])
                    if param_value:
                        flattened[f"{supkey}/{param.lower()}"] = param_value

        if self.uniform_energy_steps:
            uniform = False

            try:
                x_arr = flattened["kinetic_energy/data"]
                uniform = check_uniform_step_width(x_arr)
            except KeyError:
                x_arr = flattened["binding_energy/data"]
                uniform = check_uniform_step_width(x_arr)

            if not uniform:
                return {**metadata, **flattened}
            else:
                uniform_dict = {}
                all_arrays = {}

                for key, value in flattened.copy().items():
                    if isinstance(value, np.ndarray):
                        all_arrays[key] = value
                    else:
                        uniform_dict[key] = value

                x_arr, resampled_arrays = interpolate_arrays(
                    x_arr, list(all_arrays.values())
                )

                for i, key in enumerate(all_arrays):
                    uniform_dict[key] = resampled_arrays[i]

                return {**metadata, **uniform_dict}

        return {**metadata, **flattened}

    def _build_list_of_dicts(self, blocks):
        """
        Build list of dictionaries, with each dict containing data
        and metadata of one spectrum (block).


        Parameters
        ----------
        blocks : list
            List of data blocks containing one spectrum each.

        Returns
        -------
        spectra : list
            List of dicts with spectrum data and metadata.

        """
        spectra = []
        for block in blocks:
            parsed_data = self._parse_block_data(block)

            if "binding_energy/data" not in parsed_data:
                parsed_data["binding_energy/data"] = (
                    parsed_data["characteristic_energy"]
                    - parsed_data["kinetic_energy/data"]
                )

            if "counts_per_second/data" not in parsed_data:
                parsed_data["counts_per_second/data"] = (
                    parsed_data["counts"] / parsed_data["acquisition_time"]
                )

            if check_uniform_step_width(parsed_data["kinetic_energy/data"]):
                parsed_data["step_size"] = get_minimal_step(
                    parsed_data["kinetic_energy/data"]
                )

            parsed_data["energy_label"] = "binding"

            spectra += [parsed_data]

            plt.show()

        return spectra


class CsvMapperVamasResult(XPSMapper):
    """
    Class for restructuring .csv result files from
    Casa report export (from Vamas) into python dictionary.
    """

    config_file = "config_vms.json"

    def __init__(self):
        super().__init__()

    def _select_parser(self):
        """
        Select parser based on the structure of the text file

        Returns
        -------
        TextParser
            TextParser for CasaXPS export from Vamas files.

        """
        return CsvResultParser()

    def construct_data(self):
        self._xps_dict = self.raw_data

    def update_main_file_dict(self, main_file_dicts: List[Dict[str, Any]]):
        """
        Update the dictionaries returned by the main files with specific keys from self.data_dict.

        Args:
            main_file_dicts (List[Dict[str, Any]]): List of dictionaries to update.
        """
        pattern = re.compile(r"(component\d+/)name")
        update_with = {
            "Area/(RSF*T*MFP)",
            "atomic_concentration",
        }  # Use a set for faster lookups

        for existing_dict in main_file_dicts:
            filtered_keys = {
                key: match.group(1)
                for key in existing_dict
                if (match := pattern.search(key))
            }

            for key in filtered_keys:
                value = existing_dict[key]
                if value in self.data_dict:
                    subdict = self.data_dict[value]
                    for subkey in update_with & subdict.keys():
                        new_key = f"{key.rsplit('name', 1)[0]}{subkey}"
                        existing_dict[new_key] = subdict[subkey]


class CsvResultParser:
    def parse_file(self, file: str, **kwargs):
        """
        Parse only the first table from the input file,

        Args:
            file_path (str): Path to the .vms file.

        Returns:
            dict: Parsed data including the file path, header, and rows.
        """
        table_data: Dict[str, Any] = {}
        headers: List[str] = []
        reading_table: bool = False

        with open(file, "r") as f:
            reader = csv.reader(f, delimiter="\t")

            for row in reader:
                if not row:
                    if reading_table:
                        break
                    continue

                # Detect header row
                if row[0].startswith("Name") and not reading_table:
                    headers = get_dict_keys(row)[1:]
                    reading_table = True
                    continue

                # Process rows of the table
                if reading_table:
                    table_data[row[0]] = {}
                    for header, value in zip(headers, row[1:]):
                        if value:
                            formatted_value = _format_value(value)
                            if header == "atomic_concentration" and isinstance(
                                value, (int, float)
                            ):
                                formatted_value /= 100
                            table_data[row[0]].update({header: formatted_value})

        return table_data
