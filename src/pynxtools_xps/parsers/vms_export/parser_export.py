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

import itertools
import operator
import re
from abc import abstractmethod
from pathlib import Path
from typing import Any, ClassVar

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from pynxtools_xps.mapping import convert_pascal_to_snake
from pynxtools_xps.numerics import (
    _get_minimal_step,
    check_uniform_step_width,
    interpolate_arrays,
)
from pynxtools_xps.parsers.base import ParsedSpectrum, _construct_entry_name, _XPSParser
from pynxtools_xps.parsers.vms_export.metadata import _context, _handle_repetitions


def _select_from_list(input_list: list[str], skip: int, keep_middle: int) -> list[str]:
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


def _get_dict_keys(header_lines: list[str]) -> list[str]:
    """
    Maps a list of header strings to their corresponding keys based on a predefined mapping.

    Args:
        header_lines (List[str]): A list of header strings to be mapped.

    Returns:
        List[str]: A list of keys, where each header is replaced by its mapped value
                   or left unchanged if no mapping is found.
    """
    return [_context.normalize_key(header) for header in header_lines if header]


class _TextParser(_XPSParser):  # pylint: disable=too-few-public-methods
    """
    Internal parser for ASCII files exported from CasaXPS.

    Not intended for direct use outside this module.  Use
    ``VamasExportParser`` as the public entry point.
    """

    config_file: ClassVar[str] = "config_vms.json"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".txt",)

    def __init__(self):
        super().__init__()
        self.lines: list[str] = []
        self.n_headerlines: int = 7
        self.uniform_energy_steps: bool = True

    def _parse(self, file: Path, **kwargs) -> None:
        """
        Parse the file into a list of dictionaries.

        Parsed data stored in the attribute 'self._data'.

        Parameters
        ----------
        file : Path
            XPS data filepath.
        uniform_energy_steps : bool, optional
            If true, the spectra are interpolated to have uniform
            energy steps. The default is True.
        n_headerlines : int, optional
            Number of header lines in each data block.

        """
        self.n_headerlines = kwargs.pop("n_headerlines", self.n_headerlines)
        self.uniform_energy_steps = kwargs.pop(
            "uniform_energy_steps",
            self.uniform_energy_steps,
        )

        self._read_lines(file)
        blocks = self._parse_blocks()

        self._flat_spectra = self._build_list_of_dicts(blocks)

    def _read_lines(self, file):
        """
        Read in all lines from the file as a list of strings.

        Parameters
        ----------
        file : str
            XPS data filepath.

        """
        with open(file, encoding="utf-8") as txt_file:
            for line in txt_file:
                self.lines += [line]

    @abstractmethod
    def _parse_blocks(self) -> list[list[str]]:
        """
        Extract spectrum blocks from full data string.

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

        Parameters
        ----------
        blocks : list
            List of data blocks containing one spectrum each.

        Returns
        -------
        spectra : list
            List of dicts with spectrum data and metadata.

        """
        ...

    def _separate_header_and_data(self, block):
        """
        Separate header (with metadata) from data for one measurement
        block.

        """
        header = block[: self.n_headerlines]
        data = block[self.n_headerlines :]

        return header, data


class _TextParserRows(_TextParser):
    """
    Internal parser for ASCII files exported from CasaXPS using the
    'Rows of Tables' option.
    """

    def __init__(self):
        super().__init__()
        self.n_headerlines: int = 7

    def matches_file(self, file: Path) -> bool:
        """Return True for CasaXPS VAMAS text exports in rows layout)."""
        try:
            with open(file, encoding="utf-8", errors="ignore") as f:
                first = f.readline()
            # Rows format: "Characteristic Energy eV\t<float>\tAcquisition Time s\t..."
            return "Characteristic Energy eV" in first and "Acquisition Time s" in first
        except Exception:
            return False

    def _parse_blocks(self) -> list[list[str]]:
        """
        With the 'Rows of Tables' option, there is only one block
        with common metadata.

        """
        return [self.lines]

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
        lines = blocks[0]
        header, data_lines = self._separate_header_and_data(lines)
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
                y_units = spec_header.split(":")[-1].strip()
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
                "binding_energy/data": np.array(x_bin),
                "counts_per_second/data": np.array(intensity).squeeze(),
                "start_energy": x_bin[0],
                "stop_energy": x_bin[-1],
                "energy_type": "binding",
                "y_units": "counts_per_second",
            }

            if check_uniform_step_width(x_bin):
                spectrum["step_size"] = _get_minimal_step(x_bin)

            data += [spectrum]

        return data


class _TextParserColumns(_TextParser):
    """
    Internal parser for ASCII files exported from CasaXPS using the
    'Columns of Tables' option.
    """

    def __init__(self):
        super().__init__()
        self.n_headerlines = 8

    def matches_file(self, file: Path) -> bool:
        """Return True for CasaXPS VAMAS text exports in columns layout)."""
        try:
            with open(file, encoding="utf-8", errors="ignore") as f:
                first = f.readline()
            # Columns format: "Cycle N:GroupName:SpectrumType\t..."
            if first.startswith("Cycle ") and ":" in first:
                return True
        except Exception:
            return False

        return False

    def _parse_blocks(self) -> list[list[str]]:
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
                comp_names = _handle_repetitions(_get_dict_keys(line.split("\t")[1:]))
                for comp_index, comp_name in enumerate(comp_names):
                    data[f"component{comp_index}"] = {
                        "name": comp_name,
                        "data": [],
                        "data_cps": [],
                    }
                n_components = len(comp_names)

            elif line.startswith(("Area", "Width", "FWHM", "Position", "data")):
                line_split = line.split("\t")
                fit_data[line_split[0]] = [
                    _context._format_value(val) for val in line_split[1:]
                ]

            elif "K.E." in line and "Counts" in line:
                # Parse column headers
                all_names = _get_dict_keys(line.split("\t"))

                keep_middle = 2

                for name in ["background_intensity", "fit_sum"]:
                    if name in all_names:
                        keep_middle += 1

                names = _select_from_list(
                    all_names, skip=n_components, keep_middle=keep_middle
                )

                names = _handle_repetitions(names)

                non_comp_names = []

                for name in names:
                    if (
                        not any(
                            sub_dict.get("name") == name for sub_dict in data.values()
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
                    for key, sub_dict in data.items():
                        if sub_dict.get("name") == name:
                            matching_key = key
                            break

                    if name not in lineshape_in:
                        data[matching_key]["data"].append(_context._format_value(value))
                        lineshape_in += [name]
                    else:
                        data[matching_key]["data_cps"].append(
                            _context._format_value(value)
                        )
        flattened = {}
        for i, (sup_key, sub_dict) in enumerate(data.items()):
            for sub_key, value in sub_dict.items():
                if sup_key == value:
                    continue
                if value and any(str(val).strip() for val in value):
                    if "data" in sub_key:
                        value = np.array(value)
                    flattened[f"{sup_key}/{sub_key}"] = value

            for param in ("Area", "Width", "FWHM", "Position", "data"):
                if param in fit_data and i < len(fit_data[param]):
                    param_value = fit_data[param][i]
                    if param_value:
                        param_key, param_value, unit = _context.format(
                            param, param_value
                        )
                        flattened[f"{sup_key}/{param_key}"] = param_value
                        if unit:
                            flattened[f"{sup_key}/{param_key}/@units"] = unit

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
                parsed_data["step_size"] = _get_minimal_step(
                    parsed_data["kinetic_energy/data"]
                )

            parsed_data["energy_label"] = "binding"

            spectra += [parsed_data]

            plt.show()

        return spectra


class VamasExportParser(_XPSParser):
    """
    Parser for ASCII files exported from CasaXPS (from Vamas).

    Supports both 'Rows of Tables' and 'Columns of Tables' export formats.
    Populates ``self._data`` with one ``ParsedSpectrum`` per
    (group_name, spectrum_type) pair.
    """

    config_file: ClassVar[str] = "config_vms.json"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".txt",)

    _SUB_PARSERS: ClassVar[tuple[type[_TextParser], ...]] = (
        _TextParserRows,
        _TextParserColumns,
    )

    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset(
        {
            "binding_energy/data",
            "kinetic_energy/data",
            "counts_per_second/data",
            "counts/data",
            "counts/data_cps",
        }
    )

    def matches_file(self, file: Path) -> bool:
        """Return True if any txt export variant recognizes the file."""
        return any(cls().matches_file(file) for cls in self._SUB_PARSERS)

    def _parse(self, file: Path, **kwargs) -> None:
        """
        Parse the CasaXPS text export and populate ``self._data``.

        One ``ParsedSpectrum`` is built per (group_name, spectrum_type) pair.
        The intensity DataArray has dims ``("cycle", "scan", "energy")``.
        Fitting data (component*/data, background*/data, fit_sum/data) is
        stored in ``metadata`` as-is.

        Dispatch to the first matching strategy and merge its data

        Parameters
        ----------
        file : Path
            XPS data filepath.
        n_headerlines : int, optional
            Number of header lines per data block.
        uniform_energy_steps : bool, optional
            If True, spectra are interpolated to have uniform energy steps.

        """
        for sub_parser_cls in self._SUB_PARSERS:
            sub_parser = sub_parser_cls()
            if sub_parser.matches_file(file):
                sub_parser._parse(file, **kwargs)
                self._flat_spectra = sub_parser._flat_spectra
                self._build_parsed_spectra()
                return
        raise ValueError(
            f"{self.__class__.__name__}: no sub parser matched file '{file}'"
        )

    def _build_parsed_spectra(self) -> None:
        """Group flat spectra by entry name and build ``self._data``."""
        entries: dict[str, list[dict[str, Any]]] = {}
        for spectrum in self._flat_spectra:
            entry_parts = [
                str(spectrum.get(part, ""))
                for part in ("group_name", "spectrum_type")
                if spectrum.get(part)
            ]
            entry_name = _construct_entry_name(entry_parts)
            entries.setdefault(entry_name, []).append(spectrum)

        for entry_name, spectra in entries.items():
            self._data[entry_name] = self._assemble_entry(spectra)

    def _assemble_entry(self, spectra: list[dict[str, Any]]) -> ParsedSpectrum:
        """
        Build a ``ParsedSpectrum`` from all scan dicts for one entry.

        Parameters
        ----------
        spectra : list[dict]
            All spectrum dicts from the flat parser output that share the
            same entry name.

        Returns
        -------
        ParsedSpectrum

        """
        # Resolve energy axis from the first spectrum.
        first = spectra[0]
        if "binding_energy/data" in first:
            energy_key = "binding_energy/data"
        else:
            energy_key = "kinetic_energy/data"

        # Resolve intensity axis.
        if "counts_per_second/data" in first:
            intensity_key = "counts_per_second/data"
        else:
            intensity_key = "counts/data"

        energy_axis = np.array(first[energy_key])

        # Stack intensities: shape (1, n_scans, n_energy).
        intensities = np.stack([np.array(s[intensity_key]) for s in spectra], axis=0)
        intensities = intensities[np.newaxis, ...]  # (1, n_scans, n_energy)

        data_array = xr.DataArray(
            data=intensities,
            dims=("cycle", "scan", "energy"),
            coords={"energy": energy_axis},
        )

        # Metadata: all keys except the primary energy/intensity data arrays.
        # Fitting data (component*/data, background*/data, fit_sum/data) is
        # kept as-is in metadata.
        metadata = self._filter_metadata(first)

        return ParsedSpectrum(data=data_array, raw=None, metadata=metadata)
