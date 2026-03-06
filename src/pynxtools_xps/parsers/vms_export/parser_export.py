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
Parser for reading XPS data from TXT export of CasaXPS.
"""

import itertools
import operator
import re
from abc import abstractmethod
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import xarray as xr

from pynxtools_xps.mapping import convert_pascal_to_snake
from pynxtools_xps.numerics import (
    _get_minimal_step,
    check_uniform_step_width,
    interpolate_arrays,
)
from pynxtools_xps.parsers.base import ParsedSpectrum, _construct_entry_name, _XPSParser
from pynxtools_xps.parsers.vms_export.metadata import _context

# Column-header → flat-dict key mappings for the two halves of the data section.
# The KE-side columns produce ``{key}/data`` arrays; the BE-side produce
# ``{key}/data_cps`` arrays.  Component-name columns are handled separately
# by :func:`_map_data_headers`.
_KE_SIDE_HEADERS: dict[str, str] = {
    "K.E.": "kinetic_energy/data",
    "Counts": "counts/data",
    "Background": "background/data",
    "Envelope": "envelope/data",
}

_BE_SIDE_HEADERS: dict[str, str] = {
    "B.E.": "binding_energy/data",
    "CPS": "counts_per_second/data",
    "Background CPS": "background/data_cps",
    "Envelope CPS": "envelope/data_cps",
}


def _map_data_headers(
    headers: list[str],
    comp_names: list[str],
    special_map: dict[str, str],
    comp_suffix: str,
) -> list[str]:
    """Map raw column headers from the data section to flat-dict keys.

    Special headers are resolved via *special_map*.  Headers that match a
    component name are assigned ``component{i}{comp_suffix}`` in order of
    first appearance.  All other headers are normalized with
    ``_context.normalize_key`` and *comp_suffix* appended.

    Args:
        headers:      Ordered raw column header strings (whitespace stripped).
        comp_names:   Ordered component names from the ``Name`` row.
        special_map:  Explicit header → final key mapping.
        comp_suffix:  Appended to component keys (``"/data"`` for the KE
                      side, ``"/data_cps"`` for the BE side).

    Returns:
        Flat-dict keys in the same order as *headers*.
    """
    comp_idx = 0
    keys: list[str] = []
    for h in headers:
        if h in special_map:
            keys.append(special_map[h])
        elif h in comp_names:
            keys.append(f"component{comp_idx}{comp_suffix}")
            comp_idx += 1
        else:
            keys.append(_context.normalize_key(h) + comp_suffix)
    return keys


class _TextParser(_XPSParser):
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

    def _parse_block_data(self, block_lines: list[str]) -> dict[str, Any]:
        """Parse one ``Cycle``-headed block into a flat metadata + data dict.

        The data section header contains both ``K.E.`` and ``B.E.``, separated
        by an empty tab cell (``\\t\\t``).  The left half (KE side) produces
        ``{key}/data`` arrays; the right half (BE side) produces
        ``{key}/data_cps`` arrays.  Component columns are assigned
        ``component{N}/data`` and ``component{N}/data_cps`` keys in order.

        Header rows (``Name``, ``Position``, ``FWHM``, ``Area``, ``Width``,
        ``Lineshape``) produce scalar ``component{N}/{field}`` metadata keys.

        Args:
            block_lines: Raw text lines for one spectral block.

        Returns:
            Flat dict with scalar metadata and 1-D ``np.ndarray`` values.
        """
        _fit_rows: frozenset[str] = frozenset(
            {"Position", "FWHM", "Area", "Width", "Lineshape"}
        )

        metadata: dict[str, Any] = {}
        comp_names: list[str] = []
        comp_fit: dict[int, dict[str, Any]] = {}
        ke_keys: list[str] = []
        be_keys: list[str] = []
        ke_rows: list[list[float]] = []
        be_rows: list[list[float]] = []
        in_data = False

        for raw_line in block_lines:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("Cycle"):
                parts = stripped.split(":")
                metadata["cycle"] = parts[0].strip()
                if len(parts) > 1:
                    metadata["group_name"] = parts[1].strip()
                if len(parts) > 2:
                    metadata["spectrum_type"] = parts[2].strip()

            elif stripped.startswith("Characteristic Energy"):
                for key_raw, val_str in re.findall(r"([\w\s]+)\t([\deE+\-.]+)", line):
                    key_raw = key_raw.strip()
                    key, _, unit = key_raw.rpartition(" ")
                    key = convert_pascal_to_snake(key.strip())
                    metadata[key] = float(val_str)
                    metadata[f"{key}/@units"] = unit

            elif stripped.startswith("Name"):
                comp_names = [c for c in line.split("\t")[1:] if c.strip()]

            elif any(stripped.startswith(fr) for fr in _fit_rows):
                cols = line.split("\t")
                field = _context.normalize_key(cols[0].strip())
                values = [v for v in cols[1:] if v.strip()]
                for i, val_str in enumerate(values):
                    comp_fit.setdefault(i, {})[field] = _context._format_value(val_str)

            elif "K.E." in stripped and "B.E." in stripped:
                # Data header — two halves separated by an empty tab cell.
                halves = line.split("\t\t", 1)
                left_headers = [h.strip() for h in halves[0].split("\t") if h.strip()]
                right_headers = (
                    [h.strip() for h in halves[1].split("\t") if h.strip()]
                    if len(halves) > 1
                    else []
                )
                ke_keys = _map_data_headers(
                    left_headers, comp_names, _KE_SIDE_HEADERS, "/data"
                )
                be_keys = _map_data_headers(
                    right_headers, comp_names, _BE_SIDE_HEADERS, "/data_cps"
                )
                in_data = True

            elif in_data:
                halves = line.split("\t\t", 1)
                left_vals = [v.strip() for v in halves[0].split("\t") if v.strip()]
                right_vals = (
                    [v.strip() for v in halves[1].split("\t") if v.strip()]
                    if len(halves) > 1
                    else []
                )
                if left_vals:
                    ke_rows.append([float(v) for v in left_vals[: len(ke_keys)]])
                if right_vals:
                    be_rows.append([float(v) for v in right_vals[: len(be_keys)]])

        flat: dict[str, Any] = {**metadata}

        # Component fit metadata
        for i, name in enumerate(comp_names):
            flat[f"component{i}/name"] = name
            for field, val in comp_fit.get(i, {}).items():
                flat[f"component{i}/{field}"] = val

        # KE-side data arrays
        if ke_rows:
            ke_arr = np.array(ke_rows)
            for j, key in enumerate(ke_keys):
                flat[key] = ke_arr[:, j]

        # BE-side data arrays
        if be_rows:
            be_arr = np.array(be_rows)
            for j, key in enumerate(be_keys):
                flat[key] = be_arr[:, j]

        return flat

    def _build_list_of_dicts(self, blocks: list[list[str]]) -> list[dict[str, Any]]:
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
            parsed = self._parse_block_data(block)

            # Fallback: derive binding energy axis if the BE column was absent.
            if "binding_energy/data" not in parsed and "kinetic_energy/data" in parsed:
                parsed["binding_energy/data"] = (
                    parsed["characteristic_energy"] - parsed["kinetic_energy/data"]
                )

            # Fallback: derive CPS if the CPS column was absent.
            if "counts_per_second/data" not in parsed and "counts/data" in parsed:
                parsed["counts_per_second/data"] = (
                    parsed["counts/data"] / parsed["acquisition_time"]
                )

            # Resample all arrays to a uniform energy grid when requested.
            if self.uniform_energy_steps:
                x_key = (
                    "kinetic_energy/data"
                    if "kinetic_energy/data" in parsed
                    else "binding_energy/data"
                )
                x_arr = parsed.get(x_key)
                if x_arr is not None and not check_uniform_step_width(x_arr):
                    y_keys = [
                        k
                        for k, v in parsed.items()
                        if isinstance(v, np.ndarray) and k != x_key
                    ]
                    x_uniform, resampled = interpolate_arrays(
                        x_arr, [parsed[k] for k in y_keys]
                    )
                    parsed[x_key] = x_uniform
                    for k, v in zip(y_keys, resampled):
                        parsed[k] = v

            ke_data = parsed.get("kinetic_energy/data")
            if ke_data is not None and check_uniform_step_width(ke_data):
                parsed["step_size"] = _get_minimal_step(ke_data)

            parsed["energy_label"] = "binding"
            spectra.append(parsed)

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
