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
Parsers for reading XPS data from Scienta IBW (Igor binary wave) files.

Public API
----------
ScientaIgorParser
    Top-level parser that dispatches between the two IBW variants.
ScientaIgorParserOld
    Direct parser for IBW files produced by the legacy SES software.
ScientaIgorParserPEAK
    Direct parser for IBW files produced by Scienta PEAK software.
"""

import json
import re
from abc import abstractmethod
from pathlib import Path
from typing import Any, ClassVar, cast

import jsonschema
import numpy as np
import xarray as xr
from igor2 import binarywave

from pynxtools_xps.mapping import convert_pascal_to_snake
from pynxtools_xps.parsers.base import ParsedSpectrum, _construct_entry_name, _XPSParser
from pynxtools_xps.parsers.scienta.data_model import (
    ScientaRegion,
    scienta_igor_peak_schema,
)
from pynxtools_xps.parsers.scienta.metadata import (
    _check_valid_value,
    _construct_date_time,
    _context,
    _get_key_value_unit,
)


def _flatten_dict(
    d: dict[str, Any], parent_key: str = "", sep: str = "/"
) -> dict[str, Any]:
    """
    Flatten a nested dictionary into a single level with keys representing the hierarchy.

    Parameters
    ----------
    d : dict[str, Any]
        The dictionary to flatten.
    parent_key : str
        The base key to prepend (used for recursion).
    sep : str
        The separator to use for flattened keys.

    Returns
    -------
    dict[str, Any]
        The flattened dictionary.
    """
    items: list[tuple[str, Any]] = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(_flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


class _IgorWaveParser(_XPSParser):
    """
    Abstract base for IBW (Igor binary wave) parsers.

    Subclasses implement ``_parse_note`` and ``_parse_region_metadata``
    to handle the two note-encoding variants (SES key=value vs. PEAK JSON).
    """

    config_file: ClassVar[str] = "config_scienta.json"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".ibw",)
    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset(
        {"data", "axis_labels", "data_labels", "units"}
    )

    def matches_file(self, file: Path) -> bool:
        """Internal base: format matching is ensured by the concrete subclass."""
        return True

    def __init__(self):
        super().__init__()
        self.lines: list[str] = []
        self._spectrum: dict[str, Any] = {}

    def _parse(self, file: Path, **kwargs) -> None:
        """
        Read the Igor binary wave file and populate ``self._data``.

        Parameters
        ----------
        file : Path
            Filepath of the IBW file to be read.
        """
        ibw = binarywave.load(file)
        ibw_version, wave = ibw["version"], ibw["wave"]

        axes_labels_with_units = self._parse_unit(wave["dimension_units"])
        data_labels_with_units = self._parse_unit(wave["data_units"])

        wave_header = wave["wave_header"]
        data = wave["wData"]

        notes: dict[str, Any] = self._parse_note(wave["note"])

        self.no_of_regions = len(data.shape)

        spectrum: dict[str, Any] = {
            "data": cast(dict[str, Any], {}),
            "axis_labels": cast(list[str], []),
            "data_labels": cast(list[str], []),
            "units": cast(dict[str, Any], {}),
        }

        for i, (dim, unit) in enumerate(axes_labels_with_units):
            if dim in ("kinetic_energy", "binding_energy", "analyser_energy"):
                spectrum["energy_scale"] = _context.map_value("energy_scale", dim)
                dim = "energy"
            spectrum["data"][dim] = self.axis_for_dim(wave_header, dim=i)
            spectrum["axis_labels"].append(dim)
            if not unit:
                unit = _context.get_default_unit(dim)
            spectrum["units"][dim] = _context.map_unit(unit)

        for label, unit in data_labels_with_units:
            spectrum["data"][label] = data
            spectrum["data_labels"].append(label)
            if not unit:
                unit = _context.get_default_unit(label)
            spectrum["units"][label] = _context.map_unit(unit)

        spectrum["igor_binary_wave_format_version"] = ibw_version

        region_metadata = self._parse_region_metadata(region_id=0, notes=notes)
        spectrum |= region_metadata

        self._spectrum = spectrum
        self._build_parsed_spectra()

    def _build_parsed_spectra(self) -> None:
        """Build one ``ParsedSpectrum`` from ``self._spectrum``."""
        spectrum = self._spectrum

        entry_name = _construct_entry_name(
            [spectrum.get("region_name", ""), spectrum.get("spectrum_type", "")]
        )
        if not entry_name:
            entry_name = "region_0"

        # Extract energy axis; fall back to first axis_label
        energy_axis = spectrum["data"].get("energy")
        if energy_axis is None and spectrum["axis_labels"]:
            first_ax = spectrum["axis_labels"][0]
            energy_axis = spectrum["data"].get(first_ax)
            if energy_axis is None:
                energy_axis = np.array([])

        # Extract intensity: first data_label that is not "reduced"
        intensity_arr: np.ndarray | None = None
        for label in spectrum["data_labels"]:
            if "reduced" not in label:
                arr = spectrum["data"].get(label)
                if arr is not None:
                    intensity_arr = np.asarray(arr).squeeze()
                    break
        if intensity_arr is None:
            intensity_arr = np.array([])

        # Extract unit and add to metadata
        metadata = self._filter_metadata(spectrum)
        for axis, unit in spectrum["units"].items():
            metadata[f"{axis}/@units"] = unit

        # Shape: (cycle=1, scan=1, n_energy)
        data_arr = xr.DataArray(
            data=intensity_arr[np.newaxis, np.newaxis, :],
            coords={"energy": energy_axis},
            dims=("cycle", "scan", "energy"),
        )

        self._data[entry_name] = ParsedSpectrum(
            data=data_arr,
            raw=None,
            metadata=metadata,
        )

    @abstractmethod
    def _parse_note(self, bnote: bytes) -> dict[str, Any]:
        pass

    @abstractmethod
    def _parse_region_metadata(
        self, region_id: int, notes: dict[str, Any]
    ) -> dict[str, Any]:
        return {}

    def _parse_unit(self, bunit: bytes) -> list[tuple[str, str | None]]:
        """
        Extract labels and units from a string containing one or more label-unit pairs.

        Parameters
        ----------
        bunit : bytes
            The input bytes containing the label and unit(s).

        Returns
        -------
        list of tuple
            Each tuple contains (unit_label, unit) where unit is None if absent.
        """
        unit = bunit.decode("utf-8").replace("\r", "\n")

        pattern = r"([\w\s]+?)\s*\[([\w\s.]+)\]"

        if matches := re.findall(pattern, unit):
            return [
                (convert_pascal_to_snake(label.strip()), u.strip())
                for label, u in matches
            ]

        return [(convert_pascal_to_snake(unit.strip()), None)]

    def axis_for_dim(self, wave_header: dict[str, Any], dim: int) -> np.ndarray:
        """
        Return the axis values for a given dimension from the wave header.

        Parameters
        ----------
        wave_header : dict[str, Any]
            The wave_header of the ibw file.
        dim : int
            The dimension to return the axis for.

        Returns
        -------
        np.ndarray
            Axis values for a given dimension.
        """
        return (
            wave_header["sfA"][dim] * np.arange(wave_header["nDim"][dim])
            + wave_header["sfB"][dim]
        )

    def axis_units_for_dim(self, wave_header: dict[str, Any], dim: int) -> str:
        """
        Return the unit for a given dimension from the wave header.

        Parameters
        ----------
        wave_header : dict[str, Any]
            The wave_header of the ibw file.
        dim : int
            The dimension to return the unit for.

        Returns
        -------
        str
            The axis units.
        """
        unit_arr = wave_header["dimUnits"][dim]

        unit = ""
        for elem in unit_arr:
            unit += elem.decode("utf-8")

        return unit


class ScientaIgorParserOld(_IgorWaveParser):
    """Parser for IBW files produced by the legacy Scienta SES software."""

    def matches_file(self, file: Path) -> bool:
        """Return True for IBW files with SES-format (non-JSON) note containing [SES]."""
        try:
            ibw = binarywave.load(file)
            note = ibw["wave"]["note"].decode("utf-8", errors="ignore")
            try:
                json.loads(note)
                return False  # PEAK format
            except (json.JSONDecodeError, ValueError):
                return "[SES]" in note
        except Exception:
            return False

    def _parse_note(self, bnote: bytes) -> dict[str, Any]:
        """
        Parse the note field of the Igor binary wave file.

        Assumes that the note field contains key-value pairs of the form
        ``key=value`` separated by newlines.

        Parameters
        ----------
        bnote : bytes
            The bytes of the binarywave note field.

        Returns
        -------
        dict[str, Any]
            The dictionary of the parsed note field.
        """
        note = bnote.decode("utf-8").replace("\r", "\n")

        notes: dict[str, Any] = {}

        for line in note.split("\n"):
            key, value, unit = _get_key_value_unit(line)
            if key:
                notes[key] = value
            if unit:
                notes[f"{key}/@units"] = unit

        return notes

    def _parse_region_metadata(
        self, region_id: int, notes: dict[str, Any]
    ) -> dict[str, Any]:
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

        region.time_stamp = _construct_date_time(region.start_date, region.time)

        region.validate_types()

        region_dict = region.dict()

        for key in unused_notes_keys:
            region_dict[key] = notes[key]

        return region_dict


class ScientaIgorParserPEAK(_IgorWaveParser):
    """Parser for IBW files produced by Scienta PEAK software."""

    def matches_file(self, file: Path) -> bool:
        """Return True for IBW files with JSON note containing Version and ImageSource."""
        try:
            ibw = binarywave.load(file)
            note = ibw["wave"]["note"].decode("utf-8", errors="ignore")
            data = json.loads(note)
            return "Version" in data and "ImageSource" in data
        except Exception:
            return False

    def _parse_note(self, bnote: bytes) -> dict[str, Any]:
        """
        Parse the note field of the Igor binary wave file.

        Assumes the note field contains a JSON string.

        Parameters
        ----------
        bnote : bytes
            The bytes of the binarywave note field.

        Returns
        -------
        dict[str, Any]
            The dictionary of the parsed note field.
        """
        note_str = bnote.decode("utf-8")

        data = json.loads(note_str)

        try:
            jsonschema.validate(instance=data, schema=scienta_igor_peak_schema)
            return data
        except jsonschema.ValidationError as err:
            raise jsonschema.ValidationError(
                f"JSON with metadata is invalid: {err.message}"
            ) from err

    def _parse_region_metadata(
        self,
        region_id: int,
        notes: dict[str, Any],
    ) -> dict[str, Any]:
        region: dict[str, Any] = {"region_id": region_id}
        region |= _flatten_dict(notes)
        region["timestamp"] = _construct_date_time(region["Date"], region["Time"])

        return region


class ScientaIgorParser(_XPSParser):
    """
    Public parser for Scienta IBW files.

    Dispatches to the appropriate implementation based on the IBW note encoding:

    - ``ScientaIgorParserPEAK`` for JSON notes (Scienta PEAK software).
    - ``ScientaIgorParserOld`` for SES-format notes (legacy ``[SES]`` sections).

    Both variants use the ``.ibw`` extension; the note content is the only
    reliable discriminator.
    """

    config_file: ClassVar[str] = "config_scienta.json"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".ibw",)

    _SUB_PARSERS: ClassVar[tuple[type[_IgorWaveParser], ...]] = (
        ScientaIgorParserPEAK,
        ScientaIgorParserOld,
    )

    def matches_file(self, file: Path) -> bool:
        """Return True if any IBW variant recognizes the file."""
        return any(cls().matches_file(file) for cls in self._SUB_PARSERS)

    def _parse(self, file: Path, **kwargs) -> None:
        """Dispatch to the first matching strategy and merge its data."""
        for sub_parser_cls in self._SUB_PARSERS:
            sub_parser = sub_parser_cls()
            if sub_parser.matches_file(file):
                sub_parser._parse(file, **kwargs)
                self._data.update(sub_parser.data)
                return
        raise ValueError(
            f"{self.__class__.__name__}: no sub parser matched file '{file}'"
        )
