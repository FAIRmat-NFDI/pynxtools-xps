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
Scienta spectrometers (.ibw or .txt format), to be passed to
MPES nxdl (NeXus Definition Language) template.
"""

import re
import copy
import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union, Optional, cast
import json
from abc import ABC, abstractmethod

import h5py
import numpy as np
import jsonschema
import xarray as xr
from igor2 import binarywave

from pynxtools_xps.reader_utils import (
    XPSMapper,
    _check_valid_value,
    _re_map_single_value,
    construct_data_key,
    construct_entry_name,
    convert_pascal_to_snake,
    _format_value,
)
from pynxtools_xps.value_mappers import (
    get_units_for_key,
    convert_units,
    convert_energy_type,
)

from pynxtools_xps.scienta.scienta_data_model import (
    ScientaHeader,
    ScientaRegion,
    scienta_igor_peak_schema,
)

from pynxtools_xps.scienta.scienta_mappings import (
    UNITS,
    VALUE_MAP,
    _get_key_value_pair,
    _construct_date_time,
)
from pynxtools_xps.value_mappers import convert_units, get_units_for_key

logger = logging.getLogger(__name__)


def _flatten_dict(
    d: Dict[str, Any], parent_key: str = "", sep: str = "/"
) -> Dict[str, Any]:
    """
    Flattens a nested dictionary into a single level with keys representing the hierarchy.

    Args:
        d (Dict[str, Any]): The dictionary to flatten.
        parent_key (str): The base key to prepend (used for recursion).
        sep (str): The separator to use for flattened keys.

    Returns:
        Dict[str, Any]: The flattened dictionary.
    """
    items: List[Tuple[str, Any]] = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(_flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


class MapperScienta(XPSMapper):
    """
    Class for restructuring data from
    Scienta spectrometers into a structured python
    dictionaries.
    """

    config_file = {
        ".h5": "config_scienta_hdf5.json",
        ".hdf5": "config_scienta_hdf5.json",
        ".ibw": "config_scienta.json",
        ".txt": "config_scienta.json",
    }

    __prmt_file_ext__ = [
        ".h5",
        ".hdf5",
        ".ibw",
        ".txt",
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
            try:
                with open(str(self.file), "rb") as f:
                    data_ibw = binarywave.load(f)
                check_note = json.loads(data_ibw["wave"]["note"].decode("utf-8"))[
                    "Version"
                ]
                return ScientaIgorParserPEAK()
            except json.JSONDecodeError:
                return ScientaIgorParserOld()
            return ScientaIgorParser()
        elif str(self.file).endswith((".h5", ".hdf5")):
            return ScientaHdf5Parser()
        raise ValueError(MapperScienta.__file_err_msg__)

    def construct_data(self):
        """Map Parser data to NXmpes-ready dict."""
        # pylint: disable=duplicate-code

        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"] = cast(Dict[str, Any], {})

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: Dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """

        entry_parts = []
        for part in [
            "spectrum_type",
            "region_name",
            "Name",
            "acquisition/spectrum/name",
            "title",
        ]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            if isinstance(value, dict):
                continue
            if key.startswith("entry"):
                key = key.replace("entry/", "", 1)
            mpes_key = f"{entry_parent}/{key}"
            self._xps_dict[mpes_key] = value
            units = convert_units(get_units_for_key(key, UNITS))
            if units is not None:
                self._xps_dict[f"{mpes_key}/@units"] = units

        try:
            self._fill_with_data_txt_ibw(spectrum, entry, entry_parent)
        except (IndexError, KeyError):
            self._fill_with_data_hdf5(spectrum, entry, entry_parent)

    def _fill_with_data_txt_ibw(
        self, spectrum: Dict[str, Any], entry: str, entry_parent: str
    ):
        # If multiple spectra exist to entry, only create a new
        # xr.Dataset if the entry occurs for the first time.
        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        axis_units = spectrum.get("units")
        if axis_units:
            for axis_name, unit in axis_units.items():
                unit_key = f"{entry_parent}/{axis_name}/@units"
                self._xps_dict[unit_key] = convert_units(unit)

        # Create key for writing to data
        scan_key = construct_data_key(spectrum)

        axes = {
            key: value
            for key, value in spectrum["data"].items()
            if key in spectrum["axis_labels"]
        }
        intensities = np.array(
            [
                value
                for key, value in spectrum["data"].items()
                if key in spectrum["data_labels"] and "reduced" not in key
            ]
        ).squeeze(axis=0)

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
            averaged_scans = intensities

        self._xps_dict["data"][entry][scan_key.split("_")[0]] = xr.DataArray(
            data=averaged_scans, coords=axes, dims=list(axes.keys())
        )

        # Write scan data to 'data'.
        self._xps_dict["data"][entry][scan_key] = xr.DataArray(
            data=intensities, coords=axes
        )

        # Write channel data to 'data'.
        channel_key = f"{scan_key}_chan0"
        self._xps_dict["data"][entry][channel_key] = xr.DataArray(
            data=intensities, coords=axes
        )

    def _fill_with_data_hdf5(
        self, spectrum: Dict[str, Any], entry: str, entry_parent: str
    ):
        self._xps_dict["data"][entry] = {}

        data_keys = [
            "data/data",
            "data/sum",
            "data/x_axis",
            "data/y_axis",
            "data_reduced_1d/data",
            "data_reduced_1d/x_axis",
        ]

        for key in data_keys:
            value = spectrum.get(f"acquisition/spectrum/{key}")
            if value is not None:
                self._xps_dict["data"][entry][key] = value


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

        # Convert date and time to ISO8601 date time.
        region.time_stamp = _construct_date_time(region.start_date, region.time)

        region.validate_types()

        region_dict = {**self.header.dict(), **region.dict()}
        region_dict["data"] = {
            "energy": np.array(energies),
            "intensity": np.array(intensities),
        }
        region_dict["energy/@units"] = "eV"
        region_dict["intensity/@units"] = "counts_per_second"

        region_dict["axis_labels"] = ["energy"]
        region_dict["data_labels"] = ["intensity"]

        self.spectra.append(region_dict)


class ScientaIgorParser(ABC):
    """Parser for Scienta IBW exports."""

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

        axes_labels_with_units = self._parse_unit(wave["dimension_units"])
        data_labels_with_units = self._parse_unit(wave["data_units"])

        wave_header = wave["wave_header"]
        data = wave["wData"]

        # TODO: Add support for formulas if they are written by the
        # measurement software.
        # formula = wave["formula"]

        notes: Dict[str, Any] = self._parse_note(wave["note"])

        self.no_of_regions = len(data.shape)

        spectrum: Dict[str, Any] = {
            "data": cast(Dict[str, Any], {}),
            "axis_labels": cast(List[str], []),
            "data_labels": cast(List[str], []),
            "units": cast(Dict[str, Any], {}),
        }

        for i, (dim, unit) in enumerate(axes_labels_with_units):
            if dim in ("kinetic_energy", "binding_energy", "analyser_energy"):
                spectrum["energy_scale"] = convert_energy_type(dim)
                dim = "energy"
            spectrum["data"][dim] = self.axis_for_dim(wave_header, dim=i)
            spectrum["axis_labels"].append(dim)
            spectrum["units"][dim] = convert_units(unit)

        for label, unit in data_labels_with_units:
            spectrum["data"][label] = data
            spectrum["data_labels"].append(label)
            spectrum["units"][label] = convert_units(unit)
            # Convert date and time to ISO8601 date time.

        spectrum["igor_binary_wave_format_version"] = ibw_version

        region_metadata = self._parse_region_metadata(region_id=0, notes=notes)
        spectrum |= region_metadata

        self.spectra.append(spectrum)

        return self.spectra

    @abstractmethod
    def _parse_note(self, bnote: bytes) -> Dict[str, Any]:
        pass

    @abstractmethod
    def _parse_region_metadata(
        self, region_id: int, notes: Dict[str, Any]
    ) -> Dict[str, Any]:
        return {}

    def _parse_unit(self, bunit: bytes) -> List[Tuple[str, Optional[str]]]:
        """
        Extracts labels and units from a string containing one or more label-unit pairs.
        If no unit is present, it returns just the label.

        Parameters
        ----------
        bunit: bytes
            The input string containing the label and unit(s).

        Returns
        -------
        list of tuple
            A list of tuples, each containing:
            - unit_label : str
                The extracted label.
            - unit : str or None
                The extracted unit (None if no unit is present).
        """
        # Decode the bytes to a string
        unit = bunit.decode("utf-8").replace("\r", "\n")

        # Regex to match "Label [Unit]" patterns
        pattern = r"([\w\s]+?)\s*\[([\w\s.]+)\]"

        if matches := re.findall(pattern, unit):
            # Process each match into a tuple of (label, unit)
            return [
                (convert_pascal_to_snake(label.strip()), unit.strip())
                for label, unit in matches
            ]

        # If no matches, return the entire string as a label with no unit
        return [(convert_pascal_to_snake(unit.strip()), None)]

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


class ScientaIgorParserOld(ScientaIgorParser):
    """Parser version for the old Scienta exporter (i.e., not the one used by the PEAK software)."""

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

    def _parse_region_metadata(
        self, region_id: int, notes: Dict[str, Any]
    ) -> Dict[str, Any]:
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

        # Convert date and time to ISO8601 date time.
        region.time_stamp = _construct_date_time(region.start_date, region.time)

        region.validate_types()

        region_dict = region.dict()

        for key in unused_notes_keys:
            region_dict[key] = notes[key]

        return region_dict


class ScientaIgorParserPEAK(ScientaIgorParser):
    """Parser version for data exported by Scienta's PEAK software."""

    def _parse_note(self, bnote: bytes) -> Dict[str, Any]:
        """
        Parses the note field of the igor binarywave file.

        This is the _parse_note version for the

        It assumes that the note field contains a JSON string.

        Parameters
        ----------
        bnote : bytes
            The bytes of the binarywave note field.

        Returns
        -------
        Dict[str, Any]
            The dictionary of the parsed note field.
        """
        # Decode the byte string to UTF-8
        note_str = bnote.decode("utf-8")

        data = json.loads(note_str)

        try:
            # Validate against the defined schema
            jsonschema.validate(instance=data, schema=scienta_igor_peak_schema)
            return data
        except jsonschema.ValidationError as err:
            raise jsonschema.ValidationError(
                f"JSON with metadata is invalid: {err.message}"
            ) from err

    def _parse_region_metadata(
        self,
        region_id: int,
        notes: Dict[str, Any],
    ) -> Dict[str, Any]:
        region: Dict[str, Any] = {"region_id": region_id}
        region |= _flatten_dict(notes)
        region["timestamp"] = _construct_date_time(region["Date"], region["Time"])

        return region


class ScientaHdf5Parser:
    def __init__(self):
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

        def format_value(key: str, value_str: str) -> Any:
            """
            Formats a value string (to a corresponding key) according to a series of transformations.
            This function:
            1. Formats the numeric part of the value according to its expected type.
            2. Remaps the value to a new format if specified in `VALUE_MAP`.
            Args:
                key (str): The key associated with the value, which may need mapping and formatting.
                value_str (str): The value string to format and separate into numeric value and unit.
            Returns:
                Tuple[Any, str]:
                    - The formatted key (converted to snake_case and remapped if needed).
                    - The formatted value, with numeric value processed and remapped according to `VALUE_MAP`.
            """
            kwargs: Dict[str, Any] = {}

            if key.endswith(("start_time", "stop_time")):
                kwargs["possible_date_formats"] = [
                    "%Y-%m-%dT%H:%M:%S",
                    "%Y-%m-%dT%H:%M:%S.%f%z",
                    "%Y-%m-%dT%H:%M:%S%z",
                ]
            value = _re_map_single_value(key, value_str, VALUE_MAP, **kwargs)
            value = _format_value(value)

            return value

        def recursively_read_group(group, path=""):
            result = {}
            for key, item in group.items():
                new_path = f"{path}/{key}" if path else key
                if isinstance(item, h5py.Group):
                    # Recursively read subgroups
                    result.update(recursively_read_group(item, new_path))
                elif isinstance(item, h5py.Dataset):
                    # Read datasets
                    data = item[()]
                    if isinstance(data, bytes):
                        data = data.decode("utf-8")
                    data = format_value(key, data)
                    data = format_value(new_path, data)
                    result[new_path] = data
            return result

        # Open the HDF5 file and read its contents
        with h5py.File(file, "r") as hdf:
            hdf5_data = recursively_read_group(hdf)

        try:
            length, width = (
                hdf5_data["instrument/analyser/slit/length"],
                hdf5_data["instrument/analyser/slit/width"],
            )
            hdf5_data["instrument/analyser/slit/size"] = np.array([length, width])
        except KeyError:
            pass

        # Find all axes
        pattern = re.compile(r"^acquisition/spectrum/data/(?!x_axis$)([^/]+_axis)$")
        additional_axes = list(
            reversed(
                [match.group(1) for key in hdf5_data if (match := pattern.match(key))]
            )
        )

        if additional_axes:
            hdf5_data["data_axes"] = additional_axes + ["energy"]
        else:
            hdf5_data["data_axes"] = ["energy"]

        return [hdf5_data]
