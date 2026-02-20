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
Phi PHI VersaProbe 4 instruments (.spe or .pro format) and mapping
to the NXmpes/NXxps template.
"""

import struct
import warnings
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import xarray as xr

from pynxtools_xps.logging import _logger
from pynxtools_xps.numerics import safe_arange_with_edges
from pynxtools_xps.parsers.base import _construct_entry_name, _XPSMapper, _XPSParser
from pynxtools_xps.parsers.phi.data_model import (
    PHIMetadata,
    PHISpatialArea,
    PHISpectralRegion,
)
from pynxtools_xps.parsers.phi.metadata import _context, _convert_channel_info


class PHIMapper(_XPSMapper):
    """
    Map PHI format to NXmpes-ready dict.
    """

    config_file: ClassVar[str] = "config_phi.json"

    def __init__(self):
        super().__init__()
        self.write_channels_to_data = True

    def _select_parser(self):
        """Select the proper Phi data parser."""
        return PHIParser()

    def construct_data(self, raw_data: list[dict[str, Any]]):
        """Map Phi format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code

        for spectrum in raw_data:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []

        for part in ["group_name", "spectrum_type"]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = _construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            if key.startswith("entry"):
                entry_parent = "/ENTRY[entry]"
                key = key.replace("entry/", "", 1)
            mpes_key = f"{entry_parent}/{key}"
            self._data[mpes_key] = value

            # units = get_units_for_key(key, UNITS)
            # if units is not None:
            #     self._data[f"{mpes_key}/@units"] = units

        # Create key for writing to data
        cycle_key = "cycle0"

        energy = np.array(spectrum["energy"])

        if entry not in self._data["data"]:
            self._data["data"][entry] = xr.Dataset()

        # Write averaged data to 'data'.
        all_scan_data = np.array(
            [np.array(value) for key, value in spectrum["data"].items()]
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            averaged_scans = np.mean(all_scan_data, axis=0)

        self._data["data"][entry][cycle_key] = xr.DataArray(
            data=averaged_scans,
            coords={"energy": energy},
        )

        # Write scan data to 'data'.
        for scan_no, intensity in spectrum["data"].items():
            xarr_key = f"{cycle_key}_{scan_no}"
            self._data["data"][entry][xarr_key] = xr.DataArray(
                data=intensity,
                coords={"energy": energy},
            )


class PHIParser(_XPSParser):  # pylint: disable=too-few-public-methods
    """
    A parser for reading in PHI VersaProbe 4 data in the .spe or
    .pro format.
    Tested with Software version SS 3.3.3.2.
    """

    def __init__(self):
        """
        Construct the parser.

        """
        self.metadata = PHIMetadata()

        self.binary_header_length = 4
        self._data_header_length = 24
        self.encoding: tuple[str, int] = ("<f", 4)

        self.binary_header: np.ndarray = None
        self._data_header: np.ndarray = None

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the .spe, .pro file into a list of dictionaries.

        Parsed data is stored in the attribute 'self.data'.
        Each dictionary in the data list is a grouping of related
        attributes. The dictionaries are later re-structured into a
        nested dictionary that more closely resembles the domain logic.

        Parameters
        ----------
        file : str
            XPS data filepath.

        Returns
        -------
        list
            Flat list of dictionaries containing one spectrum each.

        """
        self.raw_data = self._read_lines(file)
        header, data = self._separate_header_and_data()
        # l = b'\x01'.join(data)

        self.parse_header_into_metadata(header)
        regions = self.parse_spectral_regions(header)
        areas = self.parse_spatial_areas(header)

        self.add_regions_and_areas_to_spectra(regions, areas)

        self.parse_binary_header(data)
        self.parse_data_into_spectra(data)

        self.add_metadata_to_each_spectrum()

        return self._data

    def _read_lines(self, file: str | Path):
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
        with open(file, mode="rb") as txt_file:  # "utf-8"
            data = txt_file.read()

        return data

    def _separate_header_and_data(self):
        """
        Separate header (with metadata) from data.

        Returns
        -------
        None.

        """
        try:  # Windows
            header, data = self.raw_data.split(b"EOFH\r\n")
            header = header.decode("ascii").split("\r")
            header = [line.strip() for line in header if line.strip()]

        except ValueError:  # Linux, MacOS
            header, data = self.raw_data.split(b"EOFH\n")
            header = header.decode("ascii").split("\n")
            header = [line.strip() for line in header if line.strip()]

        return header, data

    def parse_header_into_metadata(self, header: list[str]):
        """
        Parse header into PHiMetadata dataclass.

        Parameters
        ----------
        header : str
            Header data for one spectrum as a String.

        """

        def map_keys(key: str, channel_count: int):
            """Maps key names for better identification of fields."""
            prefix_map = {
                "neut": ("neut_", "neutral_"),
                "prof": ("prof_", "profiling_"),
                "x_ray": ("x_ray", "xray"),
                "egun_neut": ("egun_neut", "flood_gun"),
            }

            for prefix, (initial, replacement) in prefix_map.items():
                if key.startswith(prefix):
                    key = key.replace(initial, replacement)

            if key.startswith("sxi_lens"):
                key += "_voltage"

            if key.startswith("channel _info"):
                key = f"channel_{channel_count}_info"
                channel_count += 1

            key_map: dict[str, str] = {
                "desc": "description",
                "defect_pos": "defect_positioner",
                "g_c_i_b": "gcib",
            }

            for old_key, new_key in key_map.items():
                if old_key in key:
                    key = key.replace(old_key, new_key)

            key = _context.normalize_key(key)

            return key, channel_count

        channel_count = 1
        datacls_fields = list(self.metadata.__dataclass_fields__.keys())
        datacls_fields = [field for field in datacls_fields if "_units" not in field]

        for line in header:
            raw_value: str
            try:
                key, raw_value = line.split(": ")

            except ValueError:
                key = line.strip(": ")
                raw_value = ""

            # TODO: use full _context
            key = _context.normalize_key(key)
            key, channel_count = map_keys(key, channel_count)

            if key in datacls_fields:
                field_type = type(getattr(self.metadata, key))

                # TODO: use _context
                # if key in KEYS_WITH_UNITS:
                #     value, unit = self.extract_unit(key, value)
                #     setattr(self.metadata, f"{key}_units", unit)

                if not key.startswith("channel_"):
                    value = _context.map_value(key, raw_value)

                else:
                    value = _convert_channel_info(raw_value)

                setattr(self.metadata, key, value)

        self.metadata.validate_types()

    def parse_spectral_regions(self, header: list[str]):
        """
        Parse spectral regions definitions.

        Parameters
        ----------
        header : list
            List of header strings.

        Returns
        -------
        regions_full : list
            List of regions with `full` keyword.
        regions : list
            List of regions without `full` keyword.

        """
        spectral_defs = [line for line in header if line.startswith("SpectralReg")]

        regions = []
        regions_full = []
        for n in range(self.metadata.no_spectral_regions_full):
            regions_full += [PHISpectralRegion(region_id=n + 1)]
        for n in range(self.metadata.no_spectral_regions):
            regions += [PHISpectralRegion(region_id=n + 1)]

        concept_map = {
            "SpectralRegDef": "region_definition",
            "SpectralRegDef2": "region_definition2",
            "SpectralRegBackground": "region_background",
            "SpectralRegHero": "region_hero",
            "SpectralRegIR": "region_ir",
        }

        for region in regions_full:
            region.full_region = True
            for file_key, metadata_key in concept_map.items():
                if spectral_defs and spectral_defs[0].startswith(file_key):
                    setattr(region, metadata_key, spectral_defs[0].split(" ", 1)[1])
                    spectral_defs.pop(0)

        for region in regions:
            region.full_region = False

            for file_key, metadata_key in concept_map.items():
                if spectral_defs and spectral_defs[0].startswith(file_key):
                    setattr(region, metadata_key, spectral_defs[0].split(" ", 1)[1])
                    spectral_defs.pop(0)

        for region in regions:
            def_split = region.region_definition.split(" ")
            region.spectrum_type = def_split[2]
            region.n_values = int(def_split[4])
            step = -float(def_split[5])
            start = float(def_split[6])
            stop = float(def_split[7])
            region.dwell_time = float(def_split[-3])
            region.dwell_time_units = "s"
            region.pass_energy = float(def_split[-2])
            region.energy_type = "binding"
            region.energy_units = "eV"

            region.energy = np.flip(safe_arange_with_edges(stop, start, step))

            region.validate_types()

        return regions

    def parse_spatial_areas(self, header: list[str]):
        """
        Parse spatial areas definitions.

        Parameters
        ----------
        header : list
            List of header strings.

        Returns
        -------
        areas: list
            List of spatial areas and their definitions.

        """
        spatial_defs = [line for line in header if line.startswith("Spatial")]

        areas = []
        for n in range(self.metadata.no_spatial_areas):
            areas += [PHISpatialArea(area_id=n)]

        concept_map = {
            "SpatialAreaDef": "area_definition",
            "SpatialAreaDesc": "area_description",
            "SpatialHRPhotoCor": "area_hr_photo_correction",
        }

        for area in areas:
            for file_key, metadata_key in concept_map.items():
                if spatial_defs and spatial_defs[0].startswith(file_key):
                    setattr(area, metadata_key, spatial_defs[0].split(" ", 1)[1])
                    spatial_defs.pop(0)

        return areas

    def add_regions_and_areas_to_spectra(self, regions: list, areas: list):
        """
        Define each spectra by its region and area defintions.

        Parameters
        ----------
        regions : list
            List of PHISpectralRegion objects.
        areas : list
            List of PHISpatialArea objects

        Returns
        -------
        None.

        """
        for region in regions:
            for area in areas:
                concatenated = {**region.dict(), **area.dict()}

                region_and_areas: dict[str, Any] = {}

                for key, value in concatenated.items():
                    replacement_map = {
                        "_units": "/@units",
                        "energy_type": "energy/@type",
                    }

                    new_key = key

                    for suffix, replacement in replacement_map.items():
                        if suffix in key:
                            new_key = key.replace(suffix, replacement)

                    region_and_areas[new_key] = value

                self._data += [region_and_areas]

    def parse_binary_header(self, binary_data):
        """
        Read the binary headers
        Assuming the headers are 4 bytes unsigned integers
        Each spectrum gets 24 unsigned 4 byte integers

        Parameters
        ----------
        binary_data : bytes
            Binary XPS data, format is 64 bit float.

        """
        binary_header = struct.unpack("I", binary_data[: self.binary_header_length])[0]

        for i, spectrum in enumerate(self._data):
            start = (
                self.binary_header_length * self._data_header_length * i
                + self.binary_header_length
            )
            stop = start + self.binary_header_length * self._data_header_length
            spectrum_header = struct.unpack(
                "I" * self._data_header_length, binary_data[start:stop]
            )

            spectrum.update(
                {
                    "binary_header": binary_header,
                    "spectrum_header": np.array(spectrum_header),
                    "n_values": spectrum_header[8],
                    "n_scans": spectrum_header[9],
                    "binary_start": spectrum_header[23],
                    "binary_len": spectrum_header[22],
                }
            )

        self._check_encoding()

    def parse_data_into_spectra(self, binary_data):
        """
        Parse the data of all spectra.

        Parameters
        ----------
        binary_data : bytes
            Binary XPS data, format is 64 bit float.

        """
        for i, spectrum in enumerate(self._data):
            float_buffer = self.encoding[1]

            # spectrum["n_scans"] = 2
            spectrum["data"] = {}
            start = spectrum["binary_start"]
            stop = start + spectrum["binary_len"]

            binary_spectrum_data = binary_data[start:stop]
            n_values_binary = spectrum["n_values"] * float_buffer

            for scan_no in range(spectrum["n_scans"]):
                spec_start = scan_no * n_values_binary
                spec_stop = spec_start + n_values_binary
                scan_data = binary_spectrum_data[spec_start:spec_stop]

                spectrum["data"][f"scan_{scan_no}"] = self._parse_binary_scan(scan_data)

    def _parse_binary_scan(self, binary_scan_data):
        """
        For each spectrum scan, parse the XPS data by unpacking the
        64 bit floats.

        Parameters
        ----------
        binary_spectrum_data : bytes
            Binary data containing the intensity data for one
            spectrum.

        Returns
        -------
        parsed_data : np.ndarray
            One-dimensional array with spectrum intensities.

        """
        float_buffer = self.encoding[1]
        encoding = self.encoding[0]
        parsed_data = []

        n_values = int(len(binary_scan_data) / float_buffer)

        for i in range(n_values):
            start = i * float_buffer
            stop = start + float_buffer
            parsed_data.append(
                struct.unpack_from(encoding, binary_scan_data[start:stop])[0]
            )

        return parsed_data

    def _check_encoding(self):
        """Check if the binary data is single or double encoded."""
        datasize = sum(
            s["spectrum_header"][8] * s["spectrum_header"][9] for s in self._data
        )
        binary_size = sum(s["spectrum_header"][-2] for s in self._data)

        encodings_map: dict[str, tuple[str, int]] = {
            "double": ("<d", 8),
            "float": ("<f", 4),
        }
        if binary_size / datasize == 4:
            self.encoding = encodings_map["float"]
        elif binary_size / datasize == 8:
            self.encoding = encodings_map["double"]
        else:
            _logger.error("This binary encoding is not supported.")

    # TODO: This should be done through the _context
    # def extract_unit(self, key: str, value):
    #     """
    #     Extract units for the metadata containing unit information.

    #     Example:
    #         analyzer_work_function: 4.506 eV
    #         -> analyzer_work_function: 4.506,
    #            analyzer_work_function_units: eV,

    #     Parameters
    #     ----------
    #     key : str
    #         Key of the associated value.
    #     value : str
    #         Combined unit and value information.

    #     Returns
    #     -------
    #     value :
    #         value with units.
    #     unit : str
    #         Associated unit.

    #     """
    #     try:
    #         value, unit = split_value_and_unit(value)
    #     except ValueError:
    #         unit = ""

    #     unit = convert_units(unit)

    #     if key in UNIT_MISSING:
    #         unit = UNIT_MISSING[key]

    #     return value, unit

    def add_metadata_to_each_spectrum(self):
        """
        Add flattened dict with PHIMetadata fields to each spectrum.

        """
        flattened_metadata = self.flatten_metadata()

        for spectrum in self.data:
            spectrum.update(flattened_metadata)
            spectrum["intensity/@units"] = UNIT_MISSING["intensity"]

    def flatten_metadata(self):
        """
        Flatten metadata dict so that key-value pairs of nested
        dictionaries are at the top level.

        Parameters
        ----------
        metadata_dict : dict
            Metadata dict with PHIMetadata fields as keys.

        Returns
        -------
        flattened_dict : dict
            Flatted metadata_dict without any nested substructure.

        """

        def shorten_sup_key(sup_key: str):
            """Shorted the key for some nested dicts."""
            shortened_key_map = {
                "xray_source": "xray",
                "xray_settings": "xray",
                "stage_positions": "stage",
            }
            if sup_key in shortened_key_map:
                return shortened_key_map[sup_key]
            return sup_key

        flattened_dict = {}

        for key, value in self.metadata.dict().items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    sup_key = shorten_sup_key(key)
                    flattened_dict[f"{sup_key}_{sub_key}"] = sub_value
            else:
                flattened_dict[key] = value

        for key in flattened_dict.copy().keys():
            if "_units" in key:
                new_key = key.replace("_units", "/@units")
                flattened_dict[new_key] = flattened_dict.pop(key)

        for key in flattened_dict.copy().keys():
            if key in KEYS_WITH_UNITS and f"{key}/@units" not in flattened_dict:
                flattened_dict[f"{key}/@units"] = UNIT_MISSING[key]

        return flattened_dict
