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
Parser for reading XPS data from Scienta TXT exports.
"""

from pathlib import Path
from typing import ClassVar

import numpy as np
import xarray as xr

from pynxtools_xps.parsers.base import (
    ParsedSpectrum,
    VendorType,
    _construct_entry_name,
    _XPSParser,
)
from pynxtools_xps.parsers.scienta.data_model import ScientaHeader, ScientaRegion
from pynxtools_xps.parsers.scienta.metadata import (
    _check_valid_value,
    _construct_date_time,
    _get_key_value_unit,
)


class ScientaTXTParser(_XPSParser):
    """Parser for Scienta TXT exports."""

    config_file: ClassVar[str] = "config_scienta.json"
    supported_vendor: ClassVar[VendorType | None] = "scienta"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".txt",)
    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset(
        {"data", "axis_labels", "data_labels", "units"}
    )

    def matches_file(self, file: Path) -> bool:
        try:
            with open(file, encoding="utf-8", errors="ignore") as f:
                first = f.readline().strip()
                if first != "[Info]":
                    return False

                found_regions = False
                found_scienta = False

                for _ in range(200):  # scan first ~200 lines
                    line = f.readline()
                    if not line:
                        break
                    if "Number of Regions=" in line:
                        found_regions = True
                    if "scienta" in line.lower():
                        found_scienta = True
                    if found_regions and found_scienta:
                        return True

                return found_regions
        except Exception:
            return False

    def __init__(self):
        super().__init__()
        self.lines: list[str] = []
        self.header = ScientaHeader()

    def _parse(self, file: Path, **kwargs) -> None:
        """
        Parse the file's data and metadata into ``self._data``.

        Parameters
        ----------
        file : Path
            Filepath of the TXT file to be read.
        """
        self._read_lines(file)
        self._parse_header()
        self._build_parsed_spectra()

    def _build_parsed_spectra(self) -> None:
        """Iterate regions and populate ``self._data``."""
        for region_id in range(1, self.header.no_of_regions + 1):
            self._assemble_entry(region_id)

    def _read_lines(self, file: str | Path):
        """
        Read all lines from the input txt files.

        Parameters
        ----------
        file : str
            Filepath of the TXT file to be read.
        """
        with open(file, encoding="utf-8") as txt_file:
            for line in txt_file:
                self.lines += [line]

    def _parse_header(self):
        """
        Parse header with information about the software version
        and the number of spectra in the file.
        """
        n_headerlines = 4
        headerlines = self.lines[:n_headerlines]
        self.lines = self.lines[n_headerlines:]

        for line in headerlines:
            key, value, unit = _get_key_value_unit(line)
            if key:
                setattr(self.header, key, value)
        self.header.validate_types()

    def _assemble_entry(self, region_id: int) -> None:
        """
        Parse data from one region (i.e., one measured spectrum) and store as
        a ``ParsedSpectrum`` in ``self._data``.

        Parameters
        ----------
        region_id : int
            Number of the region in the file.
        """
        region = ScientaRegion(region_id=region_id)
        # Workaround as units are not defined in the ScientaRegion
        units: dict[str, str] = {}

        bool_variables = {
            "in_region": False,
            "in_region_info": False,
            "in_run_mode_info": False,
            "in_ui_info": False,
            "in_manipulator": False,
            "in_data": False,
        }

        energies: list[float] = []
        intensities: list[float] = []

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
                key, value, unit = _get_key_value_unit(line)
                setattr(region, key, value)
                if unit:
                    units[f"{key}/@units"] = unit

            if bool_variables["in_ui_info"]:
                key, value, unit = _get_key_value_unit(line)
                if _check_valid_value(value):
                    if bool_variables["in_manipulator"]:
                        key = f"manipulator_{key}"
                    setattr(region, key, value)
                    if unit:
                        units[f"{key}/@units"] = unit

            if bool_variables["in_data"]:
                try:
                    [energy, intensity] = [float(s) for s in line.split(" ") if s != ""]
                    energies.append(energy)
                    intensities.append(intensity)
                except ValueError:
                    pass

        region.time_stamp = _construct_date_time(region.start_date, region.time)
        region.validate_types()

        metadata = {**self.header.dict(), **region.dict(), **units}

        entry_name = _construct_entry_name(
            [metadata.get("spectrum_type", ""), metadata.get("region_name", "")]
        )
        if not entry_name:
            entry_name = f"region_{region_id}"

        energy_arr = np.array(energies)
        intensity_arr = np.array(intensities)

        # Shape: (cycle=1, scan=1, n_energy)
        data_arr = xr.DataArray(
            data=intensity_arr[np.newaxis, np.newaxis, :],
            coords={"energy": energy_arr},
            dims=("cycle", "scan", "energy"),
        )

        metadata = self._filter_metadata(metadata)
        metadata["energy/@units"] = "eV"
        metadata["intensity/@units"] = "counts_per_second"

        self._data[entry_name] = ParsedSpectrum(
            data=data_arr,
            raw=None,
            metadata=metadata,
        )
