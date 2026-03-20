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
Parser for reading XPS data from Scienta HDF5 files.
"""

import re
from pathlib import Path
from typing import Any, ClassVar

import h5py
import numpy as np
import xarray as xr

from pynxtools_xps.parsers.base import (
    ParsedSpectrum,
    VendorType,
    _construct_entry_name,
    _XPSParser,
)


class ScientaHDF5Parser(_XPSParser):
    """Parser for Scienta Omicron HDF5 exports."""

    config_file: ClassVar[str] = "config_scienta_hdf5.json"
    supported_vendor: ClassVar[VendorType | None] = "scienta"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".h5", ".hdf5")
    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset(
        {"acquisition/spectrum/data/data", "acquisition/spectrum/data/x_axis"}
    )
    _HDF5_MAGIC: ClassVar[bytes] = b"\x89HDF\r\n\x1a\n"

    def matches_file(self, file: Path) -> bool:
        """Return True for Scienta Omicron HDF5 files (vendor root attribute)."""
        try:
            with open(file, "rb") as f:
                if f.read(8) != self._HDF5_MAGIC:
                    return False
            with h5py.File(file, "r") as hdf:
                return hdf.attrs.get("vendor") == b"Scienta Omicron"
        except Exception:
            return False

    def _filter_metadata(self, raw: dict[str, Any]) -> dict[str, Any]:
        """Exclude data-array keys and any numpy array values from metadata."""
        return {
            k: v
            for k, v in raw.items()
            if k not in self._metadata_exclude_keys and not isinstance(v, np.ndarray)
        }

    def _parse(self, file: Path, **kwargs) -> None:
        """
        Read the HDF5 file into ``self._hdf5_data`` and delegate to
        ``_build_parsed_spectra``.

        Parameters
        ----------
        file : Path
            Filepath of the HDF5 file to be read.
        """

        # TODO: this shall be removed or refactored into context
        def format_value(key: str, value_str: str) -> Any:
            """
            Format a value string according to a series of transformations.

            Args:
                key (str): The key associated with the value.
                value_str (str): The value string to format.

            Returns:
                Any: The formatted value.
            """
            if key.endswith(("start_time", "stop_time")):
                pass  # reserved for future date-format normalization
            # NOTE: _format_value is intentionally not called here; the original
            # code referenced an undefined local _format_value which was a bug.
            # Returning the value as-is preserves the original behavior.
            return value_str

        def recursively_read_group(group, path=""):
            result = {}
            for key, item in group.items():
                new_path = f"{path}/{key}" if path else key
                if isinstance(item, h5py.Group):
                    result.update(recursively_read_group(item, new_path))
                elif isinstance(item, h5py.Dataset):
                    data = item[()]
                    if isinstance(data, bytes):
                        data = data.decode("utf-8")
                    data = format_value(key, data)
                    data = format_value(new_path, data)
                    result[new_path] = data
            return result

        with h5py.File(file, "r") as hdf:
            self._hdf5_data = recursively_read_group(hdf)

        self._build_parsed_spectra()

    def _build_parsed_spectra(self) -> None:
        """
        Build one ``ParsedSpectrum`` from ``self._hdf5_data`` and store it in
        ``self._data``.

        The HDF5 layout uses ``acquisition/spectrum/data/data`` for the
        multi-dimensional intensity array and ``acquisition/spectrum/data/x_axis``
        (plus optional ``*_axis`` siblings) for the axis values.  A single
        ``ParsedSpectrum`` is produced per file with dimensions
        ``("cycle", "scan", *data_axes)``.
        """
        hdf5_data = self._hdf5_data

        try:
            length, width = (
                hdf5_data["instrument/analyser/slit/length"],
                hdf5_data["instrument/analyser/slit/width"],
            )
            hdf5_data["instrument/analyser/slit/size"] = np.array([length, width])
        except KeyError:
            pass

        # Discover additional axes (everything matching *_axis except x_axis)
        pattern = re.compile(r"^acquisition/spectrum/data/(?!x_axis$)([^/]+_axis)$")
        additional_axes = list(
            reversed(
                [match.group(1) for key in hdf5_data if (match := pattern.match(key))]
            )
        )

        data_axes = (additional_axes + ["energy"]) if additional_axes else ["energy"]

        # Retrieve the intensity array
        raw_data = hdf5_data.get("acquisition/spectrum/data/data")
        x_axis = hdf5_data.get("acquisition/spectrum/data/x_axis")

        if raw_data is not None:
            intensity_arr = np.asarray(raw_data)
        else:
            intensity_arr = np.array([])

        if x_axis is not None:
            energy_arr = np.asarray(x_axis)
        else:
            energy_arr = (
                np.arange(intensity_arr.shape[-1])
                if intensity_arr.size
                else np.array([])
            )

        # Ensure leading (cycle, scan) dims by prepending two size-1 axes
        intensity_arr = intensity_arr[np.newaxis, np.newaxis, ...]

        # Build dims: ("cycle", "scan", *data_axes)
        n_physical = intensity_arr.ndim - 2  # subtract cycle and scan
        physical_dims = (
            data_axes[-n_physical:]
            if n_physical <= len(data_axes)
            else [f"axis_{i}" for i in range(n_physical)]
        )
        dims = ("cycle", "scan") + tuple(physical_dims)

        # Build coords: only energy (x_axis) is reliably available
        coords: dict[str, Any] = {}
        if "energy" in physical_dims:
            coords["energy"] = energy_arr

        data_arr = xr.DataArray(
            data=intensity_arr,
            coords=coords,
            dims=dims,
        )

        # Entry name from HDF5 metadata
        name = hdf5_data.get("acquisition/spectrum/name", "")
        title = hdf5_data.get("title", "")
        entry_name = _construct_entry_name([p for p in [name, title] if p])
        if not entry_name:
            entry_name = "entry"

        metadata = self._filter_metadata(hdf5_data)
        metadata["data_axes"] = data_axes

        self._data[entry_name] = ParsedSpectrum(
            data=data_arr,
            raw=None,
            metadata=metadata,
        )
