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
A generic reader for loading XPS (X-ray Photoelectron Spectroscopy) data
file into mpes nxdl (NeXus Definition Language) template.
"""

import copy
import datetime
import re
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any, cast

import numpy as np
import xarray as xr
from pynxtools.dataconverter.helpers import extract_atom_types
from pynxtools.dataconverter.readers.multi.reader import MultiFormatReader
from pynxtools.dataconverter.readers.utils import parse_yml
from pynxtools.dataconverter.template import Template

from pynxtools_xps.logging import _logger
from pynxtools_xps.numerics import check_units
from pynxtools_xps.parsers import (
    PHIParser,
    ScientaHDF5Parser,
    ScientaIgorParser,
    ScientaTXTParser,
    SpecsSLEParser,
    SpecsXMLParser,
    SpecsXYParser,
    VamasExportParser,
    VamasParser,
    VamasResultParser,
)
from pynxtools_xps.parsers.base import ParsedSpectrum, _XPSMetadataParser, _XPSParser

_PROCESS_ORDER: list[tuple[str, Any]] = [
    (
        "channels_averaging",
        lambda s: s.raw is not None and s.raw.sizes.get("channel", 1) > 1,
    ),
    ("scan_averaging", lambda s: "scan" in s.data.dims and s.data.sizes["scan"] > 1),
    ("cycle_averaging", lambda s: "cycle" in s.data.dims and s.data.sizes["cycle"] > 1),
]

_CONVERT_DICT = {
    "unit": "@units",
    "version": "@version",
    "user": "USER[user]",
    "instrument": "INSTRUMENT[instrument]",
    "source_xray": "sourceTYPE[source_xray]",
    "beam_xray": "beamTYPE[beam_xray]",
    "electronanalyzer": "ELECTRONANALYZER[electronanalyzer]",
    "collectioncolumn": "COLLECTIONCOLUMN[collectioncolumn]",
    "energydispersion": "ENERGYDISPERSION[energydispersion]",
    "detector": "ELECTRON_DETECTOR[detector]",
    "manipulator": "MANIPULATOR[manipulator]",
    "pid": "PID_CONTROLLER[pid_controller]",
    "process": "PROCESS[process]",
    "sample": "SAMPLE[sample]",
    "substance": "SUBSTANCE[substance]",
}

_REPLACE_NESTED: dict[str, str] = {}


def _concatenate_values(value1: Any, value2: Any) -> Any:
    """
    Concatenate two values of same type to be stored
    in xps_data_dict. Dicts are merged and every other object is
    appended to a list.
    """
    if isinstance(value1, dict) and isinstance(value2, dict):
        concatenated = {**value1, **value2}
    else:
        if not isinstance(value1, list):
            value1 = [value1]
        if not isinstance(value2, list):
            value2 = [value2]
        concatenated = value1 + value2

    return concatenated


# ---------------------------------------------------------------------------
# Module-level helpers for get_data dispatch
# ---------------------------------------------------------------------------


def _raw_or_data_array(spectrum: "ParsedSpectrum") -> xr.DataArray:
    """Return raw DataArray, falling back to processed data."""
    data_array = spectrum.raw if spectrum.raw is not None else spectrum.data
    if data_array is not None:
        return data_array
    else:
        raise ValueError(f"Could not retrieve raw or data array for {spectrum}.")


def _raw_active_dims(spectrum: "ParsedSpectrum") -> list[str]:
    """Dims with size > 1 or named 'energy' in the raw (or processed) array."""
    da = _raw_or_data_array(spectrum)
    return [str(d) for d in da.dims if d == "energy" or da.sizes[d] > 1]


def _raw_dim_indices(spectrum: "ParsedSpectrum", dim: str) -> np.ndarray | None:
    """Index array for *dim* in raw (or processed) DA; None if absent or singular."""
    da = _raw_or_data_array(spectrum)
    if dim not in da.dims or da.sizes[dim] <= 1:
        return None
    return np.arange(da.sizes[dim])


def _get_raw(spectrum: "ParsedSpectrum") -> np.ndarray:
    """Return raw (or processed) data, squeezing singular non-energy dims."""
    da = _raw_or_data_array(spectrum)
    squeeze = [d for d in da.dims if d != "energy" and da.sizes[d] <= 1]
    return da.squeeze(squeeze).values if squeeze else da.values


def _cycle_averaging(spectrum: "ParsedSpectrum") -> np.ndarray | None:
    """Average over all scan and cycle dims; None if only one cycle present."""
    if spectrum.data is None:
        return None
    if "cycle" not in spectrum.data.dims or spectrum.data.sizes["cycle"] <= 1:
        return None
    dims = [d for d in ["scan", "cycle"] if d in spectrum.data.dims]
    return spectrum.data.mean(dims).values


def _get_raw_energy_index(spectrum: "ParsedSpectrum") -> int:
    """Position of 'energy' in the active raw dims; last dim if not found."""
    dims = _raw_active_dims(spectrum)
    return dims.index("energy") if "energy" in dims else len(dims) - 1


def _stats_reduced(
    spectrum: "ParsedSpectrum", stats_fn: Callable[[], xr.DataArray]
) -> np.ndarray | None:
    """Apply *stats_fn* and reduce extra (non-energy) dims by summation.

    Returns None when the result is already 1D (i.e. identical to the
    unreduced variant), signalling the config to drop that field.
    """
    result = stats_fn()
    extra_dims = [d for d in result.dims if d != "energy"]
    if not extra_dims:
        return None
    return result.sum(dim=extra_dims).values


def _check_multiple_extensions(
    file_paths: Iterable[str | Path] | None = None,
) -> bool:
    """
    Determines if a list of file paths contains more than one unique file extension.

    This method accepts a list of file paths (as strings or `Path` objects) and checks
    if there are multiple unique file extensions present in the list. A file extension
    is identified as the substring after the last period (`.`) in the file name.

    Parameters:
        file_paths (tuple[str]): A tuple of file paths, which can be strings or
                                 `Path` objects. Defaults to None.

    Returns:
        bool: True if more than one unique file extension is found, False otherwise.

    Raises:
        TypeError: If `file_paths` is not a tuple of strings or `Path` objects.
    """
    if file_paths is None:
        return False
    extensions = {str(path).split(".")[-1] for path in file_paths if "." in str(path)}

    return len(extensions) > 1


def _collect_supported_extensions(
    parser_classes: Iterable[type],
) -> list[str]:
    """
    Collect unique supported_file_extensions from the given parser classes,
    preserving first occurrence order.
    """
    seen: set[str] = set()
    extensions: list[str] = []

    for parser in parser_classes:
        for ext in getattr(parser, "supported_file_extensions", []):
            if ext not in seen:
                seen.add(ext)
                extensions.append(ext)

    # ".txt" comes last because of the processing_order
    if ".txt" in extensions:
        extensions = [ext for ext in extensions if ext != ".txt"] + [".txt"]

    return extensions


class XPSReader(MultiFormatReader):
    """Reader for XPS."""

    supported_nxdls: list[str] = [
        "NXmpes",
        "NXxps",
    ]

    reader_dir: Path = Path(__file__).parent
    config_file: str | Path | None = reader_dir.joinpath("config", "template.json")

    parsers: list[type[_XPSParser]] = [
        PHIParser,
        ScientaHDF5Parser,
        ScientaIgorParser,
        ScientaTXTParser,
        SpecsSLEParser,
        SpecsXMLParser,
        SpecsXYParser,
        VamasExportParser,
        VamasParser,
    ]

    metadata_parsers: list[type[_XPSMetadataParser]] = [
        VamasResultParser,
    ]

    supported_file_extensions: list[str] = _collect_supported_extensions(parsers)
    supported_metadata_file_extensions: list[str] = _collect_supported_extensions(
        metadata_parsers
    )

    def __init__(self, config_file: str | None = None, *args, **kwargs):
        super().__init__(config_file, *args, **kwargs)

        self.parsed_data_dicts: list[dict[str, Any]] = []
        self.parsed_data: dict[str, Any] = {}
        self._active_parsers: list[_XPSParser] = []
        self._pending_metadata: list[_XPSMetadataParser] = []
        self.eln_data: dict[str, Any] = {}

        self.extensions: dict[str, Callable] = {
            ".yml": self.handle_eln_file,
            ".yaml": self.handle_eln_file,
            ".json": self.set_config_file,
        }

        self.processing_order: list[str] = (
            XPSReader.supported_file_extensions
            + XPSReader.supported_metadata_file_extensions
            + list(self.extensions.keys())
        )

        for ext in (
            XPSReader.supported_file_extensions
            + XPSReader.supported_metadata_file_extensions
        ):
            self.extensions[ext] = self.handle_data_file

    def set_config_file(
        self, file_path: str | Path | None, replace: bool = True
    ) -> dict[str, Any]:
        if not file_path:
            return {}

        if replace:
            if self.config_file is not None:
                if file_path != self.config_file:
                    _logger.info(
                        f"Config file already set. Replaced by the new file {file_path}."
                    )
            self.config_file = file_path
        else:
            if self.config_file is None:
                self.config_file = file_path

        return {}

    def handle_eln_file(self, file_path: str) -> dict[str, Any]:
        """
        Loads ELN file and handles specific cases.
        """

        def _combine_and_unique_string(string: str, elements: list[str]) -> str:
            """
            Combines a comma-separated string and a list into a single string with unique elements.

            Args:
                string (str): A comma-separated string.
                elements (list): A list of elements to combine with the string.

            Returns:
                str: A comma-separated string with unique elements.
            """
            existing_elements = [
                item.strip() for item in string.split(",") if item.strip()
            ]
            combined_elements = list(set(existing_elements + elements))
            combined_elements.sort()
            return ", ".join(combined_elements)

        eln_data = parse_yml(
            file_path,
            convert_dict=_CONVERT_DICT,
            replace_nested=_REPLACE_NESTED,
            parent_key="/ENTRY",
        )

        # replace paths for entry-specific ELN data
        pattern = re.compile(r"(/ENTRY)/ENTRY(\[[^\]]+\])")

        formula_keys = ("molecular_formula_hill", "chemical_formula")

        initial_eln_keys = list(eln_data.keys())

        for key, value in eln_data.copy().items():
            new_key = pattern.sub(r"\1\2", key)

            # Parse substance/molecular_formula_hill and chemical_formula into atom_types
            for form_key in formula_keys:
                if form_key in key:
                    atom_types = list(extract_atom_types(value))

                    if atom_types:
                        modified_key = re.sub(r"SUBSTANCE\[.*?\]/", "", key)
                        modified_key = modified_key.replace(form_key, "atom_types")

                        if modified_key not in initial_eln_keys:
                            if modified_key not in self.eln_data:
                                self.eln_data[modified_key] = ", ".join(atom_types)
                            else:
                                self.eln_data[modified_key] = (
                                    _combine_and_unique_string(
                                        self.eln_data[modified_key], atom_types
                                    )
                                )
                        else:
                            _logger.info(
                                f"{key} from ELN was not parsed to 'atom_types' because "
                                f"{modified_key} already exists."
                            )

            if isinstance(value, datetime.datetime):
                eln_data[key] = value.isoformat()

            self.eln_data[new_key] = eln_data.pop(key)

        return {}

    def handle_data_file(self, file_path: str) -> dict[str, Any]:
        """Dispatch *file_path* to the matching primary or metadata parser."""
        file = Path(file_path)

        primary_matches = [P for P in self.parsers if P.is_mainfile(file)]
        metadata_matches = [M for M in self.metadata_parsers if M.is_mainfile(file)]

        if primary_matches and metadata_matches:
            raise ValueError(
                f"File '{file.name}' matches both primary and metadata parsers."
            )

        if not primary_matches and not metadata_matches:
            _logger.warning("No parser matches file: %s", file_path)
            return {}

        if primary_matches:
            if len(primary_matches) > 1:
                raise ValueError(
                    f"Ambiguous file format: {len(primary_matches)} parsers match "
                    f"'{file.name}'."
                )
            parser = primary_matches[0]()
            parser.parse_file(file_path, **(self.kwargs or {}))

            config_file = parser.config_file
            if isinstance(config_file, dict):
                config_file = config_file.get(file.suffix)
            self.set_config_file(
                XPSReader.reader_dir.joinpath("config", config_file),
                replace=False,
            )

            data = parser.data
            self.parsed_data_dicts += [data]
            self._active_parsers.append(parser)

            # Apply any metadata parsers that were waiting for a compatible parser
            still_pending: list[_XPSMetadataParser] = []
            for m_parser in self._pending_metadata:
                if m_parser.__class__.supports_parser(parser):
                    m_parser.update_main_file_data(data)
                else:
                    still_pending.append(m_parser)
            self._pending_metadata = still_pending

        else:
            if len(metadata_matches) > 1:
                raise ValueError(
                    f"Ambiguous metadata format: multiple parsers match '{file.name}'."
                )
            m_parser = metadata_matches[0]()
            m_parser.parse_file(file_path, **(self.kwargs or {}))

            attached = False
            for parser in self._active_parsers:
                if m_parser.__class__.supports_parser(parser):
                    m_parser.update_main_file_data(parser.data)
                    attached = True

            if not attached:
                self._pending_metadata.append(m_parser)

        return {}

    def get_entry_names(self) -> list[str]:
        """
        Returns a list of entry names which should be constructed from the data.
        Defaults to creating a single entry named "entry".
        """
        entries = list(getattr(self, "parsed_data", {}).keys())
        return entries or ["entry"]

    def setup_template(self) -> dict[str, Any]:
        """
        Setups the initial data in the template.
        """
        # TODO: Set fixed information, e.g., about the reader.
        return {}

    def handle_objects(self, objects: tuple[Any]) -> dict[str, Any]:
        """
        Handles the objects passed into the reader.
        """
        return {}

    def post_process(self) -> None:
        """
        Do postprocessing after all files and the config file are read .
        """
        self._combine_datasets()

        # TODO: make processing of multiple entities robust
        # self.process_multiple_entities()

    def _combine_datasets(self) -> None:
        """
        This function (which is called after all files have been read)
        combines the different data sets from the XPS files and ensures
        that entry names are different and that no data is overwritten.
        """
        # Find entry names that appear in more than one parsed data dict
        entry_to_dicts: dict[str, set[int]] = {}
        for i, d in enumerate(self.parsed_data_dicts):
            for entry_name in list(d.keys()):
                if entry_name not in entry_to_dicts:
                    entry_to_dicts[entry_name] = set()
                entry_to_dicts[entry_name].add(i)

        common_entries = {
            e for e, indices in entry_to_dicts.items() if len(indices) > 1
        }

        if common_entries and not self.overwrite_keys:
            for entry in common_entries:
                dicts_with_entry = [
                    self.parsed_data_dicts[i] for i in sorted(entry_to_dicts[entry])
                ]
                for i, data_dict in enumerate(dicts_with_entry):
                    if entry in data_dict:
                        data_dict[f"{entry}{i}"] = data_dict.pop(entry)

        for data_dict in self.parsed_data_dicts:
            # If there are multiple input data files of the same type,
            # make sure that existing keys are not overwritten.
            existing = [
                (key, self.parsed_data[key], data_dict[key])
                for key in set(self.parsed_data).intersection(data_dict)
            ]
            self.parsed_data = {**self.parsed_data, **data_dict}

            if not self.overwrite_keys:
                for key, value1, value2 in existing:
                    self.parsed_data[key] = _concatenate_values(value1, value2)

    def _get_analyzer_names(self) -> list[str]:
        """
        Returns a list of analyzer names which should be constructed
        from the data. Defaults to creating a single analyzer named
        "analyzer".

        Currently, this is not used, but can be changed if there are
        multiple analyzers in the future.
        """
        analyzers: list[str] = []

        if not analyzers:
            analyzers += ["electronanalyzer"]

        return list(dict.fromkeys(analyzers))

    def _get_detector_names(self) -> list[str]:
        """
        Returns a list of detector names which should be constructed
        from the data. Defaults to creating a single detector named
        "detector".
        """
        detectors: list[str] = []

        try:
            for spectrum in self.parsed_data.values():
                if isinstance(spectrum, ParsedSpectrum) and spectrum.raw is not None:
                    n_channels = spectrum.raw.sizes.get("channel", 0)
                    for i in range(n_channels):
                        detectors.append(f"detector{i}")
        except (KeyError, AttributeError):
            pass

        if not detectors:
            detectors += ["detector"]

        return list(dict.fromkeys(detectors))

    def _process_multiple_entities(self) -> None:
        """
        Check if there are multiple of some class and, if so, change the
        keys and values in the config file.

        This replaces all occurrences of "detector" and "electronanalyzer"
        in the config dict by the respective names (e.g., detector0, detector1)
        and removes the generic term if there are multiple different instances.

        """
        multiples_to_check = {
            "electronanalyzer": self._get_analyzer_names,
            "detector": self._get_detector_names,
        }

        """
        Currently, it only works if the config_dict is loaded BEFORE the
        parse_json_config method".

        In principle, the same replacement should be done for the eln and
        (meta)data dicts.
        """

        for config_key, config_value in self.config_dict.copy().items():
            for original_key, search_func in multiples_to_check.items():
                entity_names = search_func()

                if len(entity_names) >= 1 and entity_names[0] is not original_key:
                    for name in entity_names:
                        modified_value = copy.deepcopy(config_value)

                        modified_key = config_key.replace(
                            f"[{original_key}]",
                            f"[{name}]",
                        )

                        if (
                            isinstance(config_value, str)
                            and f"{original_key}/" in config_value
                        ):
                            modified_value = config_value.replace(
                                f"{original_key}/", f"{name}/"
                            )

                        self.config_dict[modified_key] = modified_value
                        del self.config_dict[config_key]

    def _search_metadata(self, metadata: dict[str, Any], path: str) -> Any:
        """Suffix-match ``path`` in a flat metadata dict.

        Tries an exact key match first, then falls back to any key ending with
        ``"/" + path``.  Returns ``None`` when the value is empty/missing.
        """
        _sentinel = object()
        value = metadata.get(path, _sentinel)
        if value is _sentinel:
            value = next(
                (v for k, v in metadata.items() if k.endswith(f"/{path}")), None
            )
        if (
            value is None
            or str(value) in {"None", ""}
            or (isinstance(value, list) and all(v == "" for v in value))
        ):
            return None
        return value.isoformat() if isinstance(value, datetime.datetime) else value

    def _get_process_sequence_index(
        self, key: str, spectrum: "ParsedSpectrum"
    ) -> int | None:
        """Return the 1-based sequence index for the PROCESS group in *key*.

        Returns ``None`` when that process is not applicable for *spectrum*,
        which causes the ``!@attrs:sequence_index`` notation to drop the whole
        PROCESS group.
        """
        m = re.search(r"PROCESS\[(\w+)\]", key)
        if m is None:
            return None
        name = m.group(1)
        active = [p for p, pred in _PROCESS_ORDER if pred(spectrum)]
        if name not in active:
            return None
        return active.index(name) + 1

    def get_attr(self, key: str, path: str) -> Any:
        """Return metadata stored in the parsed spectrum for the current entry."""
        spectrum: ParsedSpectrum | None = self.parsed_data.get(
            self.callbacks.entry_name
        )
        if spectrum is None:
            return None
        if path == "sequence_index":
            return self._get_process_sequence_index(key, spectrum)
        return self._search_metadata(spectrum.metadata, path)

    def get_eln_data(self, key: str, path: str) -> Any:
        """
        Returns data from the given eln path.
        Gives preference to ELN data for a given entry before searching
        the ELN data for all entries.
        Returns None if the path does not exist.
        """
        if key in self.eln_data:
            return self.eln_data.get(key)

        else:
            # check for similar key with generic /ENTRY/
            pattern = re.compile(r"(/ENTRY)\[[^\]]+\]")
            modified_key = pattern.sub(r"\1", key)
            if modified_key in self.eln_data:
                return self.eln_data.get(modified_key)
        return

    def get_data_dims(self, key: str, path: str) -> list[str]:
        """
        Returns the dimensions of the data from the given path.
        """
        entry = self.callbacks.entry_name
        spectrum: ParsedSpectrum | None = self.parsed_data.get(entry)

        if spectrum is None:
            return []

        def get_all_keys(template_key: str) -> list[str]:
            keys = {
                k[len(template_key) :].split("/")[0]
                for k in spectrum.metadata
                if k.startswith(template_key)
            }
            return sorted(keys)

        # Peak and background wildcards (fit-related)
        if re.search(r"peak", key):
            return get_all_keys("component")
        if re.search(r"background", key):
            return get_all_keys("region")

        # External channel wildcards (e.g., SPECS-XY additional axes)
        if isinstance(path, str) and path.endswith(("*.external", "*.external_unit")):
            return get_all_keys("external_")

        # Axis name wildcard
        if isinstance(path, str) and path.endswith(".unit"):
            return [
                cast(str, d) for d in spectrum.data.dims if d not in ("cycle", "scan")
            ]

        # Axis data wildcard (e.g., AXISNAME[*] → @data:*.axes)
        if isinstance(path, str) and path.endswith(".axes"):
            return [
                cast(str, d) for d in spectrum.data.dims if d not in ("cycle", "scan")
            ]

        return []

    def get_data(self, key: str, path: str) -> Any | None:
        """Retrieve XPS spectral data for the current entry by path token."""
        data_handlers: dict[str, Callable[[ParsedSpectrum], Any]] = {
            # stats
            "average": lambda s: s.average().values,
            "average_reduced": lambda s: _stats_reduced(s, s.average),
            "errors": lambda s: s.errors().values,
            "errors_reduced": lambda s: _stats_reduced(s, s.errors),
            # NXprocess averaging
            "channels_averaging": lambda s: (
                s.data.values
                if s.raw is not None and s.raw.sizes.get("channel", 1) > 1
                else None
            ),
            "scan_averaging": lambda s: (
                s.data.mean("scan").values
                if "scan" in s.data.dims and s.data.sizes["scan"] > 1
                else None
            ),
            "cycle_averaging": _cycle_averaging,
            # raw data
            "raw": _get_raw,
            # processed axis indices
            "cycle": lambda s: (
                np.arange(s.data.sizes["cycle"]) if "cycle" in s.data.dims else None
            ),
            "scan": lambda s: (
                np.arange(s.data.sizes["scan"]) if "scan" in s.data.dims else None
            ),
            "channel": lambda s: (
                np.arange(s.raw.sizes["channel"]) if s.raw is not None else None
            ),
            # raw axis indices (fall back to processed DA, require size > 1)
            "raw_cycle": lambda s: _raw_dim_indices(s, "cycle"),
            "raw_scan": lambda s: _raw_dim_indices(s, "scan"),
            "raw_channel": lambda s: (
                np.arange(s.raw.sizes["channel"])
                if s.raw is not None and s.raw.sizes.get("channel", 0) > 1
                else None
            ),
            # raw axis metadata
            "raw_axes": _raw_active_dims,
            "raw_energy_index": _get_raw_energy_index,
            "raw_units": lambda s: (
                "counts" if s.raw is not None else "counts_per_second"
            ),
        }

        entry = self.callbacks.entry_name
        spectrum: ParsedSpectrum | None = self.parsed_data.get(entry)
        if spectrum is None:
            return None

        if path in data_handlers:
            return data_handlers[path](spectrum)

        if path.endswith(".axes"):
            axis = path[: -len(".axes")]
            coord = spectrum.data.coords.get(axis)
            return np.array(coord.values) if coord is not None else None

        if path.endswith(".unit"):
            channel = path[: -len(".unit")]
            return self._search_metadata(spectrum.metadata, f"{channel}/@units")

        if path.endswith((".external", ".external_unit")):
            channel = path.split(".external", maxsplit=1)[0]
            if path.endswith("_unit"):
                return self._search_metadata(
                    spectrum.metadata, f"external_{channel}/@units"
                )
            value = self._search_metadata(spectrum.metadata, f"external_{channel}")
            return np.array(value) if value is not None else None

        # default: coordinate lookup
        coord = spectrum.data.coords.get(path)
        return np.array(coord.values) if coord is not None else None

    def set_nxdata_defaults(self, template):
        """Set the default for automatic plotting."""
        survey_count, count = 0, 0

        def _get_unique_nxfit_names(template) -> set[str]:
            """Extract unique 'ENTRY[<some-name>]/FIT[<some-other-name>]' pairs from template keys."""
            pattern = re.compile(r"^/?ENTRY\[(?P<entry>[^]]+)\]/FIT\[(?P<fit>[^]]+)\]/")

            result = set()
            for key in template:
                m = pattern.match(key)
                if m:
                    result.add(f"{m.group('entry')}/{m.group('fit')}")
            return result

        def _get_first_matching_fit(
            entry_name: str, unique_fits: list[str]
        ) -> str | None:
            """Return the first '<fit>' name that matches the given entry name, if any."""
            for fit in unique_fits:
                if fit.startswith(f"{entry_name}/"):
                    return fit.split("/", 1)[1]  # Extract only the fit name
            return None

        # Sorting for deterministic ordering
        unique_fits = sorted(_get_unique_nxfit_names(template))

        for entry in self.get_entry_names():
            if unique_fits:
                template["/@default"] = unique_fits[0]
                match = _get_first_matching_fit(entry, unique_fits)
                if match:
                    template[f"/ENTRY[{entry}]/@default"] = match
                else:
                    template[f"/ENTRY[{entry}]/@default"] = "data"

            else:
                if "Survey" in entry and survey_count == 0:
                    survey_count += 1
                    template["/@default"] = entry

                # If no Survey, set any scan for default
                elif survey_count == 0 and count == 0:
                    count += 1
                    template["/@default"] = entry

    def read(
        self,
        template: dict = None,
        file_paths: tuple[str] = None,
        objects: tuple[Any] | None = None,
        **kwargs,
    ) -> dict:
        self.overwrite_keys = _check_multiple_extensions(file_paths)

        self.set_config_file(kwargs.get("config_file", self.config_file))

        template = super().read(template, file_paths, objects, suppress_warning=True)
        self.set_nxdata_defaults(template)

        final_template = Template()
        for key, val in template.items():
            if key.endswith("@units"):
                parent = key.replace("/@units", "")
                if parent not in template:
                    continue
            if val is not None:
                if "@units" in key:
                    check_units(key, val)
                final_template[key] = val
                if isinstance(val, dict) and "link" in val:
                    final_template[f"{key}/@target"] = val["link"]

        return final_template


READER = XPSReader
