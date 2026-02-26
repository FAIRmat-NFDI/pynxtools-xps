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
import os
import re
from collections.abc import Callable, Iterable
from pathlib import Path
from re import Pattern
from typing import Any

import numpy as np
from pynxtools.dataconverter.helpers import extract_atom_types
from pynxtools.dataconverter.readers.multi.reader import MultiFormatReader
from pynxtools.dataconverter.readers.utils import parse_yml
from pynxtools.dataconverter.template import Template

from pynxtools_xps.logging import _logger
from pynxtools_xps.numerics import check_units
from pynxtools_xps.parsers import (
    PHIMapper,
    ScientaMapper,
    SpecsSLEMapper,
    SpecsXMLMapper,
    SpecsXYMapper,
    VamasExportMapper,
    VamasMapper,
    VamasResultMapper,
)
from pynxtools_xps.parsers.base import (
    # ParsedSpectrum,
    _XPSMapper,
    _XPSMetadataParser,
    _XPSParser,
)

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
    "detector": "DETECTOR[detector]",
    "manipulator": "MANIPULATOR[manipulator]",
    "pid": "PID[pid]",
    "process": "PROCESS[process]",
    "sample": "SAMPLE[sample]",
    "substance": "SUBSTANCE[substance]",
}

_REPLACE_NESTED: dict[str, str] = {}

_CHAN_COUNT = "_chan"
_SCAN_COUNT = "_scan"


def _get_channel_vars(data_vars: list[str]) -> list[str]:
    """Get all data vars that contain _chan."""
    return [data_var for data_var in data_vars if _CHAN_COUNT in data_var]


def _get_scan_vars(data_vars: list[str]) -> list[str]:
    """Get all data vars that contain _scan, but not _chan."""
    return [
        data_var
        for data_var in data_vars
        if _SCAN_COUNT in data_var and _CHAN_COUNT not in data_var
    ]


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


# TODO: enable
# def _assemble_spectra(
#     spectra: list[ParsedSpectrum],
#     file_path: str | Path,
# ) -> dict[str, Any]:
#     """Build ENTRY-keyed dict + xarray datasets from parsed spectra.

#     This is the single, generic replacement for all per-parser
#     ``_XPSMapper.construct_data`` / ``_update_xps_dict_with_spectrum``
#     methods.
#     """
#     data_dict: dict[str, Any] = {
#         "File": file_path,
#         "file_ext": os.path.splitext(file_path)[1],
#         "data": {},
#     }

#     for spectrum in spectra:
#         spectrum.validate()

#         entry = _construct_entry_name(
#             [spectrum.group_name, spectrum.spectrum_type]
#         )
#         if not entry:
#             entry = "entry"

#         entry_parent = f"/ENTRY[{entry}]"

#         # Write metadata as ENTRY-prefixed keys
#         for key, value in spectrum.metadata.items():
#             if key.startswith("entry/"):
#                 mpes_key = f"/ENTRY[entry]/{key.replace('entry/', '', 1)}"
#             else:
#                 mpes_key = f"{entry_parent}/{key}"
#             data_dict[mpes_key] = value

#         # Create or reuse xarray Dataset for this entry
#         if entry not in data_dict["data"]:
#             data_dict["data"][entry] = xr.Dataset()
#         ds = data_dict["data"][entry]

#         base_key = _construct_data_key(spectrum)

#         # Write individual scans
#         for scan in spectrum.scans:
#             if len(spectrum.scans) > 1:
#                 s_key = f"{base_key}_scan{scan.scan_id}"
#             else:
#                 s_key = base_key

#             ds[s_key] = xr.DataArray(
#                 data=scan.intensity,
#                 coords={"energy": spectrum.energy},
#             )

#             # Channel-resolved data
#             if scan.channels is not None:
#                 for i in range(scan.channels.shape[1]):
#                     ds[f"{s_key}_chan{i}"] = xr.DataArray(
#                         data=scan.channels[:, i],
#                         coords={"energy": spectrum.energy},
#                     )

#             # Raw counts (stored as chan0 when no channel separation)
#             if scan.raw_counts is not None and scan.channels is None:
#                 ds[f"{s_key}_chan0"] = xr.DataArray(
#                     data=scan.raw_counts,
#                     coords={"energy": spectrum.energy},
#                 )

#         # Write average across scans (if multiple)
#         if len(spectrum.scans) > 1:
#             with warnings.catch_warnings():
#                 warnings.simplefilter("ignore", category=RuntimeWarning)
#                 averaged = np.mean(
#                     [s.intensity for s in spectrum.scans], axis=0
#                 )
#             ds[base_key] = xr.DataArray(
#                 data=averaged,
#                 coords={"energy": spectrum.energy},
#             )

#     return data_dict


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
    extensions = {str(path).split(".")[-1] for path in file_paths if "." in str(path)}

    return len(extensions) > 1


# pylint: disable=too-few-public-methods
class XPSReader(MultiFormatReader):
    """Reader for XPS."""

    supported_nxdls: list[str] = [
        "NXmpes",
        "NXxps",
    ]

    reader_dir: Path = Path(__file__).parent
    config_file: str | Path | None = reader_dir.joinpath("config", "template.json")

    # TODO: implement
    parsers: list[type[_XPSParser]] = [
        # KratosParser,
    ]

    metadata_parsers: list[type[_XPSMetadataParser]] = [
        # TODO: implement
        # VamasExportMapper,
    ]

    mappers: list[type[_XPSMapper]] = [
        PHIMapper,
        ScientaMapper,
        SpecsSLEMapper,
        SpecsXMLMapper,
        SpecsXYMapper,
        VamasMapper,
        VamasResultMapper,
        VamasExportMapper,
    ]

    supported_file_extensions: list[str] = [
        ".h5",
        ".hdf5",
        ".ibw",
        ".npl",
        ".pro",
        ".spe",
        ".sle",
        ".slh",
        ".vms",
        ".xml",
        ".xy",
        ".txt",  # This is last because of the processing_order
    ]

    supported_metadata_file_extensions: dict[str, str] = {".csv": ".txt"}
    supported_vendors: list[str] = ["kratos", "phi", "scienta", "specs", "unknown"]
    vendor_map: dict[str, dict[str, Any]] = {
        ".csv": {"unknown": VamasResultMapper},
        ".h5": {"scienta": ScientaMapper},
        ".hdf5": {"scienta": ScientaMapper},
        ".ibw": {"scienta": ScientaMapper},
        ".npl": {"unknown": VamasMapper},
        ".pro": {"phi": PHIMapper},
        ".spe": {"phi": PHIMapper},
        ".sle": {"specs": SpecsSLEMapper},
        ".txt": {
            "scienta": ScientaMapper,
            "unknown": VamasExportMapper,
        },
        ".vms": {"unknown": VamasMapper},
        ".xml": {"specs": SpecsXMLMapper},
        ".xy": {"specs": SpecsXYMapper},
    }

    file_err_msg: str = (
        "Need an XPS data file with one of the following extensions: "
        f"data files: {supported_file_extensions}, metadata files: {supported_metadata_file_extensions}."
    )

    vendor_err_msg: str = (
        f"Need an XPS data file from one of the following vendors: {supported_vendors}"
    )

    def __init__(self, config_file: str | None = None, *args, **kwargs):
        super().__init__(config_file, *args, **kwargs)

        self.parsed_data_dicts: list[dict[str, Any]] = []
        self.parsed_data: dict[str, Any] = {}
        self.eln_data: dict[str, Any] = {}

        self.extensions: dict[str, Callable] = {
            ".yml": self.handle_eln_file,
            ".yaml": self.handle_eln_file,
            ".json": self.set_config_file,
        }

        self.processing_order: list[str] = (
            XPSReader.supported_file_extensions
            + list(XPSReader.supported_metadata_file_extensions.keys())
            + list(self.extensions.keys())
        )

        for ext in XPSReader.supported_file_extensions + list(
            XPSReader.supported_metadata_file_extensions.keys()
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
                                f"{key} from ELN was not parsed to atom_types because {modified_key} already exists."
                            )

            if isinstance(value, datetime.datetime):
                eln_data[key] = value.isoformat()

            self.eln_data[new_key] = eln_data.pop(key)

        return {}

    def handle_data_file(self, file_path: str) -> dict[str, Any]:
        # # # TODO: this is structurally incorrect, we need this:
        # # # self.datasets: list[Dataset] = []
        # # # self.pending_metadata: list[_XPSMetadataParser] = []
        # # # class Dataset:
        # # #     def __init__(self, parser: _XPSParser, file: Path):
        # # #         self.parser = parser
        # # #         self.file = file
        # # #         self.data = parser.data

        # # # # Important Final Detail: Metadata should not be applied multiple times to the same dataset.
        # # # # # If that is a risk, track application:
        # # # # dataset.applied_metadata: set[type]

        # # # file = Path(file_path)

        # # # primary_matches = [
        # # #     P for P in self.parsers
        # # #     if P.is_mainfile(file)
        # # # ]

        # # # metadata_matches = [
        # # #     M for M in self.metadata_parsers
        # # #     if M.is_mainfile(file)
        # # # ]

        # # # # Handle ambiguity strictly
        # # # if primary_matches and metadata_matches:
        # # #     raise ValueError("File matches both primary and metadata parsers.")

        # # # if not primary_matches and not metadata_matches:
        # # #     raise ValueError("No parser supports this file.")

        # # # if primary_matches:
        # # #     ParserCls = _select_unique(primary_matches)
        # # #     parser = ParserCls()
        # # #     parser.parse_file(file, **self.kwargs)

        # # #     dataset = Dataset(parser, file)
        # # #     self.datasets.append(dataset)

        # # #     # Try to apply any metadata that arrived earlier
        # # #     self._try_attach_pending_metadata(dataset)

        # # # if metadata_matches:
        # # #     MetadataCls = _select_unique(metadata_matches)
        # # #     metadata = MetadataCls()
        # # #     metadata.parse_file(file)

        # # #     attached = False

        # # #     for dataset in self.datasets:
        # # #         if MetadataCls.supports_parser(dataset.parser):
        # # #             metadata.update_main_file_data(dataset.data)
        # # #             attached = True

        # # #     if not attached:
        # # #         # No compatible dataset yet → defer
        # # #         self.pending_metadata.append(metadata)

        # # # def _try_attach_pending_metadata(self, dataset: Dataset):
        # # #     still_pending = []

        # # #     for metadata in self.pending_metadata:
        # # #         if metadata.__class__.supports_parser(dataset.parser):
        # # #             metadata.update_main_file_data(dataset.data)
        # # #         else:
        # # #             still_pending.append(metadata)

        # # #  self.pending_metadata = still_pending

        # def _select_parser(file_path: Path) -> type[_XPSParser]:
        #     matches = [
        #         ParserCls
        #         for ParserCls in self.parsers
        #         if ParserCls.is_mainfile(file)
        #     ]

        #     if len(matches) == 0:
        #         raise ValueError("No parser supports this file.")
        #     if len(matches) > 1:
        #         raise ValueError("Ambiguous file format: multiple parsers match.")

        #     parser = matches[0]

        #     return parser

        # def _select_metadata_parsers(
        #     file: Path,
        #     parser:_XPSParser,
        # ) -> list[type[_XPSMetadataParser]]:

        #     fitting_metadata_parsers: list[_XPSMetadataParser] = []

        #     matches = [
        #         ParserCls
        #         for ParserCls in self.metadata_parsers
        #         if ParserCls.is_mainfile(file)
        #     ]

        #     if not matches:
        #         return []

        #     for metadata_parser in matches:
        #         if not metadata_parser.supports_parser(parser):
        #             _logger.warning(
        #                 f"Metadata parser {metadata_parser.__name__} does not support the main parser {parser.__class__.__name__}."
        #             )
        #             continue
        #         fitting_metadata_parsers.append(metadata_parser)

        #     return fitting_metadata_parsers

        # file = Path(file_path)
        # main_file_ext = file.suffix
        # file_ext = file.suffix

        # parser = _select_parser(file)()
        # parser.parse_file(file_path, **self.kwargs)
        # main_file_data = parser.data

        # metadata_parsers = _select_metadata_parsers(file, parser)

        # for m_parser in metadata_parsers:
        #     m_parser().update_main_file_data(main_file_data)

        # config_file = parser.config_file
        # if isinstance(config_file, dict):
        #     config_file = config_file.get(file_ext)
        # self.set_config_file(
        #     XPSReader.reader_dir.joinpath("config", config_file),
        #     replace=False,
        # )

        # self.parsed_data_dicts += [parser.data]

        # return {}

        def _check_for_vendors(file_path: str) -> str:
            """
            Check for the vendor name of the XPS data file.

            """
            _, file_ext = os.path.splitext(file_path)

            vendor_dict = XPSReader.vendor_map[file_ext]

            if len(vendor_dict) == 1:
                return list(vendor_dict.keys())[0]
            if file_ext == ".txt":
                return _check_for_vendors_txt(file_path)
            raise ValueError(XPSReader.vendor_err_msg)

        def _check_for_vendors_txt(file_path: str) -> str:
            """
            Search for a vendor names in a txt file

            Args:
                file (str): XPS txt file
            Returns:
                str: vendor (str): Vendor name if that name is in the txt file or "unknown"

            """
            vendor_dict = XPSReader.vendor_map[".txt"]

            with open(file_path, encoding="utf-8") as txt_file:
                contents = txt_file.read()

            for vendor in vendor_dict:
                vendor_options = [vendor, vendor.upper(), vendor.capitalize()]

                if any(vendor_opt in contents for vendor_opt in vendor_options):
                    return vendor
                if contents[:6] == "[Info]":
                    # This is for picking the Scienta reader if "scienta"
                    # is not in the file
                    return vendor
            return "unknown"

        _, file_ext = os.path.splitext(file_path)

        if file_ext in XPSReader.supported_file_extensions:
            vendor = _check_for_vendors(file_path)

            parser = XPSReader.vendor_map[file_ext][vendor]()
            parser.parse_file(file_path, **self.kwargs)

            # New parser path: returns list[ParsedSpectrum]
            # spectra = handler.parse_file(file_path, **self.kwargs)
            # data = _assemble_spectra(spectra, file_path)
            # self.parsed_data_dicts += [data]

            config_file = parser.config_file
            if isinstance(config_file, dict):
                config_file = config_file.get(file_ext)
            self.set_config_file(
                XPSReader.reader_dir.joinpath("config", config_file),
                replace=False,
            )

            self.parsed_data_dicts += [parser.data]

        elif file_ext in XPSReader.supported_metadata_file_extensions:
            vendor = _check_for_vendors(file_path)

            metadata_parser = XPSReader.vendor_map[file_ext][vendor]()
            metadata_parser.parse_file(file_path, **self.kwargs)

            main_file_ext = XPSReader.supported_metadata_file_extensions[file_ext]

            main_file_dicts = [
                d for d in self.parsed_data_dicts if d.get("file_ext") == main_file_ext
            ]

            metadata_parser.update_main_file_dict(main_file_dicts)

        return {}

    def get_entry_names(self) -> list[str]:
        """
        Returns a list of entry names which should be constructed from the data.
        Defaults to creating a single entry named "entry".
        """
        # Track entries for using for eln data
        entries: list[str] = []

        try:
            for entry in self.parsed_data["data"]:
                entries += [entry]
        except KeyError:
            pass

        if not entries:
            entries += ["entry"]

        return list(dict.fromkeys(entries))

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
        This function (which as called after all files have been read)
        combines the different data sets from the XPS files and ensures
        that entry names are different and that no data is overwritten.

        """

        def check_for_same_entries(
            dicts: list[dict[str, Any]],
        ) -> tuple[set[str], list[set[int]]] | None:
            """
            Checks whether the input dictionaries have the same entries and identifies which
            dictionaries contain each common entry.

            Args:
                dicts (list[dict[str, any]]): A list of dictionaries with keys potentially containing entries.

            Returns:
                Optional[tuple[Set[str], list[Set[int]]]]:
                - A tuple where:
                  - The first element is a set of common entries found across the dictionaries.
                  - The second element is a list of sets, where each set contains the indices of dictionaries
                    that have the corresponding entry.
                - None if no common entries are found.
            """
            entry_pattern = re.compile(r"/ENTRY\[(.*?)\]")

            entry_to_dicts: dict[str, set] = {}

            for i, d in enumerate(dicts):
                for key in d.keys():
                    entries = entry_pattern.findall(key)
                    for entry in entries:
                        if entry not in entry_to_dicts:
                            entry_to_dicts[entry] = set()
                        entry_to_dicts[entry].add(i)

            common_entries = {
                entry
                for entry, dict_indices in entry_to_dicts.items()
                if len(dict_indices) > 1
            }

            if not common_entries:
                return None, None

            dict_indices = [entry_to_dicts[entry] for entry in common_entries]

            return common_entries, dict_indices

        common_entries, dict_indices = check_for_same_entries(self.parsed_data_dicts)

        if common_entries and not self.overwrite_keys:
            for entry, indices in zip(common_entries, dict_indices):
                dicts_with_common_entries = [self.parsed_data_dicts[i] for i in indices]

                for i, data_dict in enumerate(dicts_with_common_entries):
                    for key, value in data_dict.copy().items():
                        new_key = key.replace(f"/ENTRY[{entry}]", f"/ENTRY[{entry}{i}]")
                        if key == "data":
                            for entry_name, xarr in value.copy().items():
                                if entry_name == entry:
                                    new_entry_name = entry_name.replace(
                                        f"{entry}", f"{entry}{i}"
                                    )
                                    value[new_entry_name] = xarr
                                    del value[entry_name]
                        if new_key != key:
                            data_dict[new_key] = value
                            del data_dict[key]

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
            for entry, entry_values in self.parsed_data["data"].items():
                for data_var in entry_values:
                    if _CHAN_COUNT in data_var:
                        detector_num = data_var.split(_CHAN_COUNT)[-1]
                        detector_nm = f"detector{detector_num}"
                        detectors += [detector_nm]
        except KeyError:
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

    def _get_metadata(
        self,
        metadata_dict: dict[str, Any],
        path: str,
        entry_name: str,
    ) -> Any:
        """
        Get metadata from the ELN or XPS data dictionaries.

        Note that the keys of metadata_dict may contain more than
        the path, i.e.,
        /ENTRY[my-entry]/instrument/analyzer/collectioncolumn/voltage.
        With the regex, the path = "collectioncolumn/voltage" would
        still yield the correct value.

        Parameters
        ----------
        metadata_dict : dict[str, Any]
            One of ELN or XPS data dictionaries .
        path : str
            Path to search in the metadata_dict
        entry_name : str
            Entry name to search.

        Yields
        ------
        value: Any
            The value in the metadata_dict.

        """
        pattern = re.compile(
            rf"^/ENTRY\[{re.escape(entry_name)}\](?:/.*/|/){re.escape(path)}$"
        )

        matching_key = next((key for key in metadata_dict if pattern.match(key)), None)

        value = metadata_dict.get(matching_key)

        if (
            value is None
            or str(value) in {"None", ""}
            or (isinstance(value, list) and all(v == "" for v in value))
        ):
            return

        if isinstance(value, datetime.datetime):
            value = value.isoformat()

        return value

    def get_attr(self, key: str, path: str) -> Any:
        """
        Get the metadata that was stored in the main file.
        """
        return self._get_metadata(self.parsed_data, path, self.callbacks.entry_name)

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
        escaped_entry = re.escape(entry)
        xr_data = self.parsed_data["data"].get(entry)

        if xr_data is None:
            return []

        def get_signals(key: str) -> list[str]:
            if key == "scans":
                data_vars = _get_scan_vars(xr_data.data_vars)
            elif key == "channels":
                data_vars = _get_channel_vars(xr_data.data_vars)
                if not data_vars:
                    data_vars = _get_scan_vars(xr_data.data_vars)
            elif key == "axes":
                data_vars = list(xr_data.coords)
            else:
                data_vars = [""]

            return list(map(str, data_vars))

        def get_all_keys(template_key: str) -> list[str]:
            pattern = re.compile(rf"^/ENTRY\[{escaped_entry}]/{template_key}([^/]+)")

            keys = {
                match[1] for key in self.parsed_data if (match := pattern.search(key))
            }

            return sorted(keys)

        if isinstance(path, str) and path.endswith(("*.external", "*.external_unit")):
            data_func = lambda: get_all_keys("external_")
        else:
            if isinstance(path, str) and path.endswith(".unit"):
                data_func = lambda: list(xr_data.coords)
            else:
                data_func = lambda: get_signals(path.rsplit(":*.", maxsplit=1)[-1])

        PatternHandler = dict[Pattern[str], Callable[[], Any]]

        patterns: dict[re.Pattern, Callable[[], Any]] = {
            re.compile(r"peak"): lambda: get_all_keys("component"),
            re.compile(r"background"): lambda: get_all_keys("region"),
            re.compile(
                r"DATA\[[^\]]+\]/(?:DATA\[[^\]]+\]|AXISNAME\[[^\]]+\])(?:/@units)?"
            ): data_func,
            re.compile(r"ELECTRON_DETECTOR\[[a-zA-Z0-9_]+\]/raw_data"): lambda: (
                get_signals("channels")
            ),
        }

        # Function to match and return the handler result
        for pattern, handler in patterns.items():
            if pattern.search(key):
                return handler()

        return get_signals(key="scans")

    def _search_first(self, data: dict[str, Any], pattern: re.Pattern) -> Any | None:
        """
        Search for the first value in a dictionary whose key matches a regex pattern.

        Parameters:
            data (dict[str, Any]): The dictionary to search.
            pattern (re.Pattern): The compiled regex pattern to search for in keys.

        Returns:
            The first matching value, or None if no match is found.
        """
        return next((value for key, value in data.items() if pattern.search(key)), None)

    def get_data(self, key: str, path: str) -> Any | None:
        """
        Retrieve XPS data based on a key and a path string. The method supports multiple
        forms of data access, such as averages, raw data, axis values, units, and external links.

        Parameters:
            key (str): The key under which the data is stored (unused in current logic).
            path (str): A string path that determines what type of data to return.

        Returns:
            Any: The requested data, or None if the path is invalid or not found.
        """
        entry = self.callbacks.entry_name
        escaped_entry = re.escape(entry)
        xr_data = self.parsed_data["data"].get(entry)

        if xr_data is None:
            return None

        # Average or errors
        if path.startswith(("average", "errors")):
            data = [xr_data[var].data for var in _get_scan_vars(xr_data.data_vars)]
            if path.endswith("reduced"):
                try:
                    data = [np.sum(arr, axis=tuple(range(1, arr.ndim))) for arr in data]
                except np.exceptions.AxisError:
                    return None
            stats_func = np.mean if path.startswith("average") else np.std
            return stats_func(data, axis=0)

        # Raw data
        if path.endswith("raw_data"):
            data_vars = _get_channel_vars(xr_data.data_vars) or _get_scan_vars(
                xr_data.data_vars
            )
            return np.array(
                [xr_data[var].data for var in data_vars if _SCAN_COUNT in var]
            )

        # Channels or scans by suffix
        for suffix, extractor in {
            ".scans": xr_data.get,
            ".channels": xr_data.get,
        }.items():
            if path.endswith(suffix):
                name = path.split(suffix, maxsplit=1)[0]
                data = extractor(name)
                return np.array(data) if data is not None else None

        # Axis or coordinate data
        if path.endswith(".axes"):
            axis = path.split(".axes", maxsplit=1)[0]
            coord = xr_data.coords.get(axis)
            return np.array(coord.values) if coord is not None else None

        # Units for internal channels
        if path.endswith(".unit"):
            channel = path.split(".unit", maxsplit=1)[0]
            pattern = re.compile(
                rf"^/ENTRY\[{escaped_entry}]/{re.escape(channel)}/@units"
            )
            return self._search_first(self.parsed_data, pattern)

        # External channels and units
        if path.endswith((".external", ".external_unit")):
            channel = path.split(".external", maxsplit=1)[0]
            if path.endswith("_unit"):
                pattern = re.compile(
                    rf"^/ENTRY\[{escaped_entry}]/external_{re.escape(channel)}/@units"
                )
                return self._search_first(self.parsed_data, pattern)
            else:
                pattern = re.compile(
                    rf"^/ENTRY\[{escaped_entry}]/external_{re.escape(channel)}$"
                )
                matches = [
                    value
                    for key, value in self.parsed_data.items()
                    if pattern.search(key)
                ]
                return np.array(matches).squeeze() if matches else None

        # Default: try direct data access
        try:
            return xr_data[path]
        except KeyError:
            try:
                return np.array(xr_data.coords[path].values)
            except KeyError:
                return None

    def set_nxdata_defaults(self, template):
        """Set the default for automatic plotting."""
        survey_count, count = 0, 0

        def get_unique_nxfit_names(template) -> set[str]:
            """Extract unique 'ENTRY[<some-name>]/FIT[<some-other-name>]' pairs from template keys."""
            pattern = re.compile(r"^/?ENTRY\[(?P<entry>[^]]+)\]/FIT\[(?P<fit>[^]]+)\]/")

            result = set()
            for key in template:
                m = pattern.match(key)
                if m:
                    result.add(f"{m.group('entry')}/{m.group('fit')}")
            return result

        def get_first_matching_fit(
            entry_name: str, unique_fits: set[str]
        ) -> str | None:
            """Return the first '<fit>' name that matches the given entry name, if any."""
            for fit in unique_fits:
                if fit.startswith(f"{entry_name}/"):
                    return fit.split("/", 1)[1]  # Extract only the fit name
            return None

        unique_fits = sorted(
            get_unique_nxfit_names(template)
        )  # Sorting for deterministic ordering

        for entry in self.get_entry_names():
            if unique_fits:
                template["/@default"] = unique_fits[0]
                match = get_first_matching_fit(entry, unique_fits)
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
        objects: tuple[Any] = None,
        **kwargs,
    ) -> dict:
        self.overwrite_keys = _check_multiple_extensions(file_paths)

        if "config_file" in kwargs:
            self.set_config_file(kwargs.get("config_file"))

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
