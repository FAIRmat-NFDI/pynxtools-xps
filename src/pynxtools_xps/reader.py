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
import logging
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union, Pattern, Callable

import numpy as np
from pynxtools.dataconverter.helpers import extract_atom_types
from pynxtools.dataconverter.readers.multi.reader import MultiFormatReader
from pynxtools.dataconverter.readers.utils import parse_yml
from pynxtools.dataconverter.template import Template

from pynxtools_xps.reader_utils import check_units

from pynxtools_xps.phi.spe_pro_phi import MapperPhi
from pynxtools_xps.scienta.scienta_reader import MapperScienta
from pynxtools_xps.specs.sle.sle_specs import SleMapperSpecs
from pynxtools_xps.specs.xml.xml_specs import XmlMapperSpecs
from pynxtools_xps.specs.xy.xy_specs import XyMapperSpecs
from pynxtools_xps.vms.vamas_export import TxtMapperVamasExport, CsvMapperVamasResult
from pynxtools_xps.vms.vamas import VamasMapper

logger = logging.getLogger("pynxtools")

CONVERT_DICT = {
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

REPLACE_NESTED: Dict[str, str] = {}

CHAN_COUNT = "_chan"
SCAN_COUNT = "_scan"


def _get_channel_vars(data_vars: List[str]):
    """Get all data vars that containt _chan."""
    return [data_var for data_var in data_vars if CHAN_COUNT in data_var]


def _get_scan_vars(data_vars: List[str]):
    """Get all data vars that contain _scan, but not _chan."""
    return [
        data_var
        for data_var in data_vars
        if SCAN_COUNT in data_var and CHAN_COUNT not in data_var
    ]


def concatenate_values(value1, value2):
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


def _check_multiple_extensions(file_paths: Tuple[str] = None) -> bool:
    """
    Determines if a list of file paths contains more than one unique file extension.

    This method accepts a list of file paths (as strings or `Path` objects) and checks
    if there are multiple unique file extensions present in the list. A file extension
    is identified as the substring after the last period (`.`) in the file name.

    Parameters:
        file_paths (Tuple[str]): A tuple of file paths, which can be strings or
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

    supported_nxdls = [
        "NXmpes",
        "NXxps",
    ]

    reader_dir = Path(__file__).parent
    config_file: Optional[Union[str, Path]] = reader_dir.joinpath(
        "config", "template.json"
    )

    __prmt_file_ext__ = [
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

    __prmt_metadata_file_ext__ = {".csv": ".txt"}

    __vendors__ = ["kratos", "phi", "scienta", "specs", "unkwown"]
    __prmt_vndr_cls: Dict[str, Dict] = {
        ".csv": {"unknown": CsvMapperVamasResult},
        ".h5": {"scienta": MapperScienta},
        ".hdf5": {"scienta": MapperScienta},
        ".ibw": {"scienta": MapperScienta},
        ".npl": {"unkwown": VamasMapper},
        ".pro": {"phi": MapperPhi},
        ".spe": {"phi": MapperPhi},
        ".sle": {"specs": SleMapperSpecs},
        ".txt": {
            "scienta": MapperScienta,
            "unknown": TxtMapperVamasExport,
        },
        ".vms": {"unkwown": VamasMapper},
        ".xml": {"specs": XmlMapperSpecs},
        ".xy": {"specs": XyMapperSpecs},
    }

    __file_err_msg__ = (
        "Need an XPS data file with one of the following extensions: "
        f"data files: {__prmt_file_ext__}, metadata files: {__prmt_metadata_file_ext__}."
    )

    __vndr_err_msg__ = (
        f"Need an XPS data file from one of the following vendors: {__vendors__}"
    )

    def __init__(self, config_file: Optional[str] = None, *args, **kwargs):
        super().__init__(config_file, *args, **kwargs)

        self.xps_data_dicts: List[Dict[str, Any]] = []
        self.xps_data: Dict[str, Any] = {}
        self.eln_data: Dict[str, Any] = {}

        self.extensions = {
            ".yml": self.handle_eln_file,
            ".yaml": self.handle_eln_file,
            ".json": self.set_config_file,
        }

        self.processing_order = (
            XPSReader.__prmt_file_ext__
            + list(XPSReader.__prmt_metadata_file_ext__.keys())
            + list(self.extensions.keys())
        )

        for ext in XPSReader.__prmt_file_ext__ + list(
            XPSReader.__prmt_metadata_file_ext__.keys()
        ):
            self.extensions[ext] = self.handle_data_file

        self.processing_order = XPSReader.__prmt_file_ext__ + [
            ".yml",
            ".yaml",
            ".json",
        ]

    def set_config_file(
        self, file_path: Optional[Union[str, Path]], replace: bool = True
    ) -> Dict[str, Any]:
        if not file_path:
            return {}

        if replace:
            if self.config_file is not None:
                if file_path != self.config_file:
                    logger.info(
                        f"Config file already set. Replaced by the new file {file_path}."
                    )
            self.config_file = file_path
        else:
            if self.config_file is None:
                self.config_file = file_path

        return {}

    def handle_eln_file(self, file_path: str) -> Dict[str, Any]:
        """
        Loads ELN file and handles specific cases.
        """

        def combine_and_unique_string(string: str, elements: List[str]) -> str:
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
            convert_dict=CONVERT_DICT,
            replace_nested=REPLACE_NESTED,
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
                                self.eln_data[modified_key] = combine_and_unique_string(
                                    self.eln_data[modified_key], atom_types
                                )
                        else:
                            logger.info(
                                f"{key} from ELN was not parsed to atom_types because {modified_key} already exists."
                            )

            if isinstance(value, datetime.datetime):
                eln_data[key] = value.isoformat()

            self.eln_data[new_key] = eln_data.pop(key)

        return {}

    def handle_data_file(self, file_path: str) -> Dict[str, Any]:
        def _check_for_vendors(file_path: str) -> str:
            """
            Check for the vendor name of the XPS data file.

            """
            _, file_ext = os.path.splitext(file_path)

            vendor_dict = XPSReader.__prmt_vndr_cls[file_ext]

            if len(vendor_dict) == 1:
                return list(vendor_dict.keys())[0]
            if file_ext == ".txt":
                return _check_for_vendors_txt(file_path)
            raise ValueError(XPSReader.__vndr_err_msg__)

        def _check_for_vendors_txt(file_path: str) -> str:
            """
            Search for a vendor names in a txt file

            Parameters
            ----------
            file : str
                XPS txt file.

            Returns
            -------
            vendor
                Vendor name if that name is in the txt file.

            """
            vendor_dict = XPSReader.__prmt_vndr_cls[".txt"]

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

        if file_ext in XPSReader.__prmt_file_ext__:
            vendor = _check_for_vendors(file_path)

            parser = XPSReader.__prmt_vndr_cls[file_ext][vendor]()
            parser.parse_file(file_path, **self.kwargs)
            data_dict = parser.data_dict

            config_file = parser.config_file

            if isinstance(config_file, dict):
                config_file = config_file.get(file_ext)

            self.set_config_file(
                XPSReader.reader_dir.joinpath("config", config_file),
                replace=False,
            )

            self.xps_data_dicts += [parser.data_dict]

        elif file_ext in XPSReader.__prmt_metadata_file_ext__:
            vendor = _check_for_vendors(file_path)

            metadata_parser = XPSReader.__prmt_vndr_cls[file_ext][vendor]()
            metadata_parser.parse_file(file_path, **self.kwargs)

            main_file_ext = XPSReader.__prmt_metadata_file_ext__[file_ext]

            main_file_dicts = [
                d for d in self.xps_data_dicts if d.get("file_ext") == main_file_ext
            ]

            metadata_parser.update_main_file_dict(main_file_dicts)

        return {}

    def get_entry_names(self) -> List[str]:
        """
        Returns a list of entry names which should be constructed from the data.
        Defaults to creating a single entry named "entry".
        """
        # Track entries for using for eln data
        entries: List[str] = []

        try:
            for entry in self.xps_data["data"]:
                entries += [entry]
        except KeyError:
            pass

        if not entries:
            entries += ["entry"]

        return list(dict.fromkeys(entries))

    def setup_template(self) -> Dict[str, Any]:
        """
        Setups the initial data in the template.
        """
        # TODO: Set fixed information, e.g., about the reader.
        return {}

    def handle_objects(self, objects: Tuple[Any]) -> Dict[str, Any]:
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

        Returns
        -------
        None

        """

        def check_for_same_entries(
            dicts: List[Dict[str, Any]],
        ) -> Optional[Tuple[Set[str], List[Set[int]]]]:
            """
            Checks whether the input dictionaries have the same entries and identifies which
            dictionaries contain each common entry.

            Args:
                dicts (List[Dict[str, any]]): A list of dictionaries with keys potentially containing entries.

            Returns:
                Optional[Tuple[Set[str], List[Set[int]]]]:
                - A tuple where:
                  - The first element is a set of common entries found across the dictionaries.
                  - The second element is a list of sets, where each set contains the indices of dictionaries
                    that have the corresponding entry.
                - None if no common entries are found.
            """
            entry_pattern = re.compile(r"/ENTRY\[(.*?)\]")

            entry_to_dicts: Dict[str, Set] = {}

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

        common_entries, dict_indices = check_for_same_entries(self.xps_data_dicts)

        if common_entries and not self.overwrite_keys:
            for entry, indices in zip(common_entries, dict_indices):
                dicts_with_common_entries = [self.xps_data_dicts[i] for i in indices]

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

        for data_dict in self.xps_data_dicts:
            # If there are multiple input data files of the same type,
            # make sure that existing keys are not overwritten.
            existing = [
                (key, self.xps_data[key], data_dict[key])
                for key in set(self.xps_data).intersection(data_dict)
            ]

            self.xps_data = {**self.xps_data, **data_dict}

            if not self.overwrite_keys:
                for key, value1, value2 in existing:
                    self.xps_data[key] = concatenate_values(value1, value2)

    def _get_analyzer_names(self) -> List[str]:
        """
        Returns a list of analyzer names which should be constructed
        from the data. Defaults to creating a single analyer named
        "analyzer".

        Currently, this is not used, but can be changed if there are
        multiple analyzers in the future.
        """
        analyzers: List[str] = []

        if not analyzers:
            analyzers += ["electronanalyzer"]

        return list(dict.fromkeys(analyzers))

    def _get_detector_names(self) -> List[str]:
        """
        Returns a list of detector names which should be constructed
        from the data. Defaults to creating a single detector named
        "detector".
        """
        detectors: List[str] = []

        try:
            for entry, entry_values in self.xps_data["data"].items():
                for data_var in entry_values:
                    if CHAN_COUNT in data_var:
                        detector_num = data_var.split(CHAN_COUNT)[-1]
                        detector_nm = f"detector{detector_num}"
                        detectors += [detector_nm]
        except KeyError:
            pass

        if not detectors:
            detectors += ["detector"]

        return list(dict.fromkeys(detectors))

    def process_multiple_entities(self) -> None:
        """
        Check if there are multiple of some class and, if so, change the
        keys and values in the config file.

        This replaces all occureces of "detector" and "electronanalyzer"
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

    def get_metadata(
        self,
        metadata_dict: Dict[str, Any],
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
        metadata_dict : Dict[str, Any]
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
        return self.get_metadata(self.xps_data, path, self.callbacks.entry_name)

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

    def get_data_dims(self, key: str, path: str) -> List[str]:
        """
        Returns the dimensions of the data from the given path.
        """
        entry = self.callbacks.entry_name
        escaped_entry = re.escape(entry)
        xr_data = self.xps_data["data"].get(entry)

        def get_signals(key: str) -> List[str]:
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

        def get_all_keys(template_key: str) -> List[str]:
            pattern = re.compile(rf"^/ENTRY\[{escaped_entry}]/{template_key}([^/]+)")

            keys = set(
                match[1] for key in self.xps_data if (match := pattern.search(key))
            )

            return sorted(keys)

        if isinstance(path, str) and path.endswith(("*.external", "*.external_unit")):
            data_func = lambda: get_all_keys("external_")
        else:
            if isinstance(path, str) and path.endswith(".unit"):
                data_func = lambda: list(xr_data.coords)
            else:
                data_func = lambda: get_signals(path.split(":*.")[-1])

        PatternHandler = Dict[Pattern[str], Callable[[], Any]]

        patterns: dict[re.Pattern, Callable[[], Any]] = {
            re.compile(r"peak"): lambda: get_all_keys("component"),
            re.compile(r"background"): lambda: get_all_keys("region"),
            re.compile(
                r"DATA\[[^\]]+\]/(?:DATA\[[^\]]+\]|AXISNAME\[[^\]]+\])(?:/@units)?"
            ): data_func,
            re.compile(
                r"ELECTRON_DETECTOR\[[a-zA-Z0-9_]+\]/raw_data"
            ): lambda: get_signals("channels"),
        }

        # Function to match and return the handler result
        for pattern, handler in patterns.items():
            if pattern.search(key):
                return handler()

        return get_signals(key="scans")

    def _search_first(self, data: Dict[str, Any], pattern: re.Pattern) -> Optional[Any]:
        """
        Search for the first value in a dictionary whose key matches a regex pattern.

        Parameters:
            data (dict[str, Any]): The dictionary to search.
            pattern (re.Pattern): The compiled regex pattern to search for in keys.

        Returns:
            The first matching value, or None if no match is found.
        """
        return next((value for key, value in data.items() if pattern.search(key)), None)

    def get_data(self, key: str, path: str) -> Optional[Any]:
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
        xr_data = self.xps_data["data"].get(entry)

        # Average or errors
        if path.startswith(("average", "errors")):
            data = [xr_data[var].data for var in _get_scan_vars(xr_data.data_vars)]
            if path.endswith("reduced"):
                try:
                    data = [np.sum(arr, axis=tuple(range(1, arr.ndim))) for arr in data]
                except np.AxisError:
                    return None
            stats_func = np.mean if path.startswith("average") else np.std
            return stats_func(data, axis=0)

        # Raw data
        if path.endswith("raw_data"):
            data_vars = _get_channel_vars(xr_data.data_vars) or _get_scan_vars(
                xr_data.data_vars
            )
            return np.array(
                [xr_data[var].data for var in data_vars if SCAN_COUNT in var]
            )

        # Channels or scans by suffix
        for suffix, extractor in {
            ".scans": lambda name: xr_data.get(name),
            ".channels": lambda name: xr_data.get(name),
        }.items():
            if path.endswith(suffix):
                name = path.split(suffix)[0]
                data = extractor(name)
                return np.array(data) if data is not None else None

        # Axis or coordinate data
        if path.endswith(".axes"):
            axis = path.split(".axes")[0]
            coord = xr_data.coords.get(axis)
            return np.array(coord.values) if coord is not None else None

        # Units for internal channels
        if path.endswith(".unit"):
            channel = path.split(".unit")[0]
            pattern = re.compile(
                rf"^/ENTRY\[{escaped_entry}]/{re.escape(channel)}/@units"
            )
            return self._search_first(self.xps_data, pattern)

        # External channels and units
        if path.endswith((".external", ".external_unit")):
            channel = path.split(".external")[0]
            if path.endswith("_unit"):
                pattern = re.compile(
                    rf"^/ENTRY\[{escaped_entry}]/external_{re.escape(channel)}/@units"
                )
                return self._search_first(self.xps_data, pattern)
            else:
                pattern = re.compile(
                    rf"^/ENTRY\[{escaped_entry}]/external_{re.escape(channel)}$"
                )
                matches = [
                    value for key, value in self.xps_data.items() if pattern.search(key)
                ]
                return np.array(matches).squeeze() if matches else None

        # Default: try direct data access
        try:
            return xr_data[path]
        except KeyError:
            try:
                print(key, path)
                return np.array(xr_data.coords[path].values)
            except KeyError:
                return None

    def set_nxdata_defaults(self, template):
        """Set the default for automatic plotting."""
        survey_count, count = 0, 0

        def get_unique_nxfit_names(template) -> Set[str]:
            """Extract unique 'ENTRY[<some-name>]/FIT[<some-other-name>]' pairs from template keys."""
            pattern = re.compile(r"^/?ENTRY\[(?P<entry>[^]]+)\]/FIT\[(?P<fit>[^]]+)\]/")

            result = set()
            for key in template:
                m = pattern.match(key)
                if m:
                    result.add(f"{m.group('entry')}/{m.group('fit')}")
            return result

        def get_first_matching_fit(
            entry_name: str, unique_fits: Set[str]
        ) -> Optional[str]:
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
        file_paths: Tuple[str] = None,
        objects: Tuple[Any] = None,
        **kwargs,
    ) -> dict:
        self.overwrite_keys = _check_multiple_extensions(file_paths)
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

        return final_template


READER = XPSReader
