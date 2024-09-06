"""
A generic reader for loading XPS (X-ray Photoelectron Spectroscopy) data
file into mpes nxdl (NeXus Definition Language) template.
"""
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

import os
import sys
import re
import datetime
import copy
import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional, Union, Set
import numpy as np

from pynxtools.dataconverter.helpers import extract_atom_types
from pynxtools.dataconverter.readers.multi.reader import (
    MultiFormatReader,
)
from pynxtools.dataconverter.readers.utils import parse_yml
from pynxtools.dataconverter.template import Template

from pynxtools_xps.phi.spe_pro_phi import MapperPhi
from pynxtools_xps.scienta.scienta_reader import MapperScienta
from pynxtools_xps.specs.sle.sle_specs import SleMapperSpecs
from pynxtools_xps.specs.xy.xy_specs import XyMapperSpecs
from pynxtools_xps.specs.xml.xml_specs import XmlMapperSpecs
from pynxtools_xps.vms.txt_vamas_export import TxtMapperVamasExport
from pynxtools_xps.vms.vamas import VamasMapper

from pynxtools_xps.reader_utils import check_units


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


np.set_printoptions(threshold=sys.maxsize)


CONVERT_DICT = {
    "unit": "@units",
    "version": "@version",
    "user": "USER[user]",
    "instrument": "INSTRUMENT[instrument]",
    "source_xray": "sourceTYPE[source_xray]",
    "beam_xray": "beamTYPE[beam_xray]",
    "analyser": "ELECTRONANALYSER[electronanalyser]",
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


# pylint: disable=too-few-public-methods
class XPSReader(MultiFormatReader):
    """Reader for XPS."""

    supported_nxdls = [
        "NXmpes",
        # "NXxps",
    ]

    reader_dir = Path(__file__).parent
    config_file: Optional[Union[str, Path]] = reader_dir.joinpath(
        "config", "template.json"
    )

    __prmt_file_ext__ = [
        ".ibw",
        ".npl",
        ".pro",
        ".spe",
        ".sle",
        ".slh",
        ".txt",
        ".vms",
        ".xml",
        ".xy",
    ]
    __vendors__ = ["kratos", "phi", "scienta", "specs", "unkwown"]
    __prmt_vndr_cls: Dict[str, Dict] = {
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
        f"{__prmt_file_ext__}"
    )

    __vndr_err_msg__ = (
        "Need an XPS data file from one of the following vendors: " f"{__vendors__}"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.xps_data_dicts: List[Dict[str, Any]] = []
        self.xps_data: Dict[str, Any] = {}
        self.eln_data: Dict[str, Any] = {}

        self.extensions = {
            ".yml": self.handle_eln_file,
            ".yaml": self.handle_eln_file,
            ".json": self.set_config_file,
        }

        for ext in XPSReader.__prmt_file_ext__:
            self.extensions[ext] = self.handle_data_file

    def set_config_file(self, file_path: str) -> Dict[str, Any]:
        if self.config_file is not None:
            logger.info(
                f"Config file already set. Replaced by the new file {file_path}."
            )
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
            return None

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
                    # This is for picking the Scienta reader is "scienta"
                    # is not in the file
                    return vendor
            return "unknown"

        _, file_ext = os.path.splitext(file_path)

        if file_ext in XPSReader.__prmt_file_ext__:
            vendor = _check_for_vendors(file_path)
            try:
                parser = XPSReader.__prmt_vndr_cls[file_ext][vendor]()

                parser.parse_file(file_path, **self.kwargs)
                self.config_file = XPSReader.reader_dir.joinpath(
                    "config", parser.config_file
                )
                data_dict = parser.data_dict

            except ValueError as val_err:
                raise ValueError(XPSReader.__vndr_err_msg__) from val_err
            except KeyError as key_err:
                raise KeyError(XPSReader.__vndr_err_msg__) from key_err
        else:
            raise ValueError(XPSReader.__file_err_msg__)

        self.xps_data_dicts += [data_dict]

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

        if common_entries:
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
            for key, value1, value2 in existing:
                self.xps_data[key] = concatenate_values(value1, value2)

    def _get_analyser_names(self) -> List[str]:
        """
        Returns a list of analyser names which should be constructed
        from the data. Defaults to creating a single analyser named
        "analyser".

        Currently, this is not used, but can be changed if there are
        multiple analysers in the future.
        """
        analysers: List[str] = []

        if not analysers:
            analysers += ["electronanalyser"]

        return list(dict.fromkeys(analysers))

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

        This replaces all occureces of "detector" and "electronanalyser"
        in the config dict by the respective names (e.g., detector0, detector1)
        and removes the generic term if there are multiple different instances.

        """
        multiples_to_check = {
            "electronanalyser": self._get_analyser_names,
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
        /ENTRY[my-entry]/instrument/analyser/collectioncolumn/voltage.
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

        if value is None or str(value) == "None":
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

        def get_signals(key: str) -> List[str]:
            xr_data = self.xps_data["data"].get(f"{self.callbacks.entry_name}")

            if key == "scans":
                data_vars = _get_scan_vars(xr_data.data_vars)
            elif key == "channels":
                data_vars = _get_channel_vars(xr_data.data_vars)
                if not data_vars:
                    data_vars = _get_scan_vars(xr_data.data_vars)
            else:
                data_vars = [""]

            return list(map(str, data_vars))

        if path.startswith("@data:*"):
            return get_signals(key=path.split(":*.")[-1])

        if any(x in path for x in ["counts", "raw/@units"]):
            return get_signals(key="channels")

        return get_signals(key="scans")

    def get_data(self, key: str, path: str) -> Any:
        """
        Returns data for a given key.
        Can either return averaged, scan, or channel data.
        Should return None if the path does not exist.
        """
        xr_data = self.xps_data["data"].get(f"{self.callbacks.entry_name}")

        if path.endswith("average"):
            return np.mean(
                [xr_data[x_arr].data for x_arr in _get_scan_vars(xr_data.data_vars)],
                axis=0,
            )

        elif path.endswith("errors"):
            return np.std(
                [xr_data[x_arr].data for x_arr in _get_scan_vars(xr_data.data_vars)],
                axis=0,
            )

        elif path.endswith("raw_data"):
            data_vars = _get_channel_vars(xr_data.data_vars)

            if not data_vars:
                # If there is no channel data, use scan data.
                data_vars = _get_scan_vars(xr_data.data_vars)

            # Skip average cycle data
            return np.array(
                [
                    xr_data[data_var].data
                    for data_var in data_vars
                    if SCAN_COUNT in data_var
                ]
            )

        elif path.endswith("scans"):
            return np.array(xr_data[path.split(".scans")[0]])

        elif path.endswith("channels"):
            return np.array(xr_data[path.split(".channels")[0]])

        elif "energy" in path:
            return np.array(xr_data.coords["energy"].values)

        else:
            try:
                return xr_data[path]
            except KeyError:
                try:
                    return np.array(xr_data.coords[path].values)
                except KeyError:
                    pass

    def set_root_default(self, template):
        """Set the default for automatic plotting."""
        survey_count_ = 0
        count = 0

        for entry in self.get_entry_names():
            if "Survey" in entry and survey_count_ == 0:
                survey_count_ += 1
                template["/@default"] = entry

            # If no Survey set any scan for default
            if survey_count_ == 0 and count == 0:
                count += 1
                template["/@default"] = entry

    def read(
        self,
        template: dict = None,
        file_paths: Tuple[str] = None,
        objects: Tuple[Any] = None,
        **kwargs,
    ) -> dict:
        template = super().read(template, file_paths, objects, suppress_warning=True)
        self.set_root_default(template)

        final_template = Template()
        for key, val in template.items():
            if val is not None:
                if "@units" in key:
                    check_units(key, val)
                final_template[key] = val

        return final_template


READER = XPSReader
