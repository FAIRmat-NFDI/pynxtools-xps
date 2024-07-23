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
from typing import Any, Dict, List, Set, Tuple, Optional, Union
import numpy as np

from pynxtools.dataconverter.helpers import extract_atom_types
from pynxtools.dataconverter.readers.multi.reader import (
    MultiFormatReader,
    ParseJsonCallbacks,
)
from pynxtools.dataconverter.readers.utils import parse_yml
from pynxtools.dataconverter.template import Template
from pynxtools.dataconverter.readers.utils import (
    FlattenSettings,
    flatten_and_replace,
    parse_flatten_json,
)

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
    config_file: Optional[str] = reader_dir.joinpath("config", "template.json")

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
        "Need an XPS data file with the following extension: " f"{__prmt_file_ext__}"
    )

    __vndr_err_msg__ = (
        "Need an XPS data file from the following vendors: " f"{__vendors__}"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.xps_data: Dict[str, Any] = {}
        self.eln_data: Dict[str, Any] = {}

        self.extensions = {
            ".yml": self.handle_eln_file,
            ".yaml": self.handle_eln_file,
            ".json": self.set_config_file,
        }

        for ext in XPSReader.__prmt_file_ext__:
            self.extensions[ext] = self.handle_data_file

        self.callbacks = ParseJsonCallbacks(
            attrs_callback=self.get_attr,
            data_callback=self.get_data,
            eln_callback=self.get_eln_data,
        )

    def set_config_file(self, file_path: str) -> Dict[str, Any]:
        if self.config_file is not None:
            logger.info(
                f"Config file already set. Replaced by the new file {file_path}."
            )
        self.config_file = file_path
        return {}

    def handle_eln_file(self, file_path: str) -> Dict[str, Any]:
        eln_data = parse_yml(
            file_path,
            convert_dict=CONVERT_DICT,
            replace_nested=REPLACE_NESTED,
        )

        for key, value in eln_data.items():
            if "molecular_formula_hill" in key and "atom_types" not in key:
                atom_types: List = []
                atom_types = list(extract_atom_types(value))

                if atom_types:
                    modified_key = key.replace("chemical_formula", "atom_types")
                    eln_data[modified_key] = ", ".join(atom_types)

        # replace paths for entry-specific ELN data
        pattern = re.compile(r"/ENTRY\[entry\]/(ENTRY\[[^\]]+\]/.*)")

        for key in eln_data.copy():
            match = pattern.search(key)
            if match:
                new_key = match.group(1)
                self.eln_data[new_key] = eln_data.pop(key)
            else:
                self.eln_data[key] = eln_data.pop(key)

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
            vendor_dict = XPSReader.__prmt_vndr_cls["txt"]

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

        # If there are multiple input data files of the same type,
        # make sure that existing keys are not overwritten.
        existing = [
            (key, self.xps_data[key], data_dict[key])
            for key in set(self.xps_data).intersection(data_dict)
        ]

        self.xps_data = {**self.xps_data, **data_dict}
        for key, value1, value2 in existing:
            self.xps_data[key] = concatenate_values(value1, value2)

        return {}

    def get_entry_names(self) -> List[str]:
        """
        Returns a list of entry names which should be constructed from the data.
        Defaults to creating a single entry named "entry".
        """
        # Track entries for using for eln data
        entry_set: Set[str] = set()

        try:
            for entry, entry_values in self.xps_data["data"].items():
                entry_set.add(entry)
        except KeyError:
            pass

        if not entry_set:
            entry_set.add("entry")

        return list(entry_set)

    def setup_template(self) -> Dict[str, Any]:
        """
        Setups the initial data in the template.
        This may be used to set fixed information, e.g., about the reader.
        """
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
        pass  # self.process_multiple_entities()

    def _get_analyser_names(self) -> List[str]:
        """
        Returns a list of analyser names which should be constructed
        from the data. Defaults to creating a single analyser named
        "analyser".

        Currently, this is not used, but can be changed if there are
        multiple analysers in the future.
        """
        analyser_set: Set[str] = set()

        if not analyser_set:
            analyser_set.add("electronanalyser")

        return list(analyser_set)

    def _get_detector_names(self) -> List[str]:
        """
        Returns a list of detector names which should be constructed
        from the data. Defaults to creating a single detector named
        "detector".
        """
        detector_set: Set[str] = set()

        try:
            for entry, entry_values in self.xps_data["data"].items():
                for data_var in entry_values:
                    chan_count = "_chan"
                    if chan_count in data_var:
                        detector_num = data_var.split(chan_count)[-1]
                        detector_nm = f"detector{detector_num}"
                        detector_set.add(detector_nm)
        except KeyError:
            pass

        if not detector_set:
            detector_set.add("detector")

        return list(set(["detector"]))

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

    def get_attr(self, path: str) -> Any:
        """
        Get the metadata that was stored in the main file.
        """
        return self.get_metadata(self.xps_data, path, self.callbacks.entry_name)

    def get_eln_data(self, path: str) -> Any:
        """
        Returns data from the given eln path.
        Gives preference to ELN data for a given entry before searching
        the ELN data for all entries.
        Returns None if the path does not exist.
        """
        for entry_name in [self.callbacks.entry_name, "entry"]:
            value = self.get_metadata(self.eln_data, path, entry_name)
            if value:
                return value
        return value

    def get_data(self, path: str) -> Any:
        """
        Returns data for a given key.
        Can either return averaged, scan, or channel data.
        Should return None if the path does not exist.
        """
        xr_data = self.xps_data["data"].get(f"{self.callbacks.entry_name}")

        chan_count = "_chan"
        scan_count = "_scan"

        if path.endswith("average"):
            return np.mean(
                [
                    xr_data[x_arr].data
                    for x_arr in xr_data.data_vars
                    if (chan_count not in x_arr and scan_count in x_arr)
                ],
                axis=0,
            )

        elif path.endswith("errors"):
            return np.std(
                [
                    xr_data[x_arr].data
                    for x_arr in xr_data.data_vars
                    if (chan_count not in x_arr and scan_count in x_arr)
                ],
                axis=0,
            )

        elif path.endswith("scans"):
            scan_data: List[np.ndarray] = []
            for data_var in xr_data.data_vars:
                # Collecting only accumulated counts
                # individual channeltron counts go to detector data section
                if chan_count not in data_var and scan_count in data_var:
                    self.callbacks.dim += [str(data_var)]
                    scan_data += [xr_data[data_var].data]

            return scan_data

        elif path.endswith("raw_data"):
            raw_data: List[np.ndarray] = []

            data_vars = [
                data_var for data_var in xr_data.data_vars if chan_count in data_var
            ]

            if not data_vars:
                # If there is no channel data, use scan data.
                data_vars = [
                    data_var for data_var in xr_data.data_vars if scan_count in data_var
                ]

            raw_data = [xr_data[data_var].data for data_var in data_vars]

            return np.array(raw_data)

        elif path.endswith("channels"):
            detector_data: List[np.ndarray] = []

            # Iteration over channels
            data_vars = [
                data_var for data_var in xr_data.data_vars if chan_count in data_var
            ]

            if data_vars:
                for data_var in data_vars:
                    self.callbacks.dim += [str(data_var)]
                    detector_data += [xr_data[data_var].data]

            else:
                # If there is no channel data, iterate over scans
                data_vars = [data_var for data_var in xr_data.data_vars]
                if scan_count in data_var:
                    self.callbacks.dim += [str(data_var)]
                    detector_data += [xr_data[data_var].data]

            return detector_data

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
        template = super().read(template, file_paths, objects, **kwargs)
        self.set_root_default(template)

        final_template = Template()
        for key, val in template.items():
            if val is not None:
                if "@units" in key:
                    pass
                    check_units(key, val)
                final_template[key] = val

        return final_template


READER = XPSReader
