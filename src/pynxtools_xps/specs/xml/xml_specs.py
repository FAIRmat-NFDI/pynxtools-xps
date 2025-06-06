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
Specs Lab Prodigy XML export, to be passed to mpes nxdl
(NeXus Definition Language) template.
"""

import re
import xml.etree.ElementTree as EmtT
from typing import Tuple, List, Dict, Any
import copy
import xarray as xr
import numpy as np

from pynxtools_xps.reader_utils import XPSMapper, align_name_part
from pynxtools_xps.value_mappers import convert_energy_scan_mode

KEY_PARTS = ["RegionGroup_", "RegionData_"]

KEY_PATTERNS = [re.compile(rf"{key_part}(.*?)(?=\/|$)") for key_part in KEY_PARTS]


def construct_entry_name(key: str) -> str:
    """Construct entry name."""
    name_parts = []

    for key_pattern in KEY_PATTERNS:
        match = re.search(key_pattern, key)
        if match:
            name_part = align_name_part(match.group(1))
            name_parts.append(name_part)
    return "__".join(name_parts)


VALUE_MAP = {
    "scan_mode/name": convert_energy_scan_mode,
}


class XmlMapperSpecs(XPSMapper):
    """
    Class for restructuring xml data file from
    specs vendor into python dictionary.
    """

    config_file = "config_specs_xml.json"

    def _select_parser(self):
        """
        Select Specs XML parser.
        Currently, there is only one parser.
        Returns
        -------
        XmlParserSpecs
            Parser for reading .xml file exported by SpecsLab2.

        """
        return XmlParserSpecs()

    def parse_file(self, file, **kwargs):
        """
        Parse the file using the parser that fits the Prodigy SLE version.
        Returns flat list of dictionaries containing one spectrum each.

        """
        self.file = file
        self.parser = self._select_parser()
        metadata_dict, self.raw_data = self.parser.parse_file(file, **kwargs)

        self._xps_dict = {**self._xps_dict, **metadata_dict}

        self.construct_data()

        return self.data_dict

    # pylint: disable=too-many-locals,too-many-statements, too-many-branches
    def construct_data(self):
        """
        Construct the Binding Energy and separate the counts for
        different detectors and finally sum up all the counts for
        to find total electron counts.
        """
        spectra = copy.deepcopy(self.raw_data)
        self._xps_dict["data"]: dict = {}

        for entry, spectrum in spectra.items():
            raw_data = spectrum["raw_data"]
            mcd_num = int(raw_data["mcd_num"])

            curves_per_scan = raw_data["curves_per_scan"]
            values_per_curve = raw_data["values_per_curve"]
            values_per_scan = int(curves_per_scan * values_per_curve)
            mcd_head = int(raw_data["mcd_head"])
            mcd_tail = int(raw_data["mcd_tail"])
            excitation_energy = raw_data["excitation_energy"]
            scan_mode = raw_data["scan_mode"]
            kinetic_energy = raw_data["kinetic_energy"]
            scan_delta = raw_data["scan_delta"]
            pass_energy = raw_data["pass_energy"]
            kinetic_energy_base = raw_data["kinetic_energy_base"]
            # Adding one unit to the binding_energy_upper is added as
            # electron comes out if energy is one unit higher
            binding_energy_upper = (
                excitation_energy - kinetic_energy + kinetic_energy_base + 1
            )

            mcd_energy_shifts = raw_data["mcd_shifts"]
            mcd_energy_offsets = []
            offset_ids = []

            # consider offset values for detector with respect to
            # position at +16 which is usually large and positive value
            for mcd_shift in mcd_energy_shifts:
                mcd_energy_offset = (mcd_energy_shifts[-1] - mcd_shift) * pass_energy
                mcd_energy_offsets.append(mcd_energy_offset)
                offset_id = round(mcd_energy_offset / scan_delta)
                offset_ids.append(int(offset_id - 1 if offset_id > 0 else offset_id))

            # Skiping entry without count data
            if not mcd_energy_offsets:
                continue
            mcd_energy_offsets = np.array(mcd_energy_offsets)
            # Putting energy of the last detector as a highest energy
            starting_eng_pnts = binding_energy_upper - mcd_energy_offsets
            ending_eng_pnts = starting_eng_pnts - values_per_scan * scan_delta

            channeltron_eng_axes = np.zeros((mcd_num, values_per_scan))
            for ind in np.arange(len(channeltron_eng_axes)):
                channeltron_eng_axes[ind, :] = np.linspace(
                    starting_eng_pnts[ind], ending_eng_pnts[ind], values_per_scan
                )

            channeltron_eng_axes = np.round(channeltron_eng_axes, decimals=8)
            # construct ultimate or incorporated energy axis from
            # lower to higher energy
            scans = list(raw_data["scans"].keys())

            # Check whether array is empty or not
            if not scans:
                continue
            if not raw_data["scans"][scans[0]].any():
                continue
            # Sorting in descending order
            binding_energy = channeltron_eng_axes[-1, :]

            self._xps_dict["data"][entry] = xr.Dataset()

            for scan_nm in scans:
                channel_counts = np.zeros((mcd_num + 1, values_per_scan))
                # values for scan_nm corresponds to the data for each
                # "scan" in individual CountsSeq
                scan_counts = raw_data["scans"][scan_nm]

                if scan_mode == "fixed_analyzer_transmission":
                    for row in np.arange(mcd_num):
                        count_on_row = scan_counts[row::mcd_num]
                        # Reverse counts from lower to higher
                        # BE as in BE_eng_axis
                        count_on_row = count_on_row[mcd_head:-mcd_tail]

                        channel_counts[row + 1, :] = count_on_row
                        channel_counts[0, :] += count_on_row

                        # Storing detector's raw counts
                        self._xps_dict["data"][entry][f"{scan_nm}_chan{row}"] = (
                            xr.DataArray(
                                data=channel_counts[row + 1, :],
                                coords={"energy": binding_energy},
                            )
                        )

                        # Storing callibrated and after accumulated each scan counts
                        if row == mcd_num - 1:
                            self._xps_dict["data"][entry][scan_nm] = xr.DataArray(
                                data=channel_counts[0, :],
                                coords={"energy": binding_energy},
                            )
                else:
                    for row in np.arange(mcd_num):
                        start_id = offset_ids[row]
                        count_on_row = scan_counts[start_id::mcd_num]
                        count_on_row = count_on_row[0:values_per_scan]
                        channel_counts[row + 1, :] = count_on_row

                        # shifting and adding all the curves.
                        channel_counts[0, :] += count_on_row

                        # Storing detector's raw counts
                        self._xps_dict["data"][entry][f"{scan_nm}_chan{row}"] = (
                            xr.DataArray(
                                data=channel_counts[row + 1, :],
                                coords={"energy": binding_energy},
                            )
                        )

                        # Storing callibrated and after accumulated each scan counts
                        if row == mcd_num - 1:
                            self._xps_dict["data"][entry][scan_nm] = xr.DataArray(
                                data=channel_counts[0, :],
                                coords={"energy": binding_energy},
                            )

            # Write averaged cycle data to 'data'.
            averaged_scans = {
                key: value
                for key, value in self._xps_dict["data"][entry].items()
                if "_chan" not in key
            }

            average_cycles = {}

            for key, value in averaged_scans.items():
                cycle_number = re.findall(r"cycle(\d+)_", key)[0]
                cycle_key = f"cycle{cycle_number}"
                if cycle_number not in average_cycles:
                    average_cycles[cycle_key] = value
                else:
                    average_cycles[cycle_key] += value

            for cycle, value in average_cycles.items():
                average_cycles[cycle] = value / len(
                    [k for k in averaged_scans.keys() if cycle in k]
                )

            for cycle, value in average_cycles.items():
                self._xps_dict["data"][entry][cycle] = xr.DataArray(
                    data=value,
                    coords={"energy": binding_energy},
                )


class XmlParserSpecs:
    """Parser for SpecsLab2 XML data"""

    def __init__(self) -> None:
        self.metadata_dict: dict = {}
        self.entry_to_data: dict = {}
        self._root_path = "/ENTRY"
        self.tail_part_frm_struct = ""
        self.tail_part_frm_othr = ""
        self.child_nm_reslvers = "__child_name_resolver__"

    def parse_file(self, file: str, **kwargs):
        """Start parsing process and parse children recursively.

        Parameters
        ----------
        """
        root_element = EmtT.parse(file).getroot()
        root_element.attrib[self.child_nm_reslvers] = []  # type: ignore[assignment]
        child_num = len(root_element)
        parent_path = self._root_path
        skip_child = -1

        child_elmt_ind = 0
        while child_num > 0:
            self.pass_child_through_parsers(
                root_element, parent_path, child_elmt_ind, skip_child
            )

            child_num -= 1
            child_elmt_ind += 1

        metadata_dict: Dict[str, Any] = {}

        for key, value in self.metadata_dict.items():
            entry_name = construct_entry_name(key)
            str_to_replace = "/"
            for i, key_pattern in enumerate(KEY_PATTERNS):
                match = re.search(key_pattern, key)
                if match:
                    if KEY_PARTS[i] == "RegionData_":
                        str_to_replace += "regions/"
                    str_to_replace += f"{KEY_PARTS[i]}{match.group(1)}/"
            modified_key = key.replace(str_to_replace, f"[{entry_name}]/")

            for suffix, func in VALUE_MAP.items():
                if modified_key.endswith(suffix):
                    value = func(value)

            metadata_dict[modified_key] = value

        self.metadata_dict = metadata_dict

        self.collect_raw_data_to_construct_data()

        return self.metadata_dict, self.entry_to_data

    def pass_child_through_parsers(
        self,
        element_: EmtT.Element,
        parent_path: str,
        child_elmt_ind: int,
        skip_child: int,
    ) -> None:
        """
        Parse the element to parser according to element tag.
        Parameters
        ----------
        element_ : xml element to parse
        parent_path : Xpath of the parent element where the element_ belongs
        child_elmt_ind : Index of child element to track the children.
        skip_child : Tack the children who will be skipped to pass to
                     the parser
        Returns
        -------
        None
        """

        name_val_elmt_tag = ["ulong", "double", "string", "boolean", "enum", "any"]

        parent_element = element_
        element = parent_element[child_elmt_ind]
        element.attrib["__parent__"] = parent_element  # type: ignore[assignment]
        element.attrib["__odr_siblings__"] = child_elmt_ind  # type: ignore[assignment]

        if self.child_nm_reslvers not in parent_element.attrib.keys():
            parent_element.attrib[self.child_nm_reslvers] = []  # type: ignore[assignment]

        elmt_tag = element.tag

        if child_elmt_ind <= skip_child:
            pass

        elif elmt_tag == "sequence":
            self.parse_sequence(element, parent_path)

        elif elmt_tag == "struct":
            self.parse_struct(element, parent_path)

        elif elmt_tag in name_val_elmt_tag:
            self.last_element_parser(element, parent_path)

        else:
            raise TypeError("Needs to parse to different type of parser")

    def parse_sequence(self, element_: EmtT.Element, parent_path: str) -> None:
        """
        Parameters
        ----------
        element_ : Element with sequence tag
        parent_path : Xpath of the parent element where the element_ belongs

        Returns
        -------
        None
        """

        child_num = len(element_)
        elmt_attr = element_.attrib

        section_nm_reslvr = ""
        key_name = "name"
        if key_name in elmt_attr.keys():
            section_nm_reslvr = f"{elmt_attr[key_name]}"
            section_nm_reslvr = self.check_for_siblings_with_same_name(
                section_nm_reslvr, element_
            )

            parent_path = f"{parent_path}/{section_nm_reslvr}"

        child_elmt_ind = 0
        while child_num > 0:
            self.pass_child_through_parsers(
                element_, parent_path, child_elmt_ind, skip_child=-1
            )
            child_num -= 1
            child_elmt_ind += 1

    def struct_fc_name_sc_value(self, element_, first_child, parent_path, skip_child):
        """Struct representing parameter with first child (fc) 'name'
        and second child(sc) 'value'."""
        section_nm_reslvr = ""
        units = ["mV", "degree", "W", "kV", "ns"]

        section_nm_reslvr = self.restructure_value(first_child.text, first_child.tag)
        section_nm_reslvr = self.check_for_siblings_with_same_name(
            section_nm_reslvr, element_
        )
        skip_child += 1
        # Separating the units
        for unit in units:
            if f"_[{unit}]" in section_nm_reslvr:
                section_nm_reslvr, _ = section_nm_reslvr.split("_")
                self.metadata_dict[f"{parent_path}/{section_nm_reslvr}/@units"] = unit

        parent_path, self.tail_part_frm_struct = self.check_last_part_repetition(
            parent_path, self.tail_part_frm_struct, section_nm_reslvr
        )

        return parent_path, skip_child

    def struct_fc_name_sc_string(self, element_, first_child, parent_path, skip_child):
        """Struct representing parameter with first child(fc) 'name'
        and first child(fc) having 'string'."""

        elmt_attr = element_.attrib
        key_type_name = "type_name"
        skip_child += 1
        child_txt = self.restructure_value(first_child.text, first_child.tag)

        section_nm_reslvr = f"{elmt_attr[key_type_name]}_{child_txt}"
        section_nm_reslvr = self.check_for_siblings_with_same_name(
            section_nm_reslvr, element_
        )

        parent_path = f"{parent_path}/{section_nm_reslvr}"

        return parent_path, skip_child

    def parse_struct(self, element_: EmtT.Element, parent_path: str) -> None:
        """
        Parameters
        ----------
        element_ : Element with struct tag
        parent_path : Xpath of the parent element where the element_ belongs

        Returns
        -------
        None
        """

        child_num = len(element_)
        elmt_attr = element_.attrib

        # Resolving struct name section is here
        skip_child = -1
        section_nm_reslvr = ""
        first_child = element_[0]
        second_child = element_[1]

        key_name = "name"
        key_value = "value"
        key_type_name = "type_name"
        if key_name in elmt_attr.keys():
            section_nm_reslvr = elmt_attr[key_name]
            section_nm_reslvr = self.check_for_siblings_with_same_name(
                section_nm_reslvr, element_
            )
            parent_path, self.tail_part_frm_struct = self.check_last_part_repetition(
                parent_path, self.tail_part_frm_struct, section_nm_reslvr
            )

        elif key_name not in elmt_attr.keys():
            if (
                key_name in first_child.attrib.values()
                and key_value in second_child.attrib.values()
            ):
                parent_path, skip_child = self.struct_fc_name_sc_value(
                    element_, first_child, parent_path, skip_child
                )
            elif (
                key_name in first_child.attrib.values() and first_child.tag == "string"
            ):
                parent_path, skip_child = self.struct_fc_name_sc_string(
                    element_, first_child, parent_path, skip_child
                )

            else:
                # Check twin siblings
                section_nm_reslvr = self.restructure_value(
                    elmt_attr[key_type_name], "string"
                )
                section_nm_reslvr = (
                    section_nm_reslvr + "_" + str(elmt_attr["__odr_siblings__"])
                )
                parent_path = f"{parent_path}/{section_nm_reslvr}"

        child_elmt_ind = 0
        while child_num > 0:
            self.pass_child_through_parsers(
                element_, parent_path, child_elmt_ind, skip_child
            )
            child_num -= 1
            child_elmt_ind += 1

    def last_element_parser(self, element_: EmtT.Element, parent_path: str) -> None:
        """

        Parameters
        ----------
        element_ : Element with a tag among 'ulong', 'double', 'string',
                   'boolean', 'enum', 'any'
        parent_path : Xpath of the parent element where the element_ belongs

        Returns
        -------
        None
        """

        child_num = len(element_)
        elmt_attr = element_.attrib

        if child_num == 0:
            if "name" in elmt_attr.keys():
                section_nm_reslvr = f"{elmt_attr['name']}"
                value = self.restructure_value(element_.text, element_.tag)

                parent_path, self.tail_part_frm_othr = self.check_last_part_repetition(
                    parent_path, self.tail_part_frm_othr, section_nm_reslvr
                )
                self.metadata_dict[f"{parent_path}"] = value
            else:
                self.metadata_dict[f"{parent_path}"] = self.restructure_value(
                    element_.text, element_.tag
                )
        elif child_num == 1 and "any" == element_.tag:
            child_elmt = element_[0]
            self.metadata_dict[f"{parent_path}"] = self.restructure_value(
                child_elmt.text, child_elmt.tag
            )

    def check_for_siblings_with_same_name(self, reslv_name, new_sblings_elmt):
        """Check for the same name in the same level. For elments with the same
        write the name _1, _2... .
        """
        child_nm_reslvr_li = new_sblings_elmt.attrib["__parent__"].attrib[
            self.child_nm_reslvers
        ]
        if reslv_name not in child_nm_reslvr_li:
            parent = new_sblings_elmt.attrib["__parent__"]
            parent.attrib[self.child_nm_reslvers].append(reslv_name)
        else:
            last_twin_sib_nm = child_nm_reslvr_li[-1]
            try:
                ind = last_twin_sib_nm.split("_")[-1]
                reslv_name = f"{reslv_name}_{int(ind) + 1}"
                parent = new_sblings_elmt.attrib["__parent__"]
                parent.attrib[self.child_nm_reslvers].append(reslv_name)
            except ValueError:
                reslv_name = f"{reslv_name}_1"
                parent = new_sblings_elmt.attrib["__parent__"]
                parent.attrib[self.child_nm_reslvers].append(reslv_name)
        return reslv_name

    def check_last_part_repetition(
        self, parent_path: str, pre_tail_part: str, new_tail_part: str
    ) -> Tuple[str, str]:
        """
         Check for the data from the same group, for example repetition of the
         experiments under the same physical circumstances, make number of
         them with ..._1, ..._2.
        Parameters
        ----------
        self : XmlSpecs object
        parent_path : Xpath of the parent element
        pre_tail_part : The tail part added in the previous step
        new_tail_part : The tail part obtained from present element

        Returns
        -------
        parent_path : Newly obtained or replaced parent Xpath
        pre_tail_part : Newly obtained or replaced tail_path
        """
        if new_tail_part == pre_tail_part:
            previous_key = f"{parent_path}/{new_tail_part}"
            previous_val = self.metadata_dict.get(previous_key, None)
            if previous_val:
                self.metadata_dict[f"{parent_path}/{new_tail_part}_0"] = previous_val
                pre_tail_part = f"{new_tail_part}_1"
                parent_path = f"{parent_path}/{pre_tail_part}"

                del self.metadata_dict[previous_key]

                return parent_path, pre_tail_part

            parent_path = f"{parent_path}/{pre_tail_part}"
            return parent_path, pre_tail_part

        if new_tail_part in pre_tail_part:
            try:
                ind_ = pre_tail_part.split("_")[-1]
                ind = int(ind_)
                pre_tail_part = f"{new_tail_part}_{ind + 1}"
                parent_path = f"{parent_path}/{pre_tail_part}"
                return parent_path, pre_tail_part
            except ValueError:
                parent_path = f"{parent_path}/{new_tail_part}"
                pre_tail_part = new_tail_part

                return parent_path, pre_tail_part
            except TypeError:
                parent_path = f"{parent_path}/{new_tail_part}"
                pre_tail_part = new_tail_part

                return parent_path, pre_tail_part

        parent_path = f"{parent_path}/{new_tail_part}"
        pre_tail_part = new_tail_part

        return parent_path, pre_tail_part

    @staticmethod
    def restructure_value(value_text: str, element_tag: str) -> Any:
        """
            Collect the value_text transform it to the data_type according
            to the type name provided by element_tag.
        Parameters
        ----------
        value_text : text data that would be 'unsigned long', 'double', 'string',
                    'boolean', 'enum/string'
        element_tag : tag name among 'unsigned long', 'double', 'string',
                    'boolean', 'enum/string'

        Returns
        -------

        """

        def double_(para):
            return np.double(para)

        def ulong_(para):
            return np.uint(para)

        def bool_(para):
            return np.bool_(para)

        data_ty = {"double": double_, "ulong": ulong_, "boolean": bool_}
        string_ty = ["string", "enum"]

        if not value_text:
            return ""

        if element_tag in string_ty:
            value_text_: Any = " ".join(value_text.split()).replace(" ", "_")
            return value_text_

        for key_, _ in data_ty.items():
            if key_ == element_tag:
                value_text_ = value_text.split()
                numpy_value = data_ty[element_tag](value_text_)[...]
                if np.shape(numpy_value) == (1,):
                    return numpy_value[0]
                return numpy_value

    def cumulate_counts_series(
        self,
        scan_seq_elem: EmtT.Element,
        counts_length: int = None,
        cumulative_counts: np.ndarray = None,
    ) -> Tuple[str, np.ndarray]:
        """
        Sum the counts over different scans. Each ScanSeaq contains
        multiple scans under the same physical environment. The
        multiple scans are usually taken to make the peaks visible and
        distinguishable.

        Parameters
        ----------
        scan_seq_elem : Element with ScanSeq tag
        counts_length : Number of count (length of 1D numpy array)
                        contain in each scan
        cumulative_counts : Cumulative counts up to last scan  from the
                            same ScanSeq

        Returns
        -------
        np.ndarray : Cumulative up to last scans from the same ScanSeq
        """

        child_num = len(scan_seq_elem)
        name = "count"

        child_elmt_ind = 0
        while child_num >= 0:
            if scan_seq_elem.attrib["type_name"] == "CountsSeq":
                num_of_counts = int(scan_seq_elem.attrib["length"])
                if not counts_length:
                    counts_length = num_of_counts
                if counts_length != num_of_counts:
                    raise ValueError("Count number from all scans must be equal!!")

            if scan_seq_elem.attrib["type_name"] == "Counts":
                counts_data = self.restructure_value(
                    scan_seq_elem.text, scan_seq_elem.tag
                )

                if cumulative_counts is None:
                    cumulative_counts = counts_data
                else:
                    cumulative_counts = cumulative_counts + counts_data

            if child_num > 0:
                child_element = scan_seq_elem[child_elmt_ind]
                name, cumulative_counts = self.cumulate_counts_series(
                    child_element, counts_length, cumulative_counts
                )

            child_num = child_num - 1
            child_elmt_ind = child_elmt_ind + 1

        return (name, cumulative_counts)

    # pylint: disable=too-many-branches
    # pylint: disable=too-many-statements
    def collect_raw_data_to_construct_data(self):
        """Collect the raw data about detectors so that the binding energy and
        and counts for the corresponding nominal.
        """

        entry_list: List = []
        raw_dict = {
            "mcd_num": 0,
            "curves_per_scan": 0,
            "values_per_curve": 0,
            "mcd_head": 0,
            "mcd_tail": 0,
            "excitation_energy": 0,
            "kinetic_energy": 0,
            "kinetic_energy_base": 0,
            "effective_workfunction": 0,
            "scan_delta": 0,
            "pass_energy": 0,
            "mcd_shifts": [],
            "mcd_poss": [],
            "mcd_gains": [],
            "time": 0,
            "scans": {},
        }
        for key, val in self.metadata_dict.items():
            if not key.endswith("parameters/File"):
                entry = construct_entry_name(key)
                match = re.search(r"/ENTRY\[([^\]]+)\]", key)
                if match:
                    entry = match.group(1)
                else:
                    entry = ""

            if entry and (entry not in entry_list):
                self.entry_to_data[entry] = {"raw_data": copy.deepcopy(raw_dict)}
                entry_list.append(entry)

            if "region/curves_per_scan" in key:
                self.entry_to_data[entry]["raw_data"]["curves_per_scan"] = val
            elif "region/values_per_curve" in key:
                self.entry_to_data[entry]["raw_data"]["values_per_curve"] = val

            elif "region/excitation_energy" in key:
                self.entry_to_data[entry]["raw_data"]["excitation_energy"] = val

            elif "region/scan_mode/name" in key:
                val = convert_energy_scan_mode(val)
                self.entry_to_data[entry]["raw_data"]["scan_mode"] = val

            elif "region/kinetic_energy" in key:
                if "region/kinetic_energy_base" not in key:
                    self.entry_to_data[entry]["raw_data"]["kinetic_energy"] = val
                    continue
                if "region/kinetic_energy_base" in key:
                    self.entry_to_data[entry]["raw_data"]["kinetic_energy_base"] = val
                    continue

            elif "region/effective_workfunction" in key:
                self.entry_to_data[entry]["raw_data"]["effective_workfunction"] = val

            elif "region/scan_delta" in key:
                self.entry_to_data[entry]["raw_data"]["scan_delta"] = val

            elif "region/pass_energy" in key:
                self.entry_to_data[entry]["raw_data"]["pass_energy"] = val

            elif "mcd_head" in key:
                self.entry_to_data[entry]["raw_data"]["mcd_head"] = val

            elif "mcd_tail" in key:
                self.entry_to_data[entry]["raw_data"]["mcd_tail"] = val

            elif "shift" in key:
                self.entry_to_data[entry]["raw_data"]["mcd_shifts"].append(val)
                self.entry_to_data[entry]["raw_data"]["mcd_num"] += 1

            elif "gain" in key:
                self.entry_to_data[entry]["raw_data"]["mcd_gains"].append(val)

            elif "position" in key:
                self.entry_to_data[entry]["raw_data"]["mcd_poss"].append(val)

            # construct scan names e.g cycles2_scan0
            if "cycles/Cycle_" in key:
                _, last_part = key.split("cycles/Cycle_")
                if "/time" in last_part:
                    self.entry_to_data[entry]["raw_data"]["time"] = val
                    continue
                if "/parameters/Loop" in last_part:
                    self.entry_to_data[entry]["raw_data"]["loop_no"] = val
                    continue
                parts = last_part.split("/")
                cycle_num, scan_num = parts[0], parts[-2].split("_")[1]
                scan_name = f"cycle{cycle_num}_scan{scan_num}"

                self.entry_to_data[entry]["raw_data"]["scans"][scan_name] = val
