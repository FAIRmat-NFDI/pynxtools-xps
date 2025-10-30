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
# pylint: disable=too-many-lines,too-many-instance-attributes
"""
Parser for reading XPS (X-ray Photoelectron Spectroscopy) data from
Specs Lab Prodigy XY exports, to be passed to mpes nxdl
(NeXus Definition Language) template.
"""

import copy
import datetime
import itertools
import logging
import re
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from pynxtools_xps.reader_utils import (
    XPSMapper,
    check_uniform_step_width,
    construct_data_key,
    construct_entry_name,
    get_minimal_step,
    re_map_keys,
    re_map_values,
)
from pynxtools_xps.value_mappers import (
    convert_energy_scan_mode,
    convert_energy_type,
    convert_measurement_method,
    convert_units,
    get_units_for_key,
    parse_datetime,
)

logger = logging.getLogger("pynxtools")

SETTINGS_MAP: dict[str, str] = {
    "Group": "group_id",
    "Scan Mode": "scan_mode",
    "Analyzer Lens Voltage": "analyzer_lens_voltage",
    "Calibration File": "calibration_file",
    "Transmission File": "transmission_file",
    "Analyzer Slit": "entrance_slit",
    "Iris Diameter": "iris_diameter",
    "Energy Axis": "x_units",
    "Source": "source_label",
    "Polar Angle": "source_polar_angle",
    "Azimuth Angle": "source_azimuth_angle",
    "ex_energy": "excitation_energy",
    "Acquisition Date": "time_stamp",
    "Analysis Method": "analysis_method",
    "Analyzer": "analyzer_name",
    "Analyzer Lens": "lens_mode",
    "Analyzer Lens Mode": "lens_mode",
    "Scan Variable": "scan_variable",
    "Curves/Scan": "curves_per_scan",
    "Values/Curve": "n_values",
    "Step Size": "step_size",
    "Dwell Time": "dwell_time",
    "Excitation Energy": "excitation_energy",
    "Kinetic Energy": "kinetic_energy",
    "Pass Energy": "pass_energy",
    "Bias Voltage": "bias_voltage_electrons",
    "Binding Energy": "start_energy",
    "Detector Voltage": "detector_voltage",
    "Eff. Workfunction": "work_function",
    "Normalized By": "normalized_by",
    "Comment": "comments",
    "Spectrum ID": "spectrum_id",
    "Note": "note",
}

VALUE_MAP: dict[str, Any] = {
    "analysis_method": convert_measurement_method,
    "scan_mode": convert_energy_scan_mode,
    "bias_voltage_electrons": float,
    "n_values": int,
    "excitation_energy": float,
    "kinetic_energy": float,
    "work_function": float,
    "dwell_time": float,
    "detector_voltage": float,
    "curves_per_scan": int,
    "pass_energy": float,
    "spectrum_id": int,
    "x_units": convert_energy_type,
}

UNITS: dict[str, str] = {
    "work_function": "eV",
    "excitation_energy": "eV",
    "pass_energy": "eV",
    "bias_voltage_electrons": "V",
    "dwell_time": "s",
    "step_size": "eV",
}


class XyMapperSpecs(XPSMapper):
    """
    Class for restructuring .xy data file from
    Specs vendor into python dictionary.
    """

    config_file = "config_specs_xy.json"

    def __init__(self):
        super().__init__()
        self.write_channels_to_data = True

    def _select_parser(self):
        return XyProdigyParser()

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the file using the Specs XY parser.

        Parameters
        ----------
        file : str
            Filepath of the XY file.
        **kwargs : dict
            write_channels_to_data: bool
                If True, the spectra of each individual channel is
                written to the entry/data field in the MPES template.

        Returns
        -------
        dict
            Flattened dictionary to be passed to MPES template.

        """
        if "write_channels_to_data" in kwargs:
            self.write_channels_to_data = kwargs["write_channels_to_data"]

        return super().parse_file(file, **kwargs)

    def construct_data(self):
        """Map XY format to NXmpes-ready dict."""
        # pylint: disable=duplicate-code
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        """
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []
        for part in ["group_name", "region_name"]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            mpes_key = f"{entry_parent}/{key}"
            if "units" in key:
                if isinstance(value, dict):
                    value = {k: convert_units(v) for k, v in value.items()}
                else:
                    value = convert_units(value)
            self._xps_dict[mpes_key] = value
            units = convert_units(get_units_for_key(key, UNITS))
            if units is not None:
                self._xps_dict[f"{mpes_key}/@units"] = units

        data = spectrum["data"]

        if self.parser.export_settings["Transmission Function"]:
            self._xps_dict["transmission_function"] = data["transmission"]
            self._xps_dict["transmission_function/units"] = "counts_per_second"

        # Create key for writing to data.
        scan_key = construct_data_key(spectrum)

        x_axis_name, x_axis = list(data.items())[0]
        y_axis_name, intensity = list(data.items())[1]
        x_axis, intensity = np.array(x_axis), np.array(intensity)

        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        if not self.parser.export_settings["Separate Channel Data"]:
            averaged_channels = intensity
        else:
            all_channel_data = [
                value
                for key, value in self._xps_dict["data"][entry].items()
                if "_chan" in key
            ]
            if not all_channel_data:
                averaged_channels = intensity
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    averaged_channels = np.mean(all_channel_data, axis=0)

        if not self.parser.export_settings["Separate Scan Data"]:
            averaged_scans = intensity
        else:
            all_scan_data = [
                value
                for key, value in self._xps_dict["data"][entry].items()
                if "_scan" in key and "_chan" not in key
            ]
            if not all_scan_data:
                averaged_scans = intensity
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    averaged_scans = np.mean(all_scan_data, axis=0)

        # Write to data in order: scan, cycle, channel

        # Write averaged cycle data to 'data'.
        self._xps_dict["data"][entry][scan_key.split("_")[0]] = xr.DataArray(
            data=averaged_scans,
            coords={x_axis_name: x_axis},
        )

        # Write average cycle data to '_scan'.
        self._xps_dict["data"][entry][scan_key] = xr.DataArray(
            data=averaged_channels,
            coords={x_axis_name: x_axis},
        )

        if (
            self.parser.export_settings["Separate Channel Data"]
            and self.write_channels_to_data
        ):
            # Write channel data to '_chan'.
            channel_no = spectrum["channel_no"]
            self._xps_dict["data"][entry][f"{scan_key}_chan{channel_no}"] = (
                xr.DataArray(
                    data=intensity,
                    coords={x_axis_name: x_axis},
                )
            )

        if self.parser.export_settings["External Channel Data"]:
            self._xps_dict[f"{entry_parent}/aux_signals"] = []
            for ext_channel, channel_data in list(spectrum["data"].items()):
                if ext_channel not in [x_axis_name, y_axis_name, "transmission"]:
                    self._xps_dict[f"{entry_parent}/external_{ext_channel}"] = (
                        channel_data
                    )
                    self._xps_dict[f"{entry_parent}/external_{ext_channel}/@units"] = (
                        spectrum["channel_units"].get(ext_channel)
                    )
                    self._xps_dict[f"{entry_parent}/aux_signals"] += [ext_channel]

    def _get_units_for_key(self, unit_key: str):
        """
        Get correct units for a given key.

        Parameters
        ----------
        unit_key : str
           Key of type <mapping>:<spectrum_key>, e.g.
           detector/detector_voltage

        Returns
        -------
        str
            Unit for that unit_key.

        """
        try:
            return re.search(r"\[([A-Za-z0-9_]+)\]", unit_key).group(1)
        except AttributeError:
            try:
                return UNITS[unit_key]
            except KeyError:
                return ""


class XyProdigyParser:  # pylint: disable=too-few-public-methods
    """
    A parser for reading in ASCII-encoded .xy data from Specs Prodigy.

    Tested with SpecsLab Prodigy v 4.64.1-r88350.
    """

    def __init__(self):
        """
        Construct the parser.

        """
        self.lines = []
        self.prefix = "#"
        self.n_headerlines = 14
        self.export_settings = {}

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the .xy file into a list of dictionaries.

        Parsed data is stored in the attribute 'self.data'.
        Each dictionary in the data list is a grouping of related
        attributes. The dictionaries are later re-structured into a
        nested dictionary that more closely resembles the domain logic.

        Parameters
        ----------
        file : str
            XPS data filepath.

        **kwargs : dict
            commentprefix : str
                Prefix for comments in xy file. The default is "#".
            n_headerlines: int
                number of header_lines in each data block.

        Returns
        -------
        list
            Flat list of dictionaries containing one spectrum each.

        """
        if "commentprefix" in kwargs:
            self.prefix = kwargs["commentprefix"]
        if "n_headerlines" in kwargs:
            self.n_headerlines = kwargs["n_headerlines"]

        self.lines = self._read_lines(file)
        header, data = self._separate_header()
        self.export_settings = self._parse_export_settings(header)

        # Recursively read XPS data from flat 'lines' list.
        groups = self._handle_groups(data)

        return self._flatten_dict(groups)

    def _read_lines(self, file: str | Path):
        """
        Read all lines from the input XY files.

        Parameters
        ----------
        file : str
            Filepath of the XY file to be read.

        Returns
        -------
        lines : list
            All lines in the XY file.

        """
        with open(file, encoding="utf-8") as xy_file:
            lines = xy_file.readlines()

        return lines

    def _separate_header(self):
        """
        Split of common export header.

        Returns
        -------
        header : list
            List of export settings.
        groups : list
            List of list containing data strings.

        """
        header = self.lines[: self.n_headerlines]
        groups = self.lines[self.n_headerlines :]

        return header, groups

    def _parse_export_settings(self, header: list[str]):
        """
        Parse the top-level Prodigy export settings into a dict.

        Parameters
        ----------
        header : list
            List of header strings.

        Returns
        -------
        export_settings : dict
            Dictionary of export settings.

        """
        bool_map = {
            "yes": True,
            "no": False,
        }

        export_settings = {}
        for line in header:
            line = line.strip(self.prefix).strip()
            if len(line) == 0:
                pass
            else:
                setting = line.split(":", 1)[1].strip()
                setting_bool = bool_map.get(setting, setting)
                export_settings[line.split(":", 1)[0].strip()] = setting_bool

        export_settings = re_map_keys(export_settings, SETTINGS_MAP)
        export_settings = re_map_values(export_settings, VALUE_MAP)

        return export_settings

    def _handle_groups(self, data):
        """
        Separate the data list into a dictionary, with each
        element containing the data list of one group.

        Parameters
        ----------
        data : list
            Full data list (list of strings).

        Returns
        -------
        groups : dict
            Dict with data organized by group
            Entries are as group_name: group_data.

        """
        grouped_list = [
            list(g)
            for _, g in itertools.groupby(
                data, lambda line: "Group:" in line.strip(self.prefix).strip()
            )
        ][1:]

        groups = OrderedDict()

        for name_header, group_data in zip(grouped_list[::2], grouped_list[1::2]):
            name = self._strip_param(name_header[0], "Group:")
            group_settings = {"group_name": name}
            group_settings = re_map_keys(group_settings, SETTINGS_MAP)
            group_settings = re_map_values(group_settings, VALUE_MAP)

            group = self._handle_regions(group_data)
            group["group_settings"] = group_settings
            groups[name] = group

        return groups

    def _handle_regions(self, group_data: list[str]):
        """
        Separate the data list of an individual group into a
        dictionary, with each element containing the data list
        of one region.

        Parameters
        ----------
        group_data : list
            Group data list (list of strings).

        Returns
        -------
        regions : dict
            Dict with data organized by group
            Entries are as region_name: region_data.

        """
        grouped_list = [
            list(g)
            for _, g in itertools.groupby(
                group_data, lambda line: "Region:" in line.strip(self.prefix).strip()
            )
        ][1:]

        regions = OrderedDict()

        for name_header, region_data in zip(grouped_list[::2], grouped_list[1::2]):
            region_settings = {}
            name = self._strip_param(name_header[0], "Region:")
            region_settings["region_name"] = name

            for i, setting in enumerate(region_data):
                if not setting.startswith(self.prefix):
                    region_data = region_data[i:]
                    break
                setting_name = setting.split(self.prefix)[-1].strip().split(":")[0]

                parts = setting.split(self.prefix)[-1].strip().split(":")

                try:
                    region_settings[setting_name] = parts[1].strip()
                except IndexError:
                    pass

            region_settings = re_map_keys(region_settings, SETTINGS_MAP)
            region_settings = re_map_values(region_settings, VALUE_MAP)

            region = self._handle_cycles(region_data)
            region["region_settings"] = region_settings
            regions[name] = region

        return regions

    def _handle_cycles(self, region_data: list[str]):
        """
        Separate the data list of an individual region into a
        dictionary, with each element containing the data list
        of one cycle.

        Parameters
        ----------
        region_data : list
            Region data list (list of strings).

        Returns
        -------
        cycles : dict
            Dict with data organized by cycle
            Entries are as cycle_name: cycle_data.

        """
        cycle_pattern = re.compile(rf"{self.prefix} Cycle: \d\n", re.IGNORECASE)

        cycles = OrderedDict()
        cycle_line_nrs = {}

        for i, line in enumerate(region_data):
            if cycle_pattern.match(line):
                cycle_line_nrs[
                    "cycle_" + str(int(self._strip_param(line, "Cycle:")))
                ] = i
            if i == len(region_data) - 1:
                cycle_line_nrs["end"] = i + 1

        for i, (line_no_a, line_no_b) in enumerate(
            zip(list(cycle_line_nrs.values()), list(cycle_line_nrs.values())[1:])
        ):
            name = f"cycle_{i}"
            cycle_settings = {"loop_no": i}
            cycle_data = region_data[line_no_a:line_no_b]

            cycle_settings = re_map_keys(cycle_settings, SETTINGS_MAP)
            cycle_settings = re_map_values(cycle_settings, VALUE_MAP)

            cycle = self._handle_individual_cycles(cycle_data)
            cycle["cycle_settings"] = cycle_settings

            cycles[name] = cycle

        return cycles

    def _handle_individual_cycles(self, cycle_data: list[str]):
        """
        Separate the data list of an individual cycle into a
        dictionary, with each element containing the data list
        of one scan.

        Parameters
        cycle_data : list
            Cycle data list (list of strings).

        Returns
        -------
        scan : dict
            Dict with data organized by cycle
            Entries are as scan_name: scan_data.

        """
        spec_pattern_str = rf"{self.prefix} Cycle: \d, Curve: \d"
        if self.export_settings["Separate Scan Data"]:
            spec_pattern_str += r", Scan: \d"
        if self.export_settings["Separate Channel Data"]:
            spec_pattern_str += r", Channel: \d"
        spec_pattern = re.compile(spec_pattern_str, re.IGNORECASE)

        scans = OrderedDict()
        scan_line_nrs = {}

        for i, line in enumerate(cycle_data):
            if spec_pattern.match(line):
                name_dict = {
                    a.strip(): int(b.strip())
                    for a, b in (
                        element.split(": ")
                        for element in line.strip(self.prefix).strip().split(", ")
                    )
                }
                name = "".join(
                    [
                        f"{key.lower()}_{val}_"
                        for key, val in name_dict.items()
                        if key != "Curve"
                    ]
                ).rstrip("_")
                scan_line_nrs[name] = i
            if i == len(cycle_data) - 1:
                scan_line_nrs["end"] = i + 1

        for i, ((name, line_no_a), line_no_b) in enumerate(
            zip(list(scan_line_nrs.items()), list(scan_line_nrs.values())[1:])
        ):
            scan_data = cycle_data[line_no_a:line_no_b]
            scan = self._handle_individual_scan(scan_data)
            scan["scan_settings"].update(self._extend_scan_settings(name))
            scans[name] = scan

        return scans

    def _handle_individual_scan(self, scan_data: list[str]):
        """
        Separate the data list of an individual scan into a
        dictionary, with each element containing the data and
        metadata of one spectrum. External channels are handled as well

        Parameters
        scan_data : list
            Scan data list (list of strings).

        Returns
        -------
        scan : dict
            Dict with scan data and metadata
            Entries are as
            scan_name: {
                "data": scan_data,
                "scan_settings": scan_settings
                }

        """

        def _normalize_ext_channel_label(label: str) -> tuple[str, str | None]:
            """
            Normalize a label by converting the main part to snake_case and separating the unit.

            - Converts PascalCase or CamelCase to snake_case.
            - Preserves a single trailing capital letter (e.g., 'D' or 'L') without inserting an underscore.
            - Extracts and standardizes the unit (e.g., converts 's-1' to '1/s').
            - Extracts context from parentheses and appends it to the name in snake_case.

            Args:
                label (str): Raw input label, e.g., 'Beam_SplitterD [A]'.

            Returns:
                tuple: (normalized_label: str, unit: Optional[str])

            Examples:
                >>> _normalize_ext_channel_label("Beam_SplitterD [A]")
                ('beam_splitterD', 'A')

                >>> _normalize_ext_channel_label("Photo_count_ext [s-1]")
                ('photo_count_ext', '1/s')

                >>> _normalize_ext_channel_label("Excitation Energy [eV] (Monochromator)")
                ('excitation_energy_monochromator', 'eV')

                >>> _normalize_ext_channel_label("Ring Current [mA] (Monochromator)")
                ('ring_current_monochromator', 'mA')

                >>> _normalize_ext_channel_label("energy")
                ('energy', None)
            """
            label = label.strip()

            # Extract context in parentheses
            context = ""
            paren_match = re.search(r"\(([^)]+)\)", label)
            if paren_match:
                context = paren_match.group(1).strip()
                label = re.sub(r"\s*\([^)]*\)", "", label)

            # Extract unit in square brackets
            unit = None
            bracket_match = re.search(r"\[([^\]]+)\]", label)
            if bracket_match:
                unit = bracket_match.group(1).strip()
                label = re.sub(r"\s*\[[^\]]*\]", "", label)
                unit = convert_units(unit)

            # Preserve trailing uppercase letter (e.g., D or L)
            trailing_upper = ""
            if label and label[-1].isupper():
                trailing_upper = label[-1]
                label = label[:-1]

            # Convert label to snake_case
            name = re.sub(
                r"(?<!^)(?<!_)([A-Z])", r"_\1", label.replace(" ", "_")
            ).lower()

            if trailing_upper:
                name += trailing_upper  # Append preserved trailing capital

            if context:
                name += f"_{context.lower()}"

            return name, unit

        scan_settings: dict[str, Any] = {}
        data_channels: OrderedDict[str, list[float]] = OrderedDict()
        channel_units: OrderedDict[str, str] = OrderedDict()

        in_remote_channel = False

        if self.export_settings["External Channel Data"]:
            remote_ch_pattern_str = rf"{re.escape(self.prefix)} External Channel Data Cycle: \d+,\s*(.+?)\s*\(Remote Out Device\)"
            remote_ch_pattern = re.compile(remote_ch_pattern_str, re.IGNORECASE)

        for line in scan_data:
            if match := remote_ch_pattern.search(line):
                in_remote_channel = True
                channel = _normalize_ext_channel_label(match.group(1))[0]

            if line.startswith(self.prefix) and line.strip(self.prefix).strip():
                key, val = (
                    item.strip()
                    for item in line.strip(self.prefix).strip().split(":", 1)
                )
                if key == "Acquisition Date":
                    scan_settings[key] = self._parse_datetime(val)

                if key == "ColumnLabels":
                    if in_remote_channel:
                        column_labels = val.strip().split(" ", 1)
                        column_labels = [
                            label.split("(Remote Out Device)", 1)[0].strip()
                            for label in column_labels
                        ]
                    else:
                        column_labels = val.strip().split(" ")
                    normalized_labels_and_units = [
                        _normalize_ext_channel_label(label) for label in column_labels
                    ]

                    for i, (label, unit) in enumerate(normalized_labels_and_units):
                        if not (in_remote_channel and i == 0):
                            data_channels[label] = []
                            if unit:
                                channel_units[label] = unit

            elif not line.startswith(self.prefix) and line.strip():
                values = [float(v) for v in line.strip().split()]
                for i, (channel_label_and_unit, v) in enumerate(
                    zip(normalized_labels_and_units, values)
                ):
                    label = channel_label_and_unit[0]
                    if not (in_remote_channel and i == 0):
                        data_channels[label].append(v)

        energy_axis = data_channels.get("energy")

        if energy_axis and check_uniform_step_width(energy_axis):
            scan_settings["step_size"] = get_minimal_step(energy_axis)
        scan_settings = re_map_keys(scan_settings, SETTINGS_MAP)
        scan_settings = re_map_values(scan_settings, VALUE_MAP)

        scan: OrderedDict[str, dict[str, Any]] = OrderedDict()
        scan["data"] = data_channels
        scan["channel_units"] = channel_units
        scan["scan_settings"] = scan_settings

        return scan

    def _extend_scan_settings(self, scan_name: str):
        """
        Split the scan name and extract the scan metadata.

        Example:
            scan_name == scan_0_channel_0
            -> settings = {'scan_no': 0, 'channel_no': 0}

        Parameters
        ----------
        scan_name : str
            String with the name of the scan, in the format
            scan_0_channel_0.

            Channel number is optional.

        Returns
        -------
        settings : dict
            Dict with scan settings.

        """
        settings = {}

        split_name = scan_name.split("_")

        for param, val in zip(split_name[::2], split_name[1::2]):
            if param != "cycle":
                settings[f"{param}_no"] = int(val)

        return settings

    def _strip_param(self, line: str, key: str):
        """
        Split the scan name and extract the scan metadata.

        Example:
            if self.prefix = "#"
            line = Cycle: # Cycle: 5, key = "Cycle:\n"
            -> return 5

        Parameters
        ----------
        line : str
            String containing the string self.prefix + " " + key"
        key : str
            Keyword to strip from line.

        Returns
        -------
        str
            Stripped line without prefix and key.

        """
        if key in line:
            return line.strip().split(self.prefix + " " + key)[-1].strip()
        return line

    def _flatten_dict(self, data_dict: dict[str, Any]):
        """
        Flatten a raw data dict into a list, with each element
        being a dictionary with data and metadata for one spectrum.

        Parameters
        ----------
        data_dict : dict
            Nested dictionary containing group, regions, cycles,
            and scans.

        Returns
        -------
        spectra : list
            Flattened list of spectra dicts.

        """
        spectra = []

        for group in data_dict.values():
            group_settings = group["group_settings"]
            for region in list(group.values())[:1]:
                region_settings = region["region_settings"]
                for cycle in list(region.values())[:1]:
                    cycle_settings = cycle["cycle_settings"]
                    for scan in list(cycle.values())[:1]:
                        scan_settings = scan["scan_settings"]
                        spectrum: dict[str, Any] = {"data": {}}
                        for settings in [
                            group_settings,
                            region_settings,
                            cycle_settings,
                            scan_settings,
                        ]:
                            spectrum.update(settings)
                        for key, val in scan.items():
                            if key != "scan_settings":
                                spectrum[key] = val
                        spectrum["x_units"] = self.export_settings["x_units"]
                        spectra.append(spectrum)

        return spectra

    def _parse_datetime(self, date: str) -> str:
        """
        Parse datetime into a datetime.datetime object and return a
        string value in ISO format.

        Parameters
        ----------
        date : str
            String representation of the date, in of these formats:
            "%m/%d/%y %H:%M:%S",
            "%m-%d-%y %H:%M:%S".

        Returns
        -------
        date_object : datetime.datetime
            Datetime in datetime.datetime format.

        """
        # 2025-04-17 13:06:40 UTC
        if date.find("UTC"):
            date = date[: date.find("UTC")].strip()
            tzinfo = datetime.timezone.utc
        else:
            date = date.strip()
            tzinfo = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo  # type: ignore[assignment]

        possible_date_formats = ["%m/%d/%y %H:%M:%S", "%Y-%m-%d %H:%M:%S"]
        return parse_datetime(date, possible_date_formats, tzinfo)
