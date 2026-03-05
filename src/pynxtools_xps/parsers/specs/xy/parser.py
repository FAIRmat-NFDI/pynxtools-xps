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

import itertools
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import xarray as xr

from pynxtools_xps.mapping import (
    _convert_bool,
    _get_measurement_method_long,
    _Value,
    convert_pascal_to_snake,
)
from pynxtools_xps.numerics import _get_minimal_step, check_uniform_step_width
from pynxtools_xps.parsers.base import ParsedSpectrum, _construct_entry_name, _XPSParser
from pynxtools_xps.parsers.specs.xy.metadata import _context


class SpecsXYParser(_XPSParser):
    """
    A parser for reading in ASCII-encoded .xy data from Specs Prodigy.

    Tested with SpecsLab Prodigy v 4.64.1-r88350.
    """

    config_file: ClassVar[str] = "config_specs_xy.json"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".xy",)
    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset(
        {
            "loop_no",
            "scan_no",
            "channel_no",
            "data",
            "channel_units",
        }
    )

    def matches_file(self, file: Path) -> bool:
        """Return True for SpecsLab Prodigy XY export files."""
        try:
            with open(file, encoding="utf-8", errors="ignore") as f:
                first = f.readline()
            return first.startswith("# Created by:") and "SpecsLab" in first
        except Exception:
            return False

    def __init__(self):
        """
        Construct the parser.

        """
        super().__init__()
        self.lines = []
        self.prefix = "#"
        self.n_headerlines = 14
        self.export_settings = {}

    def _parse(self, file: Path, **kwargs) -> None:
        """
        Parse the .xy file and populate ``self._data``.

        Each entry in ``self._data`` corresponds to one
        (group_name, region_name) pair and holds a ``ParsedSpectrum``
        with dims ``("cycle", "scan", "energy")``.

        Parameters
        ----------
        file : Path
            XPS data filepath.
        commentprefix : str, optional
            Prefix for comments in xy file. The default is ``"#"``.
        n_headerlines : int, optional
            Number of header lines in each data block.

        """
        if "commentprefix" in kwargs:
            self.prefix = kwargs["commentprefix"]
        if "n_headerlines" in kwargs:
            self.n_headerlines = kwargs["n_headerlines"]

        self.lines = self._read_lines(file)
        header, data = self._separate_header()
        self.export_settings = self._parse_export_settings(header)

        groups = self._handle_groups(data)
        self._data = self._build_parsed_spectra(groups)

    def _build_parsed_spectra(
        self, groups: dict[str, Any]
    ) -> dict[str, ParsedSpectrum]:
        """
        Walk the nested group/region/cycle/scan structure and assemble
        ``ParsedSpectrum`` objects.

        The primary channel (channel_no == 0) scans are stacked into a
        DataArray with dims ``("cycle", "scan", "energy")``.  External
        channels (all other data columns besides energy and intensity)
        are stored in ``metadata`` as ``"external_{channel_name}"``.

        Parameters
        ----------
        groups : dict
            Nested dict produced by ``_handle_groups``.

        Returns
        -------
        dict[str, ParsedSpectrum]
            Mapping from NeXus entry name to typed spectrum.

        """
        # Collect all flat scan records across the file.
        flat_scans = list(self._iter_flat_scans(groups))

        # Group by entry_name.
        entries: dict[str, list[dict[str, Any]]] = {}
        for scan in flat_scans:
            entry_name = _construct_entry_name(
                [scan.get("group_name", ""), scan.get("region_name", "")]
            )
            entries.setdefault(entry_name, []).append(scan)

        parsed: dict[str, ParsedSpectrum] = {}
        for entry_name, scans in entries.items():
            parsed[entry_name] = self._assemble_entry(entry_name, scans)

        return parsed

    def _iter_flat_scans(self, groups: dict[str, Any]):
        """
        Yield one flat dict per scan from the nested groups structure.

        Each yielded dict has:
          - ``group_name``, ``region_name``
          - ``loop_no``, ``scan_no``, ``channel_no``
          - ``data``: OrderedDict (x-axis_name → values, y_axis_name → values, ...)
          - ``channel_units``: dict
          - all scalar metadata keys from group/region/cycle/scan settings

        """
        group_settings_global = groups.get("group_settings", {})

        for group_name, group in groups.items():
            if group_name == "group_settings":
                continue

            group_settings = group.get("group_settings", {})

            for region_name, region in group.items():
                if region_name == "group_settings":
                    continue

                region_settings = region.get("region_settings", {})

                for cycle_name, cycle in region.items():
                    if cycle_name == "region_settings":
                        continue

                    cycle_settings = cycle.get("cycle_settings", {})

                    for scan_name, scan in cycle.items():
                        if scan_name == "cycle_settings":
                            continue

                        scan_settings = scan.get("scan_settings", {})
                        scan_data = scan.get("data", {})
                        channel_units = scan.get("channel_units", {})

                        record: dict[str, Any] = {}
                        record.update(group_settings_global)
                        record.update(group_settings)
                        record.update(region_settings)
                        record.update(cycle_settings)
                        record.update(scan_settings)

                        # Extend scan settings (scan_no, channel_no) from
                        # the scan_name string.
                        record.update(self._extend_scan_settings(scan_name))

                        record["group_name"] = group_settings.get(
                            "group_name", group_name
                        )
                        record["region_name"] = region_settings.get(
                            "region_name", region_name
                        )
                        record["data"] = scan_data
                        record["channel_units"] = channel_units
                        record["x_units"] = self.export_settings.get("x_units", "")

                        yield record

    def _assemble_entry(
        self, entry_name: str, scans: list[dict[str, Any]]
    ) -> ParsedSpectrum:
        """
        Build a ``ParsedSpectrum`` from all scans belonging to one entry.

        Primary scans (channel_no == 0, or no channel separation) are
        stacked by (loop_no, scan_no) into a DataArray with dims
        ``("cycle", "scan", "energy")``.

        External channel columns (neither the energy axis nor the primary
        intensity column) are stored in ``metadata`` as
        ``"external_{channel_name}"`` (np.ndarray of concatenated values).

        Parameters
        ----------
        entry_name : str
            NeXus entry name, used only for error reporting.
        scans : list[dict]
            All flat scan records belonging to this entry.

        Returns
        -------
        ParsedSpectrum

        """
        # Separate primary-channel scans from external-channel scans.
        # When separate_channel_data is False, every scan is "primary".
        if self.export_settings.get("separate_channel_data", False):
            primary_scans = [s for s in scans if s.get("channel_no", 0) == 0]
        else:
            primary_scans = scans

        # Sort by (loop_no, scan_no) for deterministic ordering.
        primary_scans = sorted(
            primary_scans, key=lambda s: (s.get("loop_no", 0), s.get("scan_no", 0))
        )

        # Determine x-axis name from the first data dict.
        first_data = primary_scans[0]["data"] if primary_scans else {}
        data_items = list(first_data.items())
        x_axis_name = data_items[0][0] if data_items else "energy"
        y_axis_name = data_items[1][0] if len(data_items) > 1 else "intensity"

        energy_axis = np.array(data_items[0][1]) if data_items else np.array([])

        # Stack intensities: shape (1, n_scans, n_energy)
        intensities = np.stack(
            [np.array(s["data"][y_axis_name]) for s in primary_scans], axis=0
        )
        # Add a cycle dimension of size 1.
        intensities = intensities[np.newaxis, ...]  # (1, n_scans, n_energy)

        data_array = xr.DataArray(
            data=intensities,
            dims=("cycle", "scan", x_axis_name),
            coords={x_axis_name: energy_axis},
        )

        # Build metadata from the first primary scan's non-identity fields.
        metadata: dict[str, Any] = {}
        if primary_scans:
            first_scan = primary_scans[0]
            for key, value in first_scan.items():
                if key not in self._metadata_exclude_keys:
                    metadata[key] = value
                    if key == "analysis_method":
                        metadata["f{key}_long_name"] = _get_measurement_method_long(
                            value
                        )

        # Collect external channel columns from all primary scans and store
        # them as concatenated arrays under "external_{channel_name}".
        if self.export_settings.get("external_channel_data", False) and primary_scans:
            primary_data = primary_scans[0]["data"]
            ext_keys = [
                k
                for k in primary_data
                if k not in (x_axis_name, y_axis_name, "transmission")
            ]
            for ext_key in ext_keys:
                ext_arrays = [
                    np.array(s["data"].get(ext_key, [])) for s in primary_scans
                ]
                metadata[f"external_{ext_key}"] = np.array(ext_arrays)

                channel_units = primary_scans[0].get("channel_units", {})
                if ext_key in channel_units:
                    metadata[f"external_{ext_key}/@units"] = channel_units[ext_key]

        return ParsedSpectrum(data=data_array, raw=None, metadata=metadata)

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
        for encoding in ("utf-8", "cp1252"):
            try:
                with open(file, encoding=encoding) as xy_file:
                    return xy_file.readlines()
            except UnicodeDecodeError:
                continue
        raise ValueError("Unable to decode file with known encodings.")

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

        export_settings = {}
        for line in header:
            line = line.strip(self.prefix).strip()
            if len(line) == 0:
                pass
            else:
                key, setting_str = line.split(":", 1)
                key, setting_str = key.strip(), setting_str.strip()
                key, setting, _ = _context.format(key, setting_str)
                export_settings[key] = _convert_bool(setting)

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

        def _parse_group_header(group_data: list[str]):
            sections: dict[str, list[str]] = OrderedDict()

            current_section = "group_settings"
            sections[current_section] = []

            region_start = 0

            for i, line in enumerate(group_data):
                stripped = line.strip(self.prefix).strip()

                if stripped.startswith("Region:"):
                    region_start = i
                    break

                # ignore empty comment lines
                if not stripped:
                    continue

                # section header: no colon
                if ":" not in stripped:
                    current_section = convert_pascal_to_snake(stripped)
                    sections[current_section] = []
                    continue

                sections[current_section].append(line)

            return sections, group_data[region_start:]

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

            sections, region_data = _parse_group_header(group_data)
            group = self._handle_regions(group_data)
            for section_name, lines in sections.items():
                section_settings = self._handle_settings_block(lines)
                if section_name == "group_settings":
                    group_settings.update(section_settings)
                elif section_name.endswith("_parameters"):
                    prefix = section_name.removesuffix("_parameters")
                    group_settings.update(
                        {
                            f"{prefix}_{key}": value
                            for key, value in section_settings.items()
                        }
                    )

            group["group_settings"] = group_settings
            groups[name] = group

        return groups

    def _handle_settings_block(self, lines: list[str]) -> dict[str, _Value]:
        settings = {}

        for line in lines:
            content = line.split(self.prefix)[-1].strip()
            key, *rest = content.split(":", 1)

            if not rest:
                continue

            raw_value: _Value = rest[0].strip()

            key, value, unit = _context.format(key, raw_value)

            settings[key] = value
            if unit:
                settings[f"{key}/@units"] = unit

        return settings

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

            for i, content in enumerate(region_data):
                key, *rest = content.split(self.prefix)[-1].strip().split(":", 1)

                if not rest:
                    continue

                raw_value: _Value = rest[0].strip()
                key, value, unit = _context.format(key, raw_value)

                region_settings[key] = value
                if unit:
                    region_settings[f"{key}/@units"] = unit

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
        if self.export_settings["separate_scan_data"]:
            spec_pattern_str += r", Scan: \d"
        if self.export_settings["separate_channel_data"]:
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
                unit = _context.map_unit(unit)

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

        if self.export_settings["external_channel_data"]:
            remote_ch_pattern_str = rf"{re.escape(self.prefix)} external_channel_data Cycle: \d+,\s*(.+?)\s*\(Remote Out Device\)"
            remote_ch_pattern = re.compile(remote_ch_pattern_str, re.IGNORECASE)

        for line in scan_data:
            if match := remote_ch_pattern.search(line):
                in_remote_channel = True
                channel = _normalize_ext_channel_label(match.group(1))[0]

            if line.startswith(self.prefix) and line.strip(self.prefix).strip():
                key, value_str = (
                    item.strip()
                    for item in line.strip(self.prefix).strip().split(":", 1)
                )
                key, value, _ = _context.format(key, value_str)
                if key == "time_stamp":
                    scan_settings[key] = value

                if key == "column_labels":
                    if in_remote_channel:
                        column_labels = str(value).strip().split(" ", 1)
                        column_labels = [
                            label.split("(Remote Out Device)", 1)[0].strip()
                            for label in column_labels
                        ]
                    else:
                        column_labels = str(value).strip().split(" ")
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
            scan_settings["step_size"] = _get_minimal_step(energy_axis)
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
