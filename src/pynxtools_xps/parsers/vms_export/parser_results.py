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
Metadata parser for reading result files from CasaXPS CSV export (from Vamas).
"""

import csv
import re
from pathlib import Path
from typing import Any, ClassVar

from pynxtools_xps.parsers.base import (
    ParsedSpectrum,
    VendorType,
    _construct_entry_name,
    _XPSMetadataParser,
    _XPSParser,
)
from pynxtools_xps.parsers.vms.parser import VamasParser
from pynxtools_xps.parsers.vms_export.metadata import _context


class VamasResultParser(_XPSMetadataParser):
    """
    Metadata-only parser for CasaXPS quantification CSV exports.

    Parses a ``.csv`` result file produced by CasaXPS and injects
    fitting/quantification results into the ``metadata`` of matching
    ``ParsedSpectrum`` objects produced by ``VamasExportParser``.

    ``_parse()`` reads Table 1 (component quantities) and the first compact
    table (Sample Identifier per component), groups components by
    ``(sample_id, peak_name)``, and stores one metadata-only
    ``ParsedSpectrum`` per group in ``self._data``.

    Entry keys are built with :func:`_construct_entry_name` and therefore
    match the keys produced by ``VamasExportParser``.  If the CSV lacks a
    Sample Identifier, the base-class suffix matching in
    ``_matches_entry`` still aligns the entries correctly.

    Compatible primary parser: ``VamasExportParser``.
    """

    compatible_primary_parser: ClassVar[type[_XPSParser]] = VamasParser
    supported_vendor: ClassVar[VendorType | None] = "various"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".csv",)

    def matches_file(self, file: Path) -> bool:
        """Return True for CasaXPS quantification CSV exports."""
        try:
            with open(file, encoding="utf-8", errors="ignore") as f:
                head = f.read(2048)
            return "Goodness of Fit" in head
        except Exception:
            return False

    def _parse(self, file: Path, **kwargs) -> None:
        """Parse quantification CSV and populate ``self._data``."""
        with open(file) as f:
            all_rows = list(csv.reader(f, delimiter="\t"))

        table1_rows = self._parse_table1(all_rows)
        compact_rows = self._parse_compact_table(all_rows, n_expected=len(table1_rows))
        self._build_entries(table1_rows, compact_rows)

    def _parse_table1(self, all_rows: list[list[str]]) -> list[dict[str, Any]]:
        """Extract ordered component data from Table 1.

        Table 1 is the first ``Name``-headed block in the CSV.  It contains
        per-component quantities: ``position``, ``fwhm``, ``raw_area``,
        ``area_over_rsf_t_mfp``, ``atomic_concentration``, and
        ``goodness_of_fit``.

        Returns a list of dicts (one per component row), each with a
        ``"_peak_name"`` key holding the raw peak name from column 0.
        """
        result: list[dict[str, Any]] = []
        headers: list[str] = []
        reading = False

        for row in all_rows:
            is_blank = not row or not any(c.strip() for c in row)

            if is_blank:
                if reading:
                    break
                continue

            if row[0].strip() == "Name" and not reading:
                headers = [_context.normalize_key(h) for h in row[1:] if h]
                reading = True
                continue

            if reading:
                row_data: dict[str, Any] = {"_peak_name": row[0].strip()}
                for header, value in zip(headers, row[1:]):
                    if value.strip():
                        formatted = _context._format_value(value)
                        if header == "atomic_concentration" and isinstance(
                            formatted, int | float
                        ):
                            formatted /= 100
                        row_data[header] = formatted
                result.append(row_data)

        return result

    def _parse_compact_table(
        self,
        all_rows: list[list[str]],
        n_expected: int,
    ) -> list[tuple[str, str]]:
        """Extract ``(peak_name, sample_id)`` pairs from the compact table.

        The compact table follows Table 1 and has a header containing
        ``"Sample Identifier"``.  Data rows alternate with St.Dev. rows
        (where ``row[0]`` is empty); only data rows are collected.

        The Sample Identifier propagates forward: once set for a row, it
        remains in effect until a new non-empty identifier appears.

        Collection stops after *n_expected* data rows are gathered, which
        prevents running into the second compact table.
        """
        result: list[tuple[str, str]] = []
        reading = False
        sample_id = ""

        for row in all_rows:
            if not row or not any(c.strip() for c in row):
                continue  # skip fully blank lines

            if row[0].strip() == "Name" and any("Sample Identifier" in c for c in row):
                reading = True
                continue

            if not reading:
                continue

            # Skip sub-header row (e.g. "\tSt.Dev.") and St.Dev. data rows
            if not row[0].strip():
                continue

            peak_name = row[0].strip()
            new_id = row[2].strip() if len(row) > 2 else ""
            if new_id:
                sample_id = new_id

            result.append((peak_name, sample_id))

            if len(result) == n_expected:
                break

        return result

    def _build_entries(
        self,
        table1_rows: list[dict[str, Any]],
        compact_rows: list[tuple[str, str]],
    ) -> None:
        """Group aligned rows into ``self._data``.

        Rows are grouped by ``(sample_id, peak_name)``; a new group starts
        whenever either value changes.  Within each group,
        :func:`_handle_repetitions` disambiguates repeated peak names
        (e.g. four ``"Fe 2p"`` peaks become ``Fe_2p_1`` … ``Fe_2p_4``).

        Each group becomes one ``ParsedSpectrum(data=None)`` entry in
        ``self._data``, keyed by
        ``_construct_entry_name([sample_id, peak_name])`` or
        ``_construct_entry_name([peak_name])`` when no sample identifier
        is present.  Component field values are stored under
        ``{comp_label}/{field}`` keys in ``metadata``.
        """
        # Align table1 and compact by index; fall back to raw peak name if shorter
        aligned: list[tuple[str, str, dict[str, Any]]] = []
        for i, row_data in enumerate(table1_rows):
            if i < len(compact_rows):
                peak_name, s_id = compact_rows[i]
            else:
                peak_name = row_data["_peak_name"]
                s_id = ""
            aligned.append((peak_name, s_id, row_data))

        # Group by (sample_id, peak_name) — new group on any change
        groups: dict[str, list[dict[str, Any]]] = {}
        group_order: list[str] = []

        for peak_name, s_id, row_data in aligned:
            peak_name = peak_name.replace(" ", "")
            parts = [s_id, peak_name] if s_id else [peak_name]
            entry_name = _construct_entry_name(parts)
            if entry_name not in groups:
                groups[entry_name] = []
                group_order.append(entry_name)
            groups[entry_name].append(row_data)

        # Build ParsedSpectrum per group
        for entry_name in group_order:
            group_rows = groups[entry_name]
            comp_labels = [r["_peak_name"] for r in group_rows]

            metadata: dict[str, Any] = {}
            for i, (comp_label, row_data) in enumerate(zip(comp_labels, group_rows)):
                for field, value in row_data.items():
                    if field == "_peak_name":
                        continue
                    metadata[f"component{i}/name"] = comp_label
                    metadata[f"component{i}/{field}"] = value

            self._data[entry_name] = ParsedSpectrum(
                data=None, raw=None, metadata=metadata
            )

    def update_main_file_data(self, main_file_data: dict[str, ParsedSpectrum]) -> None:
        """
        Merge ``self._data`` metadata into matching entries of *main_file_data*.

        For each spectrum, the method searches its ``metadata`` for keys that
        match the pattern ``component{N}/name``.  When a component name matches
        an entry in the CSV data, the quantification fields
        ``area_over_rsf*t*mfp``, ``atomic_concentration``, and
        ``goodness_of_fit`` are injected into the spectrum's metadata under
        the path ``component{N}/{field}``.

        Args:
            main_file_data: Mapping from NeXus entry name to ``ParsedSpectrum``,
                as produced by the compatible primary parser.
        """
        if not self._data:
            return

        pattern = re.compile(r"(component\d+/)")
        _inject_fields = {
            "area_over_rsf*t*mfp",
            "atomic_concentration",  # overwrites the 0.0 placeholder
            "goodness_of_fit",
            "raw_area",
        }
        for meta_entry, meta_spectrum in self._data.items():
            for main_entry, main_spectrum in main_file_data.items():
                if self._matches_entry(meta_entry, main_entry):
                    for key, value in meta_spectrum.metadata.items():
                        field = key.rsplit("/", 1)[
                            -1
                        ]  # "component0/raw_area" → "raw_area"
                        if field in _inject_fields:
                            main_spectrum.metadata[key] = value
                    break
