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
Classes for reading XPS files from TXT export of CasaXPS.
"""

import csv
import re
from pathlib import Path

from pynxtools_xps.mapping import _Value
from pynxtools_xps.parsers.base import _XPSMapper, _XPSParser
from pynxtools_xps.parsers.vms_export.metadata import _context


class VamasResultMapper(_XPSMapper):
    """
    Class for restructuring .csv result files from
    Casa report export (from Vamas) into python dictionary.
    """

    config_file = "config_vms.json"

    def __init__(self):
        super().__init__()

    def _select_parser(self):
        """
        Select parser based on the structure of the text file

        Returns
        -------
        TextParser
            TextParser for CasaXPS export from Vamas files.

        """
        return CSVResultParser()

    def construct_data(self, raw_data: list[dict[str, _Value]]):
        pass

    def update_main_file_dict(self, main_file_dicts: list[dict[str, _Value]]):
        """
        Update the dictionaries returned by the main files with specific keys from self.data.

        Args:
            main_file_dicts (List[list[str, _Value]]): List of dictionaries to update.
        """
        pattern = re.compile(r"(component\d+/)name")
        update_with = {
            "Area/(RSF*T*MFP)",
            "atomic_concentration",
        }

        for existing_dict in main_file_dicts:
            filtered_keys = {
                key: match.group(1)
                for key in existing_dict
                if (match := pattern.search(key))
            }

            for key in filtered_keys:
                value = existing_dict[key]
                if value in self.data:
                    sub_dict = self.data[value]
                    for sub_key in update_with & sub_dict.keys():
                        new_key = f"{key.rsplit('name', 1)[0]}{sub_key}"
                        existing_dict[new_key] = sub_dict[sub_key]


class CSVResultParser(_XPSParser):
    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse only the first table from the input file,

        Args:
            file_path (str): Path to the .vms file.

        Returns:
            dict: Parsed data including the file path, header, and rows.
        """
        table_data: dict[str, dict[str, _Value]] = {}
        headers: list[str] = []
        reading_table: bool = False

        with open(file) as f:
            reader = csv.reader(f, delimiter="\t")

            for row in reader:
                if not row:
                    if reading_table:
                        break
                    continue

                # Detect header row
                if row[0].startswith("Name") and not reading_table:
                    headers = [_context.normalize_key(h) for h in row[1:] if h]
                    reading_table = True
                    continue

                # Process rows of the table
                if reading_table:
                    table_data[row[0]] = {}
                    for header, value in zip(headers, row[1:]):
                        if value:
                            formatted_value = _context._format_value(value)
                            if header == "atomic_concentration" and isinstance(
                                formatted_value, int | float
                            ):
                                formatted_value /= 100
                            table_data[row[0]][header] = formatted_value

        return table_data
