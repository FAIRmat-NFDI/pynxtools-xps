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
Abstract base classes for mappers and parsers.
"""

import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

from pynxtools_xps.logging import _logger
from pynxtools_xps.mapping import _Value


@dataclass
class _XPSDataclass:
    """Generic class to hold a data model and a type validation method."""

    def validate_types(self):
        ret = True
        for field_name, field_def in self.__dataclass_fields__.items():
            actual_type = type(getattr(self, field_name))
            if actual_type != field_def.type:
                _logger.warning(
                    f"Type mismatch in dataclass {type(self).__name__}. {field_name}: '{actual_type}' instead of '{field_def.type}'"
                )
                ret = False
        return ret

    def __post_init__(self):
        if not self.validate_types():
            raise ValueError(f"Type mismatch in dataclass {type(self).__name__}")

    def dict(self):
        return self.__dict__.copy()


class _XPSMapper(ABC):
    """Abstract base class from mapping from a parser to NXmpes template"""

    def __init__(self):
        self.file: str | Path = ""
        self._data: dict[str, Any] = {}
        self._data["data"] = cast(dict[str, _Value], {})
        self.parser = None

    @abstractmethod
    def _select_parser(self):
        """
        Select the correct parser for the file extension and format.

        Should be implemented by the inheriting mapper.
        """

    @property
    def data(self) -> dict:
        """Getter property."""
        return self._data

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the file using the selected parser.

        """
        self.file = file
        self.parser = self._select_parser()
        parsed_data = self.parser.parse_file(file, **kwargs)

        self._data["File"] = file
        self._data["file_ext"] = os.path.splitext(file)[1]

        self.construct_data(parsed_data)

        return self.data

    @abstractmethod
    def construct_data(self, parsed_data: list[dict[str, Any]]) -> dict[str, _Value]:
        """
        Map from individual parser format to NXmpes-ready dict.

        Should be implemented by the inheriting mapper.

        """


class _XPSParser(ABC):
    def __init__(self):
        self.file: str | Path = ""
        self._data: list[dict[str, Any]] = []

    @property
    def data(self) -> list[dict[str, Any]]:
        """Getter property."""
        return self._data

    @abstractmethod
    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the file.

        """
        self.file = file


def _construct_data_key(spectrum: dict[str, Any]) -> str:
    """
    Construct a key for the 'data' field of the xps_dict.
    Output example: cycle0_scan0.

    """
    if "loop_no" in spectrum:
        cycle_key = f"cycle{spectrum['loop_no']}"
    else:
        cycle_key = "cycle0"

    if "scan_no" in spectrum:
        scan_key = f"scan{spectrum['scan_no']}"
    else:
        scan_key = "scan0"

    return f"{cycle_key}_{scan_key}"


def _align_name_part(name_part: str):
    """Make one part of the entry name compliant with NeXus standards."""
    translation_table = str.maketrans(
        {
            " ": "_",
            ",": "",
            ".": "_",
            "-": "_",
            ":": "_",
            "+": "_",
            "/": "_",
            "=": "",
        }
    )

    return name_part.translate(translation_table)


def _construct_entry_name(parts: list[str]) -> str:
    """Construct name for the NXentry instances."""
    if len(parts) == 1:
        return _align_name_part(parts[0])
    return "__".join([_align_name_part(part) for part in parts])
