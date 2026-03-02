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
Base classes and typed intermediate representations for XPS parsers.
"""

import os
import types
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar, Union, cast, get_args, get_origin

import xarray as xr

from pynxtools_xps.logging import _logger
from pynxtools_xps.parsers.versioning import (
    VersionRange,
    VersionTuple,
    _format_version,
    is_version_supported,
)


@dataclass
class _XPSDataclass:
    """Base class for XPS dataclasses with runtime type enforcement and validation.

    Supports:
        - Plain types (int, float, str, ...)
        - PEP 604 unions (e.g. int | None)
        - typing.Union
        - list[T]
        - Any
    """

    def __setattr__(self, name: str, value: Any) -> None:
        """Intercept attribute assignment and enforce annotated types."""
        fields = getattr(self, "__dataclass_fields__", {})

        # Non-dataclass attributes → default behavior
        if name not in fields:
            super().__setattr__(name, value)
            return

        expected_type = fields[name].type
        coerced_value = self._coerce_value(value, expected_type, field_name=name)

        super().__setattr__(name, coerced_value)

    @classmethod
    def _coerce_value(
        cls,
        value: Any,
        expected_type: Any,
        *,
        field_name: str | None = None,
    ) -> Any:
        """Coerce a value to the expected type if possible.

        Args:
            value: The value being assigned.
            expected_type: The annotated type.
            field_name: Optional field name for error reporting.

        Returns:
            The validated or coerced value.

        Raises:
            TypeError: If coercion is not possible.
        """

        # Any always allowed
        if expected_type is Any:
            return value

        # Already valid → nothing to do
        if cls._is_instance_of_type(value, expected_type):
            return value

        origin = get_origin(expected_type)

        # Union / PEP 604
        if origin in (Union, types.UnionType):
            for arg in get_args(expected_type):
                if arg is type(None) and value is None:
                    return None

                # Try validation first
                if cls._is_instance_of_type(value, arg):
                    return value

                # Try coercion
                try:
                    return cls._coerce_value(value, arg)
                except TypeError:
                    continue

            raise TypeError(
                f"Cannot assign value {value!r} to field "
                f"{field_name!r}: expected {expected_type}"
            )

        # list[T]
        if origin is list:
            if not isinstance(value, list):
                raise TypeError(
                    f"Field {field_name!r} expects list, got {type(value).__name__}"
                )

            (inner_type,) = get_args(expected_type)
            return [
                cls._coerce_value(v, inner_type, field_name=field_name) for v in value
            ]

        # Plain type
        try:
            return expected_type(value)
        except Exception as exc:
            raise TypeError(
                f"Cannot assign value {value!r} to field "
                f"{field_name!r}: expected {expected_type}"
            ) from exc

    @classmethod
    def _is_instance_of_type(cls, value: Any, expected_type: Any) -> bool:
        """Check whether a value conforms to a possibly generic type."""

        if expected_type is Any:
            return True

        origin = get_origin(expected_type)

        # Union / PEP 604
        if origin in (Union, types.UnionType):
            return any(
                cls._is_instance_of_type(value, arg) for arg in get_args(expected_type)
            )

        # list[T]
        if origin is list:
            if not isinstance(value, list):
                return False
            (inner_type,) = get_args(expected_type)
            return all(cls._is_instance_of_type(v, inner_type) for v in value)

        # Plain type
        return isinstance(value, expected_type)

    def validate_types(self) -> bool:
        """Validate all dataclass fields against their type annotations.

        Returns:
            bool: True if all fields match their annotated types,
                  False otherwise.
        """
        valid = True

        for field_name, field_def in self.__dataclass_fields__.items():
            value = getattr(self, field_name)
            expected_type = field_def.type

            if not self._is_instance_of_type(value, expected_type):
                _logger.warning(
                    "Type mismatch in %s.%s: got %s, expected %s",
                    type(self).__name__,
                    field_name,
                    type(value).__name__,
                    expected_type,
                )
                valid = False

        return valid

    def __post_init__(self):
        if not self.validate_types():
            raise ValueError(f"Type mismatch in dataclass {type(self).__name__}")

    def dict(self):
        return self.__dict__.copy()


# TODO: fill this with life
# @dataclass
# class ScanData(_XPSDataclass):
#     """One scan's measured data."""

#     scan_id: int
#     data: xr.Dataset
#     channel_data: xr.Dataset | None = None
#     raw_data: xr.Dataset | None = None
#     attrs: dict[str, Any] = field(default_factory=dict)

#     def validate(self) -> None:
#         """Check internal consistency of lengths."""
#         if len(self.intensity) != n:
#             raise ValueError(
#                 f"Scan {self.scan_id}: intensity length "
#                 f"{len(self.intensity)} != energy length {n}"
#             )
#         if scan.channels is not None and self.channels.shape[0] != n:
#             raise ValueError(
#                 f"Scan {self.scan_id}: channels shape "
#                 f"{self.channels.shape} incompatible with energy length {n}"
#             )


# @dataclass
# class ParsedSpectrum(_XPSDataclass):
#     """Typed intermediate representation produced by every XPS parser.

#     All parsers return ``list[ParsedSpectrum]``.  The reader assembles
#     these into the ENTRY-keyed dict + xarray datasets consumed by the
#     config-based template mapping.
#     """

#     group_name: str
#     spectrum_type: str
#     scans: list[ScanData]
#     metadata: dict[str, Any] = field(default_factory=dict)
#     scan_no: int | None = None
#     loop_no: int | None = None
#     energy_units: str = "eV"
#     intensity_units: str = "counts_per_second"


class _Parser(ABC):
    """
    Abstract base class for all XPS file + metadata parsers.

    Subclasses define the structural and semantic rules required to
    identify and parse a specific XPS file format variant or its metadata.

    Class Attributes:
        supported_file_extensions: Tuple of supported file extensions
            (e.g., (".sle", ".txt")).
        supported_vendors; Tuple of supported instrument vendors.
        requires_version: If True, the file must explicitly provide
            version information.
        supported_versions: Tuple of supported version ranges. If empty,
            all versions are accepted unless `requires_version` is True.

    ``supported_file_extensions`` and ``supported_vendors`` must be
    set by subclasses, ``requires_version`` and ``supported_versions`` can
    be optionally set for version-aware parsers.

    Subclasses must implement:
        - matches_file()
        - _parse()

    """

    supported_file_extensions: ClassVar[tuple[str, ...]] = ()
    supported_vendors: ClassVar[tuple[str, ...]] = ()
    # TODO: implement in subclasses
    # "kratos", "phi", "scienta", "specs", "unknown"]
    requires_version: ClassVar[bool] = False
    supported_versions: ClassVar[tuple[VersionRange, ...]] = ()

    @classmethod
    def file_ext_err_msg(cls, file: Path) -> str:
        suffix = file.suffix or "<no extension>"
        allowed = ", ".join(cls.supported_file_extensions) or "<none>"

        return (
            f"Cannot process file '{file.name}' (detected extension: '{suffix}'). "
            f"{cls.__name__} supports only the following file extensions: {allowed}."
        )

    @classmethod
    def file_version_err_msg(
        cls,
        file: str | Path,
        version: VersionTuple | None,
    ) -> str:

        path = Path(file)
        file_name = path.name
        version_str = _format_version(version) if version is not None else "<unknown>"

        if version is None and cls.requires_version:
            return (
                f"File '{file_name}' does not specify a version, "
                f"but {cls.__name__} requires explicit version information."
            )

        if not cls.supported_versions:
            return f"Unsupported file version '{version_str}' for file '{file_name}'."

        ranges: list[str] = []
        for lower, upper in cls.supported_versions:
            lower_str = _format_version(lower)
            if upper is None:
                ranges.append(f">= {lower_str}")
            else:
                upper_str = _format_version(upper)
                ranges.append(f"{lower_str} – {upper_str}")

        supported_str = ", ".join(ranges)

        return (
            f"Unsupported file version '{version_str}' for file '{file_name}'. "
            f"Supported versions: {supported_str}."
        )

    @classmethod
    def is_extension_supported(cls, file: Path) -> bool:
        suffix = file.suffix.lower()
        return any(suffix == ext.lower() for ext in cls.supported_file_extensions)

    @classmethod
    def is_version_supported(
        cls,
        version: VersionTuple | None,
    ) -> bool:
        return is_version_supported(
            version,
            cls.supported_versions,
            requires_version=cls.requires_version,
        )

    @classmethod
    def is_mainfile(cls, file: str | Path) -> bool:
        """
        Check whether this parser supports the given file.

        This performs extension and structural validation
        without raising exceptions.

        Args:
            file: Path to the candidate file.

        Returns:
            True if the file can be parsed by this class, otherwise False.
        """
        try:
            parser = cls()
            parser._is_mainfile(Path(file))
            return True
        except ValueError:
            return False

    def __init__(self) -> None:
        self.file: Path
        self._data: list[dict[str, Any]] = []

    @property
    def data(self) -> list[dict[str, Any]]:
        return self._data

    def _is_mainfile(self, file: Path) -> None:
        """
        Check whether this parser supports the given file.

        This performs extension, version, and structural validation
        and raises exceptions.

        Args:
            file: Path to the candidate file.

        """
        if not self.is_extension_supported(file):
            raise ValueError(self.file_ext_err_msg(file))

        version = self.detect_version(file)

        if not self.is_version_supported(version):
            raise ValueError(self.file_version_err_msg(file, version))

        if not self.matches_file(file):
            raise ValueError(
                f"File '{file.name}' does not match the expected "
                f"format for {self.__class__.__name__}."
            )

    def detect_version(self, file: Path) -> VersionTuple | None:
        """
        Detect the version of the given file.

        This method may be overridden by subclasses that extract
        version information from file headers or metadata.

        Args:
            file: Path to the file.

        Returns:
            A version tuple (e.g., (4, 63, 1)) if detected, otherwise None.
        """
        return None

    # TODO: this should be an abstract method
    # @abstractmethod
    def matches_file(self, file: Path) -> bool:
        """
        Determine whether the file structurally matches this parser.

        This method must perform strict structural validation beyond
        extension and version checks. It should return True only if the
        file unambiguously conforms to the expected format.

        Args:
            file: Path to the candidate file.

        Returns:
            True if the file matches this parser's format, otherwise False.
        """
        ...
        # TODO: this should be an abstract method
        return True

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the given file and populate the parser's data attribute.

        Args:
            file: Path to the file to parse.
            **kwargs: Additional parser-specific keyword arguments.

        Raises:
            ValueError: If the file is not supported by this parser.
        """
        file = Path(file)
        self.file = file
        self._is_mainfile(file)
        self._parse(file, **kwargs)

    @abstractmethod
    def _parse(self, file: Path, **kwargs) -> None:
        """
        Perform the actual parsing implementation.

        Subclasses must implement this method to extract structured
        data from the validated file and populate `self._data`.

        Args:
            file: Path to the validated file.
            **kwargs: Additional parser-specific options.
        """
        ...


class _XPSParser(_Parser):
    """
    Abstract base class for all XPS file parsers.

    Subclasses define the structural and semantic rules required to
    identify and parse a specific XPS file format variant.

    Additional class attributes:
        config_file: Path or mapping used to configure the parser.

    ``config_file`` must be set by subclasses.
    """

    config_file: ClassVar[str] = ""


class _XPSMetadataParser(_Parser):
    """
    Abstract base class for supplementary parsers that enrich
    existing spectral datasets.

    Metadata parsers do not generate spectra. Instead, they parse
    auxiliary files (e.g., quantification exports) and inject
    additional information into data dictionaries produced by
    a primary `_XPSParser`.

    Additional class attributes:
        compatible_primary_parser: Tuple of compatible _XPSParser
            objects or mapping used to configure the parser.

    ``compatible_primary_parser`` must be set by subclasses.

    """

    compatible_primary_parser: ClassVar[type[_XPSParser]]

    @classmethod
    def file_ext_err_msg(cls, file: Path) -> str:
        suffix = file.suffix or "<no extension>"
        allowed = ", ".join(cls.supported_file_extensions) or "<none>"

        return (
            f"Cannot process file '{file.name}' (detected extension: '{suffix}'). "
            f"{cls.__name__} supports only the following file extensions: {allowed}."
        )

    @classmethod
    def supports_parser(cls, parser: _XPSParser) -> bool:
        return isinstance(parser, cls.compatible_primary_parser)

    def _supports_parser(self, parser: _XPSParser):
        if not self.supports_parser(parser):
            raise ValueError(
                f"{self.__class__.__name__} is not compatible with "
                f"{parser.__class__.__name__}."
            )

    # TODO: this shall be an abstract method
    # @abstractmethod
    def update_main_file_data(
        self,
        main_file_data: list[dict[str, Any]],
    ) -> None:
        """
        Inject parsed metadata into already-parsed spectral data.

        Implementations must validate compatibility with the provided
        data before modifying it.

        Args:
            main_file_dicts: List of dictionaries produced by a primary
                `_XPSParser`.

        Raises:
            ValueError: If the metadata cannot be aligned with the
                provided spectral data.
        """
        ...


class _XPSMapper(ABC):
    """
    Deprecated: retained only for backward compatibility during migration.
    """

    def __init__(self):
        self.file: Path = ""
        self._data: dict[str, Any] = {}
        self._data["data"] = cast(dict[str, xr.Dataset], {})
        self.parser = None

    @property
    def data(self) -> dict:
        """Getter property."""
        return self._data

    @abstractmethod
    def _select_parser(self):
        """Select the correct parser for the file extension and format."""

    def parse_file(self, file: str | Path, **kwargs):
        file = Path(file)
        self.file = file
        self.parser = self._select_parser()

        self._data["File"] = file
        self._data["file_ext"] = os.path.splitext(file)[1]
        self.parser.parse_file(file, **kwargs)
        self.construct_data(self.parser.data)

    @abstractmethod
    def construct_data(self, parsed_data: list[dict[str, Any]]):
        """Map from individual parser format to NXmpes-ready dict."""


# TODO: this should be implemented when deprecating mappers
def _construct_data_key(spectrum: dict[str, Any]) -> str:
    """
    Construct a key for the 'data' field of the xps_dict.

    For ParsedSpectrum: returns ``cycle0`` when scan_no is None,
    or ``cycle0_scan3`` when scan_no is set.  This lets parsers with
    all scans inside one ParsedSpectrum (PHI, SLE) use just the cycle
    key, while parsers that emit one ParsedSpectrum per scan (VAMAS)
    include the scan index.

    For legacy dicts: always returns ``cycle{N}_scan{M}`` (backward
    compatible with old mapper code).
    """
    # if isinstance(spectrum, ParsedSpectrum):
    #     loop_no = spectrum.loop_no
    #     scan_no = spectrum.scan_no
    #     cycle_key = f"cycle{loop_no}" if loop_no is not None else "cycle0"
    #     if scan_no is not None:
    #         return f"{cycle_key}_scan{scan_no}"
    #     return cycle_key
    # else:
    loop_no = spectrum.get("loop_no")
    scan_no = spectrum.get("scan_no")
    cycle_key = f"cycle{loop_no}" if loop_no is not None else "cycle0"
    scan_key = f"scan{scan_no}" if scan_no is not None else "scan0"
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
    parts = [p for p in parts if p]
    if not parts:
        return ""
    if len(parts) == 1:
        return _align_name_part(parts[0])
    return "__".join([_align_name_part(part) for part in parts])
