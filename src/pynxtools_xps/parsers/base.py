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

import types
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar, Union, get_args, get_origin

import xarray as xr

from pynxtools_xps.logging import _logger
from pynxtools_xps.parsers.versioning import (
    VersionRange,
    VersionTuple,
    _format_version,
    _version_ranges_overlap,
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


def _indent(text: str, n: int) -> str:
    """Prepend *n* spaces to every line of *text*."""
    prefix = " " * n
    return "\n".join(prefix + line for line in text.splitlines())


@dataclass
class ParsedSpectrum:
    """Typed intermediate representation for one XPS spectrum.

    Parsers return ``dict[str, ParsedSpectrum]`` where keys are NeXus entry
    names (e.g. ``"SampleName__Survey"``). The ``XPSReader`` accesses spectra
    directly through ``get_attr`` and ``get_data`` callbacks.

    Physical hierarchy::

        DETECTOR/raw_data              (cycle, scan, channel, *axes)
        PROCESS[channels_averaging]    (cycle, scan, *axes)   ← data field
        PROCESS[scan_averaging]        (cycle, *axes)          ← computed
        PROCESS[cycle_averaging]       (*axes,)                ← computed
        DATA[data]                     (*axes,) + errors       ← computed

    Channel averaging is parser-specific (detector geometry, MCD calibration)
    and must be performed by the parser before constructing this object.
    The remaining averaging steps are computed generically by the assembly layer.

    Attributes:
        data: Channel-averaged scan data, or ``None`` for metadata-only entries
            produced by ``_XPSMetadataParser`` subclasses.
            When not ``None``, required dims are ``("cycle", "scan")`` followed
            by one or more physical axes (typically ``"energy"``).
            Use ``n_cycles=1`` for formats without an explicit loop structure.
        raw:  Optional raw per-channel data.
            Required dims: ``("cycle", "scan", "channel")``, followed by the
            same physical axes as ``data``.
        metadata: Flat key-value metadata for ``@attrs:`` lookups in config
            files.  Keys should follow the same snake_case convention as the
            rest of the parser output.
    """

    data: xr.DataArray | None = None
    raw: xr.DataArray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        _BOOKKEEPING_DIMS = frozenset({"cycle", "scan", "channel"})

        def _da_summary(da: xr.DataArray | None, label: str) -> str:
            if da is None:
                return f"  {label:<6} None"
            dims_str = "(" + ", ".join(f"{d}={da.sizes[d]}" for d in da.dims) + ")"
            coord_parts = []
            for coord in da.coords:
                if coord in _BOOKKEEPING_DIMS:
                    continue
                vals = da.coords[coord].values
                unit = self.metadata.get(f"{coord}/@units", "")
                unit_str = f" {unit}" if unit else ""
                coord_parts.append(
                    f"{coord}=[{vals.min():.4g} \u2026 {vals.max():.4g}]{unit_str}"
                )
            coord_str = ("  " + "  ".join(coord_parts)) if coord_parts else ""
            return f"  {label:<6} {dims_str}{coord_str}"

        lines = ["ParsedSpectrum"]
        lines.append(_da_summary(self.data, "data:"))
        lines.append(_da_summary(self.raw, "raw:"))
        n = len(self.metadata)
        lines.append(f"  metadata ({n} keys):")
        for k, v in sorted(self.metadata.items()):
            v_str = str(v)
            if len(v_str) > 60:
                v_str = v_str[:57] + "..."
            lines.append(f"    {k:<50}  {v_str}")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Computed aggregations (used by the reader's get_data)
    # ------------------------------------------------------------------

    def _require_data(self) -> xr.DataArray:
        """Return ``self.data``, raising if this is a metadata-only entry."""
        if self.data is None:
            raise ValueError(
                "This ParsedSpectrum is metadata-only (data=None). "
                "Aggregation methods are not available."
            )
        return self.data

    def average(self) -> xr.DataArray:
        """Mean across all cycles and scans. Shape: (*axes,)."""
        return self._require_data().mean(dim=["cycle", "scan"])

    def errors(self) -> xr.DataArray:
        """Std across all cycles and scans. Shape: (*axes,)."""
        return self._require_data().std(dim=["cycle", "scan"])

    def scan_average(self) -> xr.DataArray:
        """Mean across scans within each cycle. Shape: (cycle, *axes)."""
        return self._require_data().mean(dim="scan")

    def cycle_average(self) -> xr.DataArray:
        """Mean across cycles after scan averaging. Shape: (*axes,)."""
        return self.scan_average().mean(dim="cycle")

    def validate(self) -> None:
        """Check that required dimensions are present and shapes are consistent."""
        data = self._require_data()
        required = ("cycle", "scan")
        for dim in required:
            if dim not in data.dims:
                raise ValueError(
                    f"ParsedSpectrum.data must have a '{dim}' dimension; "
                    f"got dims: {tuple(data.dims)}"
                )

        if self.raw is not None:
            required_raw = ("cycle", "scan", "channel")
            for dim in required_raw:
                if dim not in self.raw.dims:
                    raise ValueError(
                        f"ParsedSpectrum.raw must have a '{dim}' dimension; "
                        f"got dims: {tuple(self.raw.dims)}"
                    )
            # Physical axes must match
            data_axes = [d for d in data.dims if d not in ("cycle", "scan")]
            raw_axes = [
                d for d in self.raw.dims if d not in ("cycle", "scan", "channel")
            ]
            if data_axes != raw_axes:
                raise ValueError(
                    f"ParsedSpectrum.data and .raw must share the same physical "
                    f"axes; data has {data_axes}, raw has {raw_axes}"
                )


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
        self._data: dict[str, ParsedSpectrum] = {}

    def __repr__(self) -> str:
        try:
            file_str = f"{self.file.name}"
        except AttributeError:
            file_str = ""
        n = len(self._data)
        lines = [f"{self.__class__.__name__}"]
        lines.append(f"File name: {file_str}")
        lines.append(f"Number of parsed entries: {n}")
        lines.append(f"Entries:")

        for entry_name, spectrum in self._data.items():
            lines.append(f"    '{entry_name}':")
            lines.append(_indent(repr(spectrum), 8))
        return "\n".join(lines)

    @property
    def data(self) -> dict[str, ParsedSpectrum]:
        """Parsed spectra (or metadata-only entries) keyed by NeXus entry name."""
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

    @abstractmethod
    def matches_file(self, file: Path) -> bool:
        """
        Return True if *file* structurally matches this parser's format.

        Implementations must perform positive identification — not just
        extension checks or negative exclusions. The check should be fast
        (read at most a few KB), and always catch all exceptions and return
        False rather than propagating them.

        Args:
            file: Path to the candidate file.

        Returns:
            True if the file matches this parser's format, otherwise False.
        """

    def parse_file(self, file: str | Path, **kwargs):
        """
        Parse the given file and populate the parser's data attribute.

        After parsing, stamp file provenance into every spectrum's metadata.

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
        for spectrum in self._data.values():
            spectrum.metadata["File"] = str(self.file)
            spectrum.metadata["file_ext"] = self.file.suffix

    @abstractmethod
    def _parse(self, file: Path, **kwargs) -> None:
        """
        Perform the actual parsing implementation.

        Subclasses must implement this method to extract structured
        data from the validated file and populate ``self._data``.

        Args:
            file: Path to the validated file.
            **kwargs: Additional parser-specific options.
        """
        ...


class _XPSParser(_Parser):
    """Abstract base class for all XPS file parsers.

    Subclasses define the structural and semantic rules required to identify
    and parse a specific XPS file format variant.

    Subclasses must set ``config_file`` and implement ``_parse()``, which
    should populate ``self._data`` — a ``dict[str, ParsedSpectrum]``
    mapping NeXus entry names to typed spectral data.

    The ``data`` property exposes this mapping to ``XPSReader.parsed_data_dicts``.
    File provenance (``File``, ``file_ext``) is stamped into every spectrum's
    metadata by ``parse_file()`` immediately after parsing completes.
    """

    config_file: ClassVar[str] = ""
    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset()

    def __init__(self) -> None:
        super().__init__()

    def _filter_metadata(self, raw: dict[str, Any]) -> dict[str, Any]:
        """Return *raw* with all ``_metadata_exclude_keys`` removed."""
        return {k: v for k, v in raw.items() if k not in self._metadata_exclude_keys}


class _XPSMetadataParser(_Parser):
    """
    Abstract base class for supplementary parsers that enrich
    existing spectral datasets.

    Metadata parsers parse auxiliary files (e.g., quantification exports)
    and inject additional information into spectra produced by a primary
    ``_XPSParser``.  Like ``_XPSParser``, ``_parse()`` populates
    ``self._data``, but with **metadata-only** entries:
    ``ParsedSpectrum(data=None, raw=None, metadata={...})``.

    ``update_main_file_data`` merges ``self._data`` into a caller-supplied
    ``dict[str, ParsedSpectrum]`` by matching entry names.  The default
    implementation handles the common case where the metadata parser's
    entry names are either identical to — or a suffix of — the primary
    parser's entry names (see ``_matches_entry``).

    ``compatible_primary_parser`` must be set by subclasses.
    """

    compatible_primary_parser: ClassVar[type[_XPSParser]]
    supported_primary_parser_versions: ClassVar[tuple[VersionRange, ...]] = ()

    @classmethod
    def file_ext_err_msg(cls, file: Path) -> str:
        suffix = file.suffix or "<no extension>"
        allowed = ", ".join(cls.supported_file_extensions) or "<none>"

        return (
            f"Cannot process file '{file.name}' (detected extension: '{suffix}'). "
            f"{cls.__name__} supports only the following file extensions: {allowed}."
        )

    @classmethod
    def supports_parser(cls, parser: type[_XPSParser] | _XPSParser) -> bool:
        """Return True if *parser* (class or instance) is compatible with this
        metadata parser."""
        try:
            cls._supports_parser(parser)
            return True
        except ValueError:
            return False

    @classmethod
    def _supports_parser(cls, parser: type[_XPSParser] | _XPSParser) -> None:
        """Raise ValueError if *parser* is not compatible with this metadata parser.

        Checks:

        1. *parser* is a subclass/instance of ``compatible_primary_parser``.
        2. If ``supported_primary_parser_versions`` is non-empty: the parser's
           ``supported_versions`` must overlap with it.

        Subclasses may override this classmethod to add further checks, calling
        ``super()._supports_parser(parser)`` first to retain the base checks.
        Failing further checks should also raise an ValueError.
        """
        parser_cls = parser if isinstance(parser, type) else type(parser)

        if not issubclass(parser_cls, cls.compatible_primary_parser):
            raise ValueError(
                f"{cls.__name__} is not compatible with {parser_cls.__name__}: "
                f"expected a subclass of {cls.compatible_primary_parser.__name__}."
            )

        if cls.supported_primary_parser_versions and parser_cls.supported_versions:
            if not _version_ranges_overlap(
                parser_cls.supported_versions,
                cls.supported_primary_parser_versions,
            ):
                raise ValueError(
                    f"{cls.__name__} does not support any version handled by "
                    f"{parser_cls.__name__}: no overlap between "
                    f"{parser_cls.supported_versions!r} and "
                    f"{cls.supported_primary_parser_versions!r}."
                )

    def update_main_file_data(self, main_file_data: dict[str, ParsedSpectrum]) -> None:
        """
        Merge ``self._data`` metadata into matching entries of *main_file_data*.

        For each entry in ``self._data``, finds the first key in
        *main_file_data* that satisfies ``_matches_entry`` and updates its
        ``metadata`` dict in-place.

        Args:
            main_file_data: Mapping from NeXus entry name to ``ParsedSpectrum``,
                as produced by the compatible primary parser.
        """
        for meta_entry, meta_spectrum in self._data.items():
            for main_entry, main_spectrum in main_file_data.items():
                if self._matches_entry(meta_entry, main_entry):
                    main_spectrum.metadata.update(meta_spectrum.metadata)
                    break

    @staticmethod
    def _matches_entry(meta_key: str, main_key: str) -> bool:
        """Return True if *meta_key* aligns with *main_key*.

        Handles two cases:

        - **Exact match**: ``meta_key == main_key``
          (e.g. both are ``"FeO__Fe_2p"``).
        - **Suffix match**: ``main_key`` ends with ``"__" + meta_key``
          (e.g. ``meta_key="Fe_2p"``, ``main_key="FeO__Fe_2p"``), used when
          the metadata file has no sample identifier.
        """
        return main_key == meta_key or main_key.endswith("__" + meta_key)


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


def _construct_data_key(spectrum: dict[str, Any]) -> str:
    """Construct a key for the ``data`` field of the xps_dict.

    Returns ``"cycle{N}_scan{M}"`` from ``loop_no`` and ``scan_no``
    fields in *spectrum*, defaulting to ``"cycle0_scan0"`` when absent.
    """
    loop_no = spectrum.get("loop_no")
    scan_no = spectrum.get("scan_no")
    cycle_key = f"cycle{loop_no}" if loop_no is not None else "cycle0"
    scan_key = f"scan{scan_no}" if scan_no is not None else "scan0"
    return f"{cycle_key}_{scan_key}"


def _construct_entry_name(parts: list[str]) -> str:
    """Construct name for the NXentry instances."""
    parts = [p for p in parts if p]
    if not parts:
        return ""
    if len(parts) == 1:
        return _align_name_part(parts[0])
    return "__".join([_align_name_part(part) for part in parts])
