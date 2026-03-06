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
Utilities to check for versioning
"""

import re
from collections.abc import Iterable
from typing import TypeAlias

VersionTuple: TypeAlias = tuple[int | str, ...]
VersionRange: TypeAlias = tuple[VersionTuple, VersionTuple | None]


def _format_version(version: VersionTuple) -> str:
    return ".".join(str(x) for x in version)


_TOKEN_RE = re.compile(r"\d+|[A-Za-z]+")


def normalize_version(raw: str) -> VersionTuple:
    """
    Normalize a vendor version string into a comparable tuple.

    Tokenization rules
    ------------------
    - Consecutive digits are parsed as integers.
    - Consecutive letters are upper cased and kept as strings.
    - All other characters act as separators.
    - Comparison is purely tuple-based and lexicographic.

    Examples
    --------
    "4.75"    -> (4, 75)
    "4"       -> (4,)
    "2019R2"  -> (2019, "R", 2)
    "v3_1"    -> ("V", 3, 1)
    "3-2"     -> (3, 2)
    "4.100"   -> (4, 100)
    """
    if raw is None:
        raise TypeError("Version string must not be None.")

    raw = raw.strip()
    if not raw:
        raise ValueError("Version string must not be empty.")

    tokens = _TOKEN_RE.findall(raw)
    if not tokens:
        raise ValueError(f"Could not extract version tokens from '{raw}'.")

    normalized: list[int | str] = []
    for token in tokens:
        if token.isdigit():
            normalized.append(int(token))
        else:
            normalized.append(token.upper())

    return tuple(normalized)


def is_version_supported(
    version: VersionTuple | None,
    supported_versions: Iterable[VersionRange],
) -> bool:
    """
    Determine whether a normalized version tuple is supported.

    Parameters
    ----------
    version
        Normalized version tuple, or ``None`` if the file carries no version.
    supported_versions
        Iterable of half-open version intervals:
        ``(lower_inclusive, upper_exclusive_or_None)``.
        An empty iterable means no version constraint: all files are accepted,
        including those without a version.
        A non-empty iterable implicitly requires a version — ``None`` is
        rejected because it cannot fall within any declared range.

    Returns
    -------
    bool
    """
    # No constraints declared → accept everything, including files without version
    ranges = tuple(supported_versions)
    if not ranges:
        return True

    # Ranges declared but file carries no version → reject
    if version is None:
        return False

    for lower, upper in ranges:
        if upper is None:
            if version >= lower:
                return True
        elif lower <= version < upper:
            return True

    return False


def _version_ranges_overlap(
    primary_parser_versions: Iterable[VersionRange],
    required_parser_versions: Iterable[VersionRange],
) -> bool:
    """Return True if any range in *primary_parser_versions* overlaps with any
    range in *required_parser_versions*.

    Ranges are half-open intervals ``[lower, upper)`` with ``upper=None`` meaning
    unbounded. Two intervals ``[p_lo, p_hi)`` and ``[r_lo, r_hi)`` overlap iff
    ``p_lo < r_hi AND r_lo < p_hi`` (treating ``None`` as +∞).

    Used to check whether a metadata parser's ``supported_primary_parser_versions``
    is compatible with a primary parser's ``supported_versions``.
    """
    for p_lo, p_hi in primary_parser_versions:
        for r_lo, r_hi in required_parser_versions:
            p_lo_lt_r_hi = r_hi is None or p_lo < r_hi
            r_lo_lt_p_hi = p_hi is None or r_lo < p_hi
            if p_lo_lt_r_hi and r_lo_lt_p_hi:
                return True
    return False
