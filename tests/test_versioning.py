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
Parametrized tests for:
  - parsers/versioning.py
"""

import pytest

from pynxtools_xps.parsers.versioning import (
    _format_version,
    is_version_supported,
    normalize_version,
)


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("4.75", (4, 75)),
        ("4", (4,)),
        ("2019R2", (2019, "R", 2)),
        ("v3_1", ("V", 3, 1)),
        ("3-2", (3, 2)),
        ("4.100", (4, 100)),
        ("1.2.3", (1, 2, 3)),
        ("  4.75  ", (4, 75)),  # leading/trailing whitespace stripped
    ],
)
def test_normalize_version(raw, expected):
    assert normalize_version(raw) == expected


@pytest.mark.parametrize(
    "raw, exc",
    [
        (None, TypeError),  # None is rejected explicitly
        ("", ValueError),  # empty string
        ("---", ValueError),  # no alphanumeric tokens
    ],
)
def test_normalize_version_invalid(raw, exc):
    with pytest.raises(exc):
        normalize_version(raw)


@pytest.mark.parametrize(
    "version, expected",
    [
        ((4, 75), "4.75"),
        ((2019, "R", 2), "2019.R.2"),
        ((1,), "1"),
        ((1, 2, 3), "1.2.3"),
    ],
)
def test_format_version(version, expected):
    assert _format_version(version) == expected


@pytest.mark.parametrize(
    "version, supported, requires, expected",
    [
        # No version present
        (None, [], False, True),  # no constraints + no required version → True
        (None, [], True, False),  # requires_version=True but no version → False
        # No declared ranges → all versions pass
        ((4, 1), [], False, True),
        # Half-open range: [lower, upper)
        ((4, 1), [((4, 0), (5, 0))], False, True),  # inside range
        ((4, 0), [((4, 0), (5, 0))], False, True),  # on lower bound (inclusive)
        ((5, 0), [((4, 0), (5, 0))], False, False),  # on upper bound (exclusive)
        ((3, 9), [((4, 0), (5, 0))], False, False),  # below range
        # Unbounded upper: [lower, ∞)
        ((6, 0), [((4, 0), None)], False, True),
        ((3, 0), [((4, 0), None)], False, False),
        # Multiple ranges — first match wins
        ((4, 0), [((4, 0), (5, 0)), ((6, 0), None)], False, True),
        # Gap between ranges — falls through both
        ((5, 5), [((4, 0), (5, 0)), ((6, 0), None)], False, False),
    ],
)
def test_is_version_supported(version, supported, requires, expected):
    assert (
        is_version_supported(version, supported, requires_version=requires) == expected
    )
