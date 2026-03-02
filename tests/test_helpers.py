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
  - parsers/base._XPSDataclass and name-construction helpers
  - mapping.parse_datetime
"""

import datetime
import re
from dataclasses import dataclass

import pytest

from pynxtools_xps.mapping import convert_pascal_to_snake, parse_datetime
from pynxtools_xps.parsers.base import (
    _align_name_part,
    _construct_data_key,
    _construct_entry_name,
    _XPSDataclass,
)

# ── pascal_to_snake ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("ScanMode", "scan_mode"),
        ("EnergyType", "energy_type"),
        ("XPS", "xps"),  # all-uppercase → all-lowercase
        ("PassEnergy", "pass_energy"),
        ("already_snake", "already_snake"),
        ("key[unit]", "key[unit]"),  # bracketed part preserved verbatim
        ("ABCDef", "abc_def"),  # internal uppercase run split correctly
    ],
)
def test_convert_pascal_to_snake(raw, expected):
    assert convert_pascal_to_snake(raw) == expected


# ── Datetime parsing ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "datetime_string, possible_date_formats, tzinfo, expected_iso8601",
    [
        # Non-ISO format: tzinfo IS applied via strptime path
        (
            "23/10/2023 14:30:00",
            ["%d/%m/%Y %H:%M:%S"],
            datetime.timezone.utc,
            "2023-10-23T14:30:00+00:00",
        ),
        (
            "23/10/2023 14:30",
            ["%d/%m/%Y %H:%M"],
            datetime.timezone.utc,
            "2023-10-23T14:30:00+00:00",
        ),
        (
            "October 23, 2023 14:30",
            ["%B %d, %Y %H:%M"],
            datetime.timezone.utc,
            "2023-10-23T14:30:00+00:00",
        ),
        # Already-ISO strings: fromisoformat short-circuits; tzinfo is NOT applied
        (
            "2023-10-23T14:30:00Z",
            ["%Y-%m-%dT%H:%M:%SZ"],
            datetime.timezone.utc,
            "2023-10-23T14:30:00+00:00",
        ),
        (
            "2023-10-23T14:30:00+0200",
            ["%Y-%m-%dT%H:%M:%S%z"],
            datetime.timezone(datetime.timedelta(hours=2)),
            "2023-10-23T14:30:00+02:00",
        ),
        # Non-ISO format with non-UTC timezone offset
        (
            "23.10.2023 14:30:00",
            ["%d.%m.%Y %H:%M:%S"],
            datetime.timezone(datetime.timedelta(hours=1)),
            "2023-10-23T14:30:00+01:00",
        ),
        # ISO-parseable string with tzinfo=None → no offset in result
        (
            "2023-10-23 14:30:00",
            ["%Y-%m-%d %H:%M:%S"],
            None,
            "2023-10-23T14:30:00",
        ),
    ],
)
def test_parse_datetime(
    datetime_string, possible_date_formats, tzinfo, expected_iso8601
):
    result = parse_datetime(datetime_string, possible_date_formats, tzinfo)
    assert result == expected_iso8601, f"Expected {expected_iso8601} but got {result}"


@pytest.mark.parametrize(
    "datetime_string, possible_date_formats",
    [
        ("invalid date string", ["%Y-%m-%d %H:%M:%S"]),
        ("2023-10-23 invalid", ["%Y-%m-%d %H:%M:%S"]),
        ("23-10-2023", ["%Y-%m-%d %H:%M:%S"]),
    ],
)
def test_parse_datetime_invalid(datetime_string, possible_date_formats):
    expected_error = (
        f"Datetime {datetime_string} could not be converted to ISO 8601 format, "
        f"as it does not match any of these formats: {possible_date_formats}."
    )
    with pytest.raises(ValueError, match=re.escape(expected_error)):
        parse_datetime(datetime_string, possible_date_formats)


# ── base name-construction helpers ────────────────────────────────────────────


@pytest.mark.parametrize(
    "name_part, expected",
    [
        ("C 1s", "C_1s"),
        ("Fe-2p", "Fe_2p"),
        ("Au/Ag", "Au_Ag"),
        ("a.b", "a_b"),
        ("a:b", "a_b"),
        ("a+b", "a_b"),
        ("a=b", "ab"),
        ("clean", "clean"),
    ],
)
def test_align_name_part(name_part, expected):
    assert _align_name_part(name_part) == expected


@pytest.mark.parametrize(
    "parts, expected",
    [
        (["C 1s"], "C_1s"),
        (["C 1s", "scan1"], "C_1s__scan1"),
        (["", "C 1s"], "C_1s"),  # empty parts are filtered out
        ([], ""),
    ],
)
def test_construct_entry_name(parts, expected):
    assert _construct_entry_name(parts) == expected


@pytest.mark.parametrize(
    "spectrum, expected",
    [
        ({"loop_no": 0, "scan_no": 3}, "cycle0_scan3"),
        ({"loop_no": None, "scan_no": None}, "cycle0_scan0"),
        ({"loop_no": 2, "scan_no": None}, "cycle2_scan0"),
        ({"loop_no": 1, "scan_no": 0}, "cycle1_scan0"),
    ],
)
def test_construct_data_key(spectrum, expected):
    assert _construct_data_key(spectrum) == expected


# ── _XPSDataclass ─────────────────────────────────────────────────────────────


@dataclass
class _SampleDC(_XPSDataclass):
    """Minimal concrete dataclass for testing _XPSDataclass type enforcement."""

    x: int
    y: float | None
    name: str
    tags: list[str]


@pytest.mark.parametrize(
    "x, y, name, tags",
    [
        (1, None, "a", []),
        (2, 3.14, "b", ["t1"]),
        (3, None, "c", ["x", "y"]),
        (0, 0.0, "", []),
    ],
)
def test_xps_dataclass_valid(x, y, name, tags):
    dc = _SampleDC(x=x, y=y, name=name, tags=tags)
    assert dc.x == x
    assert dc.y == y
    assert dc.name == name
    assert dc.tags == tags


def test_xps_dataclass_coercion():
    """A coercible string value is converted to the annotated type."""
    dc = _SampleDC(x="2", y=None, name="a", tags=[])
    assert dc.x == 2
    assert isinstance(dc.x, int)


@pytest.mark.parametrize(
    "x, y, name, tags",
    [
        ("abc", None, "a", []),  # "abc" cannot be coerced to int
        (1, "not_a_float", "a", []),  # "not_a_float" cannot be coerced to float|None
    ],
)
def test_xps_dataclass_invalid(x, y, name, tags):
    with pytest.raises(TypeError):
        _SampleDC(x=x, y=y, name=name, tags=tags)
