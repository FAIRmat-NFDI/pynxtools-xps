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
Tests for helper functions
"""

import datetime
import re

import pytest

from pynxtools_xps.reader_utils import extract_unit
from pynxtools_xps.value_mappers import get_units_for_key, parse_datetime


@pytest.mark.parametrize(
    "unit_key, unit_map, expected_unit",
    [
        # Test cases with direct keys in the unit_map
        ("detector/detector_voltage", {"detector/detector_voltage": "V"}, "V"),
        ("temperature", {"temperature": "K"}, "K"),
        # Test cases with keys that don't exist in the unit_map
        ("nonexistent_key", {"some_key": "m"}, None),
        # Test cases with regex pattern in the key
        (
            "detector/[detector_voltage]",
            {"detector/detector_voltage": "V"},
            "detector_voltage",
        ),
        ("sensor/[sensor_current]", {"sensor/sensor_current": "A"}, "sensor_current"),
        # Test cases with key that includes a regex but isn't mapped
        ("not_mapped/[value]", {}, "value"),
        # Key with multiple slashes
        ("foo/bar/baz", {"foo/bar/baz": "kg"}, "kg"),
        # Key with trailing slash
        ("foo/bar/", {"foo/bar/": "s"}, "s"),
        # Key with leading slash
        ("/foo/bar", {"/foo/bar": "A"}, "A"),
    ],
)
def test_get_units_for_key(unit_key: str, unit_map: dict[str, str], expected_unit: str):
    result = get_units_for_key(unit_key, unit_map)
    assert result == expected_unit, f"Expected {expected_unit} but got {result}"


@pytest.mark.parametrize(
    "key, value, unit_missing, expected_value, expected_unit",
    [
        # Test cases with explicit units
        ("x_position", "0 mm", {}, 0, "mm"),
        ("polar_rotation", "0 deg", {}, 0, "deg"),
        ("analyser_work_function", "4.506eV", {}, 4.506, "eV"),
        ("temperature", "300K", {}, 300, "K"),
        ("pressure", "1.01bar", {}, 1.01, "bar"),
        # Test cases without explicit units
        ("voltage", "5.0", {"voltage": "V"}, 5.0, "V"),
        ("current", "10", {"current": "A"}, 10, "A"),
        # Test cases with scientific notation
        ("distance", "1.23e-10m", {}, 1.23e-10, "m"),
        ("charge", "1.602e-19C", {}, 1.602e-19, "C"),
        # Test cases with missing unit in value and unit_missing dictionary
        ("energy", "1000", {"energy": "J"}, 1000, "J"),
        ("mass", "0.5", {"mass": "kg"}, 0.5, "kg"),
    ],
)
def test_extract_unit(
    key: str,
    value: str,
    unit_missing: dict[str, str],
    expected_value: str,
    expected_unit: str,
):
    result_value, result_unit = extract_unit(key, value, unit_missing)

    assert isinstance(result_value, type(expected_value)), (
        f"Expected type {type(expected_value)} but got type {type(result_value)}"
    )

    assert result_value == expected_value, (
        f"Expected {expected_value} but got {result_value}"
    )
    assert result_unit == expected_unit, (
        f"Expected {expected_unit} but got {result_unit}"
    )


@pytest.mark.parametrize(
    "datetime_string, possible_date_formats, tzinfo, expected_iso8601",
    [
        # Test cases with valid datetime strings and formats
        (
            "2023-10-23 14:30:00",
            ["%Y-%m-%d %H:%M:%S"],
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
        # Test cases with timezone information
        (
            "2023-10-23 14:30:00",
            ["%Y-%m-%d %H:%M:%S"],
            datetime.timezone(datetime.timedelta(hours=1)),
            "2023-10-23T14:30:00+01:00",
        ),
        # Test case with missing timezone should still return UTC
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
        # Test cases that should raise ValueError
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
