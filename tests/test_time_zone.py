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
"""SPECS XY timestamps never get the time zone of the converting machine"""

# cspell:ignore tzset
import logging
import time
from pathlib import Path

import pytest

from pynxtools_xps.parsers.specs.xy.metadata import _parse_datetime
from pynxtools_xps.parsers.vms.parser import VamasParser

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def machine_in_tokyo(monkeypatch):
    """Pretend the converting machine runs in Tokyo (UTC+09:00)."""
    monkeypatch.setenv("TZ", "Asia/Tokyo")
    time.tzset()
    yield
    monkeypatch.undo()
    time.tzset()


def test_utc_timestamp_stays_utc(machine_in_tokyo):
    assert _parse_datetime("08/24/23 14:19:47 UTC") == "2023-08-24T14:19:47+00:00"


def test_timestamp_with_offset_is_kept(machine_in_tokyo):
    assert _parse_datetime("2023-08-24 14:19:47+02:00") == "2023-08-24T14:19:47+02:00"


def test_timestamp_without_zone_stays_without_zone(machine_in_tokyo, caplog):
    """No zone in the file: keep the time as it is and say so."""
    with caplog.at_level(logging.WARNING):
        assert _parse_datetime("08/24/23 14:19:47") == "2023-08-24T14:19:47"
    assert "no time zone" in caplog.text


def test_iso_timestamp_without_zone_stays_without_zone(machine_in_tokyo, caplog):
    with caplog.at_level(logging.WARNING):
        assert _parse_datetime("2023-08-24 14:19:47") == "2023-08-24T14:19:47"
    assert "no time zone" in caplog.text


def _parse_regular_vms(tmp_path, month: str):
    """Parse a copy of regular.vms with the month of the first block replaced."""
    lines = (DATA_DIR / "vms_regular" / "regular.vms").read_text().splitlines()
    year = lines.index("2023")
    assert lines[year + 1] == "8"  # month of the first block
    lines[year + 1] = month
    vms_file = tmp_path / "dated.vms"
    vms_file.write_text("\n".join(lines) + "\n")

    parser = VamasParser()
    parser.parse(vms_file)
    return next(iter(parser._data.values())).metadata


def test_valid_vamas_date_keeps_its_zone(tmp_path):
    metadata = _parse_regular_vms(tmp_path, month="8")
    assert metadata["time_stamp"] == "2023-08-24T14:19:47+00:00"


def test_invalid_vamas_date_is_left_out(tmp_path, caplog):
    """Month 13 is not a date: write no time stamp instead of 0001-01-01."""
    with caplog.at_level(logging.WARNING):
        metadata = _parse_regular_vms(tmp_path, month="13")
    assert metadata.get("time_stamp") is None
    assert "not a valid date" in caplog.text


def test_vamas_date_of_zeros_is_left_out_quietly(tmp_path, caplog):
    """0-0-0 means no date was recorded: no time stamp and no warning."""
    lines = (DATA_DIR / "vms_regular" / "regular.vms").read_text().splitlines()
    year = lines.index("2023")
    lines[year : year + 3] = ["0", "0", "0"]
    vms_file = tmp_path / "undated.vms"
    vms_file.write_text("\n".join(lines) + "\n")

    parser = VamasParser()
    with caplog.at_level(logging.WARNING):
        parser.parse(vms_file)
    metadata = next(iter(parser._data.values())).metadata
    assert metadata.get("time_stamp") is None
    assert caplog.text == ""
