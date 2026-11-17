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
"""Time zones of timestamps that come without one"""

# cspell:ignore tzset
import logging
import shutil
import time
from pathlib import Path

import pytest
from pynxtools import logger as pynx_logger
from pynxtools.dataconverter.convert import get_reader
from pynxtools.dataconverter.template import Template

from pynxtools_xps.mapping import attach_time_zone
from pynxtools_xps.parsers.specs.xy.metadata import _parse_datetime
from pynxtools_xps.parsers.vms.parser import VamasParser

DATA_DIR = Path(__file__).parent / "data"
XY_DIR = DATA_DIR / "specs_xy"


@pytest.fixture
def pynx_caplog(caplog):
    """caplog that also sees the pynxtools logger, which does not propagate."""
    pynx_logger.addHandler(caplog.handler)
    try:
        yield caplog
    finally:
        pynx_logger.removeHandler(caplog.handler)


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


def test_xy_parser_leaves_a_timestamp_without_zone_to_the_reader(machine_in_tokyo):
    """The parser does not guess; the reader adds the zone later."""
    assert _parse_datetime("08/24/23 14:19:47") == "2023-08-24T14:19:47"
    assert _parse_datetime("2023-08-24 14:19:47") == "2023-08-24T14:19:47"


def test_timestamp_with_zone_is_not_changed(machine_in_tokyo, pynx_caplog):
    with pynx_caplog.at_level(logging.WARNING):
        assert (
            attach_time_zone("2023-08-24T14:19:47+02:00", "America/New_York")
            == "2023-08-24T14:19:47+02:00"
        )
    assert pynx_caplog.records == []


def test_given_zone_is_attached(machine_in_tokyo, pynx_caplog):
    with pynx_caplog.at_level(logging.WARNING):
        assert (
            attach_time_zone("2023-08-24T14:19:47", "Europe/Berlin")
            == "2023-08-24T14:19:47+02:00"
        )
    assert pynx_caplog.records == []


def test_without_zone_the_machine_zone_is_assumed_with_an_info(
    machine_in_tokyo, pynx_caplog
):
    """No zone in the file and none given: use the machine's zone, but say so."""
    with pynx_caplog.at_level(logging.INFO):
        assert attach_time_zone("2023-08-24T14:19:47") == "2023-08-24T14:19:47+09:00"
    assert "time zone of the converting machine" in pynx_caplog.text


def test_unknown_zone_falls_back_to_the_machine_zone(machine_in_tokyo, pynx_caplog):
    with pynx_caplog.at_level(logging.INFO):
        assert (
            attach_time_zone("2023-08-24T14:19:47", "Mars/Olympus")
            == "2023-08-24T14:19:47+09:00"
        )
    assert "Mars/Olympus" in pynx_caplog.text
    assert "time zone of the converting machine" in pynx_caplog.text


def test_timezone_from_the_eln_is_used(machine_in_tokyo, tmp_path):
    """`timezone` in the ELN sets the zone and is not written as data."""
    for file in XY_DIR.iterdir():
        if file.suffix != ".nxs":
            shutil.copy(file, tmp_path / file.name)
    xy_file = tmp_path / "MgFe2O4_small.xy"
    xy_file.write_text(xy_file.read_text().replace(" UTC", ""))  # no zone
    eln = tmp_path / "eln_data_specs_xy.yaml"
    eln.write_text(eln.read_text() + "\ntimezone: America/New_York\n")

    template = get_reader("xps")().read(
        template=Template(), file_paths=tuple(str(f) for f in tmp_path.iterdir())
    )

    assert (
        template["/ENTRY[1_as_loaded__Survey]/start_time"]
        == "2023-08-24T14:19:47-04:00"
    )
    assert not any("timezone" in key for key in template.keys())


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


def test_invalid_vamas_date_is_left_out(tmp_path, pynx_caplog):
    """Month 13 is not a date: write no time stamp instead of 0001-01-01."""
    with pynx_caplog.at_level(logging.WARNING):
        metadata = _parse_regular_vms(tmp_path, month="13")
    assert metadata.get("time_stamp") is None
    assert "not a valid date" in pynx_caplog.text


def test_vamas_date_of_zeros_is_left_out_quietly(tmp_path, pynx_caplog):
    """0-0-0 means no date was recorded: no time stamp and no warning."""
    lines = (DATA_DIR / "vms_regular" / "regular.vms").read_text().splitlines()
    year = lines.index("2023")
    lines[year : year + 3] = ["0", "0", "0"]
    vms_file = tmp_path / "undated.vms"
    vms_file.write_text("\n".join(lines) + "\n")

    parser = VamasParser()
    with pynx_caplog.at_level(logging.WARNING):
        parser.parse(vms_file)
    metadata = next(iter(parser._data.values())).metadata
    assert metadata.get("time_stamp") is None
    assert pynx_caplog.text == ""
