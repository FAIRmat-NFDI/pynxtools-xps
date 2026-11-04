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
Tests for SPECSMetadataSLHParser, against a synthetic SQLite database built
to the real schema confirmed live against six real SPECS Lab Prodigy .slh
files (2026-08-06), rather than a committed real instrument-log fixture --
keeps this self-contained and lets each edge case (mis-encoded unit bytes,
out-of-order stored offsets, missing ParameterInfo row) be constructed
precisely instead of hoping a real file happens to exercise it.
"""

import sqlite3
from pathlib import Path

import pytest

from pynxtools_xps.parsers.base import ParsedSpectrum
from pynxtools_xps.parsers.specs.sle.parser import SPECSSLEParser
from pynxtools_xps.parsers.specs.slh.parser import SPECSMetadataSLHParser

DATA_DIR = Path(__file__).parent / "data"


def _make_slh(
    path: Path,
    *,
    version: str = "0.6",
    parameters: list[tuple[str, str, str, str, float]] | None = None,
    base_time: str = "2025-Aug-13 11:44:07.653833",
    observations: dict[int, list[tuple[float, float]]] | None = None,
) -> None:
    """Build a synthetic .slh file with the real confirmed schema.

    `parameters` is a list of (name, device, parameter, unit, scaling)
    tuples, one per logged parameter (row order = ParameterHistory.ID,
    starting at 1, same as every real file inspected).
    `observations` maps parameter id -> list of (offset_s, observation)
    rows, deliberately allowed to be given out of chronological order --
    real files store them that way.
    """
    if parameters is None:
        parameters = [("emission", "XR 50 MF", "emission", "mA", 1000.0)]
    if observations is None:
        observations = {1: [(0.0, 0.01), (10.0, 0.02), (5.0, 0.015)]}

    con = sqlite3.connect(str(path))
    cur = con.cursor()
    cur.execute("CREATE TABLE Configuration (Key TEXT, Value TEXT)")
    cur.execute(
        "INSERT INTO Configuration VALUES ('Version', ?)",
        (version,),
    )
    cur.execute(
        "CREATE TABLE ParameterHistory (ID INTEGER, Name TEXT, UniqueDev TEXT, "
        "Device TEXT, Parameter TEXT, BaseTime TEXT)"
    )
    cur.execute(
        "CREATE TABLE ParameterInfo (ID INTEGER, DeviceType TEXT, Command TEXT, "
        "ReadablePar TEXT, Prototype BLOB, Unit TEXT, Scaling FLOAT, "
        "Representation TEXT)"
    )
    cur.execute(
        "CREATE TABLE NumericalHistoryData (ID INTEGER, Offset_s FLOAT, "
        "Observation FLOAT)"
    )
    cur.execute(
        "CREATE TABLE TextualHistoryData (ID INTEGER, Offset_s FLOAT, Observation TEXT)"
    )

    for i, (name, device, parameter, unit, scaling) in enumerate(parameters, start=1):
        cur.execute(
            "INSERT INTO ParameterHistory VALUES (?, ?, ?, ?, ?, ?)",
            (i, f"{parameter} ({device})", device, device, parameter, base_time),
        )
        unit_value: str | bytes = unit
        if unit == "\xb0C":  # mis-encoded-degree-symbol case, real gotcha
            unit_value = b"\xb0C"
        cur.execute(
            "INSERT INTO ParameterInfo VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (i, device, None, parameter, b"<Par/>", unit_value, scaling, "2:0"),
        )
        for offset_s, observation in observations.get(i, []):
            cur.execute(
                "INSERT INTO NumericalHistoryData VALUES (?, ?, ?)",
                (i, offset_s, observation),
            )

    con.commit()
    con.close()


@pytest.fixture
def slh_file(tmp_path) -> Path:
    path = tmp_path / "test.slh"
    _make_slh(path)
    return path


# ── matches_file / detect_version ──────────────────────────────────────────


def test_matches_file_accepts_a_real_shaped_slh(slh_file):
    assert SPECSMetadataSLHParser().matches_file(slh_file)


def test_matches_file_rejects_an_sle_file():
    sle_path = DATA_DIR / "specs_sle" / "EX439_S718_Au.sle"
    if not sle_path.exists():
        pytest.skip(f"test file not found: {sle_path}")
    assert not SPECSMetadataSLHParser().matches_file(sle_path)


def test_sle_parser_rejects_an_slh_file(slh_file):
    """The reader refuses a file that matches both a primary and a metadata
    parser (see reader.py:handle_data_file) -- so this has to hold, not
    just the reverse direction."""
    assert not SPECSSLEParser().matches_file(slh_file)


def test_matches_file_rejects_a_non_sqlite_file(tmp_path):
    path = tmp_path / "not_sqlite.slh"
    path.write_text("not a database")
    assert not SPECSMetadataSLHParser().matches_file(path)


def test_detect_version_reads_the_configuration_table(slh_file):
    assert SPECSMetadataSLHParser().detect_version(slh_file) == (0, 6)


# ── _parse: scaling, encoding, ordering ─────────────────────────────────────


def test_parse_applies_scaling_to_raw_observations(tmp_path):
    """Confirmed live (2026-08-06) against a real file: Observation is
    stored unscaled -- e.g. raw ~0.087-0.09 for a Scaling=1000.0, "mA"
    parameter only makes sense as Amps needing the multiply."""
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("emission", "XR 50 MF", "emission", "mA", 1000.0)],
        observations={1: [(0.0, 0.0866)]},
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path)
    row = parser._parameter_logs.iloc[0]
    assert row["values"] == pytest.approx([86.6])


def test_parse_handles_mis_encoded_unit_bytes_without_crashing(tmp_path):
    """Real gotcha: ParameterInfo.Unit can contain non-UTF-8 bytes (a
    mis-encoded degree symbol in "°C")."""
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[
            ("water_temperature", "XR 50 MF", "water_temperature", "\xb0C", 1.0)
        ],
        observations={1: [(0.0, 21.6)]},
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path)  # must not raise
    row = parser._parameter_logs.iloc[0]
    assert row["values"] == pytest.approx([21.6])


def test_parse_sorts_observations_despite_out_of_order_storage(tmp_path):
    """Real gotcha: NumericalHistoryData rows are not guaranteed to be
    stored in chronological Offset_s order."""
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("pressure", "Atmion", "pressure", "mbar", 1.0)],
        observations={1: [(10.0, 3.0), (0.0, 1.0), (5.0, 2.0)]},
        base_time="2025-Aug-13 00:00:00.000000",
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path, slh_resample_interval=None)
    row = parser._parameter_logs.iloc[0]
    assert row["values"] == [1.0, 2.0, 3.0]
    assert row["timestamps"] == sorted(row["timestamps"])


def test_parse_resample_interval_none_keeps_raw_irregular_timestamps(tmp_path):
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("pressure", "Atmion", "pressure", "mbar", 1.0)],
        observations={1: [(0.0, 1.0), (0.7, 2.0), (30.0, 3.0)]},
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path, slh_resample_interval=None)
    row = parser._parameter_logs.iloc[0]
    assert len(row["timestamps"]) == 3  # noqa: PLR2004 -- one per raw observation


def test_parse_default_resamples_to_one_second(tmp_path):
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("pressure", "Atmion", "pressure", "mbar", 1.0)],
        observations={1: [(0.0, 1.0), (0.7, 2.0), (5.0, 3.0)]},
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path)  # default slh_resample_interval="1s"
    row = parser._parameter_logs.iloc[0]
    assert len(row["timestamps"]) == 6  # noqa: PLR2004 -- 0s..5s inclusive


def test_parse_skips_a_parameter_with_no_observations(tmp_path):
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("unused_param", "Dev", "unused_param", "V", 1.0)],
        observations={},
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path)
    assert len(parser._parameter_logs) == 0


# ── update_main_file_data: per-entry slicing ────────────────────────────────


def test_update_main_file_data_slices_to_the_entry_time_window(tmp_path):
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("pressure", "Atmion", "pressure", "mbar", 1.0)],
        observations={1: [(t, float(t)) for t in range(0, 20)]},
        base_time="2025-Aug-13 00:00:00.000000",
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path, slh_resample_interval=None)

    spectrum = ParsedSpectrum(
        metadata={
            "time_stamp": "2025-Aug-13T00:00:05.000000Z",
            "dwell_time": 1.0,
            "n_values": 5,
            "total_scans": 1,
        }
    )
    parser.update_main_file_data({"entry": spectrum})

    key = "pressure_Atmion"
    assert spectrum.metadata[f"{key}/@units"] == "mbar"
    # window is [5s, 5s + 1*5*1 = 10s]
    assert spectrum.metadata[key] == [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]


def test_update_main_file_data_skips_an_entry_with_no_time_stamp(tmp_path):
    path = tmp_path / "test.slh"
    _make_slh(path)
    parser = SPECSMetadataSLHParser()
    parser.parse(path)

    spectrum = ParsedSpectrum(metadata={})
    parser.update_main_file_data({"entry": spectrum})

    assert spectrum.metadata == {}


def test_update_main_file_data_leaves_metadata_untouched_when_window_matches_nothing(
    tmp_path,
):
    path = tmp_path / "test.slh"
    _make_slh(
        path,
        parameters=[("pressure", "Atmion", "pressure", "mbar", 1.0)],
        observations={1: [(0.0, 1.0)]},
        base_time="2025-Aug-13 00:00:00.000000",
    )
    parser = SPECSMetadataSLHParser()
    parser.parse(path, slh_resample_interval=None)

    spectrum = ParsedSpectrum(
        metadata={
            "time_stamp": "2030-Jan-01T00:00:00.000000Z",
            "dwell_time": 1.0,
            "n_values": 1,
            "total_scans": 1,
        }
    )
    parser.update_main_file_data({"entry": spectrum})

    assert "pressure_Atmion" not in spectrum.metadata


def test_update_main_file_data_noop_when_no_parameter_logs():
    parser = SPECSMetadataSLHParser()
    spectrum = ParsedSpectrum(metadata={"time_stamp": "2025-Aug-13T00:00:00.000000Z"})
    parser.update_main_file_data({"entry": spectrum})
    assert spectrum.metadata == {"time_stamp": "2025-Aug-13T00:00:00.000000Z"}
