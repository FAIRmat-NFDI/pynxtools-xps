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
Tests for SPECSMetadataCSVParser, against synthetic CSV fixtures built to
the real-shaped header/row convention confirmed directly by the repo owner
(semicolon-delimited despite the '.csv' extension, bracket/paren-encoded
unit+device, independently-optional time columns) -- there is no real
committed CSV export in this repo to test against, unlike .slh, so every
fixture here is hand-authored to that confirmed shape.
"""

from pathlib import Path

import pytest

from pynxtools_xps.parsers.base import ParsedSpectrum
from pynxtools_xps.parsers.specs.sle.parser import SPECSSLEParser
from pynxtools_xps.parsers.specs.slh.parser import SPECSMetadataCSVParser

DATA_DIR = Path(__file__).parent / "data"


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Write a semicolon-delimited CSV, real-export-shaped (quoted header,
    unquoted data rows, blank cells for "unchanged since last row")."""
    lines = [";".join(f'"{col}"' for col in header)]
    lines.extend(";".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ── matches_file ─────────────────────────────────────────────────────────


def test_matches_file_accepts_universal_time_header(tmp_path):
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Relative Time [s]", "Pressure [mbar] (Ceravac_AC)"],
        [["2025-Apr-03 06:57:15.625973", "0.000135", "1.2e-06"]],
    )
    assert SPECSMetadataCSVParser().matches_file(path)


def test_matches_file_accepts_local_date_time_header_only(tmp_path):
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Local Date", "Local Time of Day", "Pressure [mbar] (Ceravac_AC)"],
        [["2025-04-03", "08:57:15.625973", "1.2e-06"]],
    )
    assert SPECSMetadataCSVParser().matches_file(path)


def test_matches_file_rejects_comma_delimited_csv(tmp_path):
    path = tmp_path / "export.csv"
    path.write_text(
        "Universal Time,Relative Time [s],Pressure [mbar] (Ceravac_AC)\n"
        "2025-Apr-03 06:57:15.625973,0.000135,1.2e-06\n",
        encoding="utf-8",
    )
    assert not SPECSMetadataCSVParser().matches_file(path)


def test_matches_file_rejects_semicolon_csv_without_any_time_column(tmp_path):
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Relative Time [s]", "Pressure [mbar] (Ceravac_AC)"],
        [["0.000135", "1.2e-06"]],
    )
    assert not SPECSMetadataCSVParser().matches_file(path)


def test_matches_file_rejects_an_slh_file():
    slh_path = DATA_DIR / "specs_slh" / "xray.slh"
    if not slh_path.exists():
        pytest.skip(f"test file not found: {slh_path}")
    assert not SPECSMetadataCSVParser().matches_file(slh_path)


def test_matches_file_rejects_an_sle_file():
    sle_path = DATA_DIR / "specs_sle" / "EX439_S718_Au.sle"
    if not sle_path.exists():
        pytest.skip(f"test file not found: {sle_path}")
    assert not SPECSMetadataCSVParser().matches_file(sle_path)


def test_sle_parser_rejects_a_csv_file(tmp_path):
    """The reader refuses a file that matches both a primary and a metadata
    parser (see reader.py:handle_data_file) -- so this has to hold too."""
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Pressure [mbar] (Ceravac_AC)"],
        [["2025-Apr-03 06:57:15.625973", "1.2e-06"]],
    )
    assert not SPECSSLEParser().matches_file(path)


# ── version gate ─────────────────────────────────────────────────────────


def test_no_version_constraint_declared():
    """CSV exports carry no version marker -- no constraint declared."""
    assert SPECSMetadataCSVParser.supported_versions == ()
    assert SPECSMetadataCSVParser().detect_version(Path("does-not-matter.csv")) is None


def test_match_failure_reason_is_none_for_a_matching_file(tmp_path):
    """A structurally-matching file must pass is_mainfile()/parse()."""
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Pressure [mbar] (Ceravac_AC)"],
        [["2025-Apr-03 06:57:15.625973", "1.2e-06"]],
    )
    assert SPECSMetadataCSVParser.match_failure_reason(path) is None


# ── _parse: timestamp column fallback ───────────────────────────────────


def test_parse_uses_universal_time_when_present(tmp_path):
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Relative Time [s]", "Pressure [mbar] (Ceravac_AC)"],
        [
            ["2025-Aug-13 00:00:00.000000", "0.0", "1.0"],
            ["2025-Aug-13 00:00:05.000000", "5.0", "2.0"],
        ],
    )
    parser = SPECSMetadataCSVParser()
    parser.parse(path, csv_resample_interval=None)

    assert len(parser._parameter_logs) == 1
    row = parser._parameter_logs.iloc[0]
    assert row["measured_param"] == "Pressure"
    assert row["unit"] == "mbar"
    assert row["device"] == "Ceravac_AC"
    assert row["name"] == "Pressure (Ceravac_AC)"
    assert row["values"] == pytest.approx([1.0, 2.0])
    assert row["timestamps"] == [
        "2025-08-13T00:00:00.000000Z",
        "2025-08-13T00:00:05.000000Z",
    ]


def test_parse_falls_back_to_local_date_and_time_when_universal_time_absent(tmp_path):
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        [
            "Local Date",
            "Local Time of Day",
            "Relative Time [s]",
            "Pressure [mbar] (Ceravac_AC)",
        ],
        [["2025-08-13", "00:00:00.000000", "0.0", "1.0"]],
    )
    parser = SPECSMetadataCSVParser()
    parser.parse(path, csv_resample_interval=None)

    assert len(parser._parameter_logs) == 1
    row = parser._parameter_logs.iloc[0]
    assert row["values"] == pytest.approx([1.0])
    assert row["timestamps"] == ["2025-08-13T00:00:00.000000Z"]


def test_parse_skips_gracefully_when_no_time_column_is_present(tmp_path, caplog):
    """No usable time column -- must not crash, _parameter_logs stays empty.

    Uses ``_parse`` directly: this file also fails ``matches_file()`` (it
    requires >=1 time column), so ``parse()`` would never reach ``_parse``
    for it in practice; this exercises the defensive handling in isolation.
    """
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Relative Time [s]", "Pressure [mbar] (Ceravac_AC)"],
        [["0.0", "1.0"]],
    )
    parser = SPECSMetadataCSVParser()
    with caplog.at_level("WARNING"):
        parser._parse(path)  # must not raise
    assert parser._parameter_logs.empty
    assert any("time column" in record.getMessage() for record in caplog.records)


def test_parse_skips_gracefully_when_only_local_date_is_present(tmp_path):
    """Local Date without Local Time of Day -- partial pair, same graceful
    no-op as no time columns at all. See ``_parse`` note above."""
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Local Date", "Pressure [mbar] (Ceravac_AC)"],
        [["2025-08-13", "1.0"]],
    )
    parser = SPECSMetadataCSVParser()
    parser._parse(path)  # must not raise
    assert parser._parameter_logs.empty


# ── header parsing: unit / device extraction ────────────────────────────


@pytest.mark.parametrize(
    "header, expected",
    [
        (
            "CoilCurrent Current [µA] (Phoibos)",
            ("CoilCurrent Current", "µA", "Phoibos"),
        ),
        ("Pressure [mbar] (Ceravac_AC)", ("Pressure", "mbar", "Ceravac_AC")),
        ("SomeParam [V]", ("SomeParam", "V", "")),  # no device group
        ("NoBracketsAtAll", ("NoBracketsAtAll", "", "")),  # defensive fallback
    ],
)
def test_parse_header_extracts_param_unit_device(header, expected):
    assert SPECSMetadataCSVParser._parse_header(header) == expected


# ── sparse / blank-cell handling ────────────────────────────────────────


def test_parse_treats_blank_cells_as_no_observation_not_zero(tmp_path):
    """Blank cells mean "unchanged since last row", not 0 -- with
    resampling off, the raw series must only contain the two rows that
    actually carried a value for this parameter, skipping the blank row
    entirely (not inserting a 0.0)."""
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Pressure [mbar] (Ceravac_AC)"],
        [
            ["2025-Aug-13 00:00:00.000000", "1.0"],
            ["2025-Aug-13 00:00:01.000000", ""],
            ["2025-Aug-13 00:00:02.000000", "2.0"],
        ],
    )
    parser = SPECSMetadataCSVParser()
    parser.parse(path, csv_resample_interval=None)

    row = parser._parameter_logs.iloc[0]
    assert row["values"] == pytest.approx([1.0, 2.0])
    assert len(row["timestamps"]) == 2  # noqa: PLR2004 -- the blank row is skipped


def test_parse_default_resample_forward_fills_the_sparse_gaps(tmp_path):
    """Same forward-fill rationale as SPECSMetadataSLHParser._resample_df,
    applied to the CSV path's non-blank raw observations, so both formats
    produce the same output shape for the same underlying data."""
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Pressure [mbar] (Ceravac_AC)"],
        [
            ["2025-Aug-13 00:00:00.000000", "1.0"],
            ["2025-Aug-13 00:00:01.000000", ""],
            ["2025-Aug-13 00:00:02.000000", "2.0"],
        ],
    )
    parser = SPECSMetadataCSVParser()
    parser.parse(path)  # default csv_resample_interval="1s"

    row = parser._parameter_logs.iloc[0]
    assert len(row["timestamps"]) == 3  # noqa: PLR2004 -- 0s, 1s, 2s inclusive
    assert row["values"] == pytest.approx([1.0, 1.0, 2.0])


def test_parse_skips_a_parameter_column_with_no_observations_at_all(tmp_path):
    path = tmp_path / "export.csv"
    _write_csv(
        path,
        ["Universal Time", "Unused [V] (Dev)"],
        [["2025-Aug-13 00:00:00.000000", ""]],
    )
    parser = SPECSMetadataCSVParser()
    parser.parse(path)
    assert len(parser._parameter_logs) == 0


# ── update_main_file_data: per-entry slicing (shared with .slh) ────────


def test_update_main_file_data_slices_to_the_entry_time_window(tmp_path):
    path = tmp_path / "export.csv"
    rows = [[f"2025-Aug-13 00:00:{t:02d}.000000", str(float(t))] for t in range(0, 20)]
    _write_csv(
        path,
        ["Universal Time", "Pressure [mbar] (Ceravac_AC)"],
        rows,
    )
    parser = SPECSMetadataCSVParser()
    parser.parse(path, csv_resample_interval=None)

    spectrum = ParsedSpectrum(
        metadata={
            "time_stamp": "2025-Aug-13T00:00:05.000000Z",
            "dwell_time": 1.0,
            "n_values": 5,
            "total_scans": 1,
        }
    )
    parser.update_main_file_data({"entry": spectrum})

    key = "Pressure_Ceravac_AC"
    assert spectrum.metadata[f"{key}/@units"] == "mbar"
    # window is [5s, 5s + 1*5*1 = 10s]
    assert spectrum.metadata[key] == pytest.approx([5.0, 6.0, 7.0, 8.0, 9.0, 10.0])


def test_update_main_file_data_noop_when_no_parameter_logs():
    parser = SPECSMetadataCSVParser()
    spectrum = ParsedSpectrum(metadata={"time_stamp": "2025-Aug-13T00:00:00.000000Z"})
    parser.update_main_file_data({"entry": spectrum})
    assert spectrum.metadata == {"time_stamp": "2025-Aug-13T00:00:00.000000Z"}
