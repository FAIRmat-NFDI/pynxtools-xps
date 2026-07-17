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
"""Tests for matches_file() on all XPS parsers."""

from pathlib import Path

import pytest

from pynxtools_xps.parsers.phi.parser import PHIParser
from pynxtools_xps.parsers.scienta.igor_parser import (
    ScientaIgorParser,
    ScientaIgorParserOld,
    ScientaIgorParserPEAK,
)
from pynxtools_xps.parsers.scienta.txt_parser import ScientaTXTParser
from pynxtools_xps.parsers.specs.sle.parser import SPECSSLEParser
from pynxtools_xps.parsers.specs.xml.parser import SPECSXMLParser
from pynxtools_xps.parsers.specs.xy.parser import SPECSXYParser
from pynxtools_xps.parsers.vms.parser import VamasParser
from pynxtools_xps.parsers.vms_export.parser_export import VamasExportParser
from pynxtools_xps.parsers.vms_export.parser_results import VamasResultParser
from pynxtools_xps.reader import XPSReader

DATA_DIR = Path(__file__).parent / "data"

# True positives: each parser must recognize its own format
_POSITIVE_CASES = [
    pytest.param(VamasParser, "vms_regular/regular.vms", id="vamas-vms"),
    pytest.param(
        VamasExportParser, "vms_txt_export/vms_txt_export.txt", id="vamas-export-txt"
    ),
    pytest.param(
        VamasResultParser,
        "vms_analysis/FeO_analyzed_results.csv",
        id="vamas-result-csv",
    ),
    pytest.param(
        SPECSXMLParser, "specs_xml/In-situ_PBTTT_XPS_SPECS.xml", id="specs-xml"
    ),
    pytest.param(SPECSXYParser, "specs_xy/MgFe2O4_small.xy", id="specs-xy"),
    pytest.param(SPECSSLEParser, "specs_sle/EX439_S718_Au.sle", id="specs-sle"),
    pytest.param(PHIParser, "phi_spe/SnO2_10nm.spe", id="phi-spe"),
    pytest.param(PHIParser, "phi_pro/SnO2_10nm_1.pro", id="phi-pro"),
    pytest.param(ScientaTXTParser, "scienta_txt/Ag_0001.txt", id="scienta-txt"),
    pytest.param(
        ScientaIgorParserOld, "scienta_ibw/Ag_0001Ag__001.ibw", id="scienta-igor-old"
    ),
    pytest.param(
        ScientaIgorParser, "scienta_ibw/Ag_0001Ag__001.ibw", id="scienta-igor-wrapper"
    ),
    # ScientaIgorParserPEAK and ScientaHDF5Parser: no PEAK/HDF5 files in tests/data/
]

# Cross-format rejections: each parser must reject files of a different format
_NEGATIVE_CASES = [
    pytest.param(VamasParser, "specs_xy/MgFe2O4_small.xy", id="vamas-rejects-xy"),
    pytest.param(
        VamasParser,
        "vms_txt_export/vms_txt_export.txt",
        id="vamas-rejects-casa-xps-export",
    ),
    pytest.param(SPECSXMLParser, "vms_regular/regular.vms", id="specs-xml-rejects-vms"),
    pytest.param(
        SPECSXMLParser, "specs_xy/MgFe2O4_small.xy", id="specs-xml-rejects-xy"
    ),
    pytest.param(PHIParser, "vms_regular/regular.vms", id="phi-rejects-vms"),
    pytest.param(
        VamasExportParser,
        "scienta_txt/Ag_0001.txt",
        id="casa-xps-export-rejects-scienta-txt",
    ),
    pytest.param(
        VamasExportParser, "vms_regular/regular.vms", id="casa-xps-export-rejects-vms"
    ),
    pytest.param(
        ScientaTXTParser,
        "vms_txt_export/vms_txt_export.txt",
        id="scienta-txt-rejects-casa-xps-export",
    ),
    pytest.param(SPECSXYParser, "vms_regular/regular.vms", id="specs-xy-rejects-vms"),
    pytest.param(
        ScientaIgorParserOld,
        "vms_regular/regular.vms",
        id="scienta-igor-old-rejects-vms",
    ),
    pytest.param(
        ScientaIgorParserPEAK,
        "scienta_ibw/Ag_0001Ag__001.ibw",
        id="scienta-igor-peak-rejects-old-ibw",
    ),
]


@pytest.mark.parametrize("parser_cls, rel_path", _POSITIVE_CASES)
def test_matches_file_positive(parser_cls, rel_path):
    """Each parser must recognize its own file format."""
    path = DATA_DIR / rel_path
    if not path.exists():
        pytest.skip(f"test file not found: {path}")
    assert parser_cls().matches_file(path)


@pytest.mark.parametrize("parser_cls, rel_path", _NEGATIVE_CASES)
def test_matches_file_negative(parser_cls, rel_path):
    """Each parser must reject files of a different format."""
    path = DATA_DIR / rel_path
    if not path.exists():
        pytest.skip(f"test file not found: {path}")
    assert not parser_cls().matches_file(path)


@pytest.mark.parametrize("parser_cls, rel_path", _NEGATIVE_CASES)
def test_is_mainfile_rejection_is_silent(parser_cls, rel_path, caplog):
    """is_mainfile() is probed speculatively against every registered parser,
    so a rejection must not itself log a warning (see match_failure_reason
    for on-demand diagnostics instead)."""
    path = DATA_DIR / rel_path
    if not path.exists():
        pytest.skip(f"test file not found: {path}")
    with caplog.at_level("WARNING"):
        assert not parser_cls.is_mainfile(path)
    assert caplog.records == []


def _expected_match_failure_reason(parser_cls, path: Path) -> str:
    """The exact diagnostic match_failure_reason() must produce for *path*,
    derived independently from the same public classification the parser
    uses (extension support), not by re-deriving matches_file()'s verdict."""
    if not parser_cls.is_extension_supported(path):
        return parser_cls.file_ext_err_msg(path)
    return parser_cls.matches_file_warning(path)


@pytest.mark.parametrize("parser_cls, rel_path", _NEGATIVE_CASES)
def test_match_failure_reason_explains_rejection_without_logging(
    parser_cls, rel_path, caplog
):
    """match_failure_reason() gives the same diagnosis as the ValueError raised
    internally, but as a plain return value, and without logging on its own."""
    path = DATA_DIR / rel_path
    if not path.exists():
        pytest.skip(f"test file not found: {path}")
    with caplog.at_level("WARNING"):
        reason = parser_cls.match_failure_reason(path)
    assert reason == _expected_match_failure_reason(parser_cls, path)
    assert caplog.records == []


@pytest.mark.parametrize("parser_cls, rel_path", _POSITIVE_CASES)
def test_match_failure_reason_none_on_match(parser_cls, rel_path):
    """match_failure_reason() returns None whenever matches_file() succeeds."""
    path = DATA_DIR / rel_path
    if not path.exists():
        pytest.skip(f"test file not found: {path}")
    assert parser_cls.match_failure_reason(path) is None


def test_reader_does_not_warn_on_extension_overlap(caplog):
    """Scienta TXT and CasaXPS/Vamas-export TXT files share the .txt extension,
    so both parsers are probed for either file. Dispatching a file that one of
    them matches must not surface the other's rejection as a warning."""
    scienta_txt = DATA_DIR / "scienta_txt/Ag_0001.txt"
    vamas_export_txt = DATA_DIR / "vms_txt_export/vms_txt_export.txt"
    for path in (scienta_txt, vamas_export_txt):
        if not path.exists():
            pytest.skip(f"test file not found: {path}")

    with caplog.at_level("WARNING"):
        XPSReader().handle_data_file(str(scienta_txt))
        XPSReader().handle_data_file(str(vamas_export_txt))

    assert caplog.records == []


def test_reader_reports_all_reasons_when_no_parser_matches(tmp_path, caplog):
    """When no registered parser matches a file, the reader logs one warning
    that aggregates every parser's individual rejection reason."""
    unmatched_file = tmp_path / "unmatched.foo"
    unmatched_file.write_bytes(b"not a real XPS file")

    with caplog.at_level("WARNING"):
        XPSReader().handle_data_file(str(unmatched_file))

    assert len(caplog.records) == 1
    message = caplog.records[0].getMessage()
    assert "No parser matches file" in message
    for parser_cls in (*XPSReader.parsers, *XPSReader.metadata_parsers):
        assert parser_cls.__name__ in message
