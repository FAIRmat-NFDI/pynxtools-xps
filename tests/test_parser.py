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
from pynxtools_xps.parsers.specs.sle.parser import SpecsSLEParser
from pynxtools_xps.parsers.specs.xml.parser import SpecsXMLParser
from pynxtools_xps.parsers.specs.xy.parser import SpecsXYParser
from pynxtools_xps.parsers.vms.parser import VamasParser
from pynxtools_xps.parsers.vms_export.parser_export import VamasExportParser
from pynxtools_xps.parsers.vms_export.parser_results import VamasResultParser

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
        SpecsXMLParser, "specs_xml/In-situ_PBTTT_XPS_SPECS.xml", id="specs-xml"
    ),
    pytest.param(SpecsXYParser, "specs_xy/MgFe2O4_small.xy", id="specs-xy"),
    pytest.param(SpecsSLEParser, "specs_sle/EX439_S718_Au.sle", id="specs-sle"),
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
    pytest.param(SpecsXMLParser, "vms_regular/regular.vms", id="specs-xml-rejects-vms"),
    pytest.param(
        SpecsXMLParser, "specs_xy/MgFe2O4_small.xy", id="specs-xml-rejects-xy"
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
    pytest.param(SpecsXYParser, "vms_regular/regular.vms", id="specs-xy-rejects-vms"),
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
