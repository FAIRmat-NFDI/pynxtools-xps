#
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
Basic example based test for the XPS reader
"""

import json
import os
from typing import Literal

import numpy as np
import pytest
from pynxtools.dataconverter.convert import get_reader
from pynxtools.testing.nexus_conversion import ReaderTest

from pynxtools_xps.vms.vamas_comment_handler import handle_comments

READER_NAME = "xps"
READER_CLASS = get_reader(READER_NAME)
# TODO: make tests for all supported application definitions possible
NXDLS = ["NXxps"]  # READER_CLASS.supported_nxdls

test_cases = [
    ("phi_spe", "phi-spe-reader", {}),
    ("phi_pro", "phi-pro-reader", {}),
    ("specs_sle", "specs-sle-reader", {}),
    ("specs_xml", "specs-xml-reader", {}),
    ("specs_xy", "specs-xy-reader", {}),
    (
        "scienta_ibw",
        "scienta-ibw-reader",
        {"FIELD (//Ag__002__VB/start_time)": ["DEBUG - value:"]},
    ),
    (
        "scienta_txt",
        "scienta-txt-reader",
        {"FIELD (//Ag__002__VB/start_time)": ["DEBUG - value:"]},
    ),
    ("vms_analysis", "vms-reader-with-data-analysis", {}),
    ("vms_irregular", "irregular-vms-reader", {}),
    ("vms_regular", "regular-vms-reader", {}),
    ("vms_txt_export", "vms-txt-export-reader", {}),
]

test_params = []

for test_case in test_cases:
    # ToDo: make tests for all supported appdefs possible
    for nxdl in NXDLS:
        test_params += [
            pytest.param(
                nxdl, test_case[0], test_case[2], id=f"{test_case[1]}-{nxdl.lower()}"
            )
        ]


@pytest.mark.parametrize(
    "nxdl, sub_reader_data_dir, ignore_sections",
    test_params,
)
def test_nexus_conversion(nxdl, sub_reader_data_dir, ignore_sections, tmp_path, caplog):
    """
    Test XPS reader

    Parameters
    ----------
    nxdl : str
        Name of the NXDL application definition that is to be tested by
        this reader plugin (e.g. NXxps, NXmpes, etc)..
    sub_reader_data_dir : str
        Test data directory that contains all the files required for running the data
        conversion through one of the sub-readers. All of these data dirs
        are placed within tests/data/...
    ignore_sections: Dict[str, List[str]]
        Subsections of the log file to ignore.
    tmp_path : pathlib.PosixPath
        Pytest fixture variable, used to clean up the files generated during
        the test.
    caplog : _pytest.logging.LogCaptureFixture
        Pytest fixture variable, used to capture the log messages during the
        test.

    Returns
    -------
    None.

    """
    caplog.clear()
    reader = READER_CLASS
    assert callable(reader.read)

    files_or_dir = os.path.join(
        *[os.path.dirname(__file__), "data", sub_reader_data_dir]
    )

    test = ReaderTest(
        nxdl=nxdl,
        reader_name=READER_NAME,
        files_or_dir=files_or_dir,
        tmp_path=tmp_path,
        caplog=caplog,
    )
    test.convert_to_nexus(caplog_level="WARNING", ignore_undocumented=True)
    test.check_reproducibility_of_nexus(ignore_sections=ignore_sections)


def read_comment_file(filepath: str):
    """Read comments from one vms comment test file."""

    no_of_comments = 0
    comment_lines = []

    with open(filepath, "rb") as vms_file:
        for i, line in enumerate(vms_file):
            if i == 0:
                no_of_comments = int(line)
            else:
                comment_lines += [line.decode("utf-8", errors="ignore").strip()]

    return no_of_comments, comment_lines


@pytest.mark.parametrize(
    "file, comment_type, expected_no_of_comments",
    [
        pytest.param("kratos.vms", "block", 51, id="Kratos metadata"),
        pytest.param("phi.vms", "block", 312, id="PHI metadata"),
        pytest.param("casa_header.vms", "header", 1, id="CasaXPS header"),
        pytest.param("casa_process.vms", "block", 151, id="CasaXPS processing"),
        pytest.param("specs_header.vms", "header", 3, id="SPECS header metadata"),
        pytest.param("specs_block.vms", "block", 11, id="SPECS header metadata"),
    ],
)
def test_vms_comment_handler(
    file: str, comment_type: Literal["header", "block"], expected_no_of_comments: bool
):
    """Test for the comment handler in VAMAS files."""
    filepath = os.path.join(os.path.dirname(__file__), "data", "vms_comments", file)
    ref_json_filepath = filepath.replace(".vms", "_ref.json")

    no_of_comments, comment_lines = read_comment_file(filepath)
    assert no_of_comments == len(comment_lines)

    comments = handle_comments(comment_lines, comment_type=comment_type)

    if file == "casa_process.vms":
        casa_process = comments["casa"]
        comments = casa_process.flatten_metadata()

    for key, val in comments.items():
        if isinstance(val, np.ndarray):
            comments[key] = val.tolist()

    with open(ref_json_filepath) as json_file:
        ref_comments = json.load(json_file)

    assert len(comments) == expected_no_of_comments, (
        f"Comments ({len(comments)}) do not have the same number of lines as the reference ({expected_no_of_comments})."
    )
    assert comments == ref_comments
