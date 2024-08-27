"""
Basic example based test for the XPS reader
"""

import os
import pytest
import json
from typing import Literal
import numpy as np

from pynxtools.dataconverter.convert import get_reader
from pynxtools.testing.nexus_conversion import ReaderTest

from pynxtools_xps.vms.vamas_comment_handler import handle_comments


READER_NAME = "xps"
READER_CLASS = get_reader(READER_NAME)

test_cases = [
    ("phi_spe", "Phi .spe reader"),
    ("phi_pro", "Phi .pro reader"),
    ("specs_sle", "SPECS .sle reader"),
    ("specs_xml", "SPECS .xml reader"),
    ("specs_xy", "SPECS .xy reader"),
    ("scienta_ibw", "Scienta .ibw reader"),
    ("scienta_txt", "Scienta .txt export reader"),
    ("vms_irregular", "Irregular VAMAS reader"),
    ("vms_regular", "Regular VAMAS reader"),
    ("vms_txt_export", "Vamas txt export"),
]

test_params = []

for test_case in test_cases:
    for nxdl in READER_CLASS.supported_nxdls:
        test_params += [pytest.param(nxdl, test_case[0], id=f"{test_case[1]}, {nxdl}")]


@pytest.mark.parametrize(
    "nxdl, sub_reader_data_dir",
    test_params,
)
def test_nexus_conversion(nxdl, sub_reader_data_dir, tmp_path, caplog):
    """
    Test XPS reader

    Parameters
    ----------
    nxdl : str
        Name of the NXDL application definition that is to be tested by
        this reader plugin (e.g. NXsts, NXmpes, etc)..
    sub_reader_data_dir : str
        Test data directory that contains all the files required for running the data
        conversion through one of the sub-readers. All of these data dirs
        are placed within tests/data/...
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
    test.check_reproducibility_of_nexus()


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


COMMENT_LENGHTS = {
    "kratos.vms": 51,
    "phi.vms": 312,
    "casa_header.vms": 1,
    "casa_process.vms": 144,
    "specs_header.vms": 3,
    "specs_block.vms": 11,
}


@pytest.mark.parametrize(
    "file, comment_type",
    [
        pytest.param("kratos.vms", "block", id="Kratos metadata"),
        pytest.param("phi.vms", "block", id="PHI metadata"),
        pytest.param("casa_header.vms", "header", id="CasaXPS header"),
        pytest.param("casa_process.vms", "block", id="CasaXPS processing"),
        pytest.param("specs_header.vms", "header", id="SPECS header metadata"),
        pytest.param("specs_block.vms", "block", id="SPECS header metadata"),
    ],
)
def test_vms_comment_handler(file: str, comment_type: Literal["header", "block"]):
    """Test for the comment handler in VAMAS files."""
    filepath = os.path.join(os.path.dirname(__file__), "data", "vms_comments", file)
    ref_json_filepath = filepath.replace(".vms", "_ref.json")

    no_of_comments, comment_lines = read_comment_file(filepath)
    assert no_of_comments == len(comment_lines)

    comments = handle_comments(comment_lines, comment_type=comment_type)

    for key, val in comments.items():
        if isinstance(val, np.ndarray):
            comments[key] = val.tolist()

    with open(ref_json_filepath, "r") as json_file:
        ref_comments = json.load(json_file)

    assert len(comments) == COMMENT_LENGHTS[file]
    assert comments == ref_comments
