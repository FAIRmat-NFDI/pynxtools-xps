#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""Tests for the pynxtools reader plugin."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import numpy as np
import pytest

from pynxtools_xps.parsers.vms.comment_handler import handle_comments


def _read_comment_file(filepath: Path) -> tuple[int, list[str]]:
    """Read the declared number of comments and the comment lines from a VMS test file.

    Parameters
    ----------
    filepath
        Path to the VMS comment test file.

    Returns
    -------
    tuple[int, list[str]]
        The declared number of comments and the list of decoded comment lines.
    """
    no_of_comments: int = 0
    comment_lines: list[str] = []

    with filepath.open("rb") as vms_file:
        for i, line in enumerate(vms_file):
            if i == 0:
                no_of_comments = int(line)
            else:
                comment_lines.append(line.decode("utf-8", errors="ignore").strip())

    return no_of_comments, comment_lines


@pytest.mark.parametrize(
    "file, comment_type, expected_no_of_comments",
    [
        pytest.param("kratos.vms", "block", 52, id="Kratos metadata"),
        pytest.param("phi.vms", "block", 312, id="PHI metadata"),
        pytest.param("casa_header.vms", "header", 1, id="CasaXPS header"),
        pytest.param("casa_process.vms", "block", 151, id="CasaXPS processing"),
        pytest.param("specs_header.vms", "header", 3, id="SPECS header metadata"),
        pytest.param("specs_block.vms", "block", 8, id="SPECS block metadata"),
    ],
)
def test_vms_comment_handler(
    file: str,
    comment_type: Literal["header", "block"],
    expected_no_of_comments: int,
) -> None:
    """Test the VAMAS comment handler against reference JSON outputs."""
    base_dir: Path = Path(__file__).parent / "data" / "vms_comments"
    filepath: Path = base_dir / file
    ref_json_filepath: Path = filepath.with_name(f"{filepath.stem}_ref.json")

    no_of_comments, comment_lines = _read_comment_file(filepath)
    assert no_of_comments == len(comment_lines)

    comments: dict[str, object] = handle_comments(
        comment_lines, comment_type=comment_type
    )

    if file == "casa_process.vms":
        casa_process = comments["casa"]
        comments = casa_process.flatten_metadata()  # type: ignore[attr-defined]

    # Convert numpy arrays to JSON-comparable lists
    for key, val in list(comments.items()):
        if isinstance(val, np.ndarray):
            comments[key] = val.tolist()

    with ref_json_filepath.open() as json_file:
        ref_comments: dict[str, object] = json.load(json_file)

    assert len(comments) == expected_no_of_comments, (
        f"Comments ({len(comments)}) do not have the same number of lines "
        f"as the reference ({expected_no_of_comments})."
    )
    assert comments == ref_comments
