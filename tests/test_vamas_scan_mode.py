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
"""VAMAS files with a scan mode the parser cannot read are rejected."""

from pathlib import Path

import pytest

from pynxtools_xps.parsers.vms.parser import VamasParser

REGULAR_VMS = Path(__file__).parent / "data" / "vms_regular" / "regular.vms"


def _with_scan_mode(tmp_path, scan_mode):
    """Copy regular.vms and replace its scan mode line (REGULAR)."""
    lines = REGULAR_VMS.read_text().splitlines()
    index = lines.index("REGULAR")
    lines[index] = scan_mode
    file = tmp_path / f"{scan_mode.lower()}.vms"
    file.write_text("\n".join(lines) + "\n")
    return file


def test_mapping_scan_mode_is_not_supported_yet(tmp_path):
    # TODO: replace with a real MAPPING test once it is supported (see Issue #206).
    file = _with_scan_mode(tmp_path, "MAPPING")

    with pytest.raises(NotImplementedError, match="'MAPPING' is not supported yet"):
        VamasParser().parse(file)


def test_invalid_scan_mode_raises(tmp_path):
    file = _with_scan_mode(tmp_path, "SOMETHING_ELSE")

    with pytest.raises(ValueError, match="'SOMETHING_ELSE' is not a valid VAMAS scan"):
        VamasParser().parse(file)


def test_regular_scan_mode_still_parses():
    VamasParser().parse(REGULAR_VMS)
