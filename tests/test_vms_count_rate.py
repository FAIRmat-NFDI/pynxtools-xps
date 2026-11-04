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
"""Tests for the VAMAS count rate (counts per second)."""

from pathlib import Path

import pytest

from pynxtools_xps.parsers.vms.parser import VamasParser

DATA_DIR = Path(__file__).parent / "data"


def test_vamas_count_rate_is_not_rounded(tmp_path):
    """The count rate must be counts / dwell time, without rounding."""
    # The parser used to round the count rate to 2 decimals, so weak signals
    # lost their precision or became 0. The test copies regular.vms
    # (dwell time 0.1 s) and sets a long dwell time of 1000 s, which gives
    # count rates between about 0.02 and 11 counts/s.
    lines = (DATA_DIR / "vms_regular" / "regular.vms").read_text().splitlines()
    signal_mode = lines.index("pulse counting")
    assert lines[signal_mode + 1] == "0.1"  # dwell time in s
    lines[signal_mode + 1] = "1000"
    vms_file = tmp_path / "long_dwell_time.vms"
    vms_file.write_text("\n".join(lines) + "\n")

    parser = VamasParser()
    parser.parse(vms_file)
    spectrum = next(iter(parser._data.values()))
    counts = spectrum.raw.values.ravel()
    count_rate = spectrum.data.values.ravel()
    assert list(count_rate) == pytest.approx(list(counts / 1000))
