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
Tests for the energy step size (Bug B21 in
https://github.com/FAIRmat-NFDI/pynxtools-plugin-audit/blob/main/pynxtools-unification/FINDINGS.md).
"""

from pathlib import Path

import numpy as np
import pytest

from pynxtools_xps.numerics import _get_minimal_step, interpolate_arrays
from pynxtools_xps.parsers.vms.parser import VamasParser

DATA_DIR = Path(__file__).parent / "data"


@pytest.mark.parametrize(
    "step, expected",
    [
        (0.005, 0.005),
        (0.002, 0.002),
        (0.025, 0.025),
        (0.013130282676, 0.0131),
    ],
)
def test_minimal_step_keeps_four_decimals(step, expected):
    """The step must keep 4 decimals, not be rounded to 2 (B21)."""
    # _get_minimal_step used to round to 2 decimals, so small steps became
    # wrong (0.025 -> 0.02) or zero (0.005 -> 0.0). 4 decimals is enough
    # precision and still hides float noise. Same rule as in pynxtools-xrd.
    axis = 280.0 + step * np.arange(100)
    found = _get_minimal_step(axis)
    assert found == expected


def test_interpolation_keeps_small_steps():
    """Resampling a non-uniform axis must use the real minimal step (B21)."""
    # With the old 2-decimal rounding the step became 0.0 here, and building
    # the new axis divided by zero.
    x = [280.0, 280.005, 280.015, 280.02]
    new_x, _ = interpolate_arrays(x, [np.array([1.0, 2.0, 4.0, 5.0])])
    assert list(new_x) == pytest.approx([280.0, 280.005, 280.01, 280.015, 280.02])


def test_vamas_regular_axis_keeps_small_steps(tmp_path):
    """A VAMAS file with a 0.005 eV step must give a 0.005 eV energy axis (B21)."""
    # No real VAMAS file with such a small step is available, so the test
    # copies regular.vms (step 1 eV) and only changes the abscissa increment.
    # The parser used to round every energy to 2 decimals, so neighboring
    # points got the same energy or jumped by 0.01.
    lines = (DATA_DIR / "vms_regular" / "regular.vms").read_text().splitlines()
    start = lines.index("kinetic energy")
    assert lines[start + 3] == "1"  # abscissa increment in eV
    lines[start + 3] = "0.005"
    vms_file = tmp_path / "small_step.vms"
    vms_file.write_text("\n".join(lines) + "\n")

    parser = VamasParser()
    parser.parse(vms_file)
    spectrum = next(iter(parser._data.values()))
    energy = spectrum.data.coords["energy"].values
    assert list(np.diff(energy)) == pytest.approx([0.005] * (len(energy) - 1))
