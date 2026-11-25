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
"""Tests for VAMAS files with a binding energy axis."""

from pathlib import Path

import numpy as np
import pytest

from pynxtools_xps.parsers.vms.parser import VamasParser

DATA_DIR = Path(__file__).parent / "data"


def _parse(vms_file: Path):
    parser = VamasParser()
    parser.parse(vms_file)
    return next(iter(parser._data.values()))


@pytest.mark.parametrize(
    "file_name, kinetic_label",
    [
        pytest.param("vms_regular/regular.vms", "kinetic energy", id="regular"),
        pytest.param("vms_irregular/irregular.vms", "Kinetic Energy", id="irregular"),
    ],
)
@pytest.mark.parametrize(
    "binding_label", ["binding energy", "Binding Energy", "BINDING ENERGY"]
)
def test_binding_energy_axis_keeps_energy_and_intensity_together(
    tmp_path, file_name, kinetic_label, binding_label
):
    """Each energy must stay paired with its own intensity."""
    # The parser used to reverse only the energy axis of a "binding energy"
    # block, not the intensities, so every energy got the intensity of
    # another point. The label was also compared case-sensitively.
    # The test relabels the energy axis of a test file and checks that the
    # energy/intensity pairs are the same as in the original file.
    original_file = DATA_DIR / file_name
    lines = original_file.read_text().splitlines()
    assert lines.count(kinetic_label) == 1
    lines[lines.index(kinetic_label)] = binding_label
    binding_file = tmp_path / "binding_energy.vms"
    binding_file.write_text("\n".join(lines) + "\n")

    original = _parse(original_file)
    relabelled = _parse(binding_file)

    np.testing.assert_array_equal(
        relabelled.raw.coords["energy"].values, original.raw.coords["energy"].values
    )
    np.testing.assert_array_equal(relabelled.raw.values, original.raw.values)
    assert relabelled.metadata["energy_label"] == "binding"
