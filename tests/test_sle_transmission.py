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
"""The SPECS SLE transmission function is only written when the file has one."""

from pathlib import Path

import numpy as np

from pynxtools_xps.parsers.specs.sle.parser import SPECSSLEParser

SLE_FILE = Path(__file__).parent / "data" / "specs_sle" / "EX439_S718_Au.sle"


def _spectrum(n_energies: int, transmission: np.ndarray | None) -> dict:
    return {
        "energy_scan_mode": "fixed_analyzer_transmission",
        "energy/@type": "kinetic",
        "data": {"energy": np.arange(n_energies, dtype=float)},
        "transmission_function/relative_intensity": transmission,
    }


def test_file_without_transmission_function_gets_none():
    """No transmission function is invented (it used to be all ones)."""
    parser = SPECSSLEParser()
    parser.parse(SLE_FILE)  # this file has no transmission data

    for spectrum in parser._data.values():
        assert spectrum.metadata.get("transmission_function/relative_intensity") is None
        assert spectrum.metadata.get("transmission_function/kinetic_energy") is None


def test_transmission_function_with_matching_length_is_kept():
    transmission = np.array([0.5, 0.6, 0.7])
    spectrum = _spectrum(3, transmission)

    SPECSSLEParser()._compute_tf_energy_axis(spectrum)

    np.testing.assert_array_equal(
        spectrum["transmission_function/relative_intensity"], transmission
    )
    np.testing.assert_array_equal(
        spectrum["transmission_function/kinetic_energy"], [0.0, 1.0, 2.0]
    )


def test_transmission_function_with_other_length_is_cut_with_a_warning(caplog):
    """A length mismatch is kept as before (cut to the shorter one), but reported."""
    spectrum = _spectrum(5, np.array([0.5, 0.6, 0.7]))

    SPECSSLEParser()._compute_tf_energy_axis(spectrum)

    np.testing.assert_array_equal(
        spectrum["transmission_function/relative_intensity"], [0.5, 0.6, 0.7]
    )
    np.testing.assert_array_equal(
        spectrum["transmission_function/kinetic_energy"], [0.0, 1.0, 2.0]
    )
    assert "has 3 values, but the energy axis has 5" in caplog.text
