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
"""Tests for peak_finding.find_xps_peaks."""

import os

import h5py
import numpy as np

from pynxtools_xps.peak_finding import find_xps_peaks

REF_FILE = os.path.join(
    os.path.dirname(__file__), "data", "vms_regular", "vms_regular_ref.nxs"
)


def test_find_xps_peaks_appends_nxfit(tmp_path):
    output_path = tmp_path / "vms_regular_with_peaks.nxs"

    with h5py.File(REF_FILE, "r") as ref_file:
        original_energy = ref_file["/1_as_loaded__Survey/data/energy"][()]
        original_intensity = ref_file["/1_as_loaded__Survey/data/data"][()]

    peaks = find_xps_peaks(REF_FILE, output_path, prominence=1.0)

    assert len(peaks) > 0

    with h5py.File(output_path, "r") as out_file:
        entry = out_file["/1_as_loaded__Survey"]

        # original data is untouched
        np.testing.assert_array_equal(entry["data/energy"][()], original_energy)
        np.testing.assert_array_equal(entry["data/data"][()], original_intensity)

        fit_group = entry["fit"]
        assert fit_group.attrs.get("NX_class") == "NXfit"

        for i, peak in enumerate(peaks):
            peak_group = fit_group[f"peak_{i}"]
            assert peak_group.attrs.get("NX_class") == "NXpeak"
            assert peak_group["label"].asstr()[()] == f"peak_{i}"
            peak_data_group = peak_group["data"]
            assert peak_data_group.attrs.get("NX_class") == "NXdata"
            assert peak_data_group["position"][()] == peak["position"]
            assert peak_data_group["position"].attrs.get("units") == "eV"
            assert peak_data_group["intensity"][()] == peak["intensity"]
            assert (
                peak_data_group["intensity"].attrs.get("units") == "counts_per_second"
            )

    # find_xps_peaks(REF_FILE) with no output_path defaults to in-place, and the
    # reference file itself must stay untouched by this test either way.
    with h5py.File(REF_FILE, "r") as ref_file:
        np.testing.assert_array_equal(
            ref_file["/1_as_loaded__Survey/data/energy"][()], original_energy
        )
