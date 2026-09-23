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
Peak finding for an already-converted XPS NeXus file.

Unlike the rest of pynxtools_xps, this module does not convert vendor data to
NeXus -- it reads a spectrum back out of an existing NeXus/HDF5 file, finds
peaks in it, and appends the result as an NXfit/NXpeak group. It only depends
on pynxtools (for NXdata resolution and NeXus writing) and scipy, so it can
also run standalone, outside a NOMAD/pynxtools-xps environment.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import h5py
from pynxtools.dataconverter.helpers import add_default_root_attributes
from pynxtools.dataconverter.template import Template
from pynxtools.dataconverter.writer import Writer
from pynxtools.nexus.nxdata import (
    find_default_nxdata,
    find_default_nxentry,
    inspect_nxdata,
)
from pynxtools.nexus.utils import decode_if_string, get_nxdl_root_and_path
from scipy.signal import find_peaks


def find_xps_peaks(
    input_path: str | Path,
    output_path: str | Path | None = None,
    **find_peaks_kwargs: Any,
) -> list[dict[str, float]]:
    """Find peaks in an XPS spectrum and append them to the NeXus file as NXfit/NXpeak.

    Reads the default NXdata group's signal/axes (energy vs. intensity) from an
    already-converted NeXus file, runs `scipy.signal.find_peaks` on the intensity,
    and appends one NXfit/peakPEAK group per detected peak with the peak's energy
    (`data/position`) and intensity (`data/intensity`).

    This writes a deliberately partial NXfit: `scipy.signal.find_peaks` locates
    peaks but does not fit a background or a parametric peak shape, so NXfit's
    own required `data` (fit input/output) and `backgroundBACKGROUND` groups are
    left out rather than filled with invented values. `pynx validate` will flag
    both as missing -- that reflects this step's actual scope (peak finding, not
    peak fitting), not a bug.

    Args:
        input_path: path to an existing NeXus file with a resolvable default NXdata
            group (energy vs. intensity).
        output_path: where to write the result. Defaults to `input_path` (in place).
            If different, `input_path` is copied first and left untouched.
        **find_peaks_kwargs: forwarded to `scipy.signal.find_peaks` (e.g. `height`,
            `prominence`, `distance`).

    Returns:
        One dict per detected peak, with `position` (energy) and `intensity`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path) if output_path is not None else input_path
    if output_path != input_path:
        shutil.copyfile(input_path, output_path)

    with h5py.File(output_path, "r") as nexus_file:
        entry = find_default_nxentry(nexus_file)
        if entry is None:
            raise ValueError(f"No NXentry found in {output_path}")
        entry_name = entry.name.lstrip("/")

        definition = (
            decode_if_string(entry["definition"][()]) if "definition" in entry else None
        )
        if not definition:
            raise ValueError(
                f"Entry '{entry_name}' in {output_path} has no /definition field; "
                "cannot determine which NXDL application definition to write against."
            )

        nxdata_group = find_default_nxdata(nexus_file)
        if nxdata_group is None:
            raise ValueError(f"No default NXdata group found in {output_path}")
        nxdata_info = inspect_nxdata(nxdata_group)
        if (
            nxdata_info.signal is None
            or not nxdata_info.axes
            or not nxdata_info.axes[0]
        ):
            raise ValueError(
                f"NXdata group '{nxdata_group.name}' in {output_path} has no "
                "resolvable signal/axes."
            )

        intensity = nxdata_info.signal[()]
        energy = nxdata_info.axes[0][0][()]
        position_units = nxdata_info.axes[0][0].attrs.get("units")
        intensity_units = nxdata_info.signal.attrs.get("units")

    peak_indices, _properties = find_peaks(intensity, **find_peaks_kwargs)

    _nxdl_root, nxdl_f_path = get_nxdl_root_and_path(definition)

    template = Template()
    add_default_root_attributes(template, filename=str(output_path), append=True)

    # NXfit is a recommended child group directly under ENTRY (a sibling of the
    # spectrum's own NXdata group), not nested under it -- see NXxps.nxdl.xml.
    fit_base = f"/ENTRY[{entry_name}]/FIT[fit]"
    template[f"{fit_base}/label"] = "Peaks found via scipy.signal.find_peaks"

    peaks: list[dict[str, float]] = []
    for i, index in enumerate(peak_indices):
        position = float(energy[index])
        peak_intensity = float(intensity[index])
        peaks.append({"position": position, "intensity": peak_intensity})
        peak_id = f"peak_{i}"
        peak_base = f"{fit_base}/peakPEAK[{peak_id}]"
        template[f"{peak_base}/label"] = peak_id
        template[f"{peak_base}/DATA[data]/position"] = position
        template[f"{peak_base}/DATA[data]/position/@units"] = position_units
        template[f"{peak_base}/DATA[data]/intensity"] = peak_intensity
        template[f"{peak_base}/DATA[data]/intensity/@units"] = intensity_units

    writer = Writer(
        data=template,
        nxdl_f_path=nxdl_f_path,
        output_path=str(output_path),
        append=True,
    )
    writer.write()

    return peaks


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=find_xps_peaks.__doc__)
    parser.add_argument(
        "input_path", help="Path to an already-converted XPS NeXus file."
    )
    parser.add_argument(
        "-o", "--output-path", default=None, help="Output path (defaults to in-place)."
    )
    parser.add_argument("--height", type=float, default=None)
    parser.add_argument("--prominence", type=float, default=None)
    parser.add_argument("--distance", type=float, default=None)
    args = parser.parse_args()

    kwargs = {
        k: v
        for k, v in {
            "height": args.height,
            "prominence": args.prominence,
            "distance": args.distance,
        }.items()
        if v is not None
    }
    found_peaks = find_xps_peaks(args.input_path, args.output_path, **kwargs)
    print(json.dumps(found_peaks, indent=2))
