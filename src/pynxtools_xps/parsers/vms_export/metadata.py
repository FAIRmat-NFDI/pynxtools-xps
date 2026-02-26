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
Metadata mapping for the VAMAS TXT export parser.
"""

from collections import Counter

from pynxtools_xps.mapping import _MetadataContext, _ValueMap


def _handle_repetitions(input_list: list[str]) -> list[str]:
    """
    Process a list of strings to handle repeated items by appending a suffix
    to each duplicate item. The suffix is in the format '_n', where 'n' is the
    occurrence number of that item in the list.

    Parameters:
    - input_list (List[str]): A list of strings where repeated items are
      identified and renamed with a suffix.

    Returns:
    - List[str]: A new list where repeated items are modified by appending
      a suffix to make them unique.
    """
    counts = Counter(input_list)
    result = []
    occurrences = {}

    for item in input_list:
        if counts[item] > 1:
            # If the item has been seen before, add a suffix
            if item not in occurrences:
                occurrences[item] = 0
            occurrences[item] += 1
            result.append(f"{item}_{occurrences[item]}")
        else:
            result.append(item)

    return result


_KEY_MAP: dict[str, str] = {
    "k.e.": "kinetic_energy",
    "b.e.": "binding_energy",
    "cps": "counts_per_second",
    "background": "background_intensity",
    "background_cps": "background_intensity_cps",
    "envelope": "fit_sum",
    "envelope_cps": "fit_sum_cps",
    "%at_conc": "atomic_concentration",
    "fwhm": "width",
    "area/(rsf*t*mfp)": "area_over_rsf*t*mfp",
}

_VALUE_MAP: _ValueMap = {}
_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {"step_size": "eV", "width": "eV", "position": "eV"}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
