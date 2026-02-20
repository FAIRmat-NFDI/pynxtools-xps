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

from pynxtools_xps.mapping import _MetadataContext, _ValueMap

_KEY_MAP: dict[str, str] = {
    "K.E.": "kinetic_energy",
    "B.E.": "binding_energy",
    "Counts": "counts",
    "CPS": "counts_per_second",
    "Background": "background_intensity",
    "Background CPS": "background_intensity_cps",
    "Envelope": "fit_sum",
    "Envelope CPS": "fit_sum_cps",
    "%At Conc": "atomic_concentration",
}

_VALUE_MAP: _ValueMap = {}
_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {
    "step_size": "eV",
}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
