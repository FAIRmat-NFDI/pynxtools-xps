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
Small util functionalities.
"""

from pynxtools_xps.mapping import _Value


# TODO: can these be inlined and the file be removed?
def update_dict_without_overwrite(d1: dict[str, _Value], d2: dict[str, _Value]):
    """Update d1 with d2, but don't overwrite existing keys."""
    d1.update({k: v for k, v in d2.items() if k not in d1})


def _drop_unused_keys(dictionary: dict[str, _Value], keys_to_drop: list[str]) -> None:
    """
    Remove unwanted keys from a dictionary.

    Args:
        dictionary (dict[str, Any]):
            Dictionary containing data and metadata for a spectrum.
        keys_to_drop (list[str]):
            List of keys that should be removed from the dictionary.

    Returns:
        None
    """
    for key in keys_to_drop:
        dictionary.pop(key, None)
