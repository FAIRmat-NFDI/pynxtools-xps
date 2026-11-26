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
"""Version helpers, re-exported from pynxtools core.

They live in ``pynxtools.dataconverter.readers.parsers.versioning``. This
module stays so that existing imports keep working.
"""

from pynxtools.dataconverter.readers.parsers.versioning import (
    VersionRange,
    VersionTuple,
    format_version,
    format_version_ranges,
    is_version_supported,
    normalize_version,
    version_ranges_overlap,
)

__all__ = [
    "VersionRange",
    "VersionTuple",
    "format_version",
    "format_version_ranges",
    "is_version_supported",
    "normalize_version",
    "version_ranges_overlap",
]
