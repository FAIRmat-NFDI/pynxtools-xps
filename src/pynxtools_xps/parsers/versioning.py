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
Deprecated: use ``pynxtools.dataconverter.readers.parsers.versioning``.

This module only re-exports the version helpers from pynxtools core so that
old imports keep working. It will be removed in a future release.
"""

import warnings

from pynxtools.dataconverter.readers.parsers.versioning import (
    VersionRange,
    VersionTuple,
    format_version,
    format_version_ranges,
    is_version_supported,
    normalize_version,
    version_ranges_overlap,
)

warnings.warn(
    "pynxtools_xps.parsers.versioning is deprecated and will be removed in a "
    "future release. Import from pynxtools.dataconverter.readers.parsers.versioning "
    "instead.",
    DeprecationWarning,
    stacklevel=2,
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
