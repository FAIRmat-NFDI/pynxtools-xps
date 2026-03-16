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
"""Tests for the NOMAD examples."""

import os

import pytest

try:
    import nomad  # noqa: F401
except ImportError:
    pytest.skip(
        "Skipping NOMAD example tests because nomad-lab is not installed",
        allow_module_level=True,
    )

from pynxtools.testing.nomad_example import (
    example_upload_entry_point_valid,
    get_file_parameter,
    parse_nomad_examples,
)

from pynxtools_xps.nomad.example_uploads import xps_example_upload_entry_point

EXAMPLE_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    "src",
    "pynxtools_xps",
    "nomad",
    "example_uploads",
    "example",
)


@pytest.mark.parametrize(
    "mainfile",
    get_file_parameter(EXAMPLE_PATH),
)
def test_parse_nomad_examples(mainfile):
    """Test if NOMAD examples work."""
    archive_dict = parse_nomad_examples(mainfile)


@pytest.mark.parametrize(
    ("entrypoint", "example_path"),
    [
        pytest.param(
            xps_example_upload_entry_point,
            EXAMPLE_PATH,
            id="xps_example_upload_entry_point",
        ),
    ],
)
def test_example_upload_entry_point_valid(entrypoint, example_path):
    """Test if NOMAD ExampleUploadEntryPoint works."""
    example_upload_entry_point_valid(
        entrypoint=entrypoint,
        example_path=example_path,
    )
