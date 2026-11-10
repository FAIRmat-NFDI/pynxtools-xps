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
"""Every config file shipped with the reader must be valid JSON."""

import json
from pathlib import Path

import pytest

import pynxtools_xps

CONFIG_DIR = Path(pynxtools_xps.__file__).parent / "config"
CONFIG_FILES = sorted(CONFIG_DIR.glob("*.json"))


def test_config_files_are_found():
    assert CONFIG_FILES


@pytest.mark.parametrize("config_file", CONFIG_FILES, ids=lambda path: path.name)
def test_config_file_is_valid_json(config_file):
    with open(config_file, encoding="utf-8") as file:
        json.load(file)
