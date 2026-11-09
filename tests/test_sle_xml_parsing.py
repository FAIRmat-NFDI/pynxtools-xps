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
"""XML stored in SPECS SLE files is parsed without expanding entities."""

from pynxtools_xps.parsers.specs.sle.parser import _parse_xml

XML_WITH_ENTITY = b"""<?xml version="1.0"?>
<!DOCTYPE root [<!ENTITY greeting "expanded">]>
<root>&greeting;</root>
"""


def test_entities_are_not_expanded():
    root = _parse_xml(XML_WITH_ENTITY)

    assert "expanded" not in "".join(root.itertext())


def test_normal_xml_is_parsed():
    root = _parse_xml(b"<root><child name='a'>text</child></root>")

    assert root.find("child").get("name") == "a"
    assert root.find("child").text == "text"
