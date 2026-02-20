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
Regenerate reference JSON files for testing the VAMAS comment extraction.
"""

import os
import json
import numpy as np
from pynxtools_xps.parsers.vms.comment_handler import handle_comments


def read_comment_file(filepath: str):
    """Read comments from one vms comment test file."""

    no_of_comments = 0
    comment_lines = []

    with open(filepath, "rb") as vms_file:
        for i, line in enumerate(vms_file):
            if i == 0:
                no_of_comments = int(line)
            else:
                comment_lines += [line.decode("utf-8", errors="ignore").strip()]

    return no_of_comments, comment_lines


def generate_ref_comment_jsons():
    """This generates reference json files for the vms comment handler."""
    files = [
        ("kratos.vms", "block"),
        ("phi.vms", "block"),
        ("casa_header.vms", "header"),
        ("casa_process.vms", "block"),
        ("specs_header.vms", "header"),
        ("specs_block.vms", "block"),
    ]

    for file, comment_type in files:
        filepath = os.path.join(
            os.path.dirname(__file__), "..", "tests", "data", "vms_comments", file
        )
        ref_json_filepath = filepath.replace(".vms", "_ref.json")

        _, comment_lines = read_comment_file(filepath)

        comments = handle_comments(comment_lines, comment_type=comment_type)

        if file == "casa_process.vms":
            casa_process = comments["casa"]
            comments = casa_process.flatten_metadata()

        for key, val in comments.items():
            if isinstance(val, np.ndarray):
                comments[key] = val.tolist()

        with open(ref_json_filepath, "w") as json_file:
            json.dump(comments, json_file, indent=4)


if __name__ == "__main__":
    generate_ref_comment_jsons()
