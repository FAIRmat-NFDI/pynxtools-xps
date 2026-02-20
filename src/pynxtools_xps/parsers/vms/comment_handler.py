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
# pylint: disable=too-many-instance-attributes
"""
Comment handler for .vms files.
"""

import re
from typing import Any, Literal

from pynxtools_xps.parsers.kratos import KratosParser
from pynxtools_xps.parsers.phi import PHIParser
from pynxtools_xps.parsers.vms.data_model_casa import CasaProcess
from pynxtools_xps.parsers.vms.metadata import _context


def handle_comments(
    comment_list: list[str],
    comment_type: Literal["header", "block"],
) -> dict[str, Any]:
    comments: dict[str, Any] = {}

    def find_index(lines: list[str], marker: str) -> int | None:
        return next((i for i, line in enumerate(lines) if marker in line), None)

    def extract_block(start: int, end: int) -> list[str]:
        block = comment_list[start : end + 1]
        del comment_list[start : end + 1]
        return block

    special_keys: tuple[str, ...] = (
        "Casa Info Follows",
        "SpecsLab Prodigy",
        "SOFH",
        "Creation",
    )

    for keyword in special_keys:
        start_index = find_index(comment_list, keyword)
        if start_index is None:
            continue

        # --- CASA ------------------------------------------------------------
        if keyword == "Casa Info Follows":
            if comment_type == "header":
                special_comments = [comment_list[start_index]]
                del comment_list[start_index : start_index + 1]
            else:
                special_comments = comment_list[start_index:]
                del comment_list[start_index:]

            casa_comments, n_lines = _handle_casa_comments(
                special_comments, comment_type
            )
            comments.update(casa_comments)
            continue

        # --- SpecsLab Prodigy --------------------------------------------------------
        if keyword == "SpecsLab Prodigy":
            special_comments = [comment_list[start_index]]
            del comment_list[start_index : start_index + 1]
            comments.update(_handle_prodigy_header_comments(special_comments))
            continue

        # --- PHI HEADER ------------------------------------------------------
        if keyword == "SOFH":
            end_index = find_index(comment_list, "EOFH")
            if end_index is None:
                raise ValueError("'EOFH' not found for Phi Metadata block.")

            block = extract_block(start_index, end_index)
            comments.update(_handle_phi_header_comments(block))
            continue

        # --- KRATOS BLOCK ----------------------------------------------------
        if keyword == "Creation":
            end_index = find_index(comment_list, "Quality")
            if end_index is None:
                end_index = find_index(comment_list, "X-ray Power")

            if end_index is None:
                raise ValueError(
                    "Neither 'Quality' nor 'X-ray Power' found in Kratos Metadata block."
                )

            block = extract_block(start_index, end_index)
            comments.update(_handle_kratos_block_comments(block))
            continue

    # Remaining comments
    comments.update(_handle_misc_comments(comment_list))
    return comments


def _handle_casa_comments(
    comment_list: list[str], comment_type: Literal["header", "block"]
) -> tuple[dict[str, Any], int]:
    """Handle comments from CasaXPS, depending on their location."""
    handle_functions = {
        "header": _handle_casa_header_comments,
        "block": _handle_casa_block_comments,
    }

    return handle_functions[comment_type](comment_list)


def _handle_casa_header_comments(
    comment_list: list[str],
) -> tuple[dict[str, Any], int]:
    """Get information about CasaXPS version."""
    comment_line = comment_list[0]

    comments = {
        "casa_version": comment_line.split("Casa Info Follows CasaXPS Version")[
            1
        ].strip()
    }
    no_of_casa_lines = 1

    return comments, no_of_casa_lines


def _handle_casa_block_comments(
    comment_list: list[str],
) -> tuple[dict[str, Any], int]:
    """Get all processing and fitting data from CasaXPS comments."""
    casa = CasaProcess()
    casa.process_comments(comment_list)

    comments = {"casa": casa}

    return comments, casa.no_of_casa_lines


def _handle_prodigy_header_comments(comment_list: list[str]) -> dict[str, str]:
    """Get information about SpecsLab Prodigy version."""
    comment_line = comment_list[0]
    return {"prodigy_version": comment_line.split("Version")[1].strip()}


def _handle_phi_header_comments(comment_list: list[str]) -> dict[str, Any]:
    """Get metadata from Phi system."""
    phi_parser = PHIParser()
    phi_parser.parse_header_into_metadata(comment_list)

    phi_comments = phi_parser.metadata.dict()

    regions = phi_parser.parse_spectral_regions(comment_list)
    areas = phi_parser.parse_spatial_areas(comment_list)

    for region in regions:
        for area in areas:
            concatenated = {**region.dict(), **area.dict()}

        phi_comments.update(concatenated)

    return phi_comments


def _handle_kratos_block_comments(comment_list: list[str]) -> dict[str, Any]:
    """Get metadata from Kratos system."""
    kratos_parser = KratosParser()
    kratos_parser.parse_header_into_metadata(comment_list)

    return kratos_parser.flatten_metadata()


def _handle_misc_comments(
    comment_list: list[str],
) -> dict[str, int | float | str | Any]:
    """Handle any other comments."""
    comments = {}
    for line in comment_list:
        parts = re.split(r"[=:]", line, maxsplit=1)
        if len(parts) != 2:
            continue
        key, value = (p.strip() for p in parts)
        key, value, unit = _context.format(key, value)

        comments[key] = value
        if unit:
            comments[f"{key}/@units"] = unit
    return comments
