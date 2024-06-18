"""
Comment handler for .vms files.
"""
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

from typing import List, Literal, Any, Dict

from pynxtools_xps.vms.casa_data_model import CasaProcess
from pynxtools_xps.phi.spe_pro_phi import PhiParser
from pynxtools_xps.kratos.metadata_kratos import KratosParser

from pynxtools_xps.reader_utils import (
    convert_pascal_to_snake,
)


def handle_comments(
    comment_list: List[str], comment_type: Literal["header", "block"]
) -> Dict[str, Any]:
    comments = {}

    special_keys: List[str] = [
        "Casa Info Follows",
        "SpecsLab Prodigy",
        "SOFH",
        "Creation",
    ]

    for keyword in special_keys:
        # if any(keyword in line for line in comment_list):
        indices = [i for i, line in enumerate(comment_list) if keyword in line]

        if indices:
            index = indices[0]
            if keyword == "Casa Info Follows":
                if comment_type == "header":
                    special_comments = comment_list[index]
                    comment_list = comment_list[index + 1 :]
                if comment_type == "block":
                    special_comments = comment_list[index:]
                    # TODO: only extract casa comments
                comments.update(_handle_casa_comments(special_comments, comment_type))

            if keyword == "SpecsLab Prodigy":
                special_comments = comment_list[index]
                comment_list = comment_list[index + 1 :]
                comments.update(_handle_prodigy_header_comments(special_comments))

            if keyword == "SOFH":
                end_index = [
                    i for i, line in enumerate(comment_list) if "EOFH" in line
                ][0]
                special_comments = comment_list[index : end_index + 1]  # type: ignore[assignment]
                del comment_list[index : end_index + 1]
                comments.update(_handle_phi_header_comments(special_comments))

            if keyword == "Creation":
                end_index = [
                    i for i, line in enumerate(comment_list) if "X-ray Power" in line
                ][0]
                special_comments = comment_list[index : end_index + 1]  # type: ignore[assignment]
                del comment_list[index : end_index + 1]
                comments.update(_handle_kratos_block_comments(special_comments))

    # Handle non-special comments.
    comments.update(_handle_misc_comments(comment_list))

    return comments


def _handle_casa_comments(comment_list: str, comment_type: Literal["header", "block"]):
    """Handle comments from CasaXPS, depending on their location."""
    handle_functions = {
        "header": _handle_casa_header_comments,
        "block": _handle_casa_block_comments,
    }

    return handle_functions[comment_type](comment_list)


def _handle_casa_header_comments(comment_line: str):
    """Get information about CasaXPS version."""
    return {
        "casa_version": comment_line.split("Casa Info Follows CasaXPS Version")[
            1
        ].strip()
    }


def _handle_casa_block_comments(comment_list: List[str]) -> Dict[str, Any]:
    """Get all processing and fitting data from CasaXPS comments."""
    comments = {}

    if "Casa Info Follows" in comment_list[0]:
        #
        casa = CasaProcess()
        casa_data = casa.process_comments(comment_list)

        comments.update(casa_data)

        no_of_casa_lines = 1

        for number in (
            "n_alignments",
            "n_unknown_processes",
            "n_regions",
            "n_components",
        ):
            occurence = getattr(casa, number)
            no_of_casa_lines += 1
            if occurence >= 1:
                no_of_casa_lines += occurence

    return comments


def _handle_prodigy_header_comments(comment_line: str) -> Dict[str, str]:
    """Get information about SpecsLab Prodigy version."""
    return {"prodigy_version": comment_line.split("Version")[1].strip()}


def _handle_phi_header_comments(comment_list: List[str]) -> Dict[str, Any]:
    """Get metadata from Phi system."""
    phi_parser = PhiParser()
    phi_parser.parse_header_into_metadata(comment_list)

    phi_comments = phi_parser.metadata.dict()

    regions = phi_parser.parse_spectral_regions(comment_list)
    areas = phi_parser.parse_spatial_areas(comment_list)

    for region in regions:
        for area in areas:
            concatenated = {**region.dict(), **area.dict()}

        phi_comments.update(concatenated)

    return phi_comments


def _handle_kratos_block_comments(comment_list: List[str]) -> Dict[str, Any]:
    """Get metadata from Kratos system."""
    kratos_parser = KratosParser()
    kratos_parser.parse_header_into_metadata(comment_list)

    kratos_parser.flatten_metadata()

    return kratos_parser.metadata.dict()


def _handle_misc_comments(comment_list: List[str]) -> Dict[str, str]:
    """Handle any other comments."""
    comments = {}
    for line in comment_list:
        for sep in ("=", ":"):
            try:
                key, value = [part.strip(" ") for part in line.split("=", 1)]
                comments[convert_pascal_to_snake(key)] = value
            except ValueError:
                continue
    return comments
