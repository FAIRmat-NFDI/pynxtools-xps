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
Flatten the internal XML representation of the experimental schedule
to a list of dictionaries, with each dictionary representing one spectrum.
"""

import copy
from typing import Any

from lxml import etree as ET

from pynxtools_xps.mapping import MEASUREMENT_METHOD_MAP

# TODO: this gets another set of maps from _KEY_MAP, should be done in one go!
from pynxtools_xps.parsers.specs.sle.metadata import _KEY_MAP, _context


def _iterate_xml_at_tag(xml_elem: ET.Element, tag: str) -> dict[str, str | float | int]:
    """
    Iterates through XML elements at the specified tag and formats their attributes.

    Parameters
    ----------
    xml_elem : ET.Element
        The XML element to search within.

    tag : str
        The tag name to find in the XML structure.

    Returns
    -------
    dict[str, Union[str, float, int]]
        A dictionary containing formatted attribute values keyed by their corresponding names.
    """

    sub_elem = xml_elem.find(tag)

    settings: dict[str, Any] = {}

    # TODO: this gets another set of maps from _KEY_MAP, should be done in one go!
    special_key_map: str | dict[str, str] = _KEY_MAP.get(tag, {})

    if sub_elem is not None and isinstance(special_key_map, dict):
        for param in sub_elem.iter():
            for key, value in param.attrib.items():
                key = special_key_map.get(key, key)
                key, value, unit = _context.format(key, value)
                settings[key] = value
                if unit:
                    settings[f"{key}/@units"] = unit

    return settings


def _extract_devices(elem: ET.Element) -> dict[str, Any]:
    """
    Extract all device settings.

    Parameters
    ----------
    elem : ET.Element
        XML element containing one DeviceCommand or DeviceInfo group.

    Returns
    -------
    device_settings: Dict[str, Any]
        Dictionary containing all metadata for
        the device.

    """
    device_settings = {}

    for key, value in elem.attrib.items():
        if key == "Name":
            key = "command"
        key, value, unit = _context.format(key, value)
        device_settings[key] = value

    for param in elem.iter("Parameter"):
        key, value, unit = _context.format(param.attrib["name"], param.text)

        device_settings[key] = value
        if unit:
            device_settings[f"{key}/@units"] = unit

    return device_settings


def _extract_device_commands(elem: ET.Element) -> dict[str, Any]:
    """
    Retrieve device commands.

    Parameters
    ----------
    elem : lxml.etree._Element
         XML element containing one DeviceCommand group.

    Returns
    -------
    device_settings: Dict[str, Any]
        Dictionary containing all information metadata for
        the device.

    """

    unique_device_name = elem.attrib["UniqueDeviceName"]

    device_settings = _extract_devices(elem)

    return {
        f"{unique_device_name}/{key}": value for key, value in device_settings.items()
    }


def _extract_device_info(elem: ET.Element) -> dict[str, Any]:
    """
    Retrieve device information.

    Parameters
    ----------
    elem : lxml.etree._Element
         XML element containing one DeviceInfo group.

    Returns
    -------
    device_settings: Dict[str, Any]
        Dictionary containing all command metadata for
        the device.

    """
    unique_name = elem.attrib["UniqueName"]

    device_settings = _extract_devices(elem)

    return {f"{unique_name}/{key}": value for key, value in device_settings.items()}


def _step_profiling(elem: ET.Element) -> dict[str, Any]:
    """
    Retrieve metadata for one StepProfiling group.

    Parameters
    ----------
    elem : lxml.etree._Element
         XML element containing one StepProfiling group.

    Returns
    -------
    profiling_settings: Dict[str, Any]
        Dictionary containing all metadata for
        the spectrum group.

    """
    profiling_settings = _iterate_xml_at_tag(xml_elem=elem, tag="ProfilingParams")

    return profiling_settings


def _get_group_metadata(spectrum_group: ET.Element) -> dict[str, Any]:
    """
    Retrieve metadata for one spectrum group.

    Parameters
    ----------
    spectrum_group: lxml.etree._Element
        XML element containing one spectrum group.

    Returns
    -------
    settings: Dict[str, Any]
        Dictionary containing all metadata for
        the spectrum group.

    """
    settings = {}
    settings["group_name"] = spectrum_group.attrib["Name"]
    settings["group_id"] = spectrum_group.attrib["ID"]

    return settings


def _extract_comm_settings(elem: ET.Element) -> dict[str, Any]:
    """
    Retrieve metadata for common settings of one spectrum group.

    Parameters
    ----------
    spectrum_group: lxml.etree._Element
        XML element containing common settings for one spectrum group.

    Returns
    -------
    common_settings: Dict[str, Any]
        Dictionary containing all common metadata for
        the spectrum group.

    """
    tags_to_search = [
        "ScanMode",
        "SlitInfo",
        "Lens",
        "EnergyChannelCalibration",
        "Transmission",
        "Iris",
    ]

    common_settings = {}

    for tag in tags_to_search:
        common_setting = _iterate_xml_at_tag(
            xml_elem=elem,
            tag=tag,
        )

        if common_setting:
            common_settings.update(common_setting)

    return common_settings


def _get_spectrum_metadata(elem: ET.Element) -> dict[str, Any]:
    """
    Retrieve metadata for one spectrum.

    Parameters
    ----------
    spectrum: lxml.etree._Element
        XML element containing one spectrum.

    Returns
    -------
    spectrum_ settings: Dict[str, Any]
        Dictionary containing all metadata for
        the spectrum.

    """
    spectrum_settings = {}

    spectrum_settings["spectrum_id"] = elem.attrib["ID"]
    spectrum_settings["spectrum_type"] = elem.attrib["Name"]
    # Workaround for transitions list
    spectrum_settings["transitions"] = [spectrum_settings["spectrum_type"]]

    spectrum_types = {
        "FixedEnergiesSettings": "Alignment",
        "FixedAnalyzerTransmissionSettings": "fixed_analyzer_transmission",
        "ConstantFinalStateSettings": "constant_final_state",
    }

    for tag, spectrum_type in spectrum_types.items():
        scan_settings = _iterate_xml_at_tag(xml_elem=elem, tag=tag)

        if scan_settings:
            spectrum_settings.update(scan_settings)
            spectrum_settings["scan_type"] = spectrum_type

    if (comment := elem.find("Comment")) is not None:
        spectrum_settings.update({"spectrum_comment": comment.text})

    return spectrum_settings


FUNC_MAP = {
    "DeviceCommand": _extract_device_commands,
    "DeviceInfo": _extract_device_info,
    "StepProfiling": _step_profiling,
    "SpectrumGroup": _get_group_metadata,
    "CommonSpectrumSettings": _extract_comm_settings,
    "Spectrum": _get_spectrum_metadata,
}


def _process_xml_element(
    elem: ET.Element, settings: dict[str, Any], collect: list[dict[str, Any]]
):
    if elem.tag in FUNC_MAP:
        elem_settings = FUNC_MAP[elem.tag](elem)
        settings.update(elem_settings)

        if elem.tag == "Spectrum":
            collect.append(copy.copy(settings))
            return

    # Recursively process each child element
    for child in elem:
        _process_xml_element(child, settings, collect)

    return settings


def flatten_schedule(xml: ET.Element) -> list[dict[str, Any]]:
    """
    Flatten the nested XML schedule, keeping only the needed metadata.

    Parameters
    ----------
    xml : lxml.etree
        XML schedule of the experiment.

    Returns
    -------
    collect : list
        List of dictionary with spectra metadata.

    """
    collect: list[dict[str, Any]] = []

    for measurement_type in MEASUREMENT_METHOD_MAP:
        for group in xml.iter(measurement_type):
            data: dict[str, Any] = {}

            try:
                analysis_method, analysis_method_long = _context.map_value(  # type: ignore[misc]
                    "measurement_type", measurement_type
                )
            except ValueError:
                analysis_method = _context.map_value(
                    "measurement_type", measurement_type
                )
                analysis_method_long = "X-ray photoelectron spectroscopy"

            data["analysis_method"] = analysis_method
            data["analysis_method_long_name"] = analysis_method_long

            data["device_group_id"] = group.attrib["ID"]

            _process_xml_element(group, data, collect)

            collect += [copy.copy(data)]

    return collect


def flatten_context(xml: ET.Element) -> dict[str, Any]:
    """
    Flatten the nested XML context, keeping only the needed metadata.

    Parameters
    ----------
    xml : lxml.etree
        XML context of the experiment.

    Returns
    -------
    Dict[str, Any]:
        Dictionary with device metadata.

    """
    if xml is None:
        return {}

    device_metadata: dict[str, Any] = {}

    for elem in xml.iter("DeviceContext"):
        _process_xml_element(elem, device_metadata, [])

    return device_metadata


def flatten_metainfo(xml: ET.Element) -> dict[str, Any]:
    """
    Flatten the nested XML metainfo, keeping only the needed metadata.

    Parameters
    ----------
    xml : lxml.etree
        XML metainfo of the experiment.

    Returns
    -------
    Dict[str, Any]:
        Dictionary with metainfo.

    """
    if xml is None:
        return {}

    metainfo: dict[str, Any] = {}

    for elem in xml.iter("Parameter"):
        _process_xml_element(elem, metainfo, [])

    return metainfo
