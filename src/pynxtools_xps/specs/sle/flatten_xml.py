# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:17:23 2024

@author: pielsticker
"""

import copy
from typing import Tuple, Dict, Any, List, Union

from lxml import etree as ET

from pynxtools_xps.reader_utils import (
    convert_pascal_to_snake,
    _re_map_single_value,
    _format_value,
    extract_unit,
)
from pynxtools_xps.value_mappers import (
    convert_units,
    MEASUREMENT_METHOD_MAP,
    convert_measurement_method,
)

from pynxtools_xps.specs.sle.specs_sle_mapping import KEY_MAP, VALUE_MAP


def format_key_and_value(key: str, value_str: str) -> Tuple[Any, str]:
    key = KEY_MAP.get(key, key)
    key = convert_pascal_to_snake(key)
    value, unit = extract_unit(key, value_str)
    value = _format_value(value)
    value = _re_map_single_value(key, value, VALUE_MAP)

    return key, value


def iterate_xml_at_tag(
    xml_elem: ET.Element, tag: str
) -> Dict[str, Union[str, float, int]]:
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
    Dict[str, Union[str, float, int]]
        A dictionary containing formatted attribute values keyed by their corresponding names.
    """

    subelem = xml_elem.find(tag)

    settings = {}

    # print(f"tag: {tag}")
    # print(subelem)

    special_key_map = KEY_MAP.get(tag, {})

    if subelem is not None:
        for param in subelem.iter():
            for key, value in param.attrib.items():
                key = special_key_map.get(key, key)
                key, value = format_key_and_value(key, value)
                settings[key] = value

    return settings


def extract_devices(elem: ET.Element) -> Dict[str, Any]:
    device_settings = {}

    for key, value in elem.attrib.items():
        if key == "Name":
            key = "command"
        key, value = format_key_and_value(key, value)
        device_settings[key] = value

    for param in elem.iter("Parameter"):
        key, value = format_key_and_value(param.attrib["name"], param.text)
        device_settings[key] = value

    return device_settings


def extract_device_commands(elem: ET.Element) -> Dict[str, Any]:
    unique_device_name = elem.attrib["UniqueDeviceName"]

    device_settings = extract_devices(elem)

    return {unique_device_name: device_settings}


def extract_device_info(elem: ET.Element) -> Dict[str, Any]:
    unique_name = elem.attrib["UniqueName"]

    device_settings = extract_devices(elem)

    return {unique_name: device_settings}


def step_profiling(elem: ET.Element) -> Dict[str, Any]:
    settings: Dict[str, Any] = {}

    profiling_settings = iterate_xml_at_tag(xml_elem=elem, tag="ProfilingParams")

    return profiling_settings

    return settings


def _get_group_metadata(spectrum_group: ET.Element) -> Dict[str, Any]:
    """
    Iteratively retrieve metadata for one spectrum group.

    Parameters
    ----------
    spectrum_group: lxml.etree._Element
        XML element containing one spectrum group.

    Returns
    -------
    settings: dict
        Dictionary containing all metadata for
        the spectrum group.

    """
    settings = {}
    settings["group_name"] = spectrum_group.attrib["Name"]
    settings["group_id"] = spectrum_group.attrib["ID"]

    return settings


def _extract_comm_settings(elem: ET.Element) -> Dict[str, Any]:
    """
    Iteratively retrieve metadata for common settings of one spectrum group.

    Parameters
    ----------
    spectrum_group: lxml.etree._Element
        XML element containing common settings for one spectrum group.

    Returns
    -------
    settings: dict
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
        common_setting = iterate_xml_at_tag(
            xml_elem=elem,
            tag=tag,
        )

        if common_setting:
            common_settings.update(common_setting)

    return common_settings


def _get_spectrum_metadata(elem: ET.Element) -> Dict[str, Any]:
    """
    Iteratively retrieve metadata for one spectrum.

    Parameters
    ----------
    spectrum: lxml.etree._Element
        XML element containing one spectrum.

    Returns
    -------
    spectrum_ settings: dict
        Dictionary containing all metadata for
        the spectrum.

    """
    spectrum_settings = {}

    spectrum_settings["spectrum_id"] = elem.attrib["ID"]
    spectrum_settings["spectrum_type"] = elem.attrib["Name"]

    spectrum_types = {
        "FixedEnergiesSettings": "Alignment",
        "FixedAnalyzerTransmissionSettings": "fixed_analyer_transmission",
        "ConstantFinalStateSettings": "constant_final_state",
    }

    for tag, spectrum_type in spectrum_types.items():
        scan_settings = iterate_xml_at_tag(xml_elem=elem, tag=tag)

        if scan_settings:
            spectrum_settings.update(scan_settings)
            spectrum_settings["scan_type"] = spectrum_type

    if (comment := elem.find("Comment")) is not None:
        spectrum_settings.update({"spectrum_comment": comment.text})

    return spectrum_settings


FUNC_MAP = {
    "DeviceCommand": extract_device_commands,
    "DeviceInfo": extract_device_info,
    "StepProfiling": step_profiling,
    "SpectrumGroup": _get_group_metadata,
    "CommonSpectrumSettings": _extract_comm_settings,
    "Spectrum": _get_spectrum_metadata,
}


def process_xml_element(elem: ET.Element, settings: Dict[str, Any]):
    if elem.tag in FUNC_MAP:
        elem_settings = FUNC_MAP[elem.tag](elem)
        settings.update(elem_settings)

    # Recursively process each child element
    for child in elem:
        process_xml_element(child, settings)

    return settings


def flatten_schedule(xml: ET.Element) -> List[Dict[str, Any]]:
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
    collect = []

    for measurement_type in MEASUREMENT_METHOD_MAP:
        for group in xml.iter(measurement_type):
            data: Dict[str, Any] = {}
            data["analysis_method"] = convert_measurement_method(measurement_type)
            process_xml_element(group, data)

            collect += [copy.copy(data)]

    return collect


def flatten_context(xml: ET.Element) -> Dict[str, Any]:
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
    if xml is not None:
        return {}

    device_metadata: Dict[str, Any] = {}

    for elem in xml.iter("DeviceContext"):
        process_xml_element(elem, device_metadata)

    return device_metadata


def flatten_metainfo(xml: ET.Element) -> List[Dict[str, Any]]:
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
    if xml is not None:
        return {}

    metainfo: Dict[str, Any] = {}

    for elem in xml.iter("Parameter"):
        process_xml_element(elem, metainfo)

    return metainfo
