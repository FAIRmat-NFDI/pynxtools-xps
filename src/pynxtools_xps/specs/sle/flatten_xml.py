# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:17:23 2024

@author: pielsticker
"""

import re
import copy
from typing import Tuple, Any

import lxml

from pynxtools_xps.reader_utils import convert_pascal_to_snake

from pynxtools_xps.value_mappers import (
    convert_measurement_method,
    convert_energy_scan_mode,
    MEASUREMENT_METHOD_MAP,
    convert_units,
)


def extract_devices(elem: lxml.etree._Element):
    settings = {}

    for key, val in elem.attrib.items():
        settings[convert_pascal_to_snake(key)] = val

    for param in elem.iter("Parameter"):
        settings[param.attrib["name"]] = param.text

    print(settings)

    return settings

    # data['devices'] += [{'device_type' : j.attrib['DeviceType'],
    #                     'settings':settings}]


# data["devices"] += [device.attrib["DeviceType"]]


def step_profiling(elem: lxml.etree._Element):
    pass


def _get_group_metadata(spectrum_group):
    """
    Iteratively retrieve metadata for one spectrum group.

    Parameters
    ----------
    spectrum_group: xml.etree.ElementTree.Element
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
    for comm_settings in spectrum_group.iter("CommonSpectrumSettings"):
        common_spectrum_settings = _extract_comm_settings(comm_settings)
        settings.update(copy.copy(common_spectrum_settings))

    for spectrum in spectrum_group.iter("Spectrum"):
        spectrum_settings = _get_spectrum_metadata(spectrum)
        settings.update(copy.copy(spectrum_settings))

    return settings


def _extract_comm_settings(comm_settings):
    """
    Iteratively retrieve metadata for common settings of one spectrum group.

    Parameters
    ----------
    spectrum_group: xml.etree.ElementTree.Element
        XML element containing common settings for one spectrum group.

    Returns
    -------
    settings: dict
        Dictionary containing all common metadata for
        the spectrum group.

    """
    common_spectrum_settings = {}
    for setting in comm_settings.iter():
        if setting.tag == "ScanMode":
            energy_scan_mode = convert_energy_scan_mode(setting.attrib["Name"])
            common_spectrum_settings[setting.tag] = energy_scan_mode
        elif setting.tag == "SlitInfo":
            for key, val in setting.attrib.items():
                common_spectrum_settings[key] = val
        elif setting.tag == "Lens":
            voltage_range = setting.attrib["VoltageRange"]
            value, unit = _extract_unit(voltage_range)
            common_spectrum_settings["voltage_range"] = float(value)
            common_spectrum_settings["voltage_range/@units"] = unit
        elif setting.tag == "EnergyChannelCalibration":
            common_spectrum_settings["calibration_file/dir"] = setting.attrib["Dir"]
            common_spectrum_settings["calibration_file/path"] = setting.attrib["File"]
        elif setting.tag == "Transmission":
            common_spectrum_settings["transmission_function/file"] = setting.attrib[
                "File"
            ]
        elif setting.tag == "Iris":
            common_spectrum_settings["iris_diameter"] = setting.attrib["Diameter"]
    return common_spectrum_settings


def _get_spectrum_metadata(spectrum):
    """
    Iteratively retrieve metadata for one spectrum.

    Parameters
    ----------
    spectrum: xml.etree.ElementTree.Element
        XML element containing one spectrum.

    Returns
    -------
    spectrum_ settings: dict
        Dictionary containing all metadata for
        the spectrum.

    """
    spectrum_settings = {}

    spectrum_settings["spectrum_id"] = spectrum.attrib["ID"]
    spectrum_settings["spectrum_type"] = spectrum.attrib["Name"]
    for setting in spectrum.iter("FixedEnergiesSettings"):
        spectrum_settings["dwell_time"] = float(setting.attrib["DwellTime"])
        spectrum_settings["start_energy"] = float(copy.copy(setting.attrib["Ebin"]))
        spectrum_settings["pass_energy"] = float(setting.attrib["Epass"])
        spectrum_settings["lens_mode"] = setting.attrib["LensMode"]
        spectrum_settings["total_scans"] = int(setting.attrib["NumScans"])
        spectrum_settings["n_values"] = int(setting.attrib["NumValues"])
        spectrum_settings["end_energy"] = float(setting.attrib["End"])
        spectrum_settings["excitation_energy"] = float(setting.attrib["Eexc"])
        spectrum_settings["step_size"] = (
            spectrum_settings["start_energy"] - spectrum_settings["end_energy"]
        ) / (spectrum_settings["n_values"] - 1)
    for setting in spectrum.iter("FixedAnalyzerTransmissionSettings"):
        spectrum_settings["dwell_time"] = float(setting.attrib["DwellTime"])
        spectrum_settings["start_energy"] = float(copy.copy(setting.attrib["Ebin"]))
        spectrum_settings["pass_energy"] = float(setting.attrib["Epass"])
        spectrum_settings["lens_mode"] = setting.attrib["LensMode"]
        spectrum_settings["total_scans"] = setting.attrib["NumScans"]
        spectrum_settings["n_values"] = int(setting.attrib["NumValues"])
        spectrum_settings["end_energy"] = float(setting.attrib["End"])
        spectrum_settings["scans"] = int(setting.attrib["NumScans"])
        spectrum_settings["excitation_energy"] = float(setting.attrib["Eexc"])
        spectrum_settings["step_size"] = (
            spectrum_settings["start_energy"] - spectrum_settings["end_energy"]
        ) / (spectrum_settings["n_values"] - 1)
    return spectrum_settings


FUNC_MAP = {
    "DeviceCommand": extract_devices,
    "StepProfiling": step_profiling,
    "CommonSpectrumSettings": _extract_comm_settings,
    "Spectrum": _get_spectrum_metadata,
}


def flatten_xml(xml: lxml.etree):
    """
    Flatten the nested XML structure, keeping only the needed metadata.

    Parameters
    ----------
    xml : xml.etree.ElementTree
        XML schedule of the experiment.

    Returns
    -------
    collect : list
        List of dictionary with spectra metadata.

    """

    def process_element(elem: lxml.etree._Element):
        settings = {}
        # Check if the element's tag is in FUNC_MAP
        if elem.tag in FUNC_MAP:
            # Apply the corresponding function to the element itself
            FUNC_MAP[elem.tag](elem)

        # Recursively process each child element
        for child in elem:
            process_element(child)

    collect = []
    process_element(xml)

    # print(list(xml.iter()))

    for measurement_type in MEASUREMENT_METHOD_MAP:
        for group in xml.iter(measurement_type):
            data = {}
            data["analysis_method"] = convert_measurement_method(measurement_type)

            data["devices"] = []

            for spectrum_group in group.iter("SpectrumGroup"):
                settings = _get_group_metadata(spectrum_group)
                data.update(copy.copy(settings))
                collect += [copy.copy(data)]
    return collect


def _extract_unit(value: str) -> Tuple[Any, str]:
    """
    Extract units for the metadata containing unit information.

    Example:
        analyser_work_function: 4.506eV
        -> analyser_work_function: 4.506,
           analyser_work_function_units: eV,

    Parameters
    ----------
    key : str
        Key of the associated value.
    value : str
        Combined unit and value information.

    Returns
    -------
    value :
        value with units.
    unit : str
        Associated unit.

    """

    pattern = re.compile(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)([a-zA-Z]+)")
    match = pattern.match(value)

    if match:
        value, unit = match.groups()
    else:
        unit = ""

    unit = convert_units(unit)

    # =============================================================================
    #     if key in UNIT_MISSING:
    #         unit = UNIT_MISSING[key]
    # =============================================================================

    return value, unit
