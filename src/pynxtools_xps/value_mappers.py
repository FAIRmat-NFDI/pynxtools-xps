# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:57:41 2024

@author: pielsticker
"""

import re
from typing import Dict, Any

ENERGY_TYPE_MAP = {
    "BE": "binding",
    "KE": "kinetic",
    "Binding": "binding",
    "Binding Energy": "binding",
    "binding energy": "binding",
    "Kinetic": "kinetic",
    "Kinetic Energy": "kinetic",
    "kinetic energy": "kinetic",
}

ENERGY_SCAN_MODE_MAP = {
    "Fixed": "fixed_energy",
    "fixed": "fixed_energy",
    "FixedEnergies": "fixed_energy",
    "Swept": "fixed_analyser_transmission",
    "FixedAnalyzerTransmission": "fixed_analyser_transmission",
    "FAT": "fixed_analyser_transmission",
    "FixedRetardationRatio": "fixed_retardation_ratio",
    "FRR": "fixed_retardation_ratio",
    "Snapshot": "snapshot",
    "SnapshotFAT": "snapshot",
}

MEASUREMENT_METHOD_MAP = {
    "XPS": "X-ray photoelectron spectroscopy (XPS)",
    "UPS": "ultraviolet photoelectron spectroscopy (UPS)",
    "ElectronSpectroscopy": "electron spectroscopy for chemical analysis (ESCA)",
    "NAPXPS": "near ambient pressure X-ray photoelectron spectroscopy (NAPXPS)",
    "ARXPS": "angle-resolved X-ray photoelectron spectroscopy (ARXPS)",
}

BOOL_MAP = {
    "yes": True,
    "Yes": True,
    "no": False,
    "No": False,
    "On": True,
    "Off": False,
}

UNIT_MAP = {
    "Counts": "counts",
    "counts/s": "counts_per_second",
    "CPS": "counts_per_second",
    "u": "um",
    "KV": "kV",
    "seconds": "s",
    "(min)": "min",
    "Percent": "",  # should be changed back to percent once pint is updated
    "atom": "dimensionless",
    "eV/atom": "eV",
    "microm": "micrometer",
    "d": "degree",
}


def _replace_from_map(value: Any, value_map: Dict[str, Any]):
    """
    For a given value, return a new value if the value is
    part of the value_map.
    """
    if value in value_map:
        return value_map[value]
    return value


def convert_energy_type(energy_type: str):
    """
    Change the strings for energy type to the allowed
    values in NXmpes.

    """
    return _replace_from_map(energy_type, ENERGY_TYPE_MAP)


def convert_energy_scan_mode(energy_scan_mode: str):
    """
    Change the strings for energy scan mode to
    the allowed values in NXmpes.

    """
    return _replace_from_map(energy_scan_mode, ENERGY_SCAN_MODE_MAP)


def convert_measurement_method(measurement_method: str):
    """
    Change the strings for measurement method to
    the allowed values in NXmpes.

    """
    return _replace_from_map(measurement_method, MEASUREMENT_METHOD_MAP)


def convert_bool(bool_like: str):
    """Convert "yes", "no" to actual boooleans."""
    return _replace_from_map(bool_like, BOOL_MAP)


def convert_units(units: str):
    """Map y_units to shortened values."""
    return _replace_from_map(units, UNIT_MAP)


def get_units_for_key(unit_key: str, unit_map: Dict[str, str]) -> str:
    """
    Get correct units for a given key from a dictionary with unit map.
    Parameters
    ----------
    unit_key : str
       Key of type <mapping>:<spectrum_key>, e.g.
       detector/detector_voltage
    Returns
    -------
    str
        Unit for that unit_key.
    """
    regex_match = re.search(r"\[([A-Za-z0-9_]+)\]", unit_key)
    if regex_match is None:
        return unit_map.get(unit_key, None)
    return regex_match.group(1)
