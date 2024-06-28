# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:54:24 2024

@author: pielsticker
"""

import re
import datetime
from typing import Dict, Optional
import numpy as np

from pynxtools_xps.reader_utils import (
    convert_pascal_to_snake,
    _re_map_single_value,
)

from pynxtools_xps.value_mappers import (
    convert_energy_type,
    convert_energy_scan_mode,
)


def _extract_energy_units(energy_units: str) -> str:
    """
    Extract energy units from the strings for energy_units.
    Binding Energy [eV] -> eV

    """
    units = re.search(r"\[(.*?)\]", energy_units)
    return units.group(1) if units is not None else "eV"


def _separate_dimension_scale(scale: str) -> np.typing.ArrayLike:
    """
    Seperate the str of the dimension scale into a numpy array

    Parameters
    ----------
    scale : str
        Str of the form "600 599.5 599 598.5 5".

    Returns
    -------
    np.ndarray
        Dimension scale as a numpy array.

    """
    return np.array([float(s) for s in scale.split(" ")])


def _construct_date_time(date_string: str, time_string: str) -> Optional[str]:
    """
    Convert the native time format to the datetime string
    in the ISO 8601 format: '%Y-%b-%dT%H:%M:%S.%fZ'.

    """

    def _parse_date(date_string: str) -> datetime.datetime:
        possible_date_formats = ["%Y-%m-%d", "%m/%d/%Y"]
        for date_fmt in possible_date_formats:
            try:
                return datetime.datetime.strptime(date_string, date_fmt)
            except ValueError as err:
                continue
        raise ValueError("Date format not recognized")

    def _parse_time(time_string: str) -> datetime.time:
        possible_time_formats = ["%H:%M:%S", "%I:%M:%S %p"]
        for time_fmt in possible_time_formats:
            try:
                return datetime.datetime.strptime(time_string, time_fmt).time()
            except ValueError as err:
                continue
        raise ValueError("Time format not recognized")

    try:
        date = _parse_date(date_string)
        time = _parse_time(time_string)
        datetime_obj = datetime.datetime.combine(date, time)

        return datetime_obj.astimezone().isoformat()

    except ValueError as err:
        raise ValueError(
            "Date and time could not be converted to ISO 8601 format."
        ) from err


KEY_MAP: Dict[str, str] = {
    "number_of_regions": "no_of_regions",
    "version": "software_version",
    "dimension_name": "energy_units",
    "dimension_size": "energy_size",
    "dimension_scale": "energy_axis",
    "number_of_sweeps": "no_of_scans",
    "energy_unit": "energy_scale_2",
    "low_energy": "start_energy",
    "high_energy": "stop_energy",
    "energy_step": "step_size",
    "step_time": "dwell_time",
    "detector_first_x-_channel": "detector_first_x_channel",
    "detector_last_x-_channel": "detector_last_x_channel",
    "detector_first_y-_channel": "detector_first_y_channel",
    "detector_last_y-_channel": "detector_last_y_channel",
    "file": "data_file",
    "sequence": "sequence_file",
    "spectrum_name": "spectrum_type",
    "instrument": "instrument_name",
    "location": "vendor",
    "user": "user_name",
    "sample": "sample_name",
    "comments": "spectrum_comment",
    "date": "start_date",
    "time": "start_time",
}

VALUE_MAP = {
    "no_of_regions": int,
    "energy_size": int,
    "pass_energy": float,
    "no_of_scans": int,
    "excitation_energy": float,
    "center_energy": float,
    "start_energy": float,
    "stop_energy": float,
    "step_size": float,
    "dwell_time": float,
    "detector_first_x_channel": int,
    "detector_last_x_channel": int,
    "detector_first_y_channel": int,
    "detector_last_y_channel": int,
    "time_per_spectrum_channel": float,
    "energy_units": _extract_energy_units,
    "energy_axis": _separate_dimension_scale,
    "energy_scale": convert_energy_type,
    "energy_scale_2": convert_energy_type,
    "acquisition_mode": convert_energy_scan_mode,
    "time_per_spectrum_channel": float,
    "manipulator_r1": float,
    "manipulator_r2": float,
}

UNITS: dict = {
    "energydispersion/pass_energy": "eV",
    "beam_xray/excitation_energy": "eV",
    "region/energy_axis": "eV",
    "region/center_energy": "eV",
    "region/start_energy": "eV",
    "region/stop_energy": "eV",
    "region/step_size": "eV",
    "detector/dwell_time": "eV",
    "region/time_per_spectrum_channel": "s",
}


def _get_key_value_pair(line: str):
    """
    Split the line at the '=' sign and return a
    key-value pair. The values are mapped according
    to the desired format.

    Parameters
    ----------
    line : str
        One line from the input file.

    Returns
    -------
    Tuple[str, object]
        A tuple containing:
        - key : str
            Anything before the '=' sign, mapped to the desired
            key format.
        - value : object
            Anything after the '=' sign, mapped to the desired
            value format and type.

    """
    try:
        key, value = line.split("=")
        key = convert_pascal_to_snake(key)
        key = KEY_MAP.get(key, key)
        if "dimension" in key:
            key_part = f"dimension_{key.rsplit('_')[-1]}"
            key = KEY_MAP.get(key_part, key_part)
        value = _re_map_single_value(key, value, VALUE_MAP)

    except ValueError:
        key, value = "", ""

    return key, value
