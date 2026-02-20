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
Metadata mapping for the SPECS XML parser.
"""

import datetime

from pynxtools_xps.mapping import (
    _convert_energy_scan_mode,
    _convert_energy_type,
    _convert_measurement_method,
    _MetadataContext,
    _ValueMap,
    parse_datetime,
)

_POSSIBLE_DATE_FORMATS: list[str] = ["%m/%d/%y %H:%M:%S", "%Y-%m-%d %H:%M:%S"]


def _parse_datetime(date: str) -> str:
    """
    Parse datetime into a datetime.datetime object and return a
    string value in ISO format.

    Parameters
    ----------
    date : str
        String representation of the date, in of these formats:
        "%m/%d/%y %H:%M:%S",
        "%m-%d-%y %H:%M:%S".

    Returns
    -------
    date_object : datetime.datetime
        Datetime in datetime.datetime format.

    """
    # 2025-04-17 13:06:40 UTC
    if date.find("UTC"):
        date = date[: date.find("UTC")].strip()
        tzinfo = datetime.timezone.utc
    else:
        date = date.strip()
        tzinfo = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo  # type: ignore[assignment]

    return parse_datetime(date, _POSSIBLE_DATE_FORMATS, tzinfo)


# TODO: is this correct or a wrong copy from the SLE parser

_KEY_MAP: dict[str, str] = {
    "Group": "group_id",
    "Scan Mode": "scan_mode",
    "Analyzer Lens Voltage": "analyzer_lens_voltage",
    "Calibration File": "calibration_file",
    "Transmission File": "transmission_file",
    "Analyzer Slit": "entrance_slit",
    "Iris Diameter": "iris_diameter",
    "Energy Axis": "x_units",
    "Source": "source_label",
    "Polar Angle": "source_polar_angle",
    "Azimuth Angle": "source_azimuth_angle",
    "ex_energy": "excitation_energy",
    "Acquisition Date": "time_stamp",
    "Analysis Method": "analysis_method",
    "Analyzer": "analyzer_name",
    "Analyzer Lens": "lens_mode",
    "Analyzer Lens Mode": "lens_mode",
    "Scan Variable": "scan_variable",
    "Curves/Scan": "curves_per_scan",
    "Values/Curve": "n_values",
    "Step Size": "step_size",
    "Dwell Time": "dwell_time",
    "Excitation Energy": "excitation_energy",
    "Kinetic Energy": "kinetic_energy",
    "Pass Energy": "pass_energy",
    "Bias Voltage": "bias_voltage_electrons",
    "Binding Energy": "start_energy",
    "Detector Voltage": "detector_voltage",
    "Eff. Workfunction": "work_function",
    "Normalized By": "normalized_by",
    "Comment": "comments",
    "Spectrum ID": "spectrum_id",
    "Note": "note",
}

_VALUE_MAP: _ValueMap = {
    "acquisition_date": _parse_datetime,
    "analysis_method": _convert_measurement_method,
    "scan_mode": _convert_energy_scan_mode,
    "bias_voltage_electrons": float,
    "n_values": int,
    "excitation_energy": float,
    "kinetic_energy": float,
    "work_function": float,
    "dwell_time": float,
    "detector_voltage": float,
    "curves_per_scan": int,
    "pass_energy": float,
    "spectrum_id": int,
    "x_units": _convert_energy_type,
}

# TODO: are these still relevant?
_UNIT_MAP: dict[str, str | None] = {
    "s-1": "1/s",  # workaround for SPECS XY reader
    "norm": None,  # workaround for SPECS XY reader
}

_DEFAULT_UNITS: dict[str, str] = {
    "work_function": "eV",
    "excitation_energy": "eV",
    "pass_energy": "eV",
    "bias_voltage_electrons": "V",
    "dwell_time": "s",
    "step_size": "eV",
}

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
