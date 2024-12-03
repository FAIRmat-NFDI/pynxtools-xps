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
Data model for data from Scienta instruments.
"""

from dataclasses import dataclass, field
import numpy as np

from pynxtools_xps.reader_utils import XpsDataclass


@dataclass
class ScientaHeader(XpsDataclass):
    no_of_regions: int = 0
    software_version: str = ""


@dataclass
class ScientaRegion(XpsDataclass):
    region_id: int = 0
    region_name: str = ""
    energy_size: int = 0
    energy_scale: str = ""
    energy_scale_2: str = ""
    energy_axis: np.ndarray = field(default_factory=lambda: np.zeros(0))
    lens_mode: str = ""
    pass_energy: float = 0.0
    no_of_scans: int = 0
    excitation_energy: float = 0.0
    acquisition_mode: str = ""
    energy_units: str = "kinetic"  # ???
    center_energy: float = 0.0
    start_energy: float = 0.0
    stop_energy: float = 0.0
    step_size: float = 0.0
    dwell_time: float = 0.0
    detector_first_x_channel: int = 0
    detector_last_x_channel: int = 0
    detector_first_y_channel: int = 0
    detector_last_y_channel: int = 0
    number_of_slices: str = ""
    data_file: str = ""
    sequence_file: str = ""
    spectrum_type: str = ""
    instrument_name: str = ""
    vendor: str = ""
    user_name: str = ""
    sample_name: str = ""
    spectrum_comment: str = ""
    start_date: str = ""
    start_time: str = ""
    time_stamp: str = ""
    time_per_spectrum_channel: float = 0.0
    detector_mode: str = ""


"""
jsonschema data model for data from the Scienta PEAK software.
"""
scienta_igor_peak_schema = {
    "type": "object",
    "properties": {
        "Version": {"type": ["integer", "null"]},
        "Name": {"type": ["string", "null"]},
        "ElementSet": {"type": ["string", "null"]},
        "ImageSource": {"type": ["string", "null"]},
        "ImageDetector": {"type": ["string", "null"]},
        "LensMode": {"type": ["string", "null"]},
        "AcquisitionMode": {"type": ["string", "null"]},
        "PassEnergy": {"type": ["number", "null"]},
        "Date": {"type": ["string", "null"], "format": "date"},
        "Time": {"type": ["string", "null"], "format": "time"},
        "Duration": {"type": ["number", "null"]},
        "AcquisitionLoopCount": {"type": ["integer", "null"]},
        "Instrument": {"type": ["string", "null"]},
        "AnalyserModel": {"type": ["string", "null"]},
        "Location": {"type": ["string", "null"]},
        "SerialNumber": {"type": ["string", "null"]},
        "AnalyserRadius": {"type": ["number", "null"]},
        "AnalysedParticle": {"type": ["string", "null"]},
        "AnalyserWorkFunction": {"type": ["number", "null"]},
        "SpectrumVersion": {"type": ["integer", "null"]},
        "AcquisitionStatus": {"type": ["string", "null"]},
        "TotalAcquisitionTime": {"type": ["number", "null"]},
        "DwellTime": {"type": ["number", "null"]},
        "IsSweep": {"type": ["boolean", "null"]},
        "EnergyConversion": {
            "type": "object",
            "properties": {
                "UseAnalyserWorkFunction": {"type": ["boolean", "null"]},
                "UseSampleWorkFunction": {"type": ["boolean", "null"]},
                "UseSampleBias": {"type": ["boolean", "null"]},
                "AnalyserMode": {"type": ["string", "null"]},
            },
        },
        "EnergyMode": {"type": ["string", "null"]},
        "ConstantParameters": {"type": ["object", "null"]},
        "Slit": {
            "type": "object",
            "properties": {
                "Id": {"type": ["string", "null"]},
                "Name": {"type": ["string", "null"]},
                "SlitType": {"type": ["string", "null"]},
                "Width": {"type": ["number", "null"]},
                "Length": {"type": ["number", "null"]},
                "Radius": {"type": ["number", "null"]},
                "IsCentered": {"type": ["boolean", "null"]},
                "HasBadQuality": {"type": ["boolean", "null"]},
                "HasAperture": {"type": ["boolean", "null"]},
                "KnobPosition": {"type": ["integer", "null"]},
                "MotorizedPosition": {"type": ["integer", "null"]},
                "MotorizedPositionTolerance": {"type": ["number", "null"]},
            },
        },
        "Sample": {
            "type": "object",
            "properties": {
                "Name": {"type": ["string", "null"]},
                "Description": {"type": ["string", "null"]},
                "Temperature": {"type": ["number", "null"]},
                "Compound": {"type": ["string", "null"]},
                "Environment": {"type": ["string", "null"]},
                "WorkFunction": {"type": ["number", "null"]},
                "Bias": {"type": ["number", "null"]},
            },
        },
        "ExcitationSource": {
            "type": "object",
            "properties": {
                "ExcitationSourceServerInformation": {"type": ["string", "null"]},
                "ExcitationSourcePreset": {
                    "type": "object",
                    "properties": {
                        "Id": {"type": ["string", "null"]},
                        "Name": {"type": ["string", "null"]},
                        "Description": {"type": ["string", "null"]},
                        "TargetEnergy": {"type": ["number", "null"]},
                        "SourceType": {"type": ["string", "null"]},
                    },
                },
                "ExcitationSourceEnergyInformation": {
                    "type": "object",
                    "properties": {
                        "SourceType": {"type": ["string", "null"]},
                        "Energy": {"type": ["number", "null"]},
                        "TargetEnergy": {"type": ["number", "null"]},
                        "EnergyStatus": {"type": ["string", "null"]},
                        "EnergyStatusMessage": {"type": ["string", "null"]},
                    },
                },
                "ExcitationSourceServerStatusInformation": {"type": ["string", "null"]},
            },
        },
        "FixedAxes": {
            "type": "object",
            "patternProperties": {
                "^[XYZ]$": {
                    "type": "object",
                    "properties": {
                        "Center": {"type": ["number", "null"]},
                        "Binning": {"type": ["integer", "null"]},
                        "Name": {"type": ["string", "null"]},
                        "Unit": {"type": ["string", "null"]},
                    },
                }
            },
            "additionalProperties": True,
        },
        "SweepAxes": {
            "type": "object",
            "properties": {
                "X": {
                    "type": "object",
                    "properties": {
                        "Offset": {"type": ["number", "null"]},
                        "Delta": {"type": ["number", "null"]},
                        "Count": {"type": ["integer", "null"]},
                        "AnalyserDelta": {"type": ["number", "null"]},
                        "Name": {"type": ["string", "null"]},
                        "Unit": {"type": ["string", "null"]},
                    },
                }
            },
        },
        "UserName": {"type": ["string", "null"]},
        "StorageName": {"type": ["string", "null"]},
        "SequenceName": {"type": ["string", "null"]},
        "SequenceCounter": {"type": ["integer", "null"]},
        "SpectrumCounter": {"type": ["integer", "null"]},
        "TimePerEnergyChannel": {"type": ["number", "null"]},
    },
}
