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
        "Version": {"type": "integer"},
        "Name": {"type": "string"},
        "ElementSet": {"type": "string"},
        "ImageSource": {"type": "string"},
        "ImageDetector": {"type": "string"},
        "LensMode": {"type": "string"},
        "AcquisitionMode": {"type": "string"},
        "PassEnergy": {"type": "number"},
        "Date": {"type": "string", "format": "date"},
        "Time": {"type": "string", "format": "time"},
        "Duration": {"type": "number"},
        "AcquisitionLoopCount": {"type": "integer"},
        "Instrument": {"type": "string"},
        "AnalyserModel": {"type": "string"},
        "Location": {"type": "string"},
        "SerialNumber": {"type": "string"},
        "AnalyserRadius": {"type": "number"},
        "AnalysedParticle": {"type": "string"},
        "AnalyserWorkFunction": {"type": "number"},
        "SpectrumVersion": {"type": "integer"},
        "AcquisitionStatus": {"type": "string"},
        "TotalAcquisitionTime": {"type": "number"},
        "DwellTime": {"type": "number"},
        "IsSweep": {"type": "boolean"},
        "EnergyConversion": {
            "type": "object",
            "properties": {
                "UseAnalyserWorkFunction": {"type": "boolean"},
                "UseSampleWorkFunction": {"type": "boolean"},
                "UseSampleBias": {"type": "boolean"},
                "AnalyserMode": {"type": "string"},
            },
        },
        "EnergyMode": {"type": "string"},
        "ConstantParameters": {"type": "object"},
        "Slit": {
            "type": "object",
            "properties": {
                "Id": {"type": "string"},
                "Name": {"type": "string"},
                "SlitType": {"type": "string"},
                "Width": {"type": "number"},
                "Length": {"type": "number"},
                "Radius": {"type": ["number", "null"]},
                "IsCentered": {"type": "boolean"},
                "HasBadQuality": {"type": "boolean"},
                "HasAperture": {"type": "boolean"},
                "KnobPosition": {"type": "integer"},
                "MotorizedPosition": {"type": "integer"},
                "MotorizedPositionTolerance": {"type": "number"},
            },
        },
        "Sample": {
            "type": "object",
            "properties": {
                "Name": {"type": "string"},
                "Description": {"type": "string"},
                "Temperature": {"type": ["number", "null"]},
                "Compound": {"type": "string"},
                "Environment": {"type": "string"},
                "WorkFunction": {"type": "number"},
                "Bias": {"type": "number"},
            },
        },
        "ExcitationSource": {
            "type": "object",
            "properties": {
                "ExcitationSourceServerInformation": {"type": ["string", "null"]},
                "ExcitationSourcePreset": {
                    "type": "object",
                    "properties": {
                        "Id": {"type": "string"},
                        "Name": {"type": "string"},
                        "Description": {"type": "string"},
                        "TargetEnergy": {"type": "number"},
                        "SourceType": {"type": "string"},
                    },
                },
                "ExcitationSourceEnergyInformation": {
                    "type": "object",
                    "properties": {
                        "SourceType": {"type": "string"},
                        "Energy": {"type": "number"},
                        "TargetEnergy": {"type": "number"},
                        "EnergyStatus": {"type": "string"},
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
                        "Center": {"type": "number"},
                        "Binning": {"type": "integer"},
                        "Name": {"type": "string"},
                        "Unit": {"type": "string"},
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
                        "Offset": {"type": "number"},
                        "Delta": {"type": "number"},
                        "Count": {"type": "integer"},
                        "AnalyserDelta": {"type": "number"},
                        "Name": {"type": "string"},
                        "Unit": {"type": "string"},
                    },
                }
            },
        },
        "UserName": {"type": "string"},
        "StorageName": {"type": "string"},
        "SequenceName": {"type": "string"},
        "SequenceCounter": {"type": "integer"},
        "SpectrumCounter": {"type": "integer"},
        "TimePerEnergyChannel": {"type": "number"},
    },
}
