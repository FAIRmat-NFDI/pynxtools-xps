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
Metadata mapping for the PHI parser.
"""

import datetime
import re
from zoneinfo import ZoneInfo

from pynxtools_xps.mapping import (
    _convert_bool,
    _convert_energy_scan_mode,
    _convert_measurement_method,
    _MetadataContext,
    _ValueMap,
)


def _map_file_type(value: str):
    """Map file_type to easily understandable values."""
    value = value.strip()
    file_type_map = {
        "SPECTRUM": "single_spectrum",
        "DEPTHPRO": "depth_profile",
    }

    if value in file_type_map:
        return file_type_map[value]
    return value


def _map_to_list(value: str):
    """Map all items in value to a list."""
    try:
        sep = ", " if ", " in value else " "
        values = [float(val) for val in value.split(sep)]
    except ValueError:
        sep = ","
        values = [float(val) for val in value.split(sep)]

    return values


def _map_to_xy(value: str):
    """Map items in value to a dictionary with keys x and y."""
    x, y = value.split(" ")

    return {"x": float(x), "y": float(y)}


def _map_to_xy_with_units(value: str):
    """Map items in value to a dictionary with keys x ."""
    x, y, unit = value.split(" ")
    return {"x": float(x), "x_units": unit, "y": float(y), "y_units": unit}


def _convert_energy_referencing(value: str):
    """Map all items in energy_referencing to a dictionary."""
    peak, energy = value.split(" ")
    return {"peak": peak, "energy": float(energy), "energy_units": "eV"}


def _convert_channel_info(value: str):
    """Split channel information into list."""
    channel_number, setting_a, setting_b = value.split(" ")
    return [int(setting_a), float(setting_b)]


def _convert_xray_source_params(value: str):
    """Map all items in xray_source_params to a dictionary."""
    label, energy, mono = value.split(" ")

    return {
        "anode_material": label,
        "energy": float(energy),
        "energy_units": "eV",
        "monochromatized": mono,
    }


def _convert_xray_source_settings(value: str):
    """Map all items in xray_source_settings to a dictionary."""
    (xray_settings) = re.split(r"(\d+)", value)

    # TODO: This should be done through the _context
    # # for i, setting in enumerate(xray_settings):
    #     xray_settings[i] = convert_units(setting)

    return {
        "spot_size": float(xray_settings[1]),
        "spot_size_units": xray_settings[2],
        "power": float(xray_settings[3]),
        "power_units": xray_settings[4],
        "high_voltage": float(xray_settings[5]),
        "high_voltage_units": xray_settings[6],
    }


def _convert_stage_positions(value: str):
    """Map all items in stage_positions to a dictionary."""
    x, y, z, azimuth, polar = value.split(" ")

    return {
        "x": float(x),
        "x_units": "mm",
        "y": float(y),
        "y_units": "mm",
        "z": float(z),
        "z_units": "mm",
        "azimuth": float(azimuth),
        "azimuth_units": "degree",
        "polar": float(polar),
        "polar_units": "degree",
    }


# ToDO: this should be replaced by generic parser from mappers
def _parse_datetime(value: str):
    """
    Parse datetime into a datetime.datetime object.

    Parameters
    ----------
    value : str
        String representation of the date in the format
        "%m/%d/%y.

    Returns
    -------
    date_object : str
        Datetime in ISO8601 format.

    """
    year, month, day = value.strip().split(" ")
    date_object = datetime.datetime(
        year=int(year),
        month=int(month),
        day=int(day),
        tzinfo=ZoneInfo("UTC"),
    )

    return date_object.isoformat()


_KEY_MAP: dict[str, str] = {
    "FileDesc": "file_description",
    "acq_filename": "acquisition_filename",
    "acq_file_date": "acquisition_file_date",
    "institution": "vendor",
    "operator": "user_name",
    "experiment_i_d": "experiment_id",
    "analyzer_work_fcn": "analyzer_work_function",
    "analyzer_retard_gain": "analyzer_retardation_gain",
    "reg_image_interval": "register_image_interval",
    "reg_image_mode": "register_image_mode",
    "reg_image_last": "register_image_last",
    "platen_i_d": "platen_id",
    "s_x_i_filename": "sxi_filename",
    "intensity_recal": "intensity_recalibration",
    "intensity_cal_coeff": "intensity_calibration_coefficients",
    "energy_recal": "energy_recalibration",
    "s_c_a_multiplier_voltage": "sca_multiplier_voltage",
    "C60IonGun": "c60_ion_gun",
    "t_f_c_parameters": "tfc_parameters",
    "image_size_xy": "image_size",
    "float_volt": "float_voltage",
    "float_enable": "float_enabled",
    "grid_volt": "grid_voltage",
    "condensor_volt": "condenser_lens_voltage",
    "objective_volt": "objective_lens_voltage",
    "bend_volt": "bend_voltage",
    "neutral_float_volt": "neutral_float_voltage",
    "neutral_float_enable": "neutral_float_enabled",
    "neutral_grid_volt": "neutral_grid_voltage",
    "neutral_condensor_volt": "neutral_condenser_lens_voltage",
    "neutral_objective_volt": "neutral_objective_lens_voltage",
    "neutral_bend_volt": "neutral_bend_voltage",
    "no_d_p_data_cyc": "no_depth_profile_cycles",
    "no_pre_sputter_cyc": "no_pre_sputter_cycles",
    "sample_rotation": "profiling_sample_rotation",
    "depth_recal": "profiling_depth_recalibration",
    "sputter_mode": "profiling_sputter_mode",
    "no_depth_reg": "profiling_no_depth_regions",
    "depth_cal_def": "depth_calibration_definition",
    "analyser_mode": "energy_scan_mode",
    "surv_num_cycles": "survey_num_of_cycles",
    "surv_time_per_step": "survey_dwell_time",
    "no_spectral_reg_full": "no_spectral_regions_full",
    "no_spectral_region": "no_spectral_regions",
    "spectral_reg_def_full": "spectral_region_definition_full",
    "spectral_reg_def2_full": "spectral_region_definition2_full",
    "spectral_reg_background_full": "spectral_region_background_full",
    "spectral_reg_hero_full": "spectral_region_hero_full",
    "spectral_reg_i_r_full": "spectral_region_ir_full",
    "no_spectral_reg": "no_spectral_regions",
    "spectral_reg_def": "spectral_region_definition",
    "spectral_reg_def2": "spectral_region_definition2",
    "spectral_reg_background": "spectral_region_background",
    "spectral_reg_hero": "spectral_region_hero",
    "spectral_reg_i_r": "spectral_region_ir",
    "no_spatial_area": "no_spatial_areas",
    "spatial_area_def": "spatial_area_definition",
    "spatial_h_r_photo_cor": "spatial_hr_photo_correction",
    "xray_offset_in_um": "xray_offset",
    "xray_mag_factor": "xray_magnification_factor",
    "xray_rotation_in_deg": "xray_rotation",
    "xray_setting": "xray_settings",
    "neutralizer_current": "flood_gun_current",
    "neutralizer_energy": "flood_gun_energy",
    "flood_gun_filament": "flood_gun_filament_current",
    "flood_gun_extractor": "flood_gun_extractor_voltage",
    "detector _acq _time": "detector_acquisition_time",
    "number _of _channels": "number_of_channels",
    "stage_position": "stage_positions",
    "defect_positioner_i_d": "defect_positioner_id",
    "defect_positioner_aligment": "defect_positioner_alignment",
    "gcib_beam": "gcib_high_voltage",
    "gcib_wien": "gcib_wien_filter_voltage",
    "gcib_bend": "gcib_bend_voltage",
    "gcib_magnet": "gcib_magnet_current",
    "auto_e_gun_neut": "auto_flood_gun",
    "auto_ion_neut": "auto_neutral_ion_source",
}

_VALUE_MAP: _ValueMap = {
    "technique": _convert_measurement_method,
    "technique_ex": _convert_measurement_method,
    "file_type": _map_file_type,
    "file_date": _parse_datetime,
    "acquisition_file_date": _parse_datetime,
    "energy_reference": _convert_energy_referencing,
    "intensity_calibration_coefficients": _map_to_list,
    "energy_recalibration": _convert_bool,
    "scan_deflection_span": _map_to_xy,
    "scan_deflection_offset": _map_to_xy,
    "tfc_parameters": _map_to_list,
    "image_size": _map_to_xy,
    "float_enabled": _convert_bool,
    "sputter_raster": _map_to_xy_with_units,
    "sputter_raster_offset": _map_to_xy_with_units,
    "neutral_raster": _map_to_xy_with_units,
    "neutral_raster_offset": _map_to_xy_with_units,
    "profiling_xray_off_during_sputter": _convert_bool,
    "profiling_source_blank_during_sputter": _convert_bool,
    "profiling_depth_recalibration": _convert_bool,
    "energy_scan_mode": _convert_energy_scan_mode,
    "xray_source": _convert_xray_source_params,
    "xray_stigmator": _map_to_xy,
    "xray_offset": _map_to_xy,
    "xray_magnification_factor": _map_to_xy,
    "xray_delay_factor": _map_to_xy,
    "xray_high_power": _convert_bool,
    "xray_emission_control": _convert_bool,
    "xray_settings": _convert_xray_source_settings,
    "sxi_auto_contrast": _convert_bool,
    "sxi_shutter_bias": _convert_bool,
    "stage_positions": _convert_stage_positions,
    "gcib_raster_size": _map_to_xy_with_units,
    "gcib_raster_offset": _map_to_xy_with_units,
    "auto_flood_gun": _convert_bool,
    "auto_neutral_ion_source": _convert_bool,
    "presputter": _convert_bool,
}

_UNIT_MAP: dict[str, str | None] = {}

_DEFAULT_UNITS: dict[str, str] = {
    "intensity": "counts_per_second",
    "grid_voltage": "V",
    "condenser_lens_voltage": "V",
    "objective_lens_voltage": "V",
    "bend_voltage": "V",
    "neutral_grid_voltage": "V",
    "neutral_condenser_lens_voltage": "V",
    "neutral_objective_lens_voltage": "V",
    "stage_current_rotation_speed": "degree/s",
    "neutral_bend_voltage": "V",
    "defect_positioner_u": "mm",
    "defect_positioner_v": "mm",
    "defect_positioner_x": "mm",
    "defect_positioner_y": "mm",
    "defect_positioner_z": "mm",
    "defect_positioner_tilt": "degree",
    "defect_positioner_rotation": "degree",
    "flood_gun_pulse_frequency": "1/s",
    "flood_gun_x_steering": "mm",
    "flood_gun_y_steering": "mm",
    "profiling_sputter_delay": "s",
    "scan_deflection_span_x": "mm",
    "scan_deflection_span_y": "mm",
    "scan_deflection_offset_x": "mm",
    "scan_deflection_offset_y": "mm",
    "survey_dwell_time": "s",
    "xray_offset_x": "um",
    "xray_offset_y": "um",
    "xray_rotation": "degree",
    "xray_stigmator_x": "mm",
    "xray_stigmator_y": "mm",
}

# TODO: this should be automatic?!
# _KEYS_WITH_UNITS: list[str] = [
#     "analyzer_work_function",
#     "source_analyzer_angle",
#     "analyzer_solid_angle",
#     "scan_deflection_span_x",
#     "scan_deflection_span_y",
#     "scan_deflection_offset_x",
#     "scan_deflection_offset_y",
#     "sca_multiplier_voltage",
#     "delay_before_acquire",
#     "sputter_current",
#     "sputter_rate",
#     "sputter_energy",
#     "float_voltage",
#     "target_sputter_time",
#     "sputter_emission",
#     "grid_voltage",
#     "condenser_lens_voltage",
#     "objective_lens_voltage",
#     "bend_voltage",
#     "deflection_bias",
#     "ion_gun_gas_pressure",
#     "sputter_emission",
#     "deflection_bias",
#     "neutral_current",
#     "neutral_rate",
#     "neutral_energy",
#     "neutral_float_voltage",
#     "neutral_grid_voltage",
#     "neutral_condenser_lens_voltage",
#     "neutral_objective_lens_voltage",
#     "neutral_bend_voltage",
#     "neutral_target_timed_on_time",
#     "neutral_emission",
#     "neutral_deflection_bias",
#     "neutral_ion_gun_gas_pressure",
#     "profiling_sputter_delay",
#     "survey_dwell_time",
#     "xray_anode_power",
#     "xray_power",
#     "xray_beam_voltage",
#     "xray_beam_diameter",
#     "xray_condenser_lens_voltage",
#     "xray_objective_coil_current",
#     "xray_blanking_voltage",
#     "xray_filament_current",
#     "xray_rotation",
#     "xray_emission_current",
#     "xray_max_filament_current",
#     "xray_stigmator_x",
#     "xray_stigmator_y",
#     "xray_offset_x",
#     "xray_offset_y",
#     "flood_gun_current",
#     "flood_gun_energy",
#     "flood_gun_extractor_voltage",
#     "flood_gun_filament_current",
#     "flood_gun_pulse_length",
#     "flood_gun_pulse_frequency",
#     "flood_gun_time_per_step",
#     "flood_gun_ramp_rate",
#     "flood_gun_x_steering",
#     "flood_gun_y_steering",
#     "sxi_binding_energy",
#     "sxi_pass_energy",
#     "sxi_lens2_voltage",
#     "sxi_lens3_voltage",
#     "sxi_lens4_voltage",
#     "sxi_lens5_voltage",
#     "sxi_rotator",
#     "sxi_lens_bias_voltage",
#     "sxi_shutter_bias_voltage",
#     "detector_acquisition_time",
#     "stage_current_rotation_speed",
#     "defect_positioner_u",
#     "defect_positioner_v",
#     "defect_positioner_x",
#     "defect_positioner_y",
#     "defect_positioner_z",
#     "defect_positioner_tilt",
#     "defect_positioner_rotation",
#     "gcib_sputter_rate",
#     "gcib_high_voltage",
#     "gcib_ionization",
#     "gcib_extractor",
#     "gcib_wien_filter_voltage",
#     "gcib_bend_voltage",
#     "gcib_emission",
#     "gcib_magnet_current",
#     "gcib_focus",
#     "gcib_objective",
#     "gcib_focus",
#     "gcib_gas_pressure",
#     "gcib_cluster_size",
#     "gcib_energy_per_atom",
#     "deconvolution_pass_energy",
# ]

_context = _MetadataContext(
    key_map=_KEY_MAP,
    value_map=_VALUE_MAP,
    unit_map=_UNIT_MAP,
    default_units=_DEFAULT_UNITS,
)
