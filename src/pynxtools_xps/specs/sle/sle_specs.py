"""
Parser for reading XPS (X-ray Photoelectron Spectroscopy) data from native
Specs Lab Prodigy SLE format, to be passed to mpes nxdl
(NeXus Definition Language) template.
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
# pylint: disable=too-many-lines,too-many-instance-attributes

"""
Parser for reading XPS (X-ray Photoelectron Spectroscopy) metadata from
SPECS Lab Prodigy, to be passed to MPES nxdl (NeXus Definition Language)
template.
"""

import copy
import logging
from typing import Dict, Any, List
import warnings
from datetime import datetime
import sqlite3
import struct
from lxml import etree as ET
from abc import ABC
import numpy as np
import xarray as xr

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_data_key,
    construct_entry_name,
)
from pynxtools_xps.value_mappers import (
    convert_energy_scan_mode,
    convert_energy_type,
    convert_units,
    get_units_for_key,
    MEASUREMENT_METHOD_MAP,
    parse_datetime,
)
from pynxtools_xps.specs.sle.flatten_xml import flatten_xml

logger = logging.getLogger(__name__)

UNITS: dict[str, str] = {
    "work_function": "eV",
    "excitation_energy": "eV",
    "iris_diameter": "mm",
    "step_size": "eV",
    "detector_voltage": "V",
    "dwell_time": "s",
    "raw_data/raw": "counts_per_second ",
    "polar_angle": "degree ",
    "azimuth_angle": "degree",
    "pass_energy": "eV",
    "start_energy": "eV",
    "emission_current": "A",
    "source_voltage": "V",
    "energy_calibration/energy": "eV",
    "transmission_function/kinetic_energy": "eV",
    "transmission_correction/transmission_function/kinetic_energy": "eV",
}


class SleMapperSpecs(XPSMapper):
    """
    Class for restructuring .sle data file from
    specs vendor into python dictionary.
    """

    config_file = "config_specs_sle.json"

    def __init__(self):
        self.parsers = [
            SleProdigyParserV1,
            SleProdigyParserV4,
        ]

        self.versions_map = {}
        for parser in self.parsers:
            supported_versions = parser.supported_versions
            for version in supported_versions:
                self.versions_map[version] = parser

        self.sql_connection = None

        super().__init__()

    def _select_parser(self):
        """
        Select the correct parser for the SLE file version.

        Returns
        -------
        Parser
            Specs SLE Parser.

        """
        version = self._get_sle_version()
        try:
            return self.versions_map[version]()
        except KeyError as exc:
            raise KeyError(
                f"Version f{version} of SPECS Prodigy is currently not supported."
            ) from exc

    def _get_sle_version(self) -> str:
        """Get the Prodigy SLE version from the file."""
        con = sqlite3.connect(self.sql_connection)
        cur = con.cursor()
        query = 'SELECT Value FROM Configuration WHERE Key=="Version"'
        cur.execute(query)
        version = cur.fetchall()[0][0]
        version = version.split(".")
        version = version[0] + "." + version[1].split("-")[0]
        return version

    def parse_file(self, file: str, **kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse the file using the parser that fits the Prodigy SLE version.

        Returns flat list of dictionaries containing one spectrum each.

        Parameters
        ----------
        file : str
            String name of the file.
        **kwargs : Dict[str, Any]
            Dict with additional keyword arguments.

        Returns
        -------
        Dict[str, Any]
            Dict with parsed data.

        """
        self.sql_connection = file
        return super().parse_file(file, **kwargs)

    def construct_data(self):
        """
        Map SLE format to NXmpes-ready dict.

        Returns
        -------
        None.

        """
        # pylint: disable=duplicate-code
        spectra = copy.deepcopy(self.raw_data)

        self._xps_dict["data"]: dict = {}

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        Parameters
        ----------
        spectrum : Dict[str, Any]
            Dictionary with data and metadata for one spectrum.
        template_key_map : Dict[str, str]
            Mapping to NXmpes terms.

        Returns
        -------
        None.

        """
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []
        for part in ["group_name", "spectrum_type"]:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            if key.startswith("entry"):
                entry_parent = "/ENTRY[entry]"
                key = key.replace("entry/", "", 1)
            mpes_key = f"{entry_parent}/{key}"
            if "units" in key:
                value = convert_units(value)
            self._xps_dict[mpes_key] = value
            units = convert_units(get_units_for_key(key, UNITS))
            if units is not None:
                self._xps_dict[f"{mpes_key}/@units"] = units

        analyzer_path = f"{entry_parent}/instrument/electronanalyzer"

        self._xps_dict[f"{entry_parent}/electronanalyzer/name"] = spectrum["devices"][0]
        self._xps_dict[f"{entry_parent}/'source/name"] = spectrum["devices"][1]

        # Create keys for writing to data
        scan_key = construct_data_key(spectrum)

        energy = np.array(spectrum["data"]["x"])

        # If multiple spectra exist to entry, only create a new
        # xr.Dataset if the entry occurs for the first time.
        if entry not in self._xps_dict["data"]:
            self._xps_dict["data"][entry] = xr.Dataset()

        # Write averaged cycle data to 'data'.
        all_scan_data = [
            value
            for key, value in self._xps_dict["data"][entry].items()
            if scan_key.split("_")[0] in key
        ]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            averaged_scans = np.mean(all_scan_data, axis=0)

        if averaged_scans.size == 1:
            # on first scan in cycle
            averaged_scans = spectrum["data"]["cps_calib"]

        if averaged_scans.shape == energy.shape:
            # TODO: fix this hotfix so that all data can be written

            self._xps_dict["data"][entry][scan_key.split("_")[0]] = xr.DataArray(
                data=averaged_scans,
                coords={"energy": energy},
            )

            # Write scan data to 'data'.
            self._xps_dict["data"][entry][scan_key] = xr.DataArray(
                data=spectrum["data"]["cps_calib"], coords={"energy": energy}
            )

            channels = [key for key in spectrum["data"] if "cps_ch_" in key]
            for channel in channels:
                ch_no = channel.rsplit("_")[-1]
                channel_key = f"{scan_key}_chan{ch_no}"
                # detector_data_key = (
                #     f"{path_map['detector']}/{detector_data_key_child}"
                #     f"_channels_Channel_{ch_no}/counts"
                # )
                cps = np.array(spectrum["data"][channel])

                # # Write raw data to detector.
                # self._xps_dict[detector_data_key] = spectrum["data"]["cps_calib"]
                # Write channel data to 'data'.
                self._xps_dict["data"][entry][channel_key] = xr.DataArray(
                    data=cps, coords={"energy": energy}
                )

            # Add unit for detector data
            detector_data_unit_key = f"{entry_parent}/detector/raw_data/raw/@units"

            detector_data_units = get_units_for_key("detector/raw_data/raw", UNITS)
            if detector_data_units is not None:
                self._xps_dict[detector_data_unit_key] = detector_data_units


KEY_MAP: Dict[str, str] = {
    "Udet": "detector_voltage",
    "Comment": "comments",
    "ElectronEnergy": "start_energy",
    "SpectrumID": "spectrum_id",
    "EpassOrRR": "pass_energy",
    "EnergyType": "energy/@type",
    "Samples": "n_values",
    "Wf": "work_function",
    "Step": "step",
    "Ubias": "electron_bias",
    "DwellTime": "dwell_time",
    "NumScans": "total_scans",
    "LensMode": "lens_mode",
    "Timestamp": "time_stamp",
    "Entrance": "entrance_slit",
    "Exit": "exit_slit",
    "ScanMode": "energy_scan_mode",
    "VoltageRange": "voltage_range",
    # spectrometer settings
    "Coil Current [mA]": "coil_current [mA]",
    "Pre Defl Y [nU]": "pre_deflector_y_current [nU]",
    "Pre Defl X [nU]": "pre_deflector_x_current [nU]",
    "L1 [nU]": "lens1_voltage [nU]",
    "L2 [nU]": "lens2_voltage [nU]",
    "Focus Displacement 1 [nu]": "focus_displacement_current [nU]",
    "Detector Voltage [V]": "detector_voltage [V]",
    "Bias Voltage Electrons [V]": "bias_voltage_electrons [V]",
    "Bias Voltage Ions [V]": "bias_voltage_ions [V]",
    # source settings
    "anode": "source_label",
    "uanode": "source_voltage",
    "iemission": "emission_current",
    "ihv": "source_high_voltage",
    "ufilament": "filament_voltage",
    "ifilament": "filament_current",
    "DeviceExcitationEnergy": "excitation_energy",
    "panode": "anode_power",
    "temperature": "source_temperature",
    # sql metadata map: Dict[str, str]
    "EnergyType": "energy/@type",
    "EpassOrRR": "pass_energy",
    "Wf": "work_function",
    "Timestamp": "time_stamp",
    "Samples": "n_values",
    "ElectronEnergy": "start_energy",
}

VALUE_MAP: Dict[str, Any] = {
    "energy/@type": convert_energy_type,
    "excitation_energy": float,
    "time_stamp": parse_datetime,
    "energy_scan_mode": convert_energy_scan_mode,
}

KEYS_TO_DROP: List[str] = [
    "Work Function",
]

POSSIBLE_DATE_FORMATS: List[str] = ["%Y-%b-%d %H:%M:%S.%f"]


class SleProdigyParser(ABC):
    """
    Generic parser without reading capabilities,
    to be used as template for implementing parsers for different versions.
    """

    def __init__(self):
        self.con = ""
        self.spectra: List[Dict[str, Any]] = []
        self.xml: ET.Element = None
        self.sum_channels: bool = False
        self.remove_align: bool = True

        self.encoding: List[str, float] = ["f", 4]

    def initiate_file_connection(self, filepath: str):
        """Set the sqllite connection of the file to be opened."""
        sql_connection = filepath
        self.con = sqlite3.connect(sql_connection)

    def parse_file(
        self, filepath: str, **kwargs: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Parse the file's data and metadata into a flat list of dictionaries.

        Parameters
        ----------
        filepath : str
            Filepath of the SLE file to be read.
        **kwargs : Dict[str, Any]
            Additional keyword arguments:
               remove_align(bool):
                   Whether or not alignment spectra shall be removed.
               sum_channels(bool):
                   Whether or not channel data shall be summed.

        Returns
        -------
        List[Dict[str, Any]]
            Flat list of dictionaries containing one spectrum each.

        """
        if "remove_align" in kwargs:
            self.remove_align = kwargs["remove_align"]

        if "sum_channels" in kwargs:
            self.sum_channels = kwargs["sum_channels"]

        # initiate connection to sql file
        self.initiate_file_connection(filepath)

        # read and parse sle file
        self._get_xml_schedule()
        self.spectra = flatten_xml(self.xml)
        self._attach_node_ids()
        self._remove_empty_nodes()
        self._attach_device_protocols()
        self._get_spectrum_metadata_from_sql()
        self._check_encoding()

        self._append_scan_data()

        self._convert_to_common_format()
        self._close_con()

        if self.remove_align:
            self._remove_fixed_energies()

        self._remove_syntax()
        self._remove_snapshot()
        self._reindex_spectra()
        self._reindex_groups()

        return self.spectra

    def _get_xml_schedule(self):
        """Parse the schedule into an XML object."""
        cur = self.con.cursor()
        query = 'SELECT Value FROM Configuration WHERE Key="Schedule"'
        cur.execute(query)
        self.xml = ET.fromstring(cur.fetchall()[0][0])

    def _append_scan_data(self):
        """
        Get the signal data, convert to counts per seconds and get scan
        metadata (Scan no., loop no. and iteration no.) from each scan and
        attach to each spectrum.

        Returns
        -------
        None.

        """
        # pylint: disable=too-many-locals
        individual_scans = []
        scan_id = 0
        for spectrum in self.spectra:
            # copy node to new instance
            node_id = self._get_sql_node_id(spectrum["spectrum_id"])
            n_channels = self._check_energy_channels(node_id)
            raw_ids = self._get_raw_ids(node_id)
            n_scans = len(raw_ids)

            transmission_data = self._get_transmission(node_id)

            for scan_no in range(n_scans):
                scan = copy.copy(spectrum)

                scan["scan_id"] = scan_id
                # get signal data for each scan
                signal_data = self._get_one_scan(raw_ids[scan_no])

                # extract the individual channel data
                signal_data = self._separate_channels(signal_data, n_channels)

                # average channels if required
                if self.sum_channels:
                    signal_data = self._sum_channels(signal_data)

                # convert to counts per second
                signal_data_cps = self._convert_to_counts_per_sec(
                    signal_data, float(scan["dwell_time"])
                )

                # attach individual channel data to scan
                for ch_no, channel_data in enumerate(signal_data_cps):
                    scan[f"cps_ch_{ch_no}"] = list(channel_data)
                # no_of_scans_avg['scans'] = 1

                # scan['cps_calib'] = self._get_calibrated_data(spectrum)
                # """ This is wrong and needs to be corrected!!!"""
                scan["cps_calib"] = copy.copy(scan["cps_ch_0"])

                # Add transmission function
                scan["transmission_function/relative_intensity"] = np.array(
                    transmission_data
                )

                # add metadata including scan, loop no and datetime
                scan_metadata = self._get_scan_metadata(raw_ids[scan_no])
                for key, values in scan_metadata.items():
                    scan[key] = values

                individual_scans += [scan]
                scan_id += 1

        # update self.spectra with the scan data
        self.spectra = individual_scans

    def _get_transmission(self, node_id: int) -> np.ndarray:
        """
        Get the transmission function data.

        Parameters
        ----------
        node_id : int
            Internal node ID of spectrum in SLE sql database.

        Returns
        -------
        transmission_data : array
            Array of TF values for the spectrum at node ID.
        """
        cur = self.con.cursor()
        query = (
            f'SELECT Data, SAMPLES, Ekin FROM TransmissionData WHERE Node="{node_id}"'
        )
        cur.execute(query)
        results = cur.fetchall()
        buffer = self.encoding[1]
        encoding = self.encoding[0]

        stream = []
        for result in results:
            length = result[1] * buffer
            data = result[0]
            for i in range(0, length, buffer):
                stream.append(struct.unpack(encoding, data[i : i + buffer])[0])

        return stream

    def _separate_channels(self, data: List[float], n_channels: int) -> np.ndarray:
        """
        Separate energy channels.

        Parameters
        ----------
        data : list
            Array of measured data
        n_channels : int
            Number of channels to be summed.

        Returns
        -------
        TYPE
            Separate data across n_channels.

        """
        n_points = int(len(data) / n_channels)
        return np.reshape(np.array(data), (n_channels, n_points))

    # """ NEED TO UPDATE THIS METHOD"""
    # def _get_calibrated_data(self, raw_data):
    #     """
    #
    #
    #     Parameters
    #     ----------
    #     raw_data : List
    #         DESCRIPTION.
    #
    #     Returns
    #     -------
    #     channel_dict : TYPE
    #         DESCRIPTION.
    #
    #     """
    #     mcd_num = int(raw_data["mcd_num"])
    #
    #     curves_per_scan = raw_data["curves_per_scan"]
    #     values_per_curve = raw_data["values_per_curve"]
    #     values_per_scan = int(curves_per_scan * values_per_curve)
    #     mcd_head = int(raw_data["mcd_head"])
    #     mcd_tail = int(raw_data["mcd_tail"])
    #     excitation_energy = raw_data["excitation_energy"]
    #     energy_scan_mode = raw_data["energy_scan_mode"]
    #     kinetic_energy = raw_data["kinetic_energy"]
    #     scan_delta = raw_data["scan_delta"]
    #     pass_energy = raw_data["pass_energy"]
    #     kinetic_energy_base = raw_data["kinetic_energy_base"]
    #     # Adding one unit to the binding_energy_upper is added as
    #     # electron comes out if energy is one unit higher
    #     binding_energy_upper = excitation_energy - \
    #         kinetic_energy + kinetic_energy_base + 1
    #
    #     mcd_energy_shifts = raw_data["mcd_shifts"]
    #     mcd_energy_offsets = []
    #     offset_ids = []
    #
    #     # consider offset values for detector with respect to
    #     # position at +16 which is usually large and positive value
    #     for mcd_shift in mcd_energy_shifts:
    #         mcd_energy_offset = (
    #             mcd_energy_shifts[-1] - mcd_shift) * pass_energy
    #         mcd_energy_offsets.append(mcd_energy_offset)
    #         offset_id = round(mcd_energy_offset / scan_delta)
    #         offset_ids.append(
    #             int(offset_id - 1 if offset_id > 0 else offset_id))
    #
    #     # Skipping entry without count data
    #     if not mcd_energy_offsets:
    #         continue
    #     mcd_energy_offsets = np.array(mcd_energy_offsets)
    #     # Putting energy of the last detector as a highest energy
    #     starting_eng_pnts = binding_energy_upper - mcd_energy_offsets
    #     ending_eng_pnts = (starting_eng_pnts
    #                        - values_per_scan * scan_delta)
    #
    #     channeltron_eng_axes = np.zeros((mcd_num, values_per_scan))
    #     for ind in np.arange(len(channeltron_eng_axes)):
    #         channeltron_eng_axes[ind, :] = \
    #             np.linspace(starting_eng_pnts[ind],
    #                         ending_eng_pnts[ind],
    #                         values_per_scan)
    #
    #     channeltron_eng_axes = np.round_(channeltron_eng_axes,
    #                                      decimals=8)
    #     # construct ultimate or incorporated energy axis from
    #     # lower to higher energy
    #     scans = list(raw_data["scans"].keys())
    #
    #     # Check whether array is empty or not
    #     if not scans:
    #         continue
    #     if not raw_data["scans"][scans[0]].any():
    #         continue
    #     # Sorting in descending order
    #     binding_energy = channeltron_eng_axes[-1, :]
    #
    #     self._xps_dict["data"][entry] = xr.Dataset()
    #
    #     for scan_nm in scans:
    #         channel_counts = np.zeros((mcd_num + 1,
    #                                    values_per_scan))
    #         # values for scan_nm corresponds to the data for each
    #         # "scan" in individual CountsSeq
    #         scan_counts = raw_data["scans"][scan_nm]
    #
    #         if energy_scan_mode == "fixed_analyzer_transmission":
    #             for row in np.arange(mcd_num):
    #
    #                 count_on_row = scan_counts[row::mcd_num]
    #                 # Reverse counts from lower to higher
    #                 # BE as in BE_eng_axis
    #                 count_on_row = \
    #                     count_on_row[mcd_head:-mcd_tail]
    #
    #                 channel_counts[row + 1, :] = count_on_row
    #                 channel_counts[0, :] += count_on_row
    #
    #                 # Storing detector's raw counts
    #                 self._xps_dict["data"][entry][f"{scan_nm}_chan_{row}"] = \
    #                     xr.DataArray(data=channel_counts[row + 1, :],
    #                                  coords={"BE": binding_energy})
    #
    #                 # Storing calibrated and after accumulated each scan counts
    #                 if row == mcd_num - 1:
    #                     self._xps_dict["data"][entry][scan_nm] = \
    #                         xr.DataArray(data=channel_counts[0, :],
    #                                      coords={"BE": binding_energy})
    #         else:
    #             for row in np.arange(mcd_num):
    #
    #                 start_id = offset_ids[row]
    #                 count_on_row = scan_counts[start_id::mcd_num]
    #                 count_on_row = count_on_row[0:values_per_scan]
    #                 channel_counts[row + 1, :] = count_on_row
    #
    #                 # shifting and adding all the curves.
    #                 channel_counts[0, :] += count_on_row
    #
    #                 # Storing detector's raw counts
    #                 self._xps_dict["data"][entry][f"{scan_nm}_chan{row}"] = \
    #                     xr.DataArray(data=channel_counts[row + 1, :],
    #                                  coords={"BE": binding_energy})
    #
    #                 # Storing calibrated and after accumulated each scan counts
    #                 if row == mcd_num - 1:
    #                     self._xps_dict["data"][entry][scan_nm] = \
    #                         xr.DataArray(data=channel_counts[0, :],
    #                                      coords={"BE": binding_energy})
    #
    #     # Skipping entry without count data
    #     if not mcd_energy_offsets:
    #         continue
    #     mcd_energy_offsets = np.array(mcd_energy_offsets)
    #     # Putting energy of the last detector as a highest energy
    #     starting_eng_pnts = binding_energy_upper - mcd_energy_offsets
    #     ending_eng_pnts = (starting_eng_pnts
    #                        - values_per_scan * scan_delta)
    #
    #     channeltron_eng_axes = np.zeros((mcd_num, values_per_scan))
    #     for ind in np.arange(len(channeltron_eng_axes)):
    #         channeltron_eng_axes[ind, :] = \
    #             np.linspace(starting_eng_pnts[ind],
    #                         ending_eng_pnts[ind],
    #                         values_per_scan)
    #
    #     channeltron_eng_axes = np.round_(channeltron_eng_axes,
    #                                      decimals=8)
    #     # construct ultimate or incorporated energy axis from
    #     # lower to higher energy
    #     scans = list(raw_data["scans"].keys())
    #
    #     # Check whether array is empty or not
    #     if not scans:
    #         continue
    #     if not raw_data["scans"][scans[0]].any():
    #         continue
    #     # Sorting in descending order
    #     binding_energy = channeltron_eng_axes[-1, :]
    #
    #     self._xps_dict["data"][entry] = xr.Dataset()
    #
    #     for scan_nm in scans:
    #         channel_counts = np.zeros((mcd_num + 1,
    #                                    values_per_scan))
    #         # values for scan_nm corresponds to the data for each
    #         # "scan" in individual CountsSeq
    #         scan_counts = raw_data["scans"][scan_nm]
    #
    #         if energy_scan_mode == "fixed_analyzer_transmission":
    #             for row in np.arange(mcd_num):
    #
    #                 count_on_row = scan_counts[row::mcd_num]
    #                 # Reverse counts from lower to higher
    #                 # BE as in BE_eng_axis
    #                 count_on_row = \
    #                     count_on_row[mcd_head:-mcd_tail]
    #
    #                 channel_counts[row + 1, :] = count_on_row
    #                 channel_counts[0, :] += count_on_row
    #
    #                 # Storing detector's raw counts
    #                 self._xps_dict["data"][entry][f"{scan_nm}_chan_{row}"] = \
    #                     xr.DataArray(data=channel_counts[row + 1, :],
    #                                  coords={"BE": binding_energy})
    #
    #                 # Storing calibrated and after accumulated each scan counts
    #                 if row == mcd_num - 1:
    #                     self._xps_dict["data"][entry][scan_nm] = \
    #                         xr.DataArray(data=channel_counts[0, :],
    #                                      coords={"BE": binding_energy})
    #         else:
    #             for row in np.arange(mcd_num):
    #
    #                 start_id = offset_ids[row]
    #                 count_on_row = scan_counts[start_id::mcd_num]
    #                 count_on_row = count_on_row[0:values_per_scan]
    #                 channel_counts[row + 1, :] = count_on_row
    #
    #                 # shifting and adding all the curves.
    #                 channel_counts[0, :] += count_on_row
    #
    #                 # Storing detector's raw counts
    #                 self._xps_dict["data"][entry][f"{scan_nm}_chan{row}"] = \
    #                     xr.DataArray(data=channel_counts[row + 1, :],
    #                                  coords={"BE": binding_energy})
    #
    #                 # Storing calibrated and after accumulated each scan counts
    #                 if row == mcd_num - 1:
    #                     self._xps_dict["data"][entry][scan_nm] = \
    #                         xr.DataArray(data=channel_counts[0, :],
    #                                      coords={"BE": binding_energy})
    #
    #     return channel_dict

    def _check_energy_channels(self, node_id: int) -> int:
        """
        Get the number of separate energy channels for the spectrum.

        This checks to see if the spectrum was saved with separated energy
        channels.

        Parameters
        ----------
        node_id : int
            Internal node ID of spectrum in SLE sql database.

        Returns
        -------
        n_channels : int
            Number of separate energy channels for the spectrum at node ID.
        """
        cur = self.con.cursor()
        query = f'SELECT EnergyChns FROM Spectrum WHERE Node="{node_id}"'
        cur.execute(query)
        result = cur.fetchall()
        if len(result) != 0:
            n_channels = result[0][0]
        return n_channels

    def _get_raw_ids(self, node_id: int) -> List[int]:
        """
        Get the raw IDs from SQL.

        There is one raw_id for each individual scan when scans were not
        already averaged in the sle file.
        To know which rows in the detector data table belong to which scans,
        one needs to first get the raw_id from the RawData table.

        Parameters
        ----------
        node_id : int
            Internal node ID of spectrum in SLE sql database.

        Returns
        -------
        List[int]
            List of raw IDs for the given note ID.

        """
        cur = self.con.cursor()
        query = f'SELECT RawId FROM RawData WHERE Node="{node_id}"'
        cur.execute(query)

        return [i[0] for i in cur.fetchall()]

    def _check_number_of_scans(self, node_id: int) -> int:
        """
        Get the number of separate scans for the spectrum.

        Parameters
        ----------
        node_id : int
            Internal node ID of spectrum in SLE sql database.

        Returns
        -------
        int
            Number of separate scans for the spectrum.

        """
        cur = self.con.cursor()
        query = f'SELECT RawId FROM RawData WHERE Node="{node_id}"'
        cur.execute(query)
        return len(cur.fetchall())

    def _get_detector_data(self, node_id: int) -> List[float]:
        """
        Get the detector data from sle file.

        The detector data is stored in the SQLite database as a blob.
        To know which blobs belong to which scans, one needs to first get the
        raw_id from the RawData table.

        Parameters
        ----------
        node_id : int
            Internal node ID of spectrum in SLE sql database.

        Returns
        -------
        detector_data : List[float]
            List of lists with measured data.

        """
        cur = self.con.cursor()
        query = f'SELECT RawID FROM RawData WHERE Node="{node_id}"'
        cur.execute(query)
        raw_ids = [i[0] for i in cur.fetchall()]
        detector_data = []
        if len(raw_ids) > 1:
            for raw_id in raw_ids:
                detector_data += [self._get_one_scan(raw_id)]
        else:
            raw_id = raw_ids[0]
            detector_data = self._get_one_scan(raw_id)

        return detector_data

    def _attach_device_protocols(self):
        """
        Get the device protocol for each node and add the parameters of
        the Phoibos to the spectra table. Occasionally these are not
        recorded, if this is the case just skip the group.

        Returns
        -------
        None.

        """
        # iterate through each spectrum
        for spectrum in self.spectra:
            # convert the xml xps id to the node ID and get the device protocol
            cur = self.con.cursor()
            protocol_node_id = self._get_sql_node_id(spectrum["device_group_id"])
            query = (
                f'SELECT Protocol FROM DeviceProtocol WHERE Node="{protocol_node_id}"'
            )
            result = cur.execute(query).fetchone()

            # if a record was accessed then parse, if not skip
            if result:
                protocol = ET.fromstring(result[0])
                protocol_params = self._get_one_device_protocol(protocol)
                spectrum.update(protocol_params)

    def _get_one_device_protocol(self, protocol: ET.Element) -> Dict[str, Any]:
        """
         Get all parameters for one device protocol

        Parameters
        ----------
        protocol : lxml.etree.Element
            One device protocol.

        Returns
        -------
        protocol_params : dict
            All parameters given in the device protocol.

        """
        protocol_params = {}
        for device in protocol.iter("Command"):
            if "Phoibos" in device.attrib["UniqueDeviceName"]:
                # iterate through the parameters and add to spectrum
                # dict
                for parameter in device.iter("Parameter"):
                    if parameter.attrib["type"] == "double":
                        param_text = float(parameter.text)
                    else:
                        param_text = parameter.text
                    protocol_params[parameter.attrib["name"]] = param_text
            elif "XRC1000" in device.attrib["UniqueDeviceName"]:
                for parameter in device.iter("Parameter"):
                    if parameter.attrib["type"] == "double":
                        param_text = float(parameter.text)
                    else:
                        param_text = parameter.text
                    protocol_params[parameter.attrib["name"]] = param_text

        return protocol_params

    def _get_one_scan(self, raw_id: int) -> List[float]:
        """
        Get the detector data for a single scan and convert it to float.

        The detector data is stored in the SQLite database as a blob.
        This function decodes the blob into python float. The blob can be
        encoded as float or double in the SQLite table.

        Parameters
        ----------
        raw_id : int
            Raw ID of the single scan.

        Returns
        -------
        List[float]
            List with measured data.

        """
        cur = self.con.cursor()
        query = f'SELECT Data, ChunkSize FROM CountRateData WHERE RawId="{raw_id}"'
        cur.execute(query)
        results = cur.fetchall()
        buffer = self.encoding[1]
        encoding = self.encoding[0]
        stream = []
        for result in results:
            length = result[1] * buffer
            data = result[0]
            for i in range(0, length, buffer):
                stream.append(struct.unpack(encoding, data[i : i + buffer])[0])
        return stream

    def _parse_external_channels(self, channel: int):
        """
        Parse additional external channels by channel number.

        Parameters
        ----------
        channel : int
            Channel number.

        Returns
        -------
        None.

        """
        if len(channel) != 0:
            pass

    def _get_spectrum_metadata_from_sql(self):
        """
        Get the metadata stored in the SQLite Spectrum table

        Returns
        -------
        None.

        """
        for spectrum in self.spectra:
            node_id = self._get_sql_node_id(spectrum["spectrum_id"])
            cur = self.con.cursor()
            query = f'SELECT * FROM Spectrum WHERE Node="{node_id}"'
            cur.execute(query)
            results = cur.fetchall()
            if len(results) != 0:
                results = results[0]

            column_names = self._get_column_names("Spectrum")
            combined = {
                k: v
                for k, v in dict(zip(column_names, results)).items()
                if k in self.sql_metadata_map
            }
            combined = copy.copy(combined)
            if "EnergyType" not in combined.keys():
                combined["EnergyType"] = "Binding"
            for key, value in combined.items():
                spectrum[key] = value

            query = f'SELECT Data FROM NodeData WHERE Node="{node_id}"'
            cur.execute(query)
            results = ET.fromstring(cur.fetchall()[0][0])
            for i in results.iter("AnalyzerSpectrumParameters"):
                spectrum["work_function"] = i.attrib["Workfunction"]
                spectrum["step_size"] = float(i.attrib["ScanDelta"])

    def _get_scan_metadata(self, raw_id: int) -> Dict[str, Any]:
        """
        Get metadata for each scan.

        Get the scan and the loop/iteration number of each spectrum scan
        and the datetime it was taken from the RawData table.

        Parameters
        ----------
        raw_id : int
            raw id of the scan in the RawData table.

        Returns
        -------
        Dict[str, Any]
            dictionary containing scan metadata.

        """
        # get string Trace from RawData
        cur = self.con.cursor()
        query = f'SELECT ScanDate, Trace FROM RawData WHERE RawID="{raw_id}"'
        result = cur.execute(query).fetchone()
        # process metadata into a dictionary
        scan_meta: Dict[str, Any] = {}
        scan_meta["time_stamp_trace"] = result[0]
        scan_meta.update(self._process_trace(result[1]))

        return scan_meta

    def _process_trace(self, trace: str) -> Dict[str, Any]:
        """
        Parse Trace string to determine the scan, loop, and iteration
        for the given trace.

        Parameters
        ----------
        trace : str
            Trace string to be parsed.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing scan loop and iteration params

        """
        trace_dict: Dict[str, Any] = {}
        loop = re.findall(r"Loop=([0-9]+)u", trace)
        if len(loop) != 0:
            trace_dict["loop_no"] = loop[0]

        scan = re.findall(r"Scan [\[Idx\] ]+=([0-9]+)u", trace)
        if len(scan) != 0:
            trace_dict["scan_no"] = scan[0]

        ramp = re.findall(r"Ramping Iteration [\[Idx\] ]+=([0-9]+)u", trace)
        if len(ramp) != 0:
            trace_dict["iteration_no"] = ramp[0]

        return trace_dict

    def _convert_to_counts_per_sec(
        self, signal_data: np.ndarray, dwell_time: float
    ) -> np.ndarray:
        """
        Convert signal data given in counts to counts per second.

        Parameters
        ----------
        signal_data : np.ndarray
            2D array of floats representing counts
            Shape: (n_channel, n_value)
        dwell_time : float
            Value of dwell_time per scan.

        Returns
        -------
        cps : np.ndarray
            2D array of values converted to counts per second.
            Shape: (n_channel, n_value)

        """
        cps = signal_data / dwell_time
        return cps

    def _get_sql_node_id(self, xml_id: int) -> int:
        """
        Get the SQL internal ID for the NodeID taken from XML.

        Sometimes the NodeID used in XML does not exactly map to the IDs for
        Spectra in the SQL tables. To fix this, there is a node mapping.

        Parameters
        ----------
        xml_id : int
            ID in the XML schedule.

        Returns
        -------
        node_id : int
            ID in the SQL tables.

        """
        cur = self.con.cursor()
        query = f'SELECT Node FROM NodeMapping WHERE InternalID="{xml_id}"'
        cur.execute(query)
        node_id = cur.fetchall()[0][0]
        return node_id

    def _attach_node_ids(self):
        """
        Attach the node_id to each spectrum in the spectra list.

        Returns
        -------
        None.

        """
        for spectrum in self.spectra:
            xml_id = spectrum["spectrum_id"]
            node_id = self._get_sql_node_id(xml_id)
            spectrum["node_id"] = node_id

    def _remove_empty_nodes(self):
        """
        Remove entries from spectra list that have no spectrum in SQLite.

        Returns
        -------
        None.
        """
        for j in reversed(list(enumerate(self.spectra))):
            idx = j[0]
            spectrum = j[1]
            node_id = spectrum["node_id"]
            cur = self.con.cursor()
            query = f'SELECT Node FROM Spectrum WHERE Node="{node_id}"'
            cur.execute(query)
            result = cur.fetchall()
            if len(result) == 0:
                del self.spectra[idx]

    def _get_energy_data(self, spectrum: Dict[str, Any]) -> np.ndarray:
        """
        Create an array of energy values.

        Parameters
        ----------
        spectrum : Dict[str, Any]
            Dictionary with spectrum data and metadata.

        Returns
        -------
        np.ndarray
            Array of uniformly separated energy values.

        """
        if spectrum["energy/@type"] == "binding":
            start = spectrum["start_energy"]
            step = spectrum["step_size"]
            points = spectrum["n_values"]
            energy = [start - i * step for i in range(points)]
        elif spectrum["energy/@type"] == "kinetic":
            start = spectrum["start_energy"]
            step = spectrum["step_size"]
            points = spectrum["n_values"]
            energy = [start + i * step for i in range(points)]
        return np.array(energy)

    def _get_table_names(self) -> List[str]:
        """
        Get a list of table names in the current database file.

        Returns
        -------
        List[str]
            List of spectrum names.

        """
        cur = self.con.cursor()
        cur.execute('SELECT name FROM sqlite_master WHERE type= "table"')
        return [i[0] for i in cur.fetchall()]

    def _get_column_names(self, table_name: str) -> List[str]:
        """
        Get the names of the columns in the table.

        Parameters
        ----------
        table_name : str
            Name of SQL table.

        Returns
        -------
        List[str]
            List of column names.

        """
        cur = self.con.cursor()
        cur.execute(f"SELECT * FROM {table_name}")
        names = [description[0] for description in cur.description]
        return names

    def _close_con(self):
        """
        Close the database connection.

        Returns
        -------
        None.

        """
        self.con.close()

    def _sum_channels(self, data: List[float]) -> np.ndarray:
        """
        Sum together energy channels.

        Parameters
        ----------
        data : List[float]
            List of measured data.

        Returns
        -------
        np.ndarray
            Summed energy channels.

        """
        summed = np.sum(data, axis=0)
        return np.reshape(summed, (1, -1))

    def _check_encoding(self):
        """
        Check whether the binary data should be decoded float or double.

        Returns
        -------
        None.

        """
        cur = self.con.cursor()
        query = "SELECT LENGTH(Data),ChunkSize FROM CountRateData LIMIT 1"
        cur.execute(query)
        data, chunksize = cur.fetchall()[0]

        encodings_map: Dict[str, List[str, float]] = {
            "double": ["d", 8],
            "float": ["f", 4],
        }

        if data / chunksize == 4:
            self.encoding = encodings_map["float"]
        elif data / chunksize == 8:
            self.encoding = encodings_map["double"]
        else:
            logger.error("This binary encoding is not supported.")

    def _reindex_spectra(self):
        """Re-number the spectrum_id."""
        for idx, spectrum in enumerate(self.spectra):
            spectrum["spectrum_id"] = idx

    def _reindex_groups(self):
        """Re-number the group_id."""
        group_ids = list({spec["group_id"] for spec in self.spectra})
        for idx, group_id in enumerate(group_ids):
            for spec in self.spectra:
                if int(spec["group_id"]) == int(group_id):
                    spec["group_id"] = copy.copy(idx)

    def _convert_to_common_format(self):
        """
        Reformat spectra into the format needed for the Mapper object
        """
        maps = {}
        for key_map in self.key_maps:
            maps.update(key_map)
        for spec in self.spectra:
            re_map_keys(spec, maps)
            re_map_values(spec, self.value_map)
            drop_unused_keys(spec, self.keys_to_drop)
            spec["data"] = {}
            spec["data"]["x"] = self._get_energy_data(spec)

            channels = [
                key
                for key in spec
                if any(name in key for name in ["cps_ch_", "cps_calib"])
            ]

            for channel_key in channels:
                spec["data"][channel_key] = np.array(spec[channel_key])
            for channel_key in channels:
                spec.pop(channel_key)

            spec["energy/@units"] = "eV"
            spec["intensity/@units"] = "counts_per_second"

            # Add energy axis for TF data.
            if spec["energy/@type"] == "binding":
                tf_energy = np.array(
                    [spec["excitation_energy"] - x for x in spec["data"]["x"]]
                )
            elif spec["energy/@type"] == "kinetic":
                tf_energy = spec["data"]["x"]

            spec["transmission_function/kinetic_energy"] = tf_energy

    def _remove_fixed_energies(self):
        """
        Remove spectra measured with the scan mode FixedEnergies.
        """
        self.spectra = [
            spec for spec in self.spectra if spec["energy_scan_mode"] != "fixed_energy"
        ]

    def _remove_syntax(self):
        """
        Remove the extra syntax in the group name.
        """
        for spectrum in self.spectra:
            new_name = spectrum["group_name"].split("#", 1)[0]
            new_name = new_name.rstrip(", ")
            spectrum["group_name"] = new_name

    def _remove_snapshot(self):
        """
        Remove spectra required in Snapshot mode.
        """
        self.spectra = [
            spec for spec in self.spectra if "snapshot" not in spec["energy_scan_mode"]
        ]

    def get_sle_version(self) -> str:
        """
        Get the Prodigy SLE version from the file.

        Returns
        -------
        version : str
            Prodigy SLE version of SLE file.

        """
        cur = self.con.cursor()
        query = 'SELECT Value FROM Configuration WHERE Key=="Version"'
        cur.execute(query)
        version = cur.fetchall()[0][0]
        return version


class SleProdigyParserV1(SleProdigyParser):
    """
    Parser for SLE version 1.
    """

    supported_versions = ["1.2", "1.8", "1.9", "1.10", "1.11", "1.12", "1.13"]


class SleProdigyParserV4(SleProdigyParser):
    """
    Parser for SLE version 4.
    """

    supported_versions = [
        "4.63",
        "4.64",
        "4.65",
        "4.66",
        "4.67",
        "4.68",
        "4.69",
        "4.70",
        "4.71",
        "4.72",
        "4.73",
        "4.100",
    ]
