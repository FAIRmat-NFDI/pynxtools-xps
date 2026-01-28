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
Parser for reading XPS (X-ray Photoelectron Spectroscopy) data from native
Specs Lab Prodigy SLE format, to be passed to MPES NXDL
(NeXus Definition Language) template.
"""

import copy
import logging
import re
import sqlite3
import warnings
import zlib
from pathlib import Path
from typing import Any, cast

import numpy as np
import xarray as xr
from lxml import etree as ET
from packaging.version import InvalidVersion, Version
from scipy.interpolate import interp1d

from pynxtools_xps.reader_utils import (
    XPSMapper,
    construct_data_key,
    construct_entry_name,
    re_map_keys,
    re_map_values,
    update_dict_without_overwrite,
)
from pynxtools_xps.specs.sle.flatten_xml import (
    flatten_context,
    flatten_metainfo,
    flatten_schedule,
)
from pynxtools_xps.specs.sle.utils import (
    KEY_MAP,
    UNITS,
    VALUE_MAP,
    format_key_value_and_unit,
    iterate_xml_at_tag,
)
from pynxtools_xps.value_mappers import get_units_for_key

logger = logging.getLogger("pynxtools")


def execute_sql_query_with_cur(cur: sqlite3.Cursor, query: str):
    """Execute a query with a sqlite3 Cursor object."""
    cur.execute(query)
    return cur.fetchall()


class SleMapperSpecs(XPSMapper):
    """
    Class for restructuring .sle data file from
    specs vendor into python dictionary.
    """

    config_file = "config_specs_sle.json"

    def __init__(self):
        self.parsers = [
            SleProdigyParser,
        ]

        self.file: str | Path = ""
        self.multiple_spectra_groups: bool = True

        super().__init__()

    def _get_sle_version(self):
        con = sqlite3.connect(self.file)
        query = 'SELECT Value FROM Configuration WHERE Key="Version"'
        return execute_sql_query_with_cur(con.cursor(), query)[0][0]

    def _select_parser(self):
        """
        Select the correct parser for the SLE file version. If multiple parsers
        support the version, the most recent parser is returned.

        Returns
        -------
        Parser
            Specs SLE Parser.

        Raises
        ------
        KeyError
            If no parser supports the given version.

        """

        def is_in_range(version_tuple, start, end):
            """
            Check if a version tuple falls within a specified range.

            Parameters
            ----------
            version_tuple : tuple
                Tuple containing the major and minor version (e.g., (1, 0)).
            start : tuple
                Starting version of the range.
            end : tuple or None
                Ending version of the range. If None, there is no upper limit.

            Returns
            -------
            bool
                True if the version is within the range, False otherwise.
            """

            if end is None:
                return version_tuple >= start
            return start <= version_tuple <= end

        def is_version_supported(version, supported_version_ranges):
            """
            Determine if a version is supported by a parser based on version ranges.

            Parameters
            ----------
            version : str
                Version string to check (e.g., 'v4.75').
            supported_version_ranges : list of tuple
                List of version range tuples, where each tuple has a start and optional end.

            Returns
            -------
            bool
                True if the version is supported, False otherwise.
            """
            try:
                parsed_version = Version(version.lstrip("v"))
            except InvalidVersion:
                raise ValueError(f"Invalid version: {version}")

            version_tuple = (parsed_version.major, parsed_version.minor)

            for start, end in supported_version_ranges:
                if is_in_range(version_tuple, start, end):
                    return True
            return False

        version = self._get_sle_version()

        supporting_parsers = []

        for parser in self.parsers:
            if is_version_supported(version, parser.supported_version_ranges):
                supporting_parsers += [parser]

        if not supporting_parsers:
            raise KeyError(
                f"Version {version} of SPECS Prodigy is currently not supported."
            )

        return supporting_parsers[-1]()  # always use newest parser

    def parse_file(self, file: str | Path, **kwargs: dict[str, Any]) -> dict[str, Any]:
        """
        Parse the file using the parser that fits the Prodigy SLE version.

        Returns flat list of dictionaries containing one spectrum each.

        Parameters
        ----------
        file : str
            String name of the file.
        **kwargs : dict[str, Any]
            Dict with additional keyword arguments.

        Returns
        -------
        dict[str, Any]
            Dict with parsed data.

        """
        self.file = file
        return super().parse_file(file, **kwargs)

    def construct_data(self):
        """
        Map SLE format to NXmpes-ready dict.

        Returns
        -------
        None.

        """
        spectra = copy.deepcopy(self.raw_data)

        if len({spectrum.get("group_name") for spectrum in spectra}) == 1:
            self.multiple_spectra_groups = False

        self._xps_dict["data"] = cast(dict[str, Any], {})

        for spectrum in spectra:
            self._update_xps_dict_with_spectrum(spectrum)

    def _update_xps_dict_with_spectrum(self, spectrum: dict[str, Any]):
        """
        Map one spectrum from raw data to NXmpes-ready dict.

        Parameters
        ----------
        spectrum : dict[str, Any]
            Dictionary with data and metadata for one spectrum.

        Returns
        -------
        None.

        """
        # pylint: disable=too-many-locals,duplicate-code
        entry_parts = []

        parts_to_use = ["group_name"] * bool(self.multiple_spectra_groups) + [
            "spectrum_type"
        ]
        for part in parts_to_use:
            val = spectrum.get(part, None)
            if val:
                entry_parts += [val]

        entry = construct_entry_name(entry_parts)
        entry_parent = f"/ENTRY[{entry}]"

        for key, value in spectrum.items():
            mpes_key = f"{entry_parent}/{key}"
            self._xps_dict[mpes_key] = value

            units = get_units_for_key(key, UNITS)
            if units is not None:
                self._xps_dict[f"{mpes_key}/@units"] = units

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

        scans = spectrum["data"].get("scans", [])
        if scans:
            # Use shortcut x or first scan axis
            energy = np.array(spectrum["data"].get("x", scans[0]["x"]))

            # Write each scan individually
            for scan in scans:
                s_key = f"{scan_key}_scan{scan['scan_id']}"
                self._xps_dict["data"][entry][s_key] = xr.DataArray(
                    data=scan["merged"],
                    coords={"energy": scan["x"]},
                    attrs={
                        k: v
                        for k, v in scan.items()
                        if k not in ["x", "raw", "channels", "merged"]
                    },
                )
                self._xps_dict[f"{entry}/{s_key}/@units"] = "counts_per_second"

                # Also write channel-resolved data
                for i in range(scan["channels"].shape[1]):
                    ch_key = f"{s_key}_chan{i}"
                    self._xps_dict["data"][entry][ch_key] = xr.DataArray(
                        data=scan["channels"][:, i], coords={"energy": scan["x"]}
                    )
                    self._xps_dict[f"{entry}/{ch_key}/@units"] = "counts_per_second"

            # Average scan
            merged_arrays = [scan["merged"] for scan in scans]
            averaged = np.mean(merged_arrays, axis=0)
            self._xps_dict["data"][entry][scan_key] = xr.DataArray(
                data=averaged, coords={"energy": energy}
            )
            self._xps_dict[f"{entry}/{scan_key}/@units"] = "counts_per_second"

            # Add unit for detector data
            detector_data_units: str = get_units_for_key("detector/raw_data/raw", UNITS)
            self._xps_dict[f"{entry_parent}/detector/raw_data/raw/@units"] = (
                detector_data_units if detector_data_units is not None else "counts"
            )


class SleProdigyParser:
    """
    Generic parser without reading capabilities,
    to be used as template for implementing parsers for different versions.
    """

    supported_version_ranges = [
        ((1, 0), None),  # Supports v1.* and higher
        ((2, 0), (2, None)),  # Supports v2.* and higher
        ((3, 0), (3, None)),  # Supports v3.* and higher
        ((4, 0), (4, None)),  # Supports 4.1 through 4.100
    ]

    def __init__(self):
        self.con: sqlite3.Conncetion = None
        self.cur = sqlite3.Cursor = None

        self.spectra: list[dict[str, Any]] = []
        self.xml_schedule: ET.Element = None
        self.xml_context: ET.Element = None
        self.xml_metainfo: ET.Element = None

        self.remove_align: bool = True

        self.encodings_dtype = {
            "short": np.int16,
            "double": np.float64,
            "float": np.float32,
        }
        self.encoding = np.float32

    def parse_file(self, file: str, **kwargs: Any) -> list[dict[str, Any]]:
        """
        Parse the file's data and metadata into a flat list of dictionaries.

        Parameters
        ----------
        filepath : str
            Filepath of the SLE file to be read.
        **kwargs : dict[str, Any]
            Additional keyword arguments:
               remove_align(bool):
                   Whether or not alignment spectra shall be removed.
               sum_channels(bool):
                   Whether or not channel data shall be summed.

        Returns
        -------
        list[dict[str, Any]]
            Flat list of dictionaries containing one spectrum each.

        """
        self.remove_align = kwargs.get("remove_align", True)
        self.sum_channels = kwargs.get("sum_channels", False)

        # initiate connection to sql file
        self.initiate_file_connection(file)

        query = "SELECT COUNT(*) FROM RawData"
        try:
            n_scans = self._execute_sql_query(query)[0][0]
            logger.info(f"RawData table has {n_scans} scans")
        except Exception as e:
            logger.warning(f"Could not count RawData rows ({e})")

        self.version = self._get_version()
        self.app_version = self._get_app_version()

        # read and parse sle file
        self._get_xml_schedule()
        self._get_xml_context()
        self._get_xml_metainfo()

        self.spectra = flatten_schedule(self.xml_schedule)

        for spectrum in self.spectra:
            update_dict_without_overwrite(spectrum, flatten_context(self.xml_context))
            update_dict_without_overwrite(spectrum, flatten_metainfo(self.xml_metainfo))

        self._attach_node_ids()
        self._get_spectrum_metadata_from_sql()
        self._remove_empty_nodes()
        # ToDO: Figure out what to do with the detector data
        # self._get_detector_data()
        self._attach_device_protocols()

        self._check_encoding()

        self._append_scan_data()

        self._convert_to_common_format()
        self._close_con()

        if self.remove_align:
            self._remove_fixed_energies()

        self._remove_syntax()
        # self._remove_snapshot()
        self._reindex_spectra()
        self._reindex_groups()

        return self.spectra

    def initiate_file_connection(self, file: str):
        """Set the SQLlite connection of the file to be opened."""
        sql_connection = file
        self.con = sqlite3.connect(sql_connection)
        self.cur = self.con.cursor()

    def _execute_sql_query(self, query: str):
        """Excute a query on the file."""
        return execute_sql_query_with_cur(self.cur, query)

    def _get_version(self):
        query = 'SELECT Value FROM Configuration WHERE Key="Version"'
        return self._execute_sql_query(query)[0][0]

    def _get_app_version(self):
        query = 'SELECT Value FROM Configuration WHERE Key="AppVersion"'
        return self._execute_sql_query(query)[0][0]

    def _get_xml_from_key(self, key: str):
        query = f"SELECT Value FROM Configuration WHERE Key='{key}'"
        try:
            return ET.fromstring(self._execute_sql_query(query)[0][0])
        except IndexError:
            return None

    def _get_xml_schedule(self):
        """Parse the schedule into an XML object."""
        self.xml_schedule = self._get_xml_from_key("Schedule")

    def _get_xml_context(self):
        """Parse the context into an XML object."""
        self.xml_context = self._get_xml_from_key("Context")

    def _get_xml_metainfo(self):
        """Parse the metainfo into an XML object."""
        self.xml_metainfo = self._get_xml_from_key("MetaInfo")

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

        total = len(self.spectra)
        for idx, spectrum in enumerate(self.spectra, start=1):
            spectrum["data"]: dict[str, Any] = {}

            # copy node to new instance
            group_node_id = self._get_sql_node_id(spectrum["group_id"])

            spectrum["detector_calib"] = self._get_detector_calibration(group_node_id)
            try:
                pass_energy = spectrum["pass_energy_or_retardation_ratio"]

                detector_shifts = [
                    item["shift"]
                    for key, item in spectrum["detector_calib"].items()
                    if key.startswith("detector")
                ]
                spectrum["detector_calib"]["shifts"] = (
                    np.array(detector_shifts) * pass_energy
                )
            except KeyError:
                pass

            n_channels = spectrum["energy_channels"]
            node_id = self._get_sql_node_id(spectrum["spectrum_id"])
            if node_id is None:
                node_id = spectrum["spectrum_id"]
            raw_ids = self._get_raw_ids(node_id)
            if not raw_ids:
                logger.warning(f"No raw_ids found for node {node_id}")

            # Add transmission function
            transmission_data = self._get_transmission(node_id)
            spectrum["transmission_function/relative_intensity"] = (
                np.array(transmission_data) if transmission_data is not None else None
            )

            spectrum["abscissa_info"] = self._get_sql_abscissa_info(node_id)

            abscissa_info = spectrum["abscissa_info"]
            if abscissa_info is not None:
                n_expected = abscissa_info.get("NumValues")
                # logger.info(
                #    f"Node {node_id} AbscissaInfo -> "
                #    f"Start={abscissa_info.get('Start')}, "
                #    f"Step={abscissa_info.get('StepSize')}, "
                #    f"NumValues={n_expected} (actual length checked after scan build)"
                # )

            for scan_id, raw_id in enumerate(raw_ids):
                scan = {"scan_id": scan_id}

                raw_data = self._get_one_scan(raw_id)
                data = self._separate_channels(raw_data, n_channels)

                raw_x = np.arange(data.shape[0]) * spectrum["step_size"]
                scan["x"] = raw_x

                if spectrum["energy_scan_mode"] == "fixed_analyzer_transmission":
                    shifts = spectrum["detector_calib"].get("shifts", None)
                    num_values = spectrum["num_values"]

                    if shifts is not None:
                        raw_spectrum = [
                            np.vstack((raw_x + shifts[i], data[:, i])).T
                            for i in range(n_channels)
                        ]

                        xmin = max(
                            raw_spectrum[n][:, 0].min() for n in range(n_channels)
                        )
                        xmax = min(
                            raw_spectrum[n][:, 0].max() for n in range(n_channels)
                        )

                        new_x = np.linspace(xmin, xmax, num_values)

                        new_spectrum = []
                        for i in range(n_channels):
                            f = interp1d(
                                raw_spectrum[i][:, 0],
                                raw_spectrum[i][:, 1],
                                kind="linear",
                            )
                            new_spectrum.append(f(new_x))

                        data = np.array(new_spectrum).T
                        scan["x"] = new_x
                    else:
                        scan["x"] = raw_x

                scan["raw"] = raw_data
                scan["channels"] = data
                scan["merged"] = np.sum(data, axis=1)

                # add metadata including scan, loop no and datetime
                scan_metadata = self._get_scan_metadata(raw_id)
                for key, values in scan_metadata.items():
                    scan[key] = values

                if "data" not in spectrum or not isinstance(spectrum["data"], dict):
                    spectrum["data"] = {}
                if "scans" not in spectrum["data"]:
                    spectrum["data"]["scans"] = []

                spectrum["data"]["scans"].append(scan)

                # keep a shortcut axis for compatibility
                if "x" not in spectrum["data"] or spectrum["data"]["x"] is None:
                    spectrum["data"]["x"] = scan["x"]

                if abscissa_info is not None:
                    actual_len = len(scan["x"])
                    expected = abscissa_info.get("NumValues")
                    if expected is not None and expected != actual_len:
                        logger.warning(
                            f"Node {node_id} axis length mismatch -> "
                            f"expected {expected}, got {actual_len}"
                        )

            if "x" not in spectrum["data"] or spectrum["data"]["x"] is None:
                logger.error(
                    f"No valid x-axis information available for node {node_id}: "
                    "missing abscissa_info and scan data"
                )
                raise ValueError("Missing x-axis / abscissa information")

            if spectrum.get("transmission_function/relative_intensity") is None:
                n_points = len(spectrum["data"]["x"])
                spectrum["transmission_function/relative_intensity"] = np.ones(n_points)

            # TODO: make this working
            # extension_channels = self._get_extension_channel_info(node_id)
            # setattr(extension_channels, self._format_name(extension_channel.detector), extension_channel)
            # if len(extension_channels)-1 == len(n_spectrum_channels):
            #     for extension_channel in extension_channel
            #         for channel in spectrum.channels:
            #             if getattr(extension_channels, k).name == channel.name:
            #                 for attr in channel.__members__():
            #                     setattr(getattr(extension_channels, k),
            #                             attr, getattr(channel, attr))
            # spectrum.channels = extension_channels

    def _get_detector_calibration(self, node_id: int):
        """Extract detector calibration for given node_id."""
        query = f'SELECT Data FROM NodeData WHERE Node="{node_id}"'
        elem = ET.fromstring(self._execute_sql_query(query)[0][0])

        detectors = {}

        detectors["info"] = iterate_xml_at_tag(elem, "DetectorCalibration")

        for detector_no, sub_elem in enumerate(elem.iter("Detector")):
            detector = {}
            for key, value in sub_elem.attrib.items():
                key, value, unit = format_key_value_and_unit(key, value)
                detector[key] = value
            detectors[f"detector{detector_no}"] = detector
        return detectors

    def _get_transmission(self, node_id: int) -> np.ndarray | None:
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
        try:
            query = f'SELECT TransmissionData, TransmissionLen FROM Spectrum WHERE Node="{node_id}"'
            result = self._execute_sql_query(query)[0]
            blob, n_values = result

            if blob is None:
                return None

            transmission_data = np.frombuffer(blob, dtype=np.float32)

            if n_values and len(transmission_data) >= n_values:
                transmission_data = transmission_data[:n_values]

            return transmission_data

        except (IndexError, sqlite3.OperationalError) as e:
            return None

    def _get_sql_abscissa_info(self, node_id: int):
        """
        Get the Abscissa Info.
        Falls back to AnalyzerSpectrumParameters in Schedule XML
        if AbscissaInfo row does not exist in SQL.
        """
        query = f'SELECT * FROM AbscissaInfo WHERE Node="{node_id}"'
        try:
            results = self._execute_sql_query(query)[0]
            abscissa_info: dict[str, Any] = {}
            for idx, key in enumerate(self._get_column_names("AbscissaInfo")):
                abscissa_info[key] = results[idx]
            return abscissa_info

        except (IndexError, sqlite3.OperationalError):
            abscissa_info = {}
            asp = self.xml_schedule.find(".//AnalyzerSpectrumParameters")
            if asp is not None:
                if "KineticEnergy" in asp.attrib:
                    abscissa_info["Start"] = float(asp.attrib["KineticEnergy"])
                elif "Ebin" in asp.attrib:
                    abscissa_info["Start"] = float(asp.attrib["Ebin"])

                if "ScanDelta" in asp.attrib:
                    abscissa_info["StepSize"] = float(asp.attrib["ScanDelta"])
                elif all(k in asp.attrib for k in ("Ebin", "End", "ValuesPerCurve")):
                    ebin = float(asp.attrib.get("Ebin", 0))
                    end = float(asp.attrib.get("End", 0))
                    n_points = int(asp.attrib.get("ValuesPerCurve", 1))
                    if n_points > 1:
                        abscissa_info["StepSize"] = (end - ebin) / (n_points - 1)

                if "ValuesPerCurve" in asp.attrib:
                    abscissa_info["NumValues"] = int(asp.attrib["ValuesPerCurve"])

            return abscissa_info if abscissa_info else None

    def _get_extension_channel_info(self, node_id) -> list[dict[str, Any]] | None:
        def _parse_channel_name(channel_name):
            device = re.findall(r"\((.*?)\)", channel_name)[0]
            unit = re.findall(r"\[(.*?)\]", channel_name)[0]
            detector = channel_name.split("(")[0].split("[")[0].strip()
            return {"device": device, "unit": unit, "detector": detector}

        query = f'SELECT Node, Channel, Name FROM ExtensionChannelInfo WHERE Node="{node_id}"'
        try:
            info = self._execute_sql_query(query)[0]
        except (IndexError, sqlite3.OperationalError):
            logger.info(f"No ExtensionChannelInfo found for node {node_id}.")
            return None

        if info:
            extension_channels: list[dict[str, Any]] = []
            # detectors = [entry[2] for entry in info]
            # extract detector names from ()
            # device = [re.findall(r'\((.*?)\)', detector)[0] for detector in detectors]
            # unit = [re.findall(r'\[(.*?)\]', detector)[0] for detector in detectors]
            # # remove content inside () and []
            # detector = [detector.split('(')[0].split('[')[0].strip() for detector in detectors]

            for entry in info:
                extension_channel: dict[str, Any] = {}
                extension_channel["node_id"] = entry[0]
                extension_channel["channel"] = entry[1]
                extension_channel["name"] = entry[2]
                name_info = _parse_channel_name(entry[2])
                for key, value in name_info.items():
                    extension_channel[key] = value

                extension_channels.append(extension_channel)

            return extension_channels

        logger.info(f"No ExtensionChannelInfo found for node {node_id}.")
        return None

    # def _add_extension_data(self):
    #     for channel in spectrum.channels.values()[1:]:
    #         # TODO: this is a temporary fix, could add __iter__ to DataSet
    #         channel.signal = []
    #         for raw_id in spectrum.scans.raw_ids:
    #             channel.signal.append(np.frombuffer(self._getExtensionData(
    #                 raw_id, channel=channel.channel)[-1], dtype=np.float64))
    #         channel.signal = np.array(channel.signal)

    def _separate_channels(self, data: np.ndarray, n_channels: int) -> np.ndarray:
        """
        Separate energy channels.

        Parameters
        ----------
        data : list[float]
            List of measured data.
        n_channels : int
            Number of channels to be summed.

        Returns
        -------
        TYPE
            Separate data across n_channels.

        """
        return data.reshape(data.size // n_channels, n_channels)

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
        query = f'SELECT EnergyChns FROM Spectrum WHERE Node="{node_id}"'
        result = self._execute_sql_query(query)[0][0]
        if len(result) != 0:
            n_channels = result[0][0]
        return n_channels

    def _get_raw_ids(self, node_id: int | str) -> list[int]:
        """
        Ensure proper RawData lookup using Node, not SpectrumID.
        Falls back to inferred Node if needed.
        """
        # Try direct lookup by Node (correct in this SLE schema)
        query = f'SELECT RawId FROM RawData WHERE Node="{node_id}"'
        rows = self._execute_sql_query(query)
        if rows:
            return [r[0] for r in rows]

        # Fallback: maybe node_id was actually a SpectrumID (string)
        try:
            node_row = self._execute_sql_query(
                f'SELECT Node FROM Spectrum WHERE SpectrumID="{node_id}"'
            )
            if node_row:
                node = node_row[0][0]
                rows = self._execute_sql_query(
                    f'SELECT RawId FROM RawData WHERE Node="{node}"'
                )
                if rows:
                    return [r[0] for r in rows]
        except Exception as e:
            logger.warning(f"could not resolve raw_ids for {node_id}: {e}")

        return []

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
        query = f'SELECT RawId FROM RawData WHERE Node="{node_id}"'
        return len(self._execute_sql_query(query))

    def _get_detector_data(self):
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
        detector_data : list[NDArray[np.float_]]
            List of numpy arrays with measured data.
        """
        for spectrum in self.spectra:
            node_id = spectrum.get("node_id")
            query = f'SELECT RawID FROM RawData WHERE Node="{node_id}"'
            raw_ids = [i[0] for i in self._execute_sql_query(query)]

            detector_data: list[np.ndarray] = []

            for raw_id in raw_ids:
                detector_data.append(self._get_one_scan(raw_id))

            spectrum["detector_data"] = detector_data

    def _attach_device_protocols(self):
        """
        Get the device protocol for each node and add the parameters to
        the spectra table. Occassionally these are not
        recorded, if this is the case just skip the group.

        Returns
        -------
        None.

        """
        # iterate through each spectrum
        for spectrum in self.spectra:
            # convert the xml xps id to the node ID and get the device protocol
            protocol_node_id = self._get_sql_node_id(spectrum["device_group_id"])
            query = (
                f'SELECT Protocol FROM DeviceProtocol WHERE Node="{protocol_node_id}"'
            )

            rows = self._execute_sql_query(query)
            if not rows:
                continue  # nothing recorded → skip safely (as docstring says)

            result = rows[0]

            protocol = ET.fromstring(result[0])
            protocol_params = self._get_one_device_protocol(protocol)
            spectrum.update(protocol_params)

    def _get_one_device_protocol(self, protocol: ET.Element) -> dict[str, Any]:
        """
         Get all parameters for one device protocol

        Parameters
        ----------
        protocol : lxml.etree.Element
            One device protocol.

        Returns
        -------
        protocol_params:  dict[str, Any]
            All parameters given in the device protocol.

        """
        protocol_params: dict[str, Any] = {}
        for elem in protocol.iter("Command"):
            unique_device_name = elem.attrib["UniqueDeviceName"]

            for parameter in elem.iter("Parameter"):
                key, value, unit = format_key_value_and_unit(
                    parameter.attrib["name"], parameter.text
                )
                protocol_params[f"{unique_device_name}/{key}"] = value
                if unit:
                    protocol_params[f"{unique_device_name}/{key}/@units"] = unit

        return protocol_params

    def _get_one_scan(self, raw_id: int) -> np.ndarray:
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
        list[float]
            List with measured data.

        """
        query = f'SELECT Data FROM CountRateData WHERE RawId="{raw_id}"'
        data = self._execute_sql_query(query)[0][0]

        data = self._decompress_data(data)

        return np.frombuffer(data, dtype=self.encoding)

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
        if channel != 0:
            pass

    def _get_spectrum_metadata_from_sql(self):
        """
        Get the metadata stored in the SQLite Spectrum table.
        Also patch in step_size from Schedule XML if missing.

        Returns
        -------
        None.
        """
        for spectrum in self.spectra:
            node_id = self._get_sql_node_id(spectrum["spectrum_id"])
            query = f'SELECT * FROM Spectrum WHERE Node="{node_id}"'
            results = self._execute_sql_query(query)

            if len(results) != 0:
                results = results[0]

            column_names = self._get_column_names("Spectrum")
            combined = {
                k: v
                for k, v in dict(zip(column_names, results)).items()
                if k in KEY_MAP
            }
            combined = copy.copy(combined)
            # Default energy type if missing
            if "EnergyType" not in combined.keys():
                combined["EnergyType"] = "Binding"

            for key, value in combined.items():
                spectrum[KEY_MAP[key]] = value

            if "step_size" not in spectrum or spectrum["step_size"] is None:
                try:
                    # 1. Prefer ScanDelta from AnalyzerSpectrumParameters
                    for spec_xml in self.xml_schedule.findall(
                        ".//AnalyzerSpectrumParameters"
                    ):
                        if "ScanDelta" in spec_xml.attrib:
                            spectrum["step_size"] = float(spec_xml.attrib["ScanDelta"])
                            break

                    # 2. Fallback: calculate from Ebin/End/NumValues if available
                    if "step_size" not in spectrum or spectrum["step_size"] is None:
                        for spec_xml in self.xml_schedule.findall(
                            ".//FixedAnalyzerTransmissionSettings"
                        ):
                            ebin = float(spec_xml.attrib.get("Ebin", 0))
                            end = float(spec_xml.attrib.get("End", 0))
                            n_points = int(spec_xml.attrib.get("NumValues", 1))
                            if n_points > 1:
                                spectrum["step_size"] = (end - ebin) / (n_points - 1)
                                break

                except Exception as e:
                    raise RuntimeError(
                        f"Failed to assign step_size for "
                        f"spectrum_id={spectrum.get('spectrum_id')}, error={e}"
                    )

    def _get_scan_metadata(self, raw_id: int) -> dict[str, Any]:
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
        dict[str, Any]
            dictionary containing scan metadata.

        """
        # get string Trace from RawData
        query = f'SELECT ScanDate, Trace FROM RawData WHERE RawID="{raw_id}"'
        result = self._execute_sql_query(query)[0]
        # process metadata into a dictionary
        scan_meta: dict[str, Any] = {}
        scan_meta["time_stamp_trace"] = result[0]
        scan_meta.update(self._process_trace(result[1]))

        return scan_meta

    def _process_trace(self, trace: str) -> dict[str, Any]:
        """
        Parse Trace string to determine the scan, loop, and iteration
        for the given trace.

        Parameters
        ----------
        trace : str
            Trace string to be parsed.

        Returns
        -------
        dict[str, Any]
            Dictionary containing scan loop and iteration params

        """
        trace_dict: dict[str, Any] = {}
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
        query = f'SELECT Node FROM NodeMapping WHERE InternalID="{xml_id}"'
        try:
            return self._execute_sql_query(query)[0][0]
        except IndexError:
            return None

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
            query = f'SELECT Node FROM Spectrum WHERE Node="{node_id}"'
            result = self._execute_sql_query(query)

            if len(result) == 0:
                del self.spectra[idx]

    def _get_energy_data(self, spectrum: dict[str, Any]) -> np.ndarray:
        """
        Create an array of energy values.

        Parameters
        ----------
        spectrum : dict[str, Any]
            Dictionary with spectrum data and metadata.

        Returns
        -------
        np.ndarray
            Array of uniformly separated energy values.
        """
        if "data" in spectrum and isinstance(spectrum["data"], dict):
            if "energy" in spectrum["data"] and spectrum["data"]["energy"] is not None:
                return np.array(spectrum["data"]["energy"])
            if "x" in spectrum["data"] and spectrum["data"]["x"] is not None:
                return np.array(spectrum["data"]["x"])
            if "scans" in spectrum["data"] and len(spectrum["data"]["scans"]) > 0:
                first_scan = spectrum["data"]["scans"][0]
                if "x" in first_scan and first_scan["x"] is not None:
                    return np.array(first_scan["x"])

        if spectrum["energy/@type"] == "binding":
            start = spectrum["binding_energy"]
            step = spectrum["step_size"]
            points = spectrum["n_values"]
            energy = [start - i * step for i in range(points)]
        elif spectrum["energy/@type"] == "kinetic":
            start = spectrum["kinetic_energy"]
            step = spectrum["step_size"]
            points = spectrum["n_values"]
            energy = [start + i * step for i in range(points)]
        else:
            logger.error(
                "Energy axis could not be constructed in _get_energy_data: "
                f"unknown energy type '{spectrum.get('energy/@type')}' or missing metadata"
            )
            raise ValueError("Missing or invalid energy axis information")
        return np.array(energy)

    def _get_table_names(self) -> list[str]:
        """
        Get a list of table names in the current database file.

        Returns
        -------
        list[str]
            List of spectrum names.

        """
        query = 'SELECT name FROM sqlite_master WHERE type= "table"'
        return [i[0] for i in self._execute_sql_query(query)]

    def _get_column_names(self, table_name: str) -> list[str]:
        """
        Get the names of the columns in the table.

        Parameters
        ----------
        table_name : str
            Name of SQL table.

        Returns
        -------
        list[str]
            List of column names.

        """
        self.cur.execute(f"SELECT * FROM {table_name}")
        names = [description[0] for description in self.cur.description]
        return names

    def _close_con(self):
        """
        Close the database connection.

        Returns
        -------
        None.

        """
        self.con.close()

    def _check_encoding(self) -> None:
        """
        Check whether the binary data should be decoded float or double.

        Returns
        -------
        None.

        """
        query = "SELECT Data, ChunkSize FROM CountRateData LIMIT 1"
        binary_data, chunksize = self._execute_sql_query(query)[0]

        binary_data = self._decompress_data(binary_data)

        length_ratio = len(binary_data) / chunksize
        if length_ratio == 2:
            self.encoding = self.encodings_dtype["double"]
        elif length_ratio == 4:
            self.encoding = self.encodings_dtype["float"]
        elif length_ratio == 8:
            self.encoding = self.encodings_dtype["double"]
        else:
            logger.error(
                "Unsupported binary encoding for length ratio: %s", length_ratio
            )

    def _decompress_data(self, binary_data: bytes | bytearray) -> bytes:
        """
        Attempts to decompress binary data using zlib. If decompression fails,
        returns the original data.

        Args:
            binary_data (bytes | bytearray): Compressed binary data to decompress.

        Returns:
            bytes: Decompressed data if successful, otherwise the original binary data.
        """
        try:
            return zlib.decompress(binary_data)
        except zlib.error:  # Catch only zlib-specific decompression errors
            return binary_data

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
        for spec in self.spectra:
            re_map_keys(spec, KEY_MAP)
            re_map_values(spec, VALUE_MAP)

            if "data" not in spec:
                spec["data"] = {}
            if "x" not in spec["data"] or spec["data"]["x"] is None:
                spec["data"]["x"] = self._get_energy_data(spec)
            if "energy" not in spec or spec["energy"] is None:
                if "x" in spec["data"] and spec["data"]["x"] is not None:
                    spec["energy"] = np.array(spec["data"]["x"])
            if "energy" not in spec["data"]:
                spec["data"]["energy"] = np.array(spec["data"]["x"])

            # Move channel arrays into data
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
            if spec["energy_scan_mode"] != "fixed_energy":
                if spec["energy/@type"] == "binding":
                    excitation_energy = spec.get("excitation_energy")
                    if excitation_energy:
                        tf_energy = np.array(
                            [excitation_energy - x for x in spec["data"]["x"]]
                        )
                    else:
                        if spec["data"].get("x") is None and "scans" in spec["data"]:
                            spec["data"]["x"] = self._get_energy_data(spec)
                        tf_energy = spec["data"]["x"]
                elif spec["energy/@type"] == "kinetic":
                    if spec["data"].get("x") is None and "scans" in spec["data"]:
                        spec["data"]["x"] = self._get_energy_data(spec)
                    tf_energy = spec["data"]["x"]
                else:
                    tf_energy = spec["data"]["x"]

                rel_intensity = spec.get("transmission_function/relative_intensity")
                if rel_intensity is not None:
                    if len(tf_energy) > len(rel_intensity):
                        tf_energy = tf_energy[: len(rel_intensity)]
                    elif len(tf_energy) < len(rel_intensity):
                        rel_intensity = rel_intensity[: len(tf_energy)]
                    spec["transmission_function/relative_intensity"] = rel_intensity

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
