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
Parser for reading XPS data from SPECS instruments in the
proprietary SPECS Lab Prodigy SLE format.
"""

import copy
import re
import sqlite3
import zlib
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import xarray as xr
from lxml import etree as ET
from scipy.interpolate import interp1d

from pynxtools_xps.logging import _logger
from pynxtools_xps.mapping import _format_dict, update_dict_without_overwrite
from pynxtools_xps.parsers.base import ParsedSpectrum, _construct_entry_name, _XPSParser
from pynxtools_xps.parsers.specs.sle.flatten_xml import (
    _iterate_xml_at_tag,
    flatten_context,
    flatten_metainfo,
    flatten_schedule,
)
from pynxtools_xps.parsers.specs.sle.metadata import _context
from pynxtools_xps.parsers.versioning import (
    VersionRange,
    VersionTuple,
    normalize_version,
)


def _execute_sql_query_with_cur(cur: sqlite3.Cursor, query: str):
    """Execute a query with a sqlite3 Cursor object."""
    cur.execute(query)
    return cur.fetchall()


class SPECSSLEParser(_XPSParser):
    """
    Generic parser without reading capabilities,
    to be used as template for implementing parsers for different versions.
    """

    config_file: ClassVar[str] = "config_specs_sle.json"
    supported_file_extensions: ClassVar[tuple[str, ...]] = (".sle",)
    _metadata_exclude_keys: ClassVar[frozenset[str]] = frozenset({"data"})
    supported_versions: ClassVar[tuple[VersionRange, ...]] = (
        ((1, 1), (4, 0)),  # 1.*, 2.*, 3.*
        ((4, 1), (4, 101)),  # 4.1 – 4.100
    )
    _SQLITE_MAGIC: ClassVar[bytes] = b"SQLite format 3\x00"

    def matches_file(self, file: Path) -> bool:
        """Return True for SPECS SLE files (SQLite 3 with SpecsLabSchedule schema)."""
        try:
            with open(file, "rb") as f:
                if f.read(16) != self._SQLITE_MAGIC:
                    return False
            conn = sqlite3.connect(str(file))
            cur = conn.cursor()
            cur.execute("SELECT Value FROM Configuration WHERE Key='Schedule' LIMIT 1")
            row = cur.fetchone()
            conn.close()
            return bool(row and row[0] and "SpecsLabSchedule" in row[0])
        except Exception:
            return False

    def __init__(self):
        super().__init__()
        self.con: sqlite3.Connection
        self.cur: sqlite3.Cursor

        self.xml_schedule: ET.Element
        self.xml_context: ET.Element
        self.xml_metainfo: ET.Element

        self.remove_align: bool = True

        self.encodings_dtype = {
            "short": np.int16,
            "double": np.float64,
            "float": np.float32,
        }
        self.encoding = np.float32

    def detect_version(self, file: Path) -> VersionTuple | None:
        self._initiate_file_connection(file)
        version = self._get_version()
        version = normalize_version(version)

        return version

    def _parse(self, file: Path, **kwargs) -> None:
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
        self._initiate_file_connection(file)

        query = "SELECT COUNT(*) FROM RawData"
        try:
            n_scans = self._execute_sql_query(query)[0][0]
            _logger.info(f"RawData table has {n_scans} scans")
        except Exception as e:
            _logger.warning(f"Could not count RawData rows ({e})")

        self.version = self._get_version()
        self.app_version = self._get_app_version()

        # read and parse sle file
        self._get_xml_schedule()
        self._get_xml_context()
        self._get_xml_metainfo()

        self._flat_spectra = flatten_schedule(self.xml_schedule)

        for spectrum in self._flat_spectra:
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

        self._close_con()

        if self.remove_align:
            self._remove_fixed_energies()

        self._remove_syntax()
        # self._remove_snapshot()
        self._reindex_spectra()
        self._reindex_groups()

        self._build_parsed_spectra()

    def _build_parsed_spectra(self) -> None:
        """Convert ``self._flat_spectra`` into ``self._data``.

        Each spectrum dict in ``self._flat_spectra`` becomes one ``ParsedSpectrum`` in ``self._data``
        entry.  Scans are stacked along the ``"scan"`` axis; a single
        synthetic ``"cycle"`` dimension is prepended so that the output
        conforms to the ``(cycle, scan, *axes)`` contract.

        Channel data (``scan["channels"]``, shape ``(n_energy, n_channels)``)
        is transposed to ``(n_channels, n_energy)`` and stored in
        ``ParsedSpectrum.raw``.
        """
        for spectrum in self._flat_spectra:
            entry_name = (
                _construct_entry_name(
                    [spectrum.get("group_name", ""), spectrum.get("spectrum_type", "")]
                )
                or "entry"
            )

            scans = spectrum.get("data", {}).get("scans", [])
            if not scans:
                continue

            energy = self._get_energy_data(spectrum)

            # Channel-averaged signal: (n_scans, n_energy) → (1, n_scans, n_energy)
            merged_stack = np.stack(
                [np.array(scan["merged"]) for scan in scans], axis=0
            )
            data_da = xr.DataArray(
                data=merged_stack[np.newaxis, ...],
                dims=("cycle", "scan", "energy"),
                coords={"energy": energy},
            )

            # Raw per-channel data: (n_scans, n_energy, n_channels) →
            #   (1, n_scans, n_channels, n_energy)
            raw_da: xr.DataArray | None = None
            if all(
                scan.get("channels") is not None
                and np.array(scan["channels"]).ndim == 2
                for scan in scans
            ):
                channels_stack = np.stack(
                    [np.array(scan["channels"]) for scan in scans], axis=0
                )  # (n_scans, n_energy, n_channels)
                channels_arr = channels_stack.transpose(0, 2, 1)
                raw_da = xr.DataArray(
                    data=channels_arr[np.newaxis, ...],
                    dims=("cycle", "scan", "channel", "energy"),
                    coords={"energy": energy},
                )

            metadata = self._filter_metadata(spectrum)

            for key in list(metadata):
                unit_key = f"{key}/@units"
                if unit_key in metadata:
                    continue
                unit = _context.get_default_unit(key)
                if unit:
                    metadata[unit_key] = unit

            self._data[entry_name] = ParsedSpectrum(
                data=data_da,
                raw=raw_da,
                metadata=metadata,
            )

    def _initiate_file_connection(self, file: str | Path):
        """Set the SQLlite connection of the file to be opened."""
        sql_connection = file
        self.con = sqlite3.connect(sql_connection)
        self.cur = self.con.cursor()

    def _execute_sql_query(self, query: str):
        """Execute a query on the file."""
        return _execute_sql_query_with_cur(self.cur, query)

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
        Fetch scan signal data and metadata from the SQLite database and
        attach them to each spectrum in ``_flat_spectra``.

        For each spectrum the method:

        1. Loads detector calibration and computes per-channel energy shifts.
        2. Loads transmission-function data.
        3. Iterates over raw scan IDs and builds one scan dict per scan via
           :math:`_build_scan`.
        4. Falls back to unit-transmission when no TF data are available.
        5. Normalizes vendor XML key names via :func:`_format_dict`.
        6. Computes the TF kinetic-energy axis via
           :math:`_compute_tf_energy_axis`.
        """
        transmission_key = "transmission_function/relative_intensity"

        for spectrum in self._flat_spectra:
            spectrum["energy/@units"] = "eV"
            spectrum["intensity/@units"] = "counts_per_second"

            spectrum["data"]: dict[str, Any] = {"scans": [], "energy": None}

            group_node_id = self._get_sql_node_id(spectrum["group_id"])
            if not group_node_id:
                continue
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
            node_id = (
                self._get_sql_node_id(spectrum["spectrum_id"])
                or spectrum["spectrum_id"]
            )
            raw_ids = self._get_raw_ids(node_id)
            if not raw_ids:
                _logger.warning(f"No raw_ids found for node {node_id}")

            spectrum[transmission_key] = self._get_transmission(node_id)
            spectrum["abscissa_info"] = self._get_sql_abscissa_info(node_id)

            for scan_id, raw_id in enumerate(raw_ids):
                scan = self._build_scan(raw_id, scan_id, spectrum, n_channels)
                spectrum["data"]["scans"].append(scan)
                if spectrum["data"]["energy"] is None:
                    spectrum["data"]["energy"] = scan["energy"]
                self._check_scan_length(
                    node_id, scan["energy"], spectrum["abscissa_info"]
                )

            if spectrum["data"]["energy"] is None:
                _logger.error(
                    f"No valid x-axis information available for node {node_id}: "
                    "missing abscissa_info and scan data"
                )
                raise ValueError("Missing x-axis / abscissa information")

            if spectrum.get(transmission_key) is None:
                spectrum[transmission_key] = np.ones(len(spectrum["data"]["energy"]))

            _format_dict(spectrum, _context)
            self._compute_tf_energy_axis(spectrum)

    def _build_scan(
        self,
        raw_id: int,
        scan_id: int,
        spectrum: dict[str, Any],
        n_channels: int,
    ) -> dict[str, Any]:
        """
        Build and return a single scan dict from raw SQLite data.

        Applies per-channel energy shifts and interpolation for FAT mode via
        :meth:`_apply_channel_shifts`.
        """
        raw_data = self._get_one_scan(raw_id)
        data = self._separate_channels(raw_data, n_channels)

        step_size = spectrum.get("step_size")
        if step_size is None:
            raise ValueError(
                f"step_size could not be determined for spectrum "
                f"{spectrum.get('spectrum_id')}. "
                "ARPES/snapshot modes may require a different energy axis source."
            )
        energy = np.arange(data.shape[0]) * step_size

        if spectrum["energy_scan_mode"] == "fixed_analyzer_transmission":
            energy, data = self._apply_channel_shifts(energy, data, spectrum)

        scan_metadata = self._get_scan_metadata(raw_id)

        return {
            "scan_id": scan_id,
            "energy": energy,
            "channels": data,
            "merged": np.sum(data, axis=1),
            **scan_metadata,
        }

    def _apply_channel_shifts(
        self,
        raw_energy: np.ndarray,
        data: np.ndarray,
        spectrum: dict[str, Any],
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Apply per-channel energy shifts and interpolate onto a common grid.

        Used for fixed-analyser-transmission (FAT) mode, where each detector
        channel is offset by a calibrated shift.  If no shift data are
        available the inputs are returned unchanged.

        Returns
        -------
        energy_calib : np.ndarray
            Common energy axis after shift and interpolation.
        data_interpolated : np.ndarray
            Channel data interpolated onto ``new_x``.
        """
        shifts = spectrum["detector_calib"].get("shifts")
        if shifts is None:
            return raw_energy, data

        n_channels = data.shape[1]
        shifted = [
            np.vstack((raw_energy + shifts[i], data[:, i])).T for i in range(n_channels)
        ]

        xmin = max(s[:, 0].min() for s in shifted)
        xmax = min(s[:, 0].max() for s in shifted)
        energy_calib = np.linspace(xmin, xmax, spectrum["num_values"])

        data_interpolated = np.array(
            [interp1d(s[:, 0], s[:, 1], kind="linear")(energy_calib) for s in shifted]
        ).T

        if spectrum.get("energy/@type") == "binding":
            energy_calib = np.flip(energy_calib)

        return energy_calib, data_interpolated

    def _check_scan_length(
        self,
        node_id: int,
        x: np.ndarray,
        abscissa_info: dict[str, Any] | None,
    ) -> None:
        """Warn when a scan's axis length differs from the expected AbscissaInfo value."""
        if abscissa_info is None:
            return

        expected = abscissa_info.get("num_values")
        if expected is not None and expected != len(x):
            _logger.warning(
                f"Node {node_id} axis length mismatch -> expected {expected}, got {len(x)}"
            )

    def _compute_tf_energy_axis(self, spectrum: dict[str, Any]) -> None:
        """
        Compute and store the kinetic-energy axis for the transmission function.

        For binding-energy axes the KE is derived as
        ``excitation_energy - BE``.  Fixed-energy mode has no meaningful TF
        energy axis and is skipped.
        """
        if spectrum.get("energy_scan_mode") == "fixed_energy":
            return

        energy = spectrum.get("data", {}).get("energy")
        energy_type = spectrum.get("energy/@type")
        excitation_energy = spectrum.get("excitation_energy")

        if energy_type == "binding" and excitation_energy:
            tf_energy = np.array([excitation_energy - v for v in energy])
        else:
            tf_energy = energy

        rel_intensity = spectrum.get("transmission_function/relative_intensity")
        if rel_intensity is not None and tf_energy is not None:
            n = min(len(tf_energy), len(rel_intensity))
            tf_energy = tf_energy[:n]
            spectrum["transmission_function/relative_intensity"] = rel_intensity[:n]

        spectrum["transmission_function/kinetic_energy"] = tf_energy

    def _get_detector_calibration(self, node_id: int):
        """Extract detector calibration for given node_id."""
        query = f'SELECT Data FROM NodeData WHERE Node="{node_id}"'
        elem = ET.fromstring(self._execute_sql_query(query)[0][0])

        detectors: dict[str, Any] = {}

        detectors["info"] = _iterate_xml_at_tag(elem, "DetectorCalibration")

        for detector_no, sub_elem in enumerate(elem.iter("Detector")):
            detector = {}
            for key, value in sub_elem.attrib.items():
                key, value, unit = _context.format(key, value)
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
                    abscissa_info["start"] = float(asp.attrib["KineticEnergy"])
                elif "Ebin" in asp.attrib:
                    abscissa_info["start"] = float(asp.attrib["Ebin"])

                if "ScanDelta" in asp.attrib:
                    abscissa_info["step_size"] = float(asp.attrib["ScanDelta"])
                elif all(k in asp.attrib for k in ("Ebin", "End", "ValuesPerCurve")):
                    ebin = float(asp.attrib.get("Ebin", 0))
                    end = float(asp.attrib.get("End", 0))
                    n_points = int(asp.attrib.get("ValuesPerCurve", 1))
                    if n_points > 1:
                        abscissa_info["StepSize"] = (end - ebin) / (n_points - 1)

                if "ValuesPerCurve" in asp.attrib:
                    abscissa_info["num_values"] = int(asp.attrib["ValuesPerCurve"])

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
            _logger.info(f"No ExtensionChannelInfo found for node {node_id}.")
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

        _logger.info(f"No ExtensionChannelInfo found for node {node_id}.")
        return None

    # def _add_extension_data(self):
    #     for channel in spectrum.channels.values()[1:]:
    #         # TODO: this is a temporary fix, could add __iter__ to spectrum
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
            _logger.warning(f"could not resolve raw_ids for {node_id}: {e}")

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
        for spectrum in self._flat_spectra:
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
        the spectra table. Occasionally these are not
        recorded, if this is the case just skip the group.

        Returns
        -------
        None.

        """
        # iterate through each spectrum
        for spectrum in self._flat_spectra:
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
                key, value, unit = _context.format(
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
        for spectrum in self._flat_spectra:
            node_id = self._get_sql_node_id(spectrum["spectrum_id"])
            query = f'SELECT * FROM Spectrum WHERE Node="{node_id}"'
            results = self._execute_sql_query(query)

            if len(results) != 0:
                results = results[0]

            column_names = self._get_column_names("Spectrum")
            combined = {k: v for k, v in dict(zip(column_names, results)).items()}

            _format_dict(combined, _context)
            spectrum.update(combined)

            if spectrum.get("step_size") is None:
                try:
                    # 1. Prefer ScanDelta from AnalyzerSpectrumParameters
                    for spec_xml in self.xml_schedule.findall(
                        ".//AnalyzerSpectrumParameters"
                    ):
                        if "ScanDelta" in spec_xml.attrib:
                            spectrum["step_size"] = float(spec_xml.attrib["ScanDelta"])
                            break

                    # 2. Fallback: calculate from Ebin/End/NumValues if available
                    if spectrum.get("step_size") is None:
                        for spec_xml in self.xml_schedule.findall(
                            ".//FixedAnalyzerTransmissionSettings"
                        ):
                            ebin = float(spec_xml.attrib.get("Ebin", 0))
                            end = float(spec_xml.attrib.get("End", 0))
                            n_points = int(spec_xml.attrib.get("NumValues", 1))
                            if n_points > 1:
                                spectrum["step_size"] = abs(end - ebin) / (n_points - 1)
                                break

                    # 3. Last resort: derive from energy range and n_values
                    if spectrum.get("step_size") is None:
                        start = spectrum.get("binding_energy") or spectrum.get(
                            "kinetic_energy"
                        )
                        end = spectrum.get("end_energy")
                        n = spectrum.get("n_values")
                        if start is not None and end is not None and n and n > 1:
                            spectrum["step_size"] = abs(end - start) / (n - 1)

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

    def _get_sql_node_id(self, xml_id: int) -> int | None:
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
        for spectrum in self._flat_spectra:
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
        for j in reversed(list(enumerate(self._flat_spectra))):
            idx = j[0]
            spectrum = j[1]
            node_id = spectrum["node_id"]
            query = f'SELECT Node FROM Spectrum WHERE Node="{node_id}"'
            result = self._execute_sql_query(query)

            if len(result) == 0:
                del self._flat_spectra[idx]

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
        if isinstance(spectrum.get("data"), dict):
            if spectrum["data"].get("energy") is not None:
                return np.array(spectrum["data"]["energy"])
            if "scans" in spectrum["data"] and len(spectrum["data"]["scans"]) > 0:
                first_scan = spectrum["data"]["scans"][0]
                if first_scan["data"].get("energy") is not None:
                    return np.array(first_scan["energy"])

        step = spectrum.get("step_size")
        if step is None:
            raise ValueError(
                f"step_size could not be determined for spectrum "
                f"{spectrum.get('spectrum_id')}."
            )
        points = spectrum["n_values"]

        energy_type = spectrum.get("energy/@type")
        if energy_type == "binding":
            start = spectrum["binding_energy"]
            energy = [start - i * step for i in range(points)]
        elif energy_type == "kinetic":
            start = spectrum["kinetic_energy"]
            energy = [start + i * step for i in range(points)]
        else:
            _logger.error(
                "Energy axis could not be constructed in _get_energy_data: "
                f"unknown energy type '{energy_type}' or missing metadata"
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
            _logger.error(
                "Unsupported binary encoding for length ratio: %s", length_ratio
            )

    def _decompress_data(self, binary_data: bytes | bytearray) -> bytes | bytearray:
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
        for idx, spectrum in enumerate(self._flat_spectra):
            spectrum["spectrum_id"] = idx

    def _reindex_groups(self):
        """Re-number the group_id."""
        group_ids = list({spec["group_id"] for spec in self._flat_spectra})
        for idx, group_id in enumerate(group_ids):
            for spec in self._flat_spectra:
                if int(spec["group_id"]) == int(group_id):
                    spec["group_id"] = copy.copy(idx)

    def _remove_fixed_energies(self):
        """
        Remove spectra measured with the scan mode FixedEnergies.
        """
        self._flat_spectra = [
            spec
            for spec in self._flat_spectra
            if spec["energy_scan_mode"] != "fixed_energy"
        ]

    def _remove_syntax(self):
        """
        Remove the extra syntax in the group name.
        """
        for spectrum in self._flat_spectra:
            new_name = spectrum["group_name"].split("#", 1)[0]
            new_name = new_name.rstrip(", ")
            spectrum["group_name"] = new_name

    def _remove_snapshot(self):
        """
        Remove spectra required in Snapshot mode.
        """
        self._flat_spectra = [
            spec
            for spec in self._flat_spectra
            if "snapshot" not in spec["energy_scan_mode"]
        ]
