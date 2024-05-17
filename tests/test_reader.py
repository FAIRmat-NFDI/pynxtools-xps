"""
Basic example based test for the XPS reader
"""

import os
import sys
import logging
import xml.etree.ElementTree as ET
from glob import glob
import pytest

from pynxtools.dataconverter.convert import get_reader
from pynxtools.dataconverter.helpers import (
    generate_template_from_nxdl,
    write_nexus_def_to_entry,
)
from pynxtools.dataconverter.validation import validate_dict_against
from pynxtools.dataconverter.template import Template
from pynxtools.definitions.dev_tools.utils.nxdl_utils import get_nexus_definitions_path

from pynxtools_xps.reader import XPSReader


READER = get_reader("xps")

test_cases = [
    ("phi_spe", "Phi .spe reader"),
    ("phi_pro", "Phi .pro reader"),
    # ("specs_sle", "SPECS .sle reader"),
    ("specs_xml", "SPECS .xml reader"),
    ("specs_xy", "SPECS .xy reader"),
    ("scienta_ibw", "Scienta .ibw reader"),
    ("scienta_txt", "Scienta .txt export reader"),
    ("vms_irregular", "Irregular VAMAS reader"),
    ("vms_regular", "Regular VAMAS reader"),
    # ("vms_txt_export", "Vamas txt export"),
]

test_params = []

for test_case in test_cases:
    for nxdl in READER.supported_nxdls:
        test_params += [pytest.param(nxdl, test_case[0], id=f"{test_case[1]}, {nxdl}")]


@pytest.mark.parametrize(
    "nxdl, sub_reader_data_dir",
    test_params,
)
def test_example_data(nxdl, sub_reader_data_dir, tmp_path, caplog) -> None:
    """
    Test the example data for the XPS reader
    """
    caplog.clear()
    reader = XPSReader
    assert callable(reader.read)

    def_dir = get_nexus_definitions_path()

    data_dir = os.path.join(os.path.dirname(__file__), "data")
    reader_dir = os.path.join(data_dir, sub_reader_data_dir)

    input_files = sorted(glob(os.path.join(reader_dir, "*")))

    nxdl_file = os.path.join(def_dir, "contributed_definitions", f"{nxdl}.nxdl.xml")

    root = ET.parse(nxdl_file).getroot()
    template = Template()
    generate_template_from_nxdl(root, template)

    read_data = reader().read(
        template=Template(template), file_paths=tuple(input_files)
    )

    entry_names = read_data.get_all_entry_names()
    for entry_name in entry_names:
        write_nexus_def_to_entry(read_data, entry_name, nxdl)

    assert isinstance(read_data, Template)

    with caplog.at_level(logging.WARNING):
        is_success = validate_dict_against(nxdl, read_data, ignore_undocumented=True)
        sys.stdout.write(caplog.text)
        assert is_success
    assert caplog.text == ""


## This will be implemented in the future.
# =============================================================================
# def test_vms_mapper():
#     mapper = VamasMapper
#     data_dir = os.path.join(os.path.dirname(__file__), "data", "vms")
#
#     files_and_keys = {
#         "regular": {
#             "/ENTRY[entry]/Group_1 as-loaded/regions/Region_Survey/scan_mode": "REGULAR",
#             "data['2 S1110, UHV, RT, Epass = 20 eV__MgKLL_1']['cycle0']": 1351,
#         },
#         "irregular": {
#             "/ENTRY[entry]/Group_1 as-loaded/regions/Region_Survey/scan_mode": "IRREGULAR",
#             "/ENTRY[entry]/Group_1 as-loaded/regions/Region_Fe2p/instrument/analyser/energydispersion/scan_mode": "FixedAnalyserTransmission",
#             "data['2 S1110, UHV, RT, Epass = 20 eV__Fe3s-Si2p-Mg2s']['cycle0']": 761,
#         },
#     }
#
#     # (C["benz"], "20240122_SBenz_102_20240122_SBenz_SnO2_10nm.vms"),
#     # (r"C:\Users\pielsticker\Lukas\MPI-CEC\Projects\deepxps\RUB MDI\example_spectra_Florian\Co", "Co 2p 0008751 M1.vms"),
#     # (C["pielst"], C["EX889"], "vms", f"{C['EX889']}_regular.vms"),
#     ((C["pielst"], C["EX889"], "vms", f"{C['EX889']}_irregular.vms"),)
#     # (C["schu"], "CleanData-alphaII VOPO4 C2 BC4316.vms"),
#     # (C["pielst"], C["EX889"], "vms", "d_reg.vms"),
#     ((C["pielst"], C["EX889"], "vms", "d_irreg.vms"),)
#
#     for vms_file in os.listdir(data_dir):
#         data = mapper.parse_file(file=file)
#
#
# for k, v in d.items():
#     if isinstance(v, str):
#         # print(k)
#         if "REG" in v:
#             print(k)
# =============================================================================
