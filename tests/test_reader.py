"""
Basic example based test for the XPS reader
"""

import os
import xml.etree.ElementTree as ET
from glob import glob
import pytest

from pynxtools.dataconverter.helpers import (
    generate_template_from_nxdl,
    validate_data_dict,
)
from pynxtools.dataconverter.template import Template
from pynxtools.definitions.dev_tools.utils.nxdl_utils import get_nexus_definitions_path

from pynxtools_xps.reader import XPSReader


@pytest.mark.parametrize(
    "sub_reader_data_dir",
    [
        pytest.param(
            "spe",
            id="Phi .spe reader",
        ),
        pytest.param(
            "pro",
            id="Phi .pro reader",
        ),
        pytest.param(
            "scienta_txt",
            id="Scienta txt export reader",
        ),
        pytest.param(
            "vms_regular",
            id="Regular VAMAS reader",
        ),
        pytest.param(
            "vms_irregular",
            id="Irregular VAMAS reader",
        ),
        pytest.param(
            "xml",
            id="Specs XML reader",
        ),
    ],
)
def test_example_data(sub_reader_data_dir):
    """
    Test the example data for the XPS reader
    """
    reader = XPSReader
    assert callable(reader.read)

    def_dir = get_nexus_definitions_path()

    data_dir = os.path.join(os.path.dirname(__file__), "data")
    reader_dir = os.path.join(data_dir, sub_reader_data_dir)

    input_files = sorted(glob(os.path.join(reader_dir, "*")))

    for supported_nxdl in reader.supported_nxdls:
        nxdl_file = os.path.join(
            def_dir, "contributed_definitions", f"{supported_nxdl}.nxdl.xml"
        )

        root = ET.parse(nxdl_file).getroot()
        template = Template()
        generate_template_from_nxdl(root, template)

        read_data = reader().read(
            template=Template(template), file_paths=tuple(input_files)
        )

        assert isinstance(read_data, Template)
        assert validate_data_dict(template, read_data, root)


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
