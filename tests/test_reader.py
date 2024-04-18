"""
Basic example based test for the XPS reader
"""

import os
import xml.etree.ElementTree as ET
from glob import glob

from pynxtools.dataconverter.helpers import (
    generate_template_from_nxdl,
    validate_data_dict,
)
from pynxtools.dataconverter.template import Template
from pynxtools.definitions.dev_tools.utils.nxdl_utils import get_nexus_definitions_path

from pynxtools_xps.reader import XPSReader


def test_example_data():
    """
    Test the example data for the XPS reader
    """
    reader = XPSReader
    assert callable(reader.read)

    def_dir = get_nexus_definitions_path()

    data_dir = os.path.join(os.path.dirname(__file__), "data")

    reader_dirs = sorted(glob(os.path.join(data_dir, "*")))

    for reader_dir in reader_dirs:
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
