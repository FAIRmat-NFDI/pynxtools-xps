#
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
"""Entry points for XPS examples."""

try:
    from nomad.config.models.plugins import ExampleUploadEntryPoint
except ImportError as exc:
    raise ImportError(
        "Could not import nomad package. Please install the package 'nomad-lab'."
    ) from exc

xps_example = ExampleUploadEntryPoint(
    title="X-ray Photoelectron Spectroscopy (XPS)",
    category="FAIRmat examples",
    description="""
        This example presents the capabilities of the NOMAD platform to store and standardize X-ray Photoelectron Spectroscopy XPS data.
        It shows the generation of a NeXus file according to the
        [NXmpes](https://manual.nexusformat.org/classes/contributed_definitions/NXmpes.html#nxmpes)
        application definition from an example measurement file and a subseqeuent analysis of that data set.
    """,
    resources=['nomad/examples/*']
)
