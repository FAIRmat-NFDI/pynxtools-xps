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
"""Entry points for XPS apps."""

try:
    from nomad.config.models.plugins import AppEntryPoint
    from nomad.config.models.ui import (
        App,
        Axis,
        Column,
        Menu,
        MenuItemHistogram,
        MenuItemPeriodicTable,
        MenuItemTerms,
        MenuSizeEnum,
        SearchQuantities,
    )
except ImportError as exc:
    raise ImportError(
        "Could not import nomad package. Please install the package 'nomad-lab'."
    ) from exc


schema = "pynxtools.nomad.schema.Root"

xps_app = AppEntryPoint(
    name="XPS App",
    description="App for X-ray photoelectron spectroscopy (XPS) data.",
    app=App(
        label="XPS",
        path="xpsapp",
        category="Experiment",
        description="A search app for NXxps and NXmpes NeXus-based XPS experiment entries.",
        readme=(
            "This app supports search and exploration of X-ray photoelectron spectroscopy (XPS) "
            "data stored in the NeXus format according to the NXxps and NXmpes application definitions."
        ),
        search_quantities=SearchQuantities(
            include=[f"*#{schema}"],
        ),
        columns=[
            Column(title="Entry ID", search_quantity="entry_id", selected=True),
            Column(
                title="File Name",
                search_quantity="mainfile",
                selected=True,
            ),
            Column(
                title="Start Time",
                search_quantity=f"data.ENTRY[*].start_time__field#{schema}",
                selected=True,
            ),
            Column(
                title="Method",
                search_quantity=f"data.ENTRY[*].method__field#{schema}",
                selected=True,
            ),
            Column(
                title="Transitions",
                search_quantity=f"data.ENTRY[*].transitions__field#{schema}",
                selected=True,
            ),
            Column(
                title="Author",
                search_quantity=f"data.ENTRY[*].USER[*].name__field#{schema}",
                selected=True,
            ),
            Column(
                title="Sample",
                search_quantity=f"data.ENTRY[*].SAMPLE[*].name__field#{schema}",
                selected=True,
            ),
            Column(
                title="Definition",
                search_quantity=f"data.ENTRY[*].definition__field#{schema}",
                selected=True,
            ),
        ],
        filters_locked={
            f"data.ENTRY.definition__field#{schema}": ["NXxps", "NXmpes"],
        },
        menu=Menu(
            size=MenuSizeEnum.MD,
            title="Menu",
            items=[
                Menu(
                    title="Material",
                    size=MenuSizeEnum.XXL,
                    items=[
                        MenuItemPeriodicTable(
                            quantity="results.material.elements",
                        ),
                        MenuItemTerms(
                            title="Chemical Formula",
                            quantity="results.material.chemical_formula_iupac",
                            width=6,
                            options=10,
                        ),
                        MenuItemTerms(
                            title="Sample Name",
                            quantity=f"data.ENTRY.SAMPLE.name__field#{schema}",
                            width=6,
                            options=10,
                        ),
                        MenuItemHistogram(
                            x="results.material.n_elements",
                        ),
                    ],
                ),
                Menu(
                    title="Experiment",
                    size=MenuSizeEnum.LG,
                    items=[
                        MenuItemTerms(
                            title="Method",
                            quantity=f"data.ENTRY.method__field#{schema}",
                            width=12,
                            options=10,
                        ),
                        MenuItemTerms(
                            title="Probed Core Levels / Transitions",
                            quantity=f"data.ENTRY.transitions__field#{schema}",
                            width=12,
                            options=15,
                        ),
                        MenuItemTerms(
                            title="Definition",
                            quantity=f"data.ENTRY.definition__field#{schema}",
                            width=12,
                            options=5,
                        ),
                    ],
                ),
                Menu(
                    title="Authors / Origin",
                    size=MenuSizeEnum.LG,
                    items=[
                        MenuItemTerms(
                            title="Entry Author",
                            search_quantity=f"data.ENTRY.USER.name__field#{schema}",
                            width=12,
                            options=5,
                        ),
                        MenuItemTerms(
                            title="Upload Author",
                            search_quantity="authors.name",
                            width=12,
                            options=5,
                        ),
                        MenuItemTerms(
                            title="Affiliation",
                            search_quantity=f"data.ENTRY.USER.affiliation__field#{schema}",
                            width=12,
                            options=5,
                        ),
                    ],
                ),
                Menu(
                    title="Instrument",
                    size=MenuSizeEnum.LG,
                    items=[
                        MenuItemTerms(
                            title="Instrument Name",
                            quantity=f"data.ENTRY.INSTRUMENT.name__field#{schema}",
                            width=12,
                            options=10,
                        ),
                        MenuItemTerms(
                            title="X-ray Source Type",
                            quantity=f"data.ENTRY.INSTRUMENT.source_probe.probe__field#{schema}",
                            width=6,
                            options=5,
                        ),
                        MenuItemTerms(
                            title="X-ray Source Name",
                            quantity=f"data.ENTRY.INSTRUMENT.source_probe.name__field#{schema}",
                            width=6,
                            options=10,
                        ),
                        MenuItemHistogram(
                            title="Photon Energy",
                            x=Axis(
                                title="Photon Energy (eV)",
                                search_quantity=f"data.ENTRY.INSTRUMENT.beam_probe.incident_energy__field#{schema}#float",
                            ),
                        ),
                        MenuItemHistogram(
                            title="Energy Resolution",
                            x=Axis(
                                title="Energy Resolution (eV)",
                                search_quantity=f"data.ENTRY.INSTRUMENT.energy_resolution.resolution__field#{schema}#float",
                            ),
                        ),
                        MenuItemTerms(
                            title="Energy Scan Mode",
                            quantity=f"data.ENTRY.INSTRUMENT.ELECTRONANALYZER.ENERGYDISPERSION.energy_scan_mode__field#{schema}",
                            width=12,
                            options=5,
                        ),
                        MenuItemHistogram(
                            title="Pass Energy",
                            x=Axis(
                                title="Pass Energy (eV)",
                                search_quantity=f"data.ENTRY.INSTRUMENT.ELECTRONANALYZER.ENERGYDISPERSION.pass_energy__field#{schema}#float",
                            ),
                        ),
                        MenuItemHistogram(
                            title="Work Function",
                            x=Axis(
                                title="Work Function (eV)",
                                search_quantity=f"data.ENTRY.INSTRUMENT.ELECTRONANALYZER.work_function__field#{schema}#float",
                            ),
                        ),
                    ],
                ),
                Menu(
                    title="Sample",
                    size=MenuSizeEnum.LG,
                    items=[
                        MenuItemTerms(
                            title="Situation",
                            quantity=f"data.ENTRY.SAMPLE.situation__field#{schema}",
                            width=12,
                            options=5,
                        ),
                        MenuItemHistogram(
                            title="Sample Temperature",
                            x=Axis(
                                title="Sample Temperature",
                                search_quantity=f"data.ENTRY.SAMPLE.temperature_env.temperature_sensor.value__field#{schema}#float",
                            ),
                        ),
                    ],
                ),
                Menu(
                    title="Data Range",
                    size=MenuSizeEnum.LG,
                    items=[
                        MenuItemHistogram(
                            title="Min. Binding Energy",
                            x=Axis(
                                title="Min. Binding Energy (eV)",
                                search_quantity=f"data.ENTRY.DATA.energy__min#{schema}#float",
                            ),
                            width=6,
                        ),
                        MenuItemHistogram(
                            title="Max. Binding Energy",
                            x=Axis(
                                title="Max. Binding Energy (eV)",
                                search_quantity=f"data.ENTRY.DATA.energy__max#{schema}#float",
                            ),
                            width=6,
                        ),
                    ],
                ),
                Menu(
                    title="Peak Fitting",
                    size=MenuSizeEnum.LG,
                    items=[
                        MenuItemTerms(
                            title="Fit Region Label",
                            quantity=f"data.ENTRY.FIT.label__field#{schema}",
                            width=12,
                            options=10,
                        ),
                        MenuItemTerms(
                            title="Peak Label",
                            quantity=f"data.ENTRY.FIT.PEAK.label__field#{schema}",
                            width=12,
                            options=15,
                        ),
                        MenuItemTerms(
                            title="Peak Function Type",
                            quantity=f"data.ENTRY.FIT.PEAK.function.function_type__field#{schema}",
                            width=12,
                            options=10,
                        ),
                        MenuItemHistogram(
                            title="Peak Position (Binding Energy)",
                            x=Axis(
                                title="Peak Position (eV)",
                                search_quantity=f"data.ENTRY.FIT.PEAK.data.position__field#{schema}#float",
                            ),
                        ),
                    ],
                ),
                MenuItemHistogram(
                    title="Start Time",
                    x=f"data.ENTRY.start_time__field#{schema}",
                    autorange=True,
                ),
                MenuItemHistogram(
                    title="Upload Creation Time",
                    x="upload_create_time",
                    autorange=True,
                ),
            ],
        ),
        dashboard={
            "widgets": [
                {
                    "type": "periodic_table",
                    "scale": "linear",
                    "quantity": "results.material.elements",
                    "title": "Periodic Table",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 5, "w": 8, "y": 0, "x": 0},
                        "md": {"minH": 3, "minW": 3, "h": 7, "w": 12, "y": 0, "x": 0},
                        "lg": {"minH": 3, "minW": 3, "h": 7, "w": 14, "y": 0, "x": 0},
                        "xl": {"minH": 3, "minW": 3, "h": 7, "w": 14, "y": 0, "x": 0},
                        "xxl": {"minH": 3, "minW": 3, "h": 8, "w": 14, "y": 0, "x": 0},
                    },
                },
                {
                    "type": "terms",
                    "show_input": False,
                    "scale": "linear",
                    "quantity": f"data.ENTRY.method__field#{schema}",
                    "title": "Method",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 5, "w": 4, "y": 0, "x": 8},
                        "md": {"minH": 3, "minW": 3, "h": 7, "w": 6, "y": 0, "x": 12},
                        "lg": {"minH": 3, "minW": 3, "h": 7, "w": 5, "y": 0, "x": 14},
                        "xl": {"minH": 3, "minW": 3, "h": 7, "w": 5, "y": 0, "x": 14},
                        "xxl": {"minH": 3, "minW": 3, "h": 8, "w": 5, "y": 0, "x": 14},
                    },
                },
                {
                    "type": "terms",
                    "show_input": False,
                    "scale": "linear",
                    "quantity": f"data.ENTRY.transitions__field#{schema}",
                    "title": "Core Levels / Transitions",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 5, "w": 4, "y": 0, "x": 8},
                        "md": {"minH": 3, "minW": 3, "h": 7, "w": 6, "y": 0, "x": 12},
                        "lg": {"minH": 3, "minW": 3, "h": 7, "w": 5, "y": 0, "x": 19},
                        "xl": {"minH": 3, "minW": 3, "h": 7, "w": 5, "y": 0, "x": 19},
                        "xxl": {"minH": 3, "minW": 3, "h": 8, "w": 5, "y": 0, "x": 19},
                    },
                },
                {
                    "type": "histogram",
                    "show_input": False,
                    "autorange": True,
                    "nbins": 30,
                    "scale": "linear",
                    "x": Axis(
                        title="Photon Energy (eV)",
                        search_quantity=f"data.ENTRY.INSTRUMENT.beam_probe.incident_energy__field#{schema}#float",
                    ),
                    "title": "Photon Energy",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 4, "w": 6, "y": 5, "x": 0},
                        "md": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 7, "x": 0},
                        "lg": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 7, "x": 0},
                        "xl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 7, "x": 0},
                        "xxl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 8, "x": 0},
                    },
                },
                {
                    "type": "histogram",
                    "show_input": False,
                    "autorange": True,
                    "nbins": 30,
                    "scale": "linear",
                    "x": Axis(
                        title="Pass Energy (eV)",
                        search_quantity=f"data.ENTRY.INSTRUMENT.ELECTRONANALYZER.ENERGYDISPERSION.pass_energy__field#{schema}#float",
                    ),
                    "title": "Pass Energy",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 4, "w": 6, "y": 5, "x": 6},
                        "md": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 7, "x": 9},
                        "lg": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 7, "x": 9},
                        "xl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 7, "x": 9},
                        "xxl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 8, "x": 9},
                    },
                },
                {
                    "type": "terms",
                    "show_input": False,
                    "scale": "linear",
                    "quantity": f"data.ENTRY.FIT.PEAK.function.function_type__field#{schema}",
                    "title": "Peak Function Type",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 4, "w": 6, "y": 9, "x": 0},
                        "md": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 11, "x": 0},
                        "lg": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 11, "x": 0},
                        "xl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 11, "x": 0},
                        "xxl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 12, "x": 0},
                    },
                },
                {
                    "type": "terms",
                    "show_input": False,
                    "scale": "linear",
                    "quantity": f"data.ENTRY.FIT.PEAK.label__field#{schema}",
                    "title": "Peak Labels",
                    "layout": {
                        "sm": {"minH": 3, "minW": 3, "h": 4, "w": 6, "y": 9, "x": 6},
                        "md": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 11, "x": 9},
                        "lg": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 11, "x": 9},
                        "xl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 11, "x": 9},
                        "xxl": {"minH": 3, "minW": 3, "h": 4, "w": 9, "y": 12, "x": 9},
                    },
                },
            ]
        },
    ),
)
