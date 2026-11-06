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
"""Tests for the NOMAD app."""

import re

import pytest

try:
    import nomad  # noqa: F401
except ImportError:
    pytest.skip(
        "Skipping NOMAD app tests because nomad-lab is not installed",
        allow_module_level=True,
    )

# this will raise an exception if pydantic model validation fails for the app
from pynxtools_xps.nomad.apps import schema, xps_app  # noqa: PLC0415

# ---------------------------------------------------------------------------
# Basic app metadata
# ---------------------------------------------------------------------------


def test_xps_app_basic_properties():
    """Verify basic metadata of the XPS app."""
    app = xps_app.app

    assert app.label == "XPS"
    assert app.path == "xpsapp"
    assert app.category == "Experiment"
    assert app.description
    assert app.readme


def test_xps_app_schema():
    """App must reference the correct Xps class."""
    assert schema == "pynxtools.nomad.metainfo.applications.xps.Xps"
    filters = xps_app.app.filters_locked
    assert "section_defs.definition_qualified_name" in filters
    assert filters["section_defs.definition_qualified_name"] == [schema]


def test_xps_app_locked_filters():
    """Ensure required locked filters are defined and well-formed."""
    app = xps_app.app

    assert "section_defs.definition_qualified_name" in app.filters_locked
    assert isinstance(
        app.filters_locked["section_defs.definition_qualified_name"], list
    )
    assert len(app.filters_locked["section_defs.definition_qualified_name"]) == 1


def test_xps_app_search_quantities_include_schema():
    """Ensure the app loads search quantities from its own schema."""
    app = xps_app.app

    assert app.search_quantities.include == [f"*#{schema}"]


# ---------------------------------------------------------------------------
# Columns
# ---------------------------------------------------------------------------


def test_xps_app_columns():
    """Check every column is present, titled, and points at a plausible quantity."""
    app = xps_app.app

    expected_titles = {
        "Entry ID",
        "File Name",
        "Start Time",
        "Method",
        "Transitions",
        "Author",
        "Sample",
        "Sample ID",
        "Definition",
    }
    titles = {col.title for col in app.columns}
    assert expected_titles == titles

    definition_column = next(col for col in app.columns if col.title == "Definition")
    assert definition_column.selected is True
    assert "data.definition" in definition_column.search_quantity

    method_column = next(col for col in app.columns if col.title == "Method")
    assert "data.method" in method_column.search_quantity

    sample_id_column = next(col for col in app.columns if col.title == "Sample ID")
    assert sample_id_column.selected is False
    assert "data.sample" in sample_id_column.search_quantity
    assert "identifier" in sample_id_column.search_quantity

    # Sample ID is intentionally not selected by default; everything else is.
    for col in app.columns:
        if col.title == "Sample ID":
            continue
        assert col.selected is True, f"{col.title} should be selected by default"


# ---------------------------------------------------------------------------
# Menu structure
# ---------------------------------------------------------------------------


def test_xps_app_top_level_menu():
    """Every top-level named menu section from the source must be present."""
    app = xps_app.app

    assert app.menu.title == "Menu"
    section_titles = [item.title for item in app.menu.items if hasattr(item, "items")]
    assert section_titles == [
        "Material",
        "Experiment",
        "Authors / Origin",
        "Instrument",
        "Sample",
        "Data Range",
        "Peak Fitting",
    ]


def test_xps_app_menu_contains_material_section():
    """Validate presence and structure of the Material menu section."""
    app = xps_app.app

    material_menu = next(item for item in app.menu.items if item.title == "Material")

    assert material_menu.size.name == "XXL"
    assert any(
        item.__class__.__name__ == "MenuItemPeriodicTable"
        for item in material_menu.items
    )
    assert any(
        getattr(item, "title", None) == "Chemical Formula"
        for item in material_menu.items
    )


def test_xps_app_menu_contains_experiment_section():
    """Validate presence and structure of the Experiment menu section."""
    app = xps_app.app

    experiment_menu = next(
        item for item in app.menu.items if item.title == "Experiment"
    )

    assert experiment_menu.size.name == "LG"
    titles = [item.title for item in experiment_menu.items]
    assert "Method" in titles
    assert "Probed Core Levels / Transitions" in titles


def test_xps_app_menu_contains_authors_origin_section():
    """Validate presence and structure of the Authors / Origin menu section."""
    app = xps_app.app

    authors_menu = next(
        item for item in app.menu.items if item.title == "Authors / Origin"
    )

    titles = [item.title for item in authors_menu.items]
    assert titles == ["Entry Author", "Upload Author", "Affiliation"]


def test_xps_app_menu_contains_instrument_section():
    """Validate presence and structure of the Instrument menu section."""
    app = xps_app.app

    instrument_menu = next(
        item for item in app.menu.items if item.title == "Instrument"
    )

    assert instrument_menu.size.name == "LG"
    titles = [item.title for item in instrument_menu.items]
    assert titles == [
        "Instrument Name",
        "X-ray Source Type",
        "X-ray Source Name",
        "Energy Resolution",
        "Energy Scan Mode",
        "Pass Energy",
        "Work Function",
    ]


def test_xps_app_menu_contains_sample_section():
    """Validate presence and structure of the Sample menu section."""
    app = xps_app.app

    sample_menu = next(item for item in app.menu.items if item.title == "Sample")

    titles = [item.title for item in sample_menu.items]
    assert titles == ["Situation", "Sample Temperature"]


def test_xps_app_menu_contains_data_range_section():
    """Validate presence of the Data Range menu section (FIELD_STATISTICS-derived)."""
    app = xps_app.app

    data_range_menu = next(
        item for item in app.menu.items if item.title == "Data Range"
    )

    titles = [item.title for item in data_range_menu.items]
    assert titles == ["Min. Binding Energy", "Max. Binding Energy"]
    for item in data_range_menu.items:
        assert "energy__m" in item.x.search_quantity  # __min / __max


def test_xps_app_menu_contains_peak_fitting_section():
    """Validate presence and structure of the Peak Fitting menu section."""
    app = xps_app.app

    peak_fitting_menu = next(
        item for item in app.menu.items if item.title == "Peak Fitting"
    )

    assert peak_fitting_menu.size.name == "LG"
    titles = [item.title for item in peak_fitting_menu.items]
    assert titles == [
        "Fit Region Label",
        "Peak Label",
        "Peak Function Type",
        "Peak Position (Binding Energy)",
    ]


def test_xps_app_menu_contains_top_level_histograms():
    """The Start Time / Upload Creation Time histograms sit directly under the
    top-level menu, not inside a named section."""
    app = xps_app.app

    top_level_histograms = [
        item
        for item in app.menu.items
        if item.__class__.__name__ == "MenuItemHistogram"
    ]
    titles = {item.title for item in top_level_histograms}
    assert titles == {"Start Time", "Upload Creation Time"}


# ---------------------------------------------------------------------------
# Dashboard
# ---------------------------------------------------------------------------


def test_xps_app_dashboard_widgets():
    """Check every dashboard widget is present with the expected type/title."""
    dashboard = xps_app.app.dashboard

    expected = {
        "Periodic Table": "periodic_table",
        "Method": "terms",
        "Core Levels / Transitions": "terms",
        "Photon Energy": "histogram",
        "Pass Energy": "histogram",
        "Peak Function Type": "terms",
        "Peak Labels": "terms",
    }
    actual = {w.title: w.type for w in dashboard.widgets}
    assert actual == expected

    for widget in dashboard.widgets:
        assert widget.layout, f"widget {widget.title!r} has no layout"

    periodic_table = next(w for w in dashboard.widgets if w.type == "periodic_table")
    assert periodic_table.search_quantity == "results.material.elements"

    photon_energy = next(w for w in dashboard.widgets if w.title == "Photon Energy")
    assert "data.data.photon_energy" in photon_energy.x.search_quantity


# ---------------------------------------------------------------------------
# Deep check: every search_quantity path referencing the Xps schema must
# actually resolve against the generated metainfo. This is the check that
# would have caught, e.g., the `name_quantity`-vs-`name` and
# FIELD_STATISTICS `__min`/`__max` gaps from Phase 4 automatically, instead
# of only being caught by manual inspection.
# ---------------------------------------------------------------------------


def _iter_search_quantities(node):
    """Recursively yield every populated search_quantity string on a Menu/menu-item
    tree or a dashboard widget list."""
    items = getattr(node, "items", None)
    if items is not None:
        for item in items:
            yield from _iter_search_quantities(item)
        return

    search_quantity = getattr(node, "search_quantity", None)
    if search_quantity:
        yield search_quantity

    axis = getattr(node, "x", None)
    if axis is not None and hasattr(axis, "search_quantity") and axis.search_quantity:
        yield axis.search_quantity


def _schema_relative_path(search_quantity: str) -> str | None:
    """Strip the '#schema[#type]' suffix and the leading 'data.' prefix from a
    search_quantity string, returning the dotted path relative to the Xps
    section, or None if the string isn't a `data.*` reference to our own
    schema (e.g. 'entry_id', 'results.material.elements', 'authors.name')."""
    path, *rest = search_quantity.split("#")
    if not rest or schema not in rest[0]:
        return None
    if not path.startswith("data."):
        return None
    path = path[len("data.") :]
    return re.sub(r"\[\*\]", "", path)


def _resolves_on_xps(dotted_path: str) -> bool:
    """Walk a dotted attribute path against the generated Xps section
    definition, descending through SubSections and terminating on a
    Quantity. Returns False for any segment that doesn't exist — this is
    exactly the class of bug a stale/guessed search_quantity string produces
    silently at the NOMAD UI layer (no error, just an empty facet)."""
    from pynxtools.nomad.metainfo.applications.xps import Xps  # noqa: PLC0415

    current = Xps.m_def
    parts = dotted_path.split(".")
    for i, part in enumerate(parts):
        is_last = i == len(parts) - 1
        if part in current.all_sub_sections:
            current = current.all_sub_sections[part].sub_section
        elif part in current.all_quantities:
            if not is_last:
                return False
        else:
            return False
    return True


def test_xps_app_all_schema_quantities_resolve():
    """Every `data.*#{schema}` search_quantity referenced anywhere in the app
    (columns, menu, dashboard) must be a real, resolvable attribute path on
    the generated Xps class."""
    app = xps_app.app

    all_strings = [col.search_quantity for col in app.columns if col.search_quantity]
    all_strings += list(_iter_search_quantities(app.menu))
    all_strings += list(_iter_search_quantities(app.dashboard))

    checked = 0
    for search_quantity in all_strings:
        relative_path = _schema_relative_path(search_quantity)
        if relative_path is None:
            continue
        assert _resolves_on_xps(relative_path), (
            f"{search_quantity!r} does not resolve on the generated Xps class "
            f"(resolved path: {relative_path!r})"
        )
        checked += 1

    # Sanity check that the walk actually inspected a meaningful number of
    # paths, so this test can't silently pass by finding nothing.
    assert checked >= 20
