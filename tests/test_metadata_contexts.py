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
Parametrized tests for:
  - mapping._MetadataContext and helper functions
  - per-vendor _context instances (kratos, vms, specs/xy)
"""

import pytest

from pynxtools_xps.mapping import (
    _format_dict,
    _MetadataContext,
    _split_key_and_unit,
    _split_value_and_unit,
    convert_pascal_to_snake,
)
from pynxtools_xps.parsers.kratos.metadata import _context as kratos_context
from pynxtools_xps.parsers.specs.xy.metadata import _context as specs_xy_context
from pynxtools_xps.parsers.vms.metadata import _context as vms_context

# ── mapping helpers ───────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "key, expected_key, expected_unit",
    [
        ("temperature [K]", "temperature", "K"),
        ("voltage [mV]", "voltage", "mV"),
        ("no_unit", "no_unit", None),
        ("key_ [eV]", "key", "eV"),  # trailing underscore stripped from key
        ("key [nU]", "key", "nU"),  # ambiguous unit passed through unchanged
        ("energy [counts/s]", "energy", "counts/s"),
    ],
)
def test_split_key_and_unit(key, expected_key, expected_unit):
    key_out, unit_out = _split_key_and_unit(key)
    assert key_out == expected_key
    assert unit_out == expected_unit


@pytest.mark.parametrize(
    "value, expected_value, expected_unit",
    [
        ("5.0 eV", "5.0", "eV"),
        ("300K", "300", "K"),
        ("1.5e-3 A", "1.5e-3", "A"),
        ("-1.2 mm", "-1.2", "mm"),
        ("hello", "hello", None),  # no numeric prefix → passthrough
        ("42", "42", None),  # integer string without unit → passthrough
    ],
)
def test_split_value_and_unit(value, expected_value, expected_unit):
    value_out, unit_out = _split_value_and_unit(value)
    assert value_out == expected_value
    assert unit_out == expected_unit


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("ScanMode", "scan_mode"),
        ("EnergyType", "energy_type"),
        ("XPS", "xps"),  # all-uppercase → all-lowercase
        ("PassEnergy", "pass_energy"),
        ("already_snake", "already_snake"),
        ("key[unit]", "key[unit]"),  # bracketed part preserved verbatim
        ("ABCDef", "abc_def"),  # internal uppercase run split correctly
    ],
)
def test_convert_pascal_to_snake(raw, expected):
    assert convert_pascal_to_snake(raw) == expected


# ── _MetadataContext — unit tests per method ──────────────────────────────────


@pytest.fixture()
def simple_context():
    """Minimal context used to test individual _MetadataContext methods."""
    return _MetadataContext(
        key_map={"old_key": "new_key", "tilt": "sample_tilt"},
        value_map={"energy_type": lambda v: v.lower()},
        unit_map={"eV": "eV", "s-1": "1/s", "norm": None},
        default_units={"sample_tilt": "degree", "pass_energy": "eV"},
    )


@pytest.mark.parametrize(
    "key, expected",
    [
        ("old_key", "new_key"),  # direct key_map hit
        ("tilt", "sample_tilt"),  # key_map hit after snake conversion
        ("ScanMode", "scan_mode"),  # PascalCase → snake, not in key_map
        ("pass_energy", "pass_energy"),  # not in key_map, already snake
    ],
)
def test_metadata_context_normalize_key(simple_context, key, expected):
    assert simple_context.normalize_key(key) == expected


@pytest.mark.parametrize(
    "value, expected_val, expected_unit",
    [
        ("5.0 eV", "5.0", "eV"),
        ("300", "300", None),  # no unit suffix
        ("text", "text", None),  # non-numeric string
        (None, None, None),  # None passthrough
        (42.0, 42.0, None),  # non-string → no unit extraction attempted
        ("", "", None),  # empty string
    ],
)
def test_metadata_context_parse_value_and_unit(
    simple_context, value, expected_val, expected_unit
):
    v, u = simple_context.parse_value_and_unit(value)
    assert v == expected_val
    assert u == expected_unit


@pytest.mark.parametrize(
    "key, unit, expected_key, expected_unit",
    [
        # Unit already present → returned as-is, key unchanged
        ("temperature", "K", "temperature", "K"),
        # No unit in value → extracted from bracketed key suffix
        ("temperature [K]", None, "temperature", "K"),
        # Nothing to extract → both unchanged
        ("no_unit", None, "no_unit", None),
        # Ambiguous "nU" resolved by key suffix
        ("bias_voltage [nU]", None, "bias_voltage", "V"),
        ("drain_current [nU]", None, "drain_current", "A"),
    ],
)
def test_metadata_context_resolve_unit_from_key(
    simple_context, key, unit, expected_key, expected_unit
):
    k, u = simple_context.resolve_unit_from_key(key, unit)
    assert k == expected_key
    assert u == expected_unit


@pytest.mark.parametrize(
    "key, expected",
    [
        ("sample_tilt", "degree"),
        ("pass_energy", "eV"),
        ("unknown_key", None),  # not in default_units
    ],
)
def test_metadata_context_get_default_unit(simple_context, key, expected):
    assert simple_context.get_default_unit(key) == expected


@pytest.mark.parametrize(
    "unit, expected",
    [
        ("eV", "eV"),  # identity mapping
        ("s-1", "1/s"),  # unit alias
        ("norm", None),  # maps to None (dimensionless / to be dropped)
        ("mm", "mm"),  # not in unit_map → passthrough
        (None, None),  # None input always returns None
    ],
)
def test_metadata_context_map_unit(simple_context, unit, expected):
    assert simple_context.map_unit(unit) == expected


@pytest.mark.parametrize(
    "key, value, expected",
    [
        ("energy_type", "Binding", "binding"),  # value_map function applied
        ("unmapped", "anything", "anything"),  # key absent → passthrough
    ],
)
def test_metadata_context_map_value(simple_context, key, value, expected):
    assert simple_context.map_value(key, value) == expected


@pytest.mark.parametrize(
    "value, expected",
    [
        ("42", 42),  # integer string → int
        ("3.14", 3.14),  # decimal string → float
        ("1e-3", 1e-3),  # scientific notation → float
        ("-5", -5),  # negative integer string → int
        ("hello", "hello"),  # non-numeric → unchanged
        (None, None),  # None → unchanged
        (42, 42),  # already int → unchanged
    ],
)
def test_metadata_context_format_value(simple_context, value, expected):
    assert simple_context._format_value(value) == expected


@pytest.mark.parametrize(
    "key, value, exp_key, exp_val, exp_unit",
    [
        # key_map + value unit
        ("old_key", "5.0 eV", "new_key", 5.0, "eV"),
        # key_map + default unit applied when no unit in value
        ("tilt", "30", "sample_tilt", 30, "degree"),
        # no key_map, default unit
        ("pass_energy", "50", "pass_energy", 50, "eV"),
        # unit embedded in key, extracted and kept
        ("temperature [K]", "300", "temperature", 300, "K"),
        # PascalCase key → snake_case, no value_map entry
        ("ScanMode", "FAT", "scan_mode", "FAT", None),
        # unit-in-key also goes through unit_map ("s-1" → "1/s")
        ("rate [s-1]", "3", "rate", 3, "1/s"),
    ],
)
def test_metadata_context_format(
    simple_context, key, value, exp_key, exp_val, exp_unit
):
    k, v, u = simple_context.format(key, value)
    assert k == exp_key
    assert v == exp_val
    assert u == exp_unit


def test_format_dict(simple_context):
    """`_format_dict` applies context.format to every key/value pair."""
    raw = {"tilt": "30", "old_key": "5.0 eV"}
    result = _format_dict(raw, simple_context)

    assert result["sample_tilt"] == 30
    assert result["sample_tilt/@units"] == "degree"
    assert result["new_key"] == 5.0
    assert result["new_key/@units"] == "eV"


# ── vendor _context integration tests ─────────────────────────────────────────


@pytest.mark.parametrize(
    "key, value, exp_key, exp_val, exp_unit",
    [
        ("tilt", "5", "sample_tilt", 5, "degree"),
        ("lens", "Wide", "lens_mode", "Wide", None),
        ("start", "100.0", "energy_start", 100.0, None),
        ("resolution", "0.5", "resolution", 0.5, "eV"),  # default unit
        ("charge_neutraliser", "yes", "charge_neutraliser", True, None),
        ("anode_library", "Al", "anode_material", "Al", None),
    ],
)
def test_kratos_context_format(key, value, exp_key, exp_val, exp_unit):
    k, v, u = kratos_context.format(key, value)
    assert k == exp_key
    assert v == exp_val
    assert u == exp_unit


@pytest.mark.parametrize(
    "key, value, exp_key, exp_val, exp_unit",
    [
        ("block_id", "C 1s", "region", "C 1s", None),
        (
            "technique",
            "XPS",
            "analysis_method",
            ("XPS", "X-ray photoelectron spectroscopy"),
            None,
        ),
        ("source_energy", "1486.6", "excitation_energy", 1486.6, "eV"),
        (
            "analyzer_mode",
            "FAT",
            "energy_scan_mode",
            "fixed_analyzer_transmission",
            None,
        ),
        ("abscissa_label", "KE", "energy_label", "kinetic", None),
    ],
)
def test_vms_context_format(key, value, exp_key, exp_val, exp_unit):
    k, v, u = vms_context.format(key, value)
    assert k == exp_key
    assert v == exp_val
    assert u == exp_unit


@pytest.mark.parametrize(
    "key, value, exp_key, exp_val, exp_unit",
    [
        # key_map: "analyzer_slit" → "entrance_slit"; no value_map entry
        ("analyzer_slit", "0.3 mm", "entrance_slit", 0.3, "mm"),
        # PascalCase key → snake; value_map "scan_mode" → _convert_energy_scan_mode
        ("ScanMode", "FAT", "scan_mode", "fixed_analyzer_transmission", None),
        # key_map: "values/curve" → "n_values"; value_map "n_values" → int
        ("values/curve", "100", "n_values", 100, None),
        # key_map: "energy_axis" → "x_units"; value_map "x_units" → _convert_energy_type
        ("energy_axis", "Binding Energy", "x_units", "binding", None),
        # unit_map: "s-1" → "1/s" (via unit embedded in key)
        ("rate [s-1]", "3", "rate", 3, "1/s"),
        # default unit
        ("pass_energy", "50", "pass_energy", 50.0, "eV"),
    ],
)
def test_specs_xy_context_format(key, value, exp_key, exp_val, exp_unit):
    k, v, u = specs_xy_context.format(key, value)
    assert k == exp_key
    assert v == exp_val
    assert u == exp_unit
