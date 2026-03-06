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
"""
MKdocs macros for the documentation
"""

import importlib


def _format_version_range(lower, upper) -> str:
    """Format a half-open version interval as a human-readable string."""
    lower_str = ".".join(str(x) for x in lower)
    if upper is None:
        return f">= {lower_str}"
    upper_str = ".".join(str(x) for x in upper)
    return f"{lower_str} – {upper_str} (exclusive)"


def define_env(env):
    """
    This is the hook for defining variables, macros and filters

    - variables: the dictionary that contains the environment variables
    - macro: a decorator function, to declare a macro.
    - filter: a function with one of more arguments,
        used to perform a transformation
    """

    # add to the dictionary of variables available to markdown pages:
    env.variables["version"] = "2023.10"  # Figure out from setuptools-scm eventually

    @env.macro
    def supported_formats_table() -> str:
        """
        Generate a Markdown list of all supported XPS file formats, grouped by vendor.

        Reads ``XPSReader.vendor_map`` and groups extensions by vendor,
        linking each item to the corresponding reference page.

        Example usage in a doc page::

            {{ supported_formats_table() }}
        """
        _VENDOR_META: dict[str, tuple[str, str]] = {
            "scienta": ("Scienta Omicron", "reference/scienta.md"),
            "specs": ("SPECS", "reference/specs.md"),
            "phi": ("PHI Electronics", "reference/phi.md"),
            "unknown": ("VAMAS (ISO 14976)", "reference/vms.md"),
        }
        # Vendor display order
        _VENDOR_ORDER = ["scienta", "specs", "phi", "unknown"]

        try:
            from pynxtools_xps.reader import XPSReader  # noqa: PLC0415
        except ImportError as exc:
            return f"*Could not load `XPSReader`: {exc}*"

        vendor_map: dict = XPSReader.vendor_map
        # Build vendor → sorted list of extensions
        vendor_extensions: dict[str, list[str]] = {}
        for ext, vendors in vendor_map.items():
            for vendor in vendors:
                if vendor not in _VENDOR_META:
                    continue
                vendor_extensions.setdefault(vendor, []).append(ext)

        rows = []
        for vendor in _VENDOR_ORDER:
            if vendor not in vendor_extensions:
                continue
            label, page = _VENDOR_META[vendor]
            exts = ", ".join(f"`{e}`" for e in sorted(vendor_extensions[vendor]))
            rows.append(f"- [{label}]({page}) — {exts}")

        return "\n".join(rows)

    @env.macro
    def parser_version_table(module_path: str, class_name: str) -> str:
        """
        Generate a Markdown table of supported version ranges for a parser class.

        Parameters
        ----------
        module_path:
            Dotted module path relative to pynxtools_xps.parsers,
            e.g. "specs.sle.parser".
        class_name:
            Name of the parser class, e.g. "SPECSSLEParser".

        Example usage in a doc page::

            {{ parser_version_table("specs.sle.parser", "SPECSSLEParser") }}
        """
        try:
            mod = importlib.import_module(f"pynxtools_xps.parsers.{module_path}")
            cls = getattr(mod, class_name)
        except (ImportError, AttributeError) as exc:
            return f"*Could not load `{class_name}` from `{module_path}`: {exc}*"

        ranges = getattr(cls, "supported_versions", ())

        if not ranges:
            return "All file versions are accepted, including files without a version."

        rows = "\n".join(
            f"| {_format_version_range(lower, upper)} |" for lower, upper in ranges
        )
        return f"| Supported versions |\n| ------------------ |\n{rows}"
