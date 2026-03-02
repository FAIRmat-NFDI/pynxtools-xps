---
render_macros: false
---

# Build a parser

Your current data is not supported yet? Don't worry, the following how-to will guide you how to write a reader your own data.

This guide covers two scenarios:

- **Your format is already partially supported** — you only need to adjust how fields are mapped to NeXus.
- **Your format is not supported at all** — you need to write a new parser from scratch.

---

## Your format is already supported, but some fields differ

Good! The basic functionality to read your data is already in place.

The simplest fix is to supply a custom ELN YAML file with the missing or corrected values.
For structural mismatches, consider one of the following before writing a new parser:

1. **Modify a config file.** Each vendor has a JSON config in
   [`src/pynxtools_xps/config/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/config)
   that maps flat dict keys to NeXus paths. Adjusting the mapping there is often enough.
2. **Open a pull request.** If your change benefits all users of that format, propose it on the
   [GitHub repository](https://github.com/FAIRmat-NFDI/pynxtools-xps/pulls).

---

## Adding a completely new format

Adding a new format requires creating a vendor subpackage, registering the parser, and
providing test data.

### 1. Set up a development environment

Follow the instructions in the [development guide](../tutorial/contributing.md).

### 2. Create the vendor subpackage

Add a new directory under
[`src/pynxtools_xps/parsers/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers)
with the following layout:

```text
parsers/<vendor>/
├── __init__.py
├── data_model.py   # typed intermediate representation
├── metadata.py     # _MetadataContext instance
└── parser.py       # _XPSParser subclass
```

Use an existing parser (for example,
[`parsers/vms/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/vms)
or
[`parsers/phi/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/phi))
as a reference for the expected structure.

### 3. Define typed data models (`data_model.py`)

Create dataclasses for each logical record in your format (header, spectrum, region, …)
by subclassing
[`_XPSDataclass`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/base.py).
`_XPSDataclass` enforces type annotations at assignment time and coerces compatible
types automatically. See the
[`phi/data_model.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/phi/data_model.py) for a worked example.

!!! Note
    Defining a data model is recommended, but remains optional. It is only applicable if the data model for a vendor format is always the same and can thus be explicitly written down.

### 4. Define the normalization context (`metadata.py`)

Create a module-level `_context` instance of
[`_MetadataContext`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/mapping.py)
with four mappings:

- `key_map` — vendor-specific key names → canonical names
- `value_map` — canonical key → converter function (use the shared converters from
  [`mapping.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/mapping.py)
  where possible: `_convert_measurement_method`, `_convert_energy_scan_mode`, etc.)
- `unit_map` — vendor unit strings → standard units (`None` = dimensionless/drop)
- `default_units` — canonical key → unit when the value carries none

See
[`parsers/vms/metadata.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/vms/metadata.py)
for a typical example.

### 5. Implement the parser (`parser.py`)

Start by subclassing
[`_XPSParser`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/base.py). Implement the following required class variables and methods:

- `supported_file_extensions` — a tuple of file extensions that your parser supports
- `matches_file(file_path)` — return `True` only if the file unambiguously conforms to
  your format (check headers, magic bytes, or required keywords).
- `_parse(file_path)` — extract all spectra and populate `self._data` as a list of dictionaries, where each dict holds the raw key-value pairs for one spectrum.

If your format embeds a version string, override `detect_version(file_path)` to return
a `VersionTuple` (produced by
[`normalize_version`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/versioning.py))
and declare `supported_versions` and `requires_version` as class variables.

Set `config_file` to the name of the JSON config you will create in step 7.

### 6. Register the parser

- Import and re-export the new mapper class from
  [`parsers/__init__.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/__init__.py).
- Add the parser and mapper to the lists in
  [`reader.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/reader.py).

### 7. Add a config file

Create a `config/config_<vendor>.json` config file by starting from
[`config/template.json`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/config/template.json). Each entry maps a canonical dict key to the NeXus path in `NXxps` or `NXmpes`.

### 8. Add tests

- Place at least one representative example file in `tests/data/<vendor>/`.
- Add your test case to the `test_cases` list in
  [`tests/test_reader.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/tests/test_reader.py).
- Generate the reference `.nxs` file by running the conversion once and inspecting the output. You may also modify [`scripts/generate_reference_files.sh`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/scripts/generate_reference_files.sh) to automate this.
- Run the full test suite to confirm nothing is broken:

```console
pytest -sv tests
```

### 9. Add documentation

Add a page to the References section of the documentation to describe what your parser does, which files and versions are supported, and how it integrates in the parser landscape.

Note that there exists a special mkdocs macro (in [`src/pynxtools_xps/mkdocs.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/mkdocs.py)) to automatically generate the supported file formats and versions of a parser.

It can be called within your reference markdown file like this:

```text
{{ parser_version_table("path.to.your.parser", "YourParser") }}
```

For example, for the `SpecsSLEParser` located in  [`src/pynxtools_xps/parsers/specs/sle/parser.py`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers), the call is:

```md
{{ parser_version_table("specs.sle.parser", "SpecsSLEParser") }}
```

---

## Further reading

- [Explanation > Parser architecture](../explanation/parser_architecture.md) — the three-layer pipeline
- [`NXxps` application definition](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXxps.html)
- [`NXmpes` application definition](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXmpes.html)
