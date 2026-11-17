# Data from SPECS instruments

The reader supports [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) and SpecsLab 2 files from [SPECS GmbH](https://www.specs-group.com/specs/).
The parsers are in
[`src/pynxtools_xps/parsers/specs/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/specs).

## Supported formats and versions

| Format | Extension | Software | Supported versions |
| ------ | --------- | -------- | ------------------ |
| SpecsLabProdigy binary | `.sle` | SpecsLabProdigy | see below |
| SpecsLab 2 XML | `.xml` | SpecsLab 2 | ≥ 4.63 (other versions likely work) |
| SpecsLabProdigy XY export | `.xy` | SpecsLabProdigy | any |
| SpecsLabProdigy parameter-history log | `.slh` | SpecsLabProdigy | see below |
| SpecsLabProdigy parameter-history CSV export | `.csv` | SpecsLabProdigy | any |

Supported `.sle` version ranges (derived from
[`SPECSSLEParser.supported_versions`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/specs/sle/parser.py)):

{{ parser_version_table("specs.sle.parser", "SPECSSLEParser") }}

If your file is rejected with a version error, check the SpecsLabProdigy version listed
in your SLE file against the ranges above.

Supported `.slh` version ranges (derived from
[`SPECSMetadataSLHParser.supported_versions`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/specs/slh/parser.py)):

{{ parser_version_table("specs.slh.parser", "SPECSMetadataSLHParser") }}

## .sle data

Example data is available in the
[`examples/specs/sle/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/sle).

```console
pynx convert --params-file params.yaml
```

The `params.yaml` file supports a `remove_align` keyword specific to the SLE parser.
Setting it to `true` removes alignment spectra acquired during the experiment, which can
considerably speed up conversion for large files.

## .xml data

Example data is available in the
[`examples/specs/xml/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/xml).

```console
pynx convert In-situ_PBTTT_XPS_SPECS.xml eln_data_xml.yaml --reader xps --nxdl NXxps --output In-situ_PBTTT.nxs
```

## .xy data

Example data is available in the
[`examples/specs/xy/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/xy).

```console
pynx convert MgFe2O4.xy eln_data_xy.yaml --reader xps --nxdl NXxps --output MgFe2O4.nxs
```

## .slh / .csv parameter-history logs

`.slh` and `.csv` files are metadata-only: they add device parameter readings (pressure,
temperature, voltages, currents, ...) to entries produced by an `.sle` file, rather than
producing entries of their own. See
[Explanation > SPECS parameter-history logs](../explanation/slh-log-mapping-and-nxxps-metadata.md)
for the file formats, the entry-merging algorithm, and full usage examples. Test fixtures
are available under
[`tests/data/specs_slh/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/tests/data/specs_slh)
and
[`tests/data/specs_csv/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/tests/data/specs_csv).

```console
pynx convert EX1559_S1710.sle xray.slh nap_parameters.slh phoibos_voltages.slh \
    eln_data.yaml --reader xps --nxdl NXxps --output EX1559_S1710.nxs
```

## Further reading

- [Explanation > Parser architecture](../explanation/parser_architecture.md)
- [Explanation > SPECS parameter-history logs](../explanation/slh-log-mapping-and-nxxps-metadata.md)
