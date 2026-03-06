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

Supported `.sle` version ranges (derived from
[`SPECSSLEParser.supported_versions`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/specs/sle/parser.py)):

{{ parser_version_table("specs.sle.parser", "SPECSSLEParser") }}

If your file is rejected with a version error, check the SpecsLabProdigy version listed
in your SLE file against the ranges above.

## .sle data

Example data is available in the
[`examples/specs/sle/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/sle).

```console
user@box:~$ dataconverter --params-file params.yaml
```

The `params.yaml` file supports a `remove_align` keyword specific to the SLE parser.
Setting it to `true` removes alignment spectra acquired during the experiment, which can
considerably speed up conversion for large files.

## .xml data

Example data is available in the
[`examples/specs/xml/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/xml).

```console
user@box:~$ dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_xml.yaml --reader xps --nxdl NXxps --output In-situ_PBTTT.nxs
```

## .xy data

Example data is available in the
[`examples/specs/xy/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/xy).

```console
user@box:~$ dataconverter MgFe2O4.xy eln_data_xy.yaml --reader xps --nxdl NXxps --output MgFe2O4.nxs
```

## Further reading

- [Explanation > Parser architecture](../explanation/parser_architecture.md)
