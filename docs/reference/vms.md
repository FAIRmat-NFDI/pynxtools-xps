# VAMAS ISO standard (VMS)

The reader supports VAMAS (`.vms`, `.npl`) files, the ISO standard data transfer format
([ISO 14976](https://www.iso.org/standard/24269.html)) for surface chemical analysis.
Data can be stored in both REGULAR (equally spaced energy axis) and IRREGULAR mode.
`.npl` files are structurally identical to `.vms` files and are handled by the same parser.

The parser is in
[`src/pynxtools_xps/parsers/vms/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/vms).

## Supported versions

The VAMAS standard does not carry explicit version information. All conforming VAMAS
files are accepted. The format is defined by
[ISO 14976](https://www.iso.org/standard/24269.html).

## Vendor metadata in comment blocks

Most XPS instruments and analysis tools write vendor-specific metadata into the VAMAS
comment lines. The reader parses this additional metadata automatically for the following
sources:

- [Kratos Analytical](https://www.kratos.com/) — see the [Kratos reference page](kratos.md)
- [SPECS GmbH](https://www.specs-group.com/specs/)
- [PHI Electronics](https://www.phi.com/) — same metadata fields as the [PHI reader](phi.md)
- [CasaXPS](http://www.casaxps.com/) — calibrations and peak fitting results

## Standard .vms data

Example data is available in the
[`examples/vms/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms).
The data was measured with and exported from
[SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/).

### REGULAR file format

In the ISO format, there are two types of scan modes: regular and irregular. In regular mode,
the energy axis is uniformly spaced. The axis itself is not saved, but only its start (“abscissa start”) and step energy (“abscissa increment”). The “abscissa label" determines whether the scan is in binding or kinetic energy.

```console
user@box:~$ dataconverter regular.vms eln_data_vms.yaml --reader xps --nxdl NXxps --output regular.vms.nxs
```

### IRREGULAR file format

In irregular mode, the energy axis is stored explicitly as a list of energy values. These values are not required to be uniformly spaced and may vary arbitrarily.

```console
user@box:~$ dataconverter irregular.vms eln_data_vms.yaml --reader xps --nxdl NXxps --output irregular.vms.nxs
```

## Data analysis and peak fitting

`pynxtools-xps` can also extract CasaXPS peak fitting data alongside the raw spectra.
Three files are required:

1. The VAMAS `.vms` file containing the original data and the peak fitting definition in
   the VAMAS comments.
2. The lineshapes exported from CasaXPS as a TXT file.
3. The analysis results (including atomic concentrations) exported from CasaXPS as a CSV file.

Example data is available in the
[`examples/vms/data_analysis/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms/data_analysis).

```console
user@box:~$ dataconverter FeO* eln.yaml --reader xps --nxdl NXxps --output vms_analysis_ref.nxs
```

For details on how to prepare CasaXPS data for NeXus conversion, see
[Explanation > Data processing with CasaXPS](../explanation/data_processing.md).

## Standalone CasaXPS export

`pynxtools-xps` also supports data exported directly from CasaXPS as a TXT file,
without the original VAMAS file.

Example data is available in the
[`examples/vms/txt_export/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms/txt_export).

```console
user@box:~$ dataconverter vms_txt_export.txt eln_data_vms_txt_export.yaml --reader xps --nxdl NXxps --output vms_txt_export.nxs
```

## Further reading

- [Explanation > Parser architecture](../explanation/parser_architecture.md)
