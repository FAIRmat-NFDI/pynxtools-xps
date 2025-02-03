# VAMAS ISO standard (VMS)

The reader supports VAMAS (.vms, .npl) files, the ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)) for X-ray photoelectron spectroscopy. The data can be stored both in REGULAR (i.e, with an equally spaced energy axis) as well as IRREGULAR mode. The reader also allows for .npl files which are structured in the same way as .vms files.

Note that most vendors and analysis software tend to write metadata from their instruments into the comment lines of the VAMAS format. Currently, the VAMAS reader supports parsing of metadata from VAMAS format for the following vendors and software solutions:

- [Kratos Analytical Ltd](https://www.kratos.com/)
- [Specs GmbH](https://www.specs-group.com/specs/)
- [Phi Electronics](https://www.phi.com/): same metadata as in the PHI reader
- [CasaXPS](http://www.casaxps.com/): calibrations and peak fitting

The reader for the VAMAS format can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/vms).

## Standard .vms data

Example data is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms). The data was measured with and exported from [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/).

### REGULAR file format

<!-- How is this data structured -->

The example conversion for the REGULAR VAMAS file can be run with the following command:

```console
user@box:~$ dataconverter regular.vms eln_data_vms.yaml --reader xps --nxdl NXxps --output regular.vms.nxs 
```

### IRREGULAR file format

<!-- How is this data structured -->

The example conversion for the IRREGULAR VAMAS file can be run with the following command:

```console
user@box:~$ dataconverter irregular.vms eln_data_vms.yaml --reader xps --nxdl NXxps --output irregular.vms.nxs
```

## Data analysis and peak fitting

```pynxtools-xps``` also supports extracting data and the description of the data analysis (i.e., peak fitting)
by the [CasaXPS data analysis software](http://www.casaxps.com/). Three files are needed for the example conversion:

1) The VAMAS (.vms) file containing the original (meta)data and the definition of the peak fitting in the VAMAS
comments
2) The lineshapes of the measurement data as well as the peak fitting, exported from CasaXPS as a TXT file.
3) The analysis results (incl. the atomic concentrations), exported from CasaXPS as a CSV file.

Example data is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms/vms_analysis).

The example conversion for the .txt export file can be run with the following command:

```console
user@box:~$ dataconverter FeO* eln.yaml --reader $READER --nxdl $NXDL --output vms_analysis_ref.nxs
```

You can learn much more about how to prepare the data in CasaXPS for NeXus conversion [here](../explanation/data_processing.md).

## Standalone export from CasaXPS

```pynxtools-xps``` also supports data exported from CasaXPS as TXT file by itself.

Example data is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms/txt_export).

The example conversion for the .txt export file can be run with the following command:

```console
user@box:~$ dataconverter vms_txt_export.txt eln_data_vms_txt_export.yaml --reader xps --nxdl NXxps --output vms_txt_export.nxs
```
