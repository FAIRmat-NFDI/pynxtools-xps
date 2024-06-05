# VAMAS ISO standard (VMS)

The reader supports VAMAS (.vms, .npl) files, the ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)) for X-ray photoelectron spectroscopy. The data can be stored both in REGULAR (i.e, with an equally spaced energy axis) as well as IRREGULAR mode. The reader also allows for .npl files which are structured in the same way as .vms files.

The reader for the VAMAS format can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/pynxtools_xps/vms).

Example data is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms). The data was measured with and exported from [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/).


## Standard .vms data

### REGULAR file format

<!-- How is this data structured -->

The example conversion for the REGULAR VAMAS file can be run with the following command:

```console
user@box:~$ dataconverter regular.vms eln_data_vms.yaml --reader xps --nxdl NXmpes --output regular.vms.nxs 
```

### IRREGULAR file format

<!-- How is this data structured -->

The example conversion for the IRREGULAR VAMAS file can be run with the following command:

```console
user@box:~$ dataconverter irregular.vms eln_data_vms.yaml --reader xps --nxdl NXmpes --output irregular.vms.nxs
```

### TXT export from CasaXPS

```pynxtools-xps``` also supports data exported from VAMAS by the [CasaXPS data analysis software](http://www.casaxps.com/) as TXT file. The example conversion for the .txt export file can be run with the following command:

```console
user@box:~$ dataconverter vms_txt_export.txt eln_data_vms_txt_export.yaml --reader xps --nxdl NXmpes --output vms_txt_export.nxs
```

<!-- ## Additional metadata in comment lines (Phi) -->
<!-- ## Data analysis in CasaXPS -->