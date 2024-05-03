# VAMAS ISO standard (VMS)

## Basic .vms data

The reader supports VAMAS (.vms) files, the ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)) for X-ray photoelectron spectroscopy. The data can be stored both in REGULAR (i.e, with an equally spaced energy axis) as well as IRREGULAR mode. The data was measured with and exported from [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/).

The reader for the VAMAS format can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/pynxtools_xps/vms).

Example data is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms).

The example conversion for the REGULAR VAMAS file can be run with the following command:

```sh
! dataconverter \
regular.vms \
eln_data_vms.yaml \
--reader xps \
--nxdl NXmpes \
--output regular.vms.nxs \
```

The example conversion for the IRREGULAR VAMAS file can be run with the following command:

```sh
! dataconverter \
irregular.vms \
eln_data_vms.yaml \
--reader xps \
--nxdl NXmpes \
--output irregular.vms.nxs \
```

The example conversion for the .txt export file can be run with the following command:

```sh
! dataconverter \
vms_txt_export.txt \
eln_data_vms_txt_export.yaml \
--reader xps \
--nxdl NXmpes \
--output vms_txt_export.nxs \
```

<!-- ## Data analysis in CasaXPS -->