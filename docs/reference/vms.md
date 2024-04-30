# VAMAS ISO standard (VMS)

## Basic .vms data

The reader supports VAMAS (.vms) files, the ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)) for X-ray photoelectron spectroscopy. The data can be stored both in REGULAR (i.e, with an equally spaced energy axis) as well as IRREGULAR mode. The data was measured with and exported from [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/).

Example data for the SLE reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/vms).

The example conversion for the REGULAR file can be run with the following command:

```sh
! dataconverter \
regular.vms \
eln_data_vms.yaml \
--reader xps \
--nxdl NXmpes \
--output vms_regular_example.nxs \
```

The example conversion for the IRREGULAR file can be run with the following command:

```sh
! dataconverter \
irregular.vms \
eln_data_vms.yaml \
--reader xps \
--nxdl NXmpes \
--output vms_irregular_example.nxs \
```

<!-- ## Data analysis in CasaXPS -->