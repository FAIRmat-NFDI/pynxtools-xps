# Example for .vms data

This is an example for VAMAS (.vms) files, the ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)) for X-ray photoelectron spectroscopy. The data has been stored both in REGULAR (i.e, with an equally spaced energy axis) as well as IRREGULAR mode. The data was measured with and exported from[SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/).

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

## Contact person in FAIRmat for this example
Lukas Pielsticker