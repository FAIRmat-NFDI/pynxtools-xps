# Example for .sle data

This is an example for [Phi MultiPak](https://www.phi.com/surface-analysis-equipment/genesis.html#software:multi-pak-data-reduction-software/) .spe (single spectra) and .pro (single spectra) files, which is the propietary format of PHI Electronics used for their VersaProbe 4 instrumens. The Phi MultiPak software version that was used to measure this data is SS 3.3.3.2.1. The example conversion can be run with the following commands.

### For the .spe data (single spectrum):
```console_
user@box:~$ dataconverter SnO2_10nm.spe eln_data_phi.yaml --reader xps --nxdl NXmpes --output SnO2_10nm.spe.nxs
```
### For the .pro data (profiling):
```console_
user@box:~$ dataconverter SnO2_10nm_1.pro eln_data_phi.yaml --reader xps --nxdl NXmpes --output SnO2_10nm_1.pro.nxs
```

## Contact person in FAIRmat for this example
Lukas Pielsticker