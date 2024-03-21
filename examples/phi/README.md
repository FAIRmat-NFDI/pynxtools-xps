# Example for .sle data

This is an example for [Phi MultiPak](https://www.phi.com/surface-analysis-equipment/genesis.html#software:multi-pak-data-reduction-software/) .spe (single spectra) and .pro (single spectra) files, which is the propietary format of PHI Electronics. The Phi MultiPak software version that was used to measure this data is SS 3.3.3.2.1. The example conversion can be run with the following commands.

### For the .spe data:
```console_
user@box:~$ dataconverter test_data.spe eln_data.yaml --reader xps --nxdl NXmpes
```
### For the .pro data:
```console_
user@box:~$ dataconverter test_data.pro eln_data.yaml --reader xps --nxdl NXmpes
```

## Contact person in FAIRmat for this example
Lukas Pielsticker