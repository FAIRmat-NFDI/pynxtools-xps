# Data from PHI instruments

The reader supports `.spe` (single spectra) and `.pro` (sputter profile / external
parameter scan) files from [PHI Electronics](https://www.phi.com/) VersaProbe 4
instruments, written by [PHI MultiPak](https://www.phi.com/surface-analysis-equipment/genesis.html#software:multi-pak-data-reduction-software/)
software.

The parser is in
[`src/pynxtools_xps/parsers/phi/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/phi).

## Supported versions

| Format | Extension | Software | Supported versions |
| ------ | --------- | -------- | ------------------ |
| PHI MultiPak single spectrum | `.spe` | PHI MultiPak | any (SS 3.3.3.2.1 tested) |
| PHI MultiPak profile | `.pro` | PHI MultiPak | any (SS 3.3.3.2.1 tested) |

## .spe data (single spectrum)

Example data is available in the
[`examples/phi/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/phi).

```console
user@box:~$ dataconverter SnO2_10nm.spe eln_data_phi.yaml --reader xps --nxdl NXxps --output SnO2_10nm.spe.nxs
```

## .pro data (profiling)

Example data is available in the
[`examples/phi/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/phi).

```console
user@box:~$ dataconverter SnO2_10nm_1.pro eln_data_phi.yaml --reader xps --nxdl NXxps --output SnO2_10nm_1.pro.nxs
```

## Acknowledgements

We thank Sebastian Benz and Dr. Joachim Sann from
[Justus-Liebig-Universität Gießen](https://www.uni-giessen.de/de)
for providing these example data sets.

## Further reading

- [Explanation > Parser architecture](../explanation/parser_architecture.md)
