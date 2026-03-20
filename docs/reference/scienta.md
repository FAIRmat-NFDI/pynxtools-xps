# Data from Scienta Omicron instruments

The reader supports data exported from
[Scienta Omicron](https://www.scientaomicron.com/en/) instruments in three formats:
plain text (`.txt`),
[Igor Binary Wave Format](https://www.wavemetrics.com/) (`.ibw`), and
HDF5 (`.h5` / `.hdf5`).

Three dedicated parsers live in
[`src/pynxtools_xps/parsers/scienta/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/scienta):

| Parser class | Source file | Format |
| ------------ | ----------- | ------ |
| `ScientaTXTParser` | `txt_parser.py` | Plain-text SES export |
| `ScientaIgorParser` | `igor_parser.py` | Igor Binary Wave (IBW) |
| `ScientaHDF5Parser` | `hdf5_parser.py` | HDF5 written by Scienta/PEAK |

## Supported versions

| Format | Extension | Software | Supported versions |
| ------ | --------- | -------- | ------------------ |
| Scienta plain text export | `.txt` | Scienta software | any |
| Igor Binary Wave | `.ibw` | WaveMetrics Igor Pro | any |
| HDF5 | `.h5`, `.hdf5` | Scienta / PEAK software | any |

## .txt data

Example data is available in the
[`examples/scienta/txt/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta/txt).

```console
user@box:~$ dataconverter Cu-HHTP_*.txt eln_data_scienta_txt.yaml --reader xps --nxdl NXxps --output Cu-HHTP.txt.nxs
```

## .ibw data

Example data is available in the
[`examples/scienta/ibw/` directory](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta/ibw).

```console
user@box:~$ dataconverter Cu-HHTP_*.ibw eln_data_scienta_ibw.yaml --reader xps --nxdl NXxps --output Cu-HHTP.ibw.nxs
```

## .h5 / .hdf5 data

HDF5 files written by Scienta or PEAK software are supported via extensions `.h5` and `.hdf5`.

```console
user@box:~$ dataconverter measurement.h5 eln_data_scienta_hdf5.yaml --reader xps --nxdl NXxps --output measurement.h5.nxs
```

## Acknowledgements

We thank Dr. Alexei Nefedov from [KIT](https://www.ifg.kit.edu/21_1296.php) for providing the example IBW and TXT data sets.

## Further reading

- [Explanation > Parser architecture](../explanation/parser_architecture.md)
