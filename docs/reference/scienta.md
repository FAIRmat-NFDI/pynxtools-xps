# Data from Scienta Omicron instruments

The reader supports reading data exported as .txt from Scienta Omicron [Scienta Omicron](https://www.scientaomicron.com/en/) instruments, both as .txt as well as [Igor Binary Wave Format](https://www.wavemetrics.com/) files.

The reader for the Scienta data can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/scienta).

## .txt data

<!-- How is this data structured --> 

Example data for the Scienta .txt reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta/txt).

The example conversion for the .txt exports can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.txt eln_data_scienta_txt.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.txt.nxs
```

## .ibw data

<!-- How is this data structured --> 

Example data for the Scienta .ibw reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta/ibw).

The example conversion for the .ibw (Igor binary wave format) data can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.ibw eln_data_scienta_ibw.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.ibw.nxs
```

## Acknowledgments
We thank Dr. Alexei Nefedov from [KIT](https://www.ifg.kit.edu/21_1296.php) for providing the example data sets.