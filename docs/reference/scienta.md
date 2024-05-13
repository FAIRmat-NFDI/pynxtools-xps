# Data from Scienta Omicron instruments

The reader supports reading data exported as .txt from Scienta Omicron [Scienta Omicron](https://www.scientaomicron.com/en/) instruments. <!--,  both exported as .txt as well as ibw ([Igor Pro](https://www.wavemetrics.com/products/igorpro) binary wave) files.--> 

The reader for the Scienta data can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/pynxtools_xps/scienta).

## .txt data

<!-- How is this data structured --> 

Example data for the Scienta TXT reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta/txt).

The example conversion can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.txt eln_data.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.txt.nxs
```

<!-- 
## .ibw ([Igor Pro](https://www.wavemetrics.com/products/igorpro) binary wave) data

How is this data structured

Example data for the Scienta Igor reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta/ibw).

The example conversion can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.ibw eln_data.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.ibw.nxs
```

## Acknowledgments
We thank Dr. Alexei Nefedov from [KIT](https://www.ifg.kit.edu/21_1296.php) for providing the example data set.--> 