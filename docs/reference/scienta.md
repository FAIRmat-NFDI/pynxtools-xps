# Data from Scienta Omicron instruments

The reader supports reading data exported as .txt from Scienta Omicron [Scienta Omicron](https://www.scientaomicron.com/en/) instruments. <!--,  both as .txt as well as IGOR files.--> 

<!-- How is this data structured --> 

The reader for the Scienta data can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/pynxtools_xps/scienta).

Example data is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/scienta).

The example conversion can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.txt eln_data.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.nxs
```

## Acknowledgments
We thank Dr. Alexei Nefedov from [KIT](https://www.ifg.kit.edu/21_1296.php) for providing the example data set.