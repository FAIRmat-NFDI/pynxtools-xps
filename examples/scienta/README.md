# Example for .txt data exported from Scienta spectrometer

This is an example for parsing the data exported by Scienta Omicron
[Scienta Omicron](https://www.scientaomicron.com/en/) instruments. We thank Dr. Alexei Nefedov from [KIT](https://www.ifg.kit.edu/21_1296.php) for providing this example data set.

The example conversion for the .txt exports can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.txt eln_data.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.txt.nxs
```

The example conversion for the .ibw (Igor binary wave format) data can be run with the following command.

```console
user@box:~$ dataconverter Cu-HHTP_*.ibw eln_data.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.ibw.nxs
```

## Contact person in FAIRmat for this example
Lukas Pielsticker
