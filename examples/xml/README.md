# Example for .sle data

This is an example for SpecsLab 2files .xml files, the XML-based data format from SPECS GmbH. The SpecsLab2 software version that was used to measure this data is v1.6. The example conversion can be run with the following command.

```console
user@box:~$ dataconverter 
```

```console
user@box:~$ dataconverter --params-file params.yaml
```

```sh
! dataconverter \
In-situ_PBTTT_XPS_SPECS_shortened.xml \
eln_data.yaml \
--reader xps \
--nxdl NXmpes \
--output xps_xlm_example.nxs \
```

## Contact person in FAIRmat for this example
Lukas Pielsticker