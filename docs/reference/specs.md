# Data from SPECS instruments

The reader supports [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) files, which is the propietary format of SPECS GmbH. Currently, the following file extensions are supported:

- .sle: [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) file software version: v1.6, >v4)
- .xml: SpecsLab 2files, XML format from SPECS GmbH (software version: v4.63 tested, other versions also work)
- .xy: SpecsLabProdigy export format in XY format (including all export settings)

The readers for the SPECS data can be found [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/specs).

## .sle data

<!-- How is this data structured --> 

Example data for the SLE reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/sle).

The example conversion can be run with the following command.
```console
user@box:~$ dataconverter --params-file params.yaml
```

Note that the `params.yaml` file contains the `remove_align` keyword which is special for the SLE parser. It allows removal of alignment spectra that were taken during the experiment. For this example, it considerably speeds up the conversion

## .xml data

<!-- How is this data structured --> 

Example data for the SPECS XML reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/xml).

The example conversion can be run with the following command.
```console
user@box:~$ dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_xml.yaml --reader xps --nxdl NXmpes --output In-situ_PBTTT.nxs
```

## .xy data

<!-- How is this data structured --> 

Example data for the SPECS XY reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/xy).

The example conversion can be run with the following command.
```console
user@box:~$ dataconverter MgFe2O4.xy eln_data_xy.yaml --reader xps --nxdl NXmpes --output MgFe2O4.nxs
``` 