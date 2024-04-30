# Data from Scienta Omicron instruments

The reader supports [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) files, which is the propietary format of SPECS GmbH. Currently, the following file extensions are supported:
- .sle: [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) file (software version: v1.6, >v4)
- .xml: SpecsLab 2files, XML format from SPECS GmbH (software version: v4.63 tested, other versions also work)
- .xy: SpecsLabProdigy export format in XY format (including all export settings)

<!-- How is this data structured --> 

Example data for the SLE reader is available [here](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs).

The example conversion can be run with the following command.
```console
user@box:~$ dataconverter --params-file params.yaml
```

Note that the `params.yaml` file contains the `remove_align` keyword which is special for the SLE parser. It allows removal of alignment spectra that were taken during the experiment. For this example, it considerably speeds up the conversion.