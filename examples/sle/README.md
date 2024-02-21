# Example for .sle data

This is an example for the [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) .sle files, which is the propietary format of SPECS GmbH. The SpecsLabProdigy software version is v4.63.1. The example conversion can be run with the following command.

```console
user@box:~$ dataconverter --params-file params.yaml
```

Not that the `params.yaml` file contains the `remove_align` keyword which is special for the SLE parser. It allows removal of alignment spectra that were taken during the experiment. For this example, it considerably speeds up the conversion.

## Contact person in FAIRmat for this example
Lukas Pielsticker