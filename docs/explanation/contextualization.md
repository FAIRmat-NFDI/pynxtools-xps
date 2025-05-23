# How to map pieces of information to NeXus

Conceptually, mapping between representations of concepts and instance data is a key tasks in information science. The plugin pynxtools-xps implements this specifically for the file and serialization formats used within the research field of photoelectron spectroscopy (PES).

In pynxtools-xps, the mapping from the vendor format is a two-step process:

1) First, each information piece is parsed from the experiment- and vendor-specific and assigned a name that describes what the reader developer thinks it semantically means. This naming can come from documentation of the original data, existing key-value infrastructure in the data file, or from domain knowledge of the reader developer. All data and metadata items are internally stored as a flat list of dictionaries, with each dictionary containing all information about a single XP spectrum.

2) This list of dicts is then mapped onto either the ([NXmpes](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html) NeXus application definition or its specialization [NXxps](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html)). For this, a JSON config file is used that provides a concept map from the originally assigned keys towards the groups, fields, and attributes in the NeXus standard. Such transformations are configured via the respective files in the [*config*](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/pynxtools_xps/config) directory of pynxtools-xps.

Upon parsing, the XPS reader uses the config file to map the (meta-)data to a *template* which follows the NeXus application definitions. It also takes metadata provided through additional means (i.e., an electronic lab notebook (ELN) file) to fill in missing required and recommended fields and attributes in the application definition that were not provided in the raw data files. It is this *template* variable from which core functions like *convert.py* of the pynxtools write the actual NeXus/HDF5 file. The latter tool is also referred to as the dataconverter of [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools).
