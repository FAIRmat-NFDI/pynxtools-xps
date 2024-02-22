[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/pytest.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/pylint.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/publish.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools_xps/actions/workflows/pytest.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools_xps/actions/workflows/pylint.yml/badge.svg)
![](https://coveralls.io/repos/github/FAIRmat-NFDI/pynxtools_xps/badge.svg?branch=master)

# A reader for XPS data

# Installation

It is recommended to use python 3.11 with a dedicated virtual environment for this package.
Learn how to manage [python versions](https://github.com/pyenv/pyenv) and
[virtual environments](https://realpython.com/python-virtual-environments-a-primer/).

Install this package with

```shell
pip install git+https://github.com/FAIRmat-NFDI/pynxtools-xps.git
```

for the latest development version.


# Purpose
This reader plugin for [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools) is used to translate diverse file formats from the scientific community and technology partners
within the field of X-ray photoelectron spectroscopy into a standardized representation using the
[NeXus](https://www.nexusformat.org/) application definition [NXmpes](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html#nxmpes).

## Supported file formats
The reader decides which parser to use based on the file extension of the files provided. For the main XPS files, the following file extensions are supported:
- .sle: [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) files, propietary format of SPECS GmbH (1 and v4)
- .xml: SpecsLab 2files, XML format from SPECS GmbH (v1.6)
- .vms: VAMAS files, ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)), both in regular and irregular format
- .xy: SpecsLabProdigy export format in XY format (including all export settings)
- .txt:
  - exported by [Scienta Omicron](https://scientaomicron.com/en) instruments
  - exported by [CasaXPS](https://www.casaxps.com/) analysis software

We are continously working on adding parsers for other data formats and technology partners. If you would like to implement a parser for your data, feel free to get in contact.

# Getting started
An example script to run the XPS reader in pynxtools:
```sh
 ! dataconverter \
--reader xps \
--nxdl NXmpes \
--input-file $<xps-file path> \
--input-file $<eln-file path> \
--output <output-file path>.nxs
```
Note that none of the supported file format have data/values for all required and recommended fields and attributes in NXmpes. In order for the validation step of the XPS reader to pass, you need to provide an ELN file that contains the missing values. Example raw and converted data can be found in  [*pynxtools_xps/examples*](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples).


# Contributing

## Development install

Install the package with its dependencies:

```shell
git clone https://github.com/FAIRmat-NFDI/pynxtools-xps.git \\
    --branch master \\
    --recursive pynxtools_xps
cd pynxtools_xps
python -m pip install --upgrade pip
python -m pip install -e .
python -m pip install -e ".[dev]"
```

There is also a [pre-commit hook](https://pre-commit.com/#intro) available
which formats the code and checks the linting before actually commiting.
It can be installed with
```shell
pre-commit install
```
from the root of this repository.

## Development Notes
The development process is modular so that new parsers can be added. The design logic is the following:
1. First, [`XpsDataFileParser`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/pynxtools_xps/file_parser.py#L36) selects the proper parser based on the file extensions of the provided files. It then calls a sub-parser that can read files with such extensions and calls the `parse_file` function of that reader. In addition, it selects a proper config file from
the `config` subfolder.
2. Afterwards, the NXmpes nxdl template is filled with the data in `XpsDataFileParser` using the [`config`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/pynxtools_xps/config) file. Data that is not in the given main files can be added through the ELN file (and must be added for required fields in NXmpes).

## Test this software

Especially relevant for developers, there exists a basic test framework written in
[pytest](https://docs.pytest.org/en/stable/) which can be used as follows:

```shell
python -m pytest -sv tests
```

## Contact person in FAIRmat for this reader
Lukas Pielsticker
