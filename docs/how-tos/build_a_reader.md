# How to build your own reader

Your current data is not supported yet? Don't worry, the following how-to will guide you how to write a reader your own data.

## pynxtools-xps supports your format, but some groups and fields are different

Good! The basic functionality to read your data is already in place. Before you start writing your own reader, consider two options:
1) You can modify the default [config files](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/config).
2) Consider opening a [pull request on the GitHub repository](https://github.com/FAIRmat-NFDI/pynxtools-xps/pulls) modifying the existing reader.

## You have a completely new data format

You will have to write a new sub-reader inside pynxtools-xps. There are multiple steps to get started:

### Development install

You should start with an devlopment install of the package with its dependencies:

```shell
git clone https://github.com/FAIRmat-NFDI/pynxtools-xps.git \\
    --branch main \\
    --recursive pynxtools_xps
cd pynxtools_xps
python -m pip install --upgrade pip
python -m pip install -e .
python -m pip install -e ".[dev,consistency_with_pynxtools]"
```

There is also a [pre-commit hook](https://pre-commit.com/#intro) available
which formats the code and checks the linting before actually commiting.
It can be installed with
```shell
pre-commit install
```
from the root of this repository.

### Design strategy
The development process is modular so that new parsers can be added. The design logic is the following:
1. First, [`XpsDataFileParser`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/file_parser.py#L36) selects the proper parser based on the file extensions of the provided files. It then calls a sub-parser that can read files with such extensions and calls the `parse_file` function of that reader. In addition, it selects a proper config file from
the `config` subfolder.
2. Afterwards, the NXmpes NXDL template is filled with the data in `XpsDataFileParser` using the [`config`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/config) file. Data that is not in the given main files can be added through the ELN file (and must be added for required fields in NXmpes).

### Write your reader
TODO!

### Test the software
There exists a basic test framework written in [pytest](https://docs.pytest.org/en/stable/) which can be used as follows:
```shell
python -m pytest -sv tests
```
You should add test data and add your reader to the `test_params` in the `test_reader.py` script.

# Further details

[NXmpes](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html)

<!-- [NXxps](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html) -->