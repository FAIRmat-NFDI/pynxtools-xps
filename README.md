[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/pytest.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/pylint.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/publish.yml/badge.svg)
![](https://img.shields.io/pypi/pyversions/pynxtools-xps)
![](https://img.shields.io/pypi/l/pynxtools-xps)
![](https://img.shields.io/pypi/v/pynxtools-xps)
![](https://coveralls.io/repos/github/FAIRmat-NFDI/pynxtools_xps/badge.svg?branch=main)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1323437.svg)]() -->

# A reader for XPS data

## Installation

It is recommended to use python 3.11 with a dedicated virtual environment for this package.
Learn how to manage [python versions](https://github.com/pyenv/pyenv) and
[virtual environments](https://realpython.com/python-virtual-environments-a-primer/).

This package is a reader plugin for [`pynxtools`](https://github.com/FAIRmat-NFDI/pynxtools) and thus should be installed together with `pynxtools`:


```shell
pip install pynxtools[xps]
```

for the latest development version.

## Purpose
This reader plugin for [`pynxtools`](https://github.com/FAIRmat-NFDI/pynxtools) is used to translate diverse file formats from the scientific community and technology partners
within the field of X-ray photoelectron spectroscopy into a standardized representation using the
[NeXus](https://www.nexusformat.org/) application definition [NXmpes](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html#nxmpes).

## Docs
Extensive documentation of this pynxtools plugin is available [here](https://fairmat-nfdi.github.io/pynxtools-xps/). You can find information about getting started, how-to guides, the supported file formats, how to get involved, and much more there.

## Contact person in FAIRmat for this reader
Lukas Pielsticker