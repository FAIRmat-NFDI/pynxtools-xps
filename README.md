[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/pytest.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/pylint.yml/badge.svg)
![](https://github.com/FAIRmat-NFDI/pynxtools-xps/actions/workflows/publish.yml/badge.svg)
![](https://img.shields.io/pypi/pyversions/pynxtools-xps)
![](https://img.shields.io/pypi/l/pynxtools-xps)
![](https://img.shields.io/pypi/v/pynxtools-xps)
![](https://coveralls.io/repos/github/FAIRmat-NFDI/pynxtools-xps/badge.svg?branch=main)
[![DOI](https://zenodo.org/badge/759916501.svg)](https://doi.org/10.5281/zenodo.13951684)

# `pynxtools-xps`: A `pynxtools` reader for XPS data

`pynxtools-xps` is a `pynxtools` reader plugin for X-ray photoelectron spectroscopy (XPS) data.

This `pynxtools` plugin was generated with [`cookiecutter`](https://github.com/cookiecutter/cookiecutter) using the [`pynxtools-plugin-template`](https://github.com/FAIRmat-NFDI/pynxtools-plugin-template) template.

## Installation

It is recommended to use python 3.12 with a dedicated virtual environment for this package.
Learn how to manage [python versions](https://github.com/pyenv/pyenv) and
[virtual environments](https://realpython.com/python-virtual-environments-a-primer/).

This package is a reader plugin for [`pynxtools`](https://github.com/FAIRmat-NFDI/pynxtools) and can be installed together with `pynxtools`:

```shell
uv pip install pynxtools[xps]
```

for the latest released version.

## Purpose

This reader plugin for [`pynxtools`](https://github.com/FAIRmat-NFDI/pynxtools) is used to translate diverse file formats from the scientific community and technology partners
within the field of X-ray photoelectron spectroscopy into a standardized representation using the
[NeXus](https://www.nexusformat.org/) application definition [NXxps](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXxps.html#nxxps).

## Docs

More information about this pynxtools plugin is available in the [documentation](https://fairmat-nfdi.github.io/pynxtools-xps/). You will find information about getting started, how-to guides, the supported file formats, how to get involved, and much more there.

## Contact person in FAIRmat for this reader

Lukas Pielsticker

## How to cite this work

Pielsticker, L., Mozumder, R., Dobener, F., Hetaba, W., Rettig, L., & Brockhauser, S. (2025). pynxtools-xps: A pynxtools reader plugin for X-ray photoelectron spectroscopy (XPS) data. Zenodo. https://doi.org/10.5281/zenodo.13951684