# Installation guide

## What should you should know before this tutorial?

You will know

- how to install `pynxtools-xps`
- how to install `pynxtools-xps` together with NOMAD

## Setup

It is recommended to use python 3.12 with a dedicated virtual environment for this package. Learn how to manage [python versions](https://github.com/pyenv/pyenv) and [virtual environments](https://realpython.com/python-virtual-environments-a-primer/).

There are many alternatives to managing virtual environments and package dependencies (requirements). We recommend using [`uv`](https://github.com/astral-sh/uv), an extremely fast manager Python package and project manager. In this tutorial, you will find paralleled descriptions, using either `uv` or a more classical approach using `venv` and `pip`.

Start by creating a virtual environment:

=== "uv"
    `uv` is capable of creating a virtual environment and install the required Python version at the same time.

    ```bash
    uv venv --python 3.12
    ```

=== "venv"

    Note that you will need to install the Python version manually beforehand.

    ```bash
    python -m venv .venv
    ```
That command creates a new virtual environment in a directory called `.venv`.

## Installation

Install the latest stable version of this package from PyPI with

=== "uv"

    ```bash
    uv pip install pynxtools-xps
    ```

=== "pip"

    ```bash
    pip install pynxtools-xps
    ```

You can also install the latest _development_ version with

=== "uv"

    ```bash
    uv pip install git+https://github.com/FAIRmat-NFDI/pynxtools-xps.git
    ```

=== "pip"

    ```bash
    pip install git+https://github.com/FAIRmat-NFDI/pynxtools-xps.git
    ```

### How to install `pynxtools-xps` with NOMAD

To use `pynxtools-xps` with NOMAD, simply install it in the same environment as the `nomad-lab` package. NOMAD will recognize `pynxtools-xps` as a plugin automatically. In addition, NOMAD will install a schema for NeXus application definitions.

## Start using `pynxtools-xps`

That's it! You can now use `pynxtools-xps`!
