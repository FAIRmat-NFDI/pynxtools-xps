# Development guide

This tutorial will guide you through on how to set up a working environment for developing `pynxtools-xps`.

## What should you should know before this tutorial?

- You should read the [guide on getting started with `pynxtools`](https://fairmat-nfdi.github.io/pynxtools/getting-started.html).
- You should read the [installation tutorial](installation.md).

## What you will know at the end of this tutorial?

You will know

- how to setup your environment for developing `pynxtools-xps`
- how to make changes to the software
- how to test the software
- how to contribute on GitHub
- how to use `pynxtools-xps` as a NOMAD plugin

## Contributing

??? info "Structure of the `pynxtools-xps` repository"
    The software tools are located inside `src/pynxtools_xps`. They are shipped with unit tests located in `tests`.

### Setup

It is recommended to use python 3.12 with a dedicated virtual environment for this package. Learn how to manage [python versions](https://github.com/pyenv/pyenv) and [virtual environments](https://realpython.com/python-virtual-environments-a-primer/). We recommend using [`uv`](https://github.com/astral-sh/uv), an extremely fast manager Python package and project manager. In this tutorial, you will find paralleled descriptions, using either `uv` or a more classical approach using `venv` and `pip`.

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

That command creates a new virtual environment in a directory called .venv.

### Development installation

We start by cloning the repository:

```console
git clone https://github.com/FAIRmat-NFDI/pynxtools-xps.git \\
    --branch main \\
    --recursive pynxtools-xps
cd pynxtools-xps
```

Next, we install the package in editable mode (together with its dependencies):

=== "uv"

    ```bash
    uv pip install -e ".[dev]"
    ```

=== "pip"

    Note that you will need to install the Python version manually beforehand.

    ```bash
    pip install --upgrade pip
    pip install -e ".[dev]"
    ```

### Linting and formatting

We are using ruff and mypy for linting, formatting, and type checking. It is recommended to use the [pre-commit hook](https://pre-commit.com/#intro) available for ruff which formats the code and checks the linting before actually making an actual Git commit.

Install the precommit by running

```console
pre-commit install
```

from the root of this repository.

### Testing

There exist unit tests for the software written in [pytest](https://docs.pytest.org/en/stable/) which can be used as follows:

```console
pytest -sv tests
```

### Editing the documentation

We are using [`mkdocs](https://www.mkdocs.org/) for the documentation. If you edit the documentation, you can build it locally. For this, you need to install an additional set of dependencies:

=== "uv"

    ```bash
    uv pip install -e ".[docs]"
    ```

=== "pip"

    ```bash
    pip install -e ".[docs]"
    ```
You can then serve the documentation locally by running

```console
mkdocs serve
```

### Contributing to the package on GitHub

Once you are happy with the changes, please commit them on a separate branch and create a pull request on GitHub. We run a number of GitHub actions that check the correct linting, run the tests in an isolated environment, and build the documentation. Once these pass and a peer review of the code has occurred, your code will be accepted.

## Developing `pynxtools-xps` as a NOMAD plugin

If you plan to contribute to the NOMAD plugin functionality of pynxtools-xps, it often makes sense to use the NOMAD development environment called `nomad-distro-dev`. You can learn more in the [NOMAD documentation](https://nomad-lab.eu/prod/v1/staging/docs/howto/develop/setup.html#nomad-distro-dev-development-environment-for-the-core-nomad-package-and-nomad-plugins).

## Troubleshooting

If you face any issues with the tool or when setting up the development environment, please create a new [Github Issue](https://github.com/FAIRmat-NFDI/pynxtools-xps/issues/new?template=bug.yaml).
