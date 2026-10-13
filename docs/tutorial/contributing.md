# Development guide

This tutorial will guide you through on how to set up a working environment for developing `pynxtools-xps`.

## What should you know before this tutorial?

- You should read the [guide on getting started with `pynxtools`](https://fairmat-nfdi.github.io/pynxtools/getting-started.html).
- You should read the [installation tutorial](installation.md).

## What you will know at the end of this tutorial?

You will know

- how to setup your environment for developing `pynxtools-xps`
- which branch to use for which kind of change
- how to make changes to the software
- how to test the software
- how to contribute on GitHub
- how a new version is released
- how to use `pynxtools-xps` as a NOMAD plugin

## Contributing

??? info "Structure of the `pynxtools-xps` repository"
    The software tools are located inside `src/pynxtools_xps`. They are shipped with unit tests located in `tests`.

### Branching model

`pynxtools-xps` follows the same branching model as the core NOMAD packages.

| Branch | Role |
|---|---|
| `develop` | Default branch. All regular development PRs target this branch. |
| `main` | Released versions only. Only `release/*` branches are merged here, so its history is the release history. |
| `release/vX.Y.Z` | Temporary release branch, created from `develop` and merged into `main`. The version may carry a pre-release suffix, e.g. `release/vX.Y.Z-beta` or `release/vX.Y.Z-rc1`. |

| Activity | Branch |
|---|---|
| Feature development | `develop` |
| PR review and integration | `develop` |
| Selecting changes for a release | `release/vX.Y.Z` |
| Release-specific bug fixes | `release/vX.Y.Z` |
| Published release history | `main` |
| Porting release fixes back to ongoing development | Merge `release/vX.Y.Z` back into `develop` |

The `Check if source is release/*` workflow
([`.github/workflows/check-source-is-release.yml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/develop/.github/workflows/check-source-is-release.yml))
allows PRs into `main` only from `release/*` branches. PRs from any other branch, **including `develop`**,
fail this check. It is enforced by marking it as a required status check in the branch protection rule
for `main`.

### Setup

It is recommended to use python 3.12 with a dedicated virtual environment for this package. Learn how to manage [python versions](https://github.com/pyenv/pyenv) and [virtual environments](https://realpython.com/python-virtual-environments-a-primer/). We recommend using [`uv`](https://github.com/astral-sh/uv), an extremely fast Python package and project manager. In this tutorial, you will find parallel descriptions, using either `uv` or a more classical approach using `venv` and `pip`.

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
    --branch develop \\
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

We are using ruff and mypy for linting, formatting, and type checking. It is recommended to use the [pre-commit hook](https://pre-commit.com/#intro) available for ruff which formats the code and checks the linting before making a Git commit.

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

Once you are happy with the changes, commit them on a separate branch created from `develop` and open a pull request against `develop`:

```console
git switch develop && git pull
git switch -c <short-descriptive-name>
# commit your changes, then
git push -u origin <short-descriptive-name>
```

Never open a pull request against `main`: it receives only `release/*` branches. We run a number of GitHub actions that check the correct linting, run the tests in an isolated environment, and build the documentation. Once these pass and a peer review of the code has occurred, your code will be accepted.

## Release workflow

1. Create the release branch from `develop`:

    ```console
    git switch develop && git pull
    git switch -c release/vX.Y.Z
    git push -u origin release/vX.Y.Z
    ```

1. Commit only release-specific changes on this branch (version bumps, changelog, last-minute bug fixes). New features keep going into `develop`.
1. Bump `version` in `CITATION.cff` to `X.Y.Z`, including any pre-release suffix (e.g. `0.1.2-rc1` for the tag `v0.1.2-rc1`). The publish workflow fails if it differs from the tag.
1. Open a PR from `release/vX.Y.Z` into `main`. Do not merge yet!
1. Create a GitHub release with the tag `vX.Y.Z` and select `release/vX.Y.Z` as its target branch. The publish workflow refuses to upload to PyPI if the release does not target `release/vX.Y.Z` or if the tagged commit is not on that branch.
1. Once CI passes and the release is published to PyPI, merge the PR.
1. If you fixed anything directly on the release branch, open a PR from `release/vX.Y.Z` into `develop` so the fixes are not lost in ongoing development.
1. Delete the `release/vX.Y.Z` branch.

## Developing `pynxtools-xps` as a NOMAD plugin

If you plan to contribute to the NOMAD plugin functionality of pynxtools-xps, it often makes sense to use the NOMAD development environment called `nomad-distro-dev`. You can learn more in the [NOMAD documentation](https://nomad-lab.eu/prod/v1/staging/docs/howto/develop/setup.html#nomad-distro-dev-development-environment-for-the-core-nomad-package-and-nomad-plugins).

## Troubleshooting

If you face any issues with the tool or when setting up the development environment, please create a new [Github Issue](https://github.com/FAIRmat-NFDI/pynxtools-xps/issues/new?template=bug.yaml).
