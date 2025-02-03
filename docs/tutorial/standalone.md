# Convert X-ray spectroscopy data and metadata to NeXus

## Who is this tutorial for?

This document is for people who want to use this reader as a standalone standardize their research data by converting these
into a NeXus standardized format.

## What should you should know before this tutorial?

- You should have a basic understanding of [FAIRmat NeXus](https://github.com/FAIRmat/nexus_definitions) and [pynxtools](https://github.com/FAIRmat/pynxtools)
- You should have a basic understanding of using Python and Jupyter notebooks via [JupyterLab](https://jupyter.org)

## What you will know at the end of this tutorial?

You will have a basic understanding how to use pynxtools-xps for converting your XPS data to a NeXus/HDF5 file.

## Steps

### Installation

See here for how to install pynxtools together with the XPS reader plugin.

### Running the reader from the command line

An example script to run the XPS reader in `pynxtools`:

```console
user@box:~$ dataconverter $<xps-file path> $<xps-file path> $<eln-file path> --reader xps --nxdl NXxps --output <output-file path>.nxs
```

Note that none of the supported file format have data/values for all required and recommended fields and attributes in ``NXxps``. In order for the validation step of the XPS reader to pass, you need to provide an ELN file that contains the missing values.

### Examples

You can find examples how to use `pynxtools-xps` for your XPS research data pipeline in [`src/pynxtools-xps/nomad/examples`](../../src/pynxtools_xps/nomad/examples/). These are designed for working with [`NOMAD`](https://nomad-lab.eu/) and its [`NOMAD Remote Tools Hub (NORTH)`](https://nomad-lab.eu/prod/v1/gui/analyze/north). Feel invited to try out the respective tutorial [here](tutorial/nomad.md).

There are also small example files with raw and converted data for using the `pynxtools` dataconverter with the `mpes` reader and the `NXmpes` application definition in the [`examples`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/) folder.

For this tutorial, we will work with the example data for the VAMAS reader (see [here](../reference/vms.md)). You can run the conversion as

```shell
dataconverter \\
    --reader xps \\
    --nxdl NXmpes \\
    regular.vms \\
    eln_data_vms.yaml \\
    -c  config_file.json \\
    --output regular.vms.nxs 
```

TODO: add more steps! <!--[The Jupyter notebook is available here](https://github.com/FAIRmat-NFDI/pynxtools-em/blob/main/examples/HowToUseTutorial.ipynb) TODO!-->

**Congrats! You now have a FAIR NeXus file!**
