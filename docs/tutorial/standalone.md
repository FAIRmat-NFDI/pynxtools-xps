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
```sh
 ! dataconverter \
--reader xps \
--nxdl NXmpes \
--input-file $<xps-file path> \
--input-file $<eln-file path> \
--output <output-file path>.nxs
```
Note that none of the supported file format have data/values for all required and recommended fields and attributes in NXmpes. In order for the validation step of the XPS reader to pass, you need to provide an ELN file that contains the missing values. Example raw and converted data can be found in  [*pynxtools_xps/examples*](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples).

TODO: add more steps! <!--[The Jupyter notebook is available here](https://github.com/FAIRmat-NFDI/pynxtools-em/blob/main/examples/HowToUseTutorial.ipynb) TODO!-->

**Congrats! You now have a FAIR NeXus file!**

The above-mentioned parsing is also integrated into the NOMAD research data management system.
Feel invited to try out the respective tutorial [here]((tutorial/nomad.md)
