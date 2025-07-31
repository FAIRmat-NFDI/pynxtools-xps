# XPS example

## Introduction

This example demonstrates how NOMAD can convert, standardize, and store XPS data. It shows the generation of a NeXus file according to the [`NXxps`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html#nxxps) application definition and a successive analysis of an example data set.

## Viewing uploaded data

Below, you find an overview of your uploaded data.
Click on the `> /` button to get a list of your file or select **FILES** from the top menu of this upload.

The general structure of this example is the following:

- The starting point is the `xps.scheme.archive.yaml` file. This file defines a schema for the electronic lab notebook (ELN) for XPS data. The ELN follows the general structure of NOMAD ELN templates. You can learn about NOMAD ELNs in the [documentation](https://nomad-lab.eu/prod/v1/staging/docs/howto/manage/eln.html). This schema contains all the NeXus concepts that are defined in the `NXxps` application definition, but are not provided through the raw experimental file. Hence, these are the information that the user should enter.
- The file `eln_conversion.archive.yaml` has two functionalities: 1) it implements the schema defined in `xps.scheme.archive.yaml` and prefills it with some values for this particular example, and 2) it uses the `ElnYamlConverter` from the NOMAD plugin [`pynxtools`](https://github.com/FAIRmat-NFDI/pynxtools) to convert the NOMAD ELN into a YAML file that can be used to populate the NeXus file.
- In the next step, the `nexus_conversion.archive.json` is used to generate a NeXus file according to the `NXxps` application definition by calling the `NexusDataConverter` from `pynxtools`. It uses as input file the raw measurement file (`EX439_S718_Au_in_25_mbar_O2.sle`) as well as the newly generated `eln_data.yaml` to create the NeXus file `Au_25_mbar_O2_no_align.nxs` (by calling the XPS reader from the package `pynxtools-xps`). You may view your supplied or generated NeXus files here with the H5Web viewer. To do so open the **FILES** tab and just select the `.nxs` file.
- Finally, the `NexusParser` from `pynxtools` is called to read in the data in the NeXus file into the NOMAD metainfo. You can explore the parsed data in the **DATA** tab.

You may add your own files to the upload or experiment with the pre-existing electronic lab notebook (ELN) example.

## Analyzing the data

The example works through the use of NOMAD remote tools hub (NORTH) containers. In addition to converting and using the uploaded XPS data in NOMAD, an analysis container can be started.  To start an analysis, note your upload id (which you find above this explanation) and select **ANALYZE** from the top menu, then **NOMAD Remote Tools Hub**. In the appearing list you'll find the `xps` container, click on it and click **LAUNCH**. After a few moments a new tab will open which displays a Jupyter environment providing the required analysis tools. To find the examples navigate to uploads inside the JupyterHub browser and select the folder with your noted upload id. There you'll find the example `ipynb` notebook for data analysis (E2).
Double-clicking the notebook will open the example in the Jupyter main window. Execute notebook to perform peak fitting analysis on the data contained in this container using the Python tool `lmfit`.

## Where to go from here?

If you're interested in using this pipeline and NOMAD in general you'll find support at [FAIRmat](https://www.fairmat-nfdi.eu/fairmat/).

For questions regarding the experiment or this specific example contact [Lukas Pielsticker](https://www.fairmat-nfdi.eu/fairmat/about-fairmat/team-fairmat) or the rest of the [FAIRmat team](https://www.fairmat-nfdi.eu/fairmat/about-fairmat/team-fairmat).
