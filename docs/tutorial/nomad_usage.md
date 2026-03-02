# Usage in NOMAD

## Who is this tutorial for?

This tutorial is for researchers who want to use `pynxtools-xps` inside
[NOMAD](https://nomad-lab.eu/){:target="_blank" rel="noopener"} to convert, store,
search, and analyze their XPS data in a full research data management system.

## What should you know before this tutorial?

- You should have `pynxtools-xps` installed with the NOMAD extras — see
  [Tutorial > Installation guide](installation.md#how-to-install-pynxtools-xps-with-nomad).
- Basic familiarity with NOMAD uploads is helpful but not required.

## What you will know at the end of this tutorial?

You will understand:

- How to create an XPS example upload in NOMAD
- What files the example contains and what each one does
- How to inspect the converted NeXus file and the parsed data
- How to run the analysis notebooks in NORTH

## Steps

### 1. Open NOMAD central

Go to [https://nomad-lab.eu/prod/v1/gui/](https://nomad-lab.eu/prod/v1/gui/){:target="_blank" rel="noopener"}
and sign in (or register for a free account).

### 2. Create the XPS example upload

Navigate to **Publish → Uploads** and click **+ Upload**.
In the dialog, select **Use an example upload** and choose
**X-ray Photoelectron Spectroscopy (XPS)** from the *NeXus Experiment Examples* category.

NOMAD creates a new upload pre-populated with all example files.

### 3. Explore the example files

The upload contains:

| File | Role |
| ---- | ---- |
| `EX439_S718_Au_in_25_mbar_O2.sle` | Raw SPECS SLE measurement file (Au in 25 mbar O₂, NAP-XPS) |
| `xps.scheme.archive.yaml` | ELN schema — defines all NXxps fields not present in the raw file |
| `eln_conversion.archive.json` | Converts the ELN schema to `eln_data.yaml` |
| `nexus_conversion.archive.json` | Converts SLE + ELN to an `NXxps` NeXus file |
| `E1 XPS data conversion to NeXus.ipynb` | Notebook: conversion pipeline walkthrough |
| `E2 XPS data analysis and fitting.ipynb` | Notebook: peak fitting with `lmfit`, Shirley background |

Open the **FILES** tab to inspect the raw and generated files.
NeXus files (`.nxs`) can be explored interactively with the built-in
[H5Web](https://h5web.panosc.eu/){:target="_blank" rel="noopener"} HDF5 viewer.

### 4. Review the ELN and trigger conversion

Click on `eln_conversion.archive.json` to view and edit the prefilled ELN metadata
(instrument, sample, user information).
NOMAD runs the `ElnYamlConverter` automatically, producing `eln_data.yaml`.

Then open `nexus_conversion.archive.json`.
NOMAD invokes the XPS reader and produces the standardized NeXus file
(`Au_25_mbar_O2_no_align.nxs`) according to the `NXxps` application definition.

The full pipeline in the example:

```text
EX439_S718_Au_in_25_mbar_O2.sle  +  eln_data.yaml
        │
        ▼  NexusDataConverter (reader: xps, nxdl: NXxps)
        │
  Au_25_mbar_O2_no_align.nxs
        │
        ├──► NexusParser → NOMAD Metainfo (DATA tab, search)
        ├──► H5Web viewer (FILES tab)
        └──► NORTH Jupyter analysis (E1 · E2 notebooks)
```

### 5. Inspect the parsed data

Open the **DATA** tab of the upload.
The `NexusParser` has extracted the NXxps structure into NOMAD Metainfo, making
instrument parameters, sample details, and spectral data browsable and searchable.

### 6. Run the analysis notebooks in NORTH

If NOMAD is deployed with NORTH, open **Analyze → NORTH** and launch a JupyterLab
container.
Navigate to the upload directory and open:

- **E1 XPS data conversion to NeXus.ipynb** — step-by-step notebook showing how the
  conversion is orchestrated.
- **E2 XPS data analysis and fitting.ipynb** — loads the NeXus file, performs Shirley
  background subtraction, and fits peaks with `lmfit`.

## Using your own data

The example uses a SPECS SLE file, but the same workflow applies to any format supported
by `pynxtools-xps` (VMS, TXT, IBW, XML, …).
Replace the raw measurement file in the upload, update the `input_files` field in
`nexus_conversion.archive.json`, and adjust the ELN metadata to match your instrument.

## Further reading

- [Explanation > NOMAD integration](../explanation/nomad_integration.md) — architectural
  overview of how `pynxtools-xps` integrates into NOMAD
- [Tutorial > Standalone data conversion](standalone.md) — using the reader without NOMAD
- [Explanation > NXmpes and NXxps](../explanation/appdefs.md) — the application definitions
  used in the conversion
