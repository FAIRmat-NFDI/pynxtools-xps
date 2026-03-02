# NOMAD integration

`pynxtools-xps` functions as both a standalone data converter and a
[NOMAD](https://nomad-lab.eu/){:target="_blank" rel="noopener"} plugin.
When installed alongside `nomad-lab`, it extends NOMAD with XPS-specific capabilities at every layer of the platform: data ingestion, standardization, search, and interactive analysis.

## Plugin architecture

The integration is two-layered.
[`pynxtools`](https://github.com/FAIRmat-NFDI/pynxtools){:target="_blank" rel="noopener"}
— a required dependency of `pynxtools-xps` — provides the generic NeXus-to-NOMAD bridge. `pynxtools-xps` then plugs into that bridge via Python entry points:

```text
NOMAD
├── pynxtools  (Schema Package · NexusParser · Normalization · Search)
│   └── pynxtools-xps  (`XPSReader` — registered via pynxtools reader entry point)
└── pynxtools-xps  (Example Upload — registered via nomad plugin entry point)
```

**`pynxtools` provides:**

- **Schema Package** — The MPES-related NeXus application definitions (including `NXmpes` and `NXxps`) are  integrated into NOMAD's data schema (Metainfo), making NeXus quantities searchable and interoperable with other NOMAD data.
- **Dataconverter in the GUI** — the `NexusDataConverter` exposes `pynxtools`' conversion pipeline inside the NOMAD interface.
- **NexusParser** — reads NeXus HDF5 files and populates NOMAD Metainfo object instances.
- **Normalization** — automatic unit handling, reference linking (sample and instrument identifiers), and population of derived quantities for search and visualization.
- **Search Application** — an Elasticsearch-powered dashboard for filtering uploaded data by technique, timestamp, and domain-specific quantities.

**`pynxtools-xps` adds:**

- The `xps` reader, registered under the `pynxtools.reader` entry point, so NOMAD's dataconverter GUI can convert raw XPS files (SLE, VMS, IBW, TXT, …) to `NXxps` or   `NXmpes` without any command-line interaction.
- An example upload registered under the `nomad.plugin` entry point, which ships a complete end-to-end XPS workflow as a ready-to-run NOMAD upload.

## Installation

Install `pynxtools-xps` in the same Python environment as `nomad-lab`:

```bash
pip install pynxtools-xps[nomad]
```

NOMAD discovers plugins automatically via Python entry points — no additional configuration is needed. See [Tutorial > Installation guide](../tutorial/installation.md#how-to-install-pynxtools-xps-with-nomad)
for the full instructions.

## What becomes available in NOMAD

### ELN-based metadata entry

Not all metadata required by `NXxps` is present in raw vendor files.
`pynxtools-xps` ships an exemplary Electronic Lab Notebook (ELN) schema
(`xps.scheme.archive.yaml`) that lets users enter the missing instrument, sample, and user metadata directly in the NOMAD GUI, without editing YAML by hand.

The schema covers:

- **User** — name, affiliation, email, ORCID
- **Instrument** — X-ray source, analyzer model, collection column, energy dispersion,
  detector, manipulator, flood gun, pressure gauge
- **Sample** — name, atom types, physical form, preparation history, CAS number

Once filled in, the ELN entry is automatically converted to an `eln_data.yaml` file (via `ElnYamlConverter`) that the XPS reader uses as supplementary input during conversion.

### Data conversion via the NOMAD GUI

The XPS reader is exposed through the `NexusDataConverter` inside NOMAD.
Users select their raw measurement file and the ELN entry; NOMAD invokes the `xps` reader and produces a standardized `NXxps` NeXus file — all from within the browser interface.

The same conversion can also be triggered from the command line (see
[Tutorial > Standalone data conversion](../tutorial/standalone.md)):

```console
user@box:~$ dataconverter measurement.sle eln_data.yaml --reader xps --nxdl NXxps --output measurement.nxs
```

### Parsing, search, and visualization

Once the NeXus file exists in a NOMAD upload:

- **DATA tab** — the `NexusParser` extracts structured data from the HDF5 file into NOMAD Metainfo, making individual fields browsable and queryable.
- **FILES tab** — NeXus files are rendered interactively with the
  [H5Web](https://h5web.panosc.eu/){:target="_blank" rel="noopener"} HDF5 viewer,
  allowing direct inspection of datasets, attributes, and groups.
- **Search** — all parsed XPS quantities (technique, energy scale, analyzer model, …) are
  indexed by Elasticsearch and available in NOMAD's search interface.

### Interactive analysis with NORTH

If NOMAD is deployed with the [NOMAD Remote Tools Hub (NORTH)](https://nomad-lab.eu/prod/v1/docs/howto/north/index.html){:target="_blank" rel="noopener"}, users can open the NeXus files in pre-configured JupyterLab containers and run peak-fitting analysis without any local Python installation.

## The XPS example upload

`pynxtools-xps` ships a complete example upload registered as a NOMAD plugin.
It appears in the NOMAD example gallery as **"X-ray Photoelectron Spectroscopy (XPS)"**
under the category *NeXus Experiment Examples* and demonstrates the full FAIR data
pipeline: raw vendor file + ELN metadata → NXxps NeXus file → parsed NOMAD entry →
interactive analysis.

The example uses a SPECS SLE file (Au in 25 mbar O₂, NAP-XPS) and includes an ELN
schema, two conversion configs, and two Jupyter analysis notebooks covering the conversion
walkthrough and lmfit peak fitting.

See [Tutorial > NOMAD usage](../tutorial/nomad_usage.md) for a step-by-step walkthrough
of how to access and run it.

## Further reading

- [Tutorial > Installation guide](../tutorial/installation.md) — setting up `pynxtools-xps` with NOMAD
- [Tutorial > NOMAD usage](../tutorial/nomad_usage.md) — step-by-step walkthrough of the NOMAD workflow
- [pynxtools NOMAD integration](https://fairmat-nfdi.github.io/pynxtools/learn/pynxtools/dataconverter-and-readers.html){:target="_blank" rel="noopener"} — the generic NeXus plugin layer
- [NOMAD example upload plugin docs](https://nomad-lab.eu/prod/v1/docs/howto/plugins/types/example_uploads.html){:target="_blank" rel="noopener"} — how example upload entry points work
- [Explanation > Parser architecture](parser_architecture.md) — how the XPS reader converts raw files to NeXus
- [Explanation > NXmpes and NXxps](appdefs.md) — the application definitions targeted by the conversion
