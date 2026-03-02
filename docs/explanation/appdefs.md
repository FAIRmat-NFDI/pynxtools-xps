# NXmpes and NXxps

`pynxtools-xps` converts heterogeneous XPS data to two NeXus application definitions:

- [**`NXmpes`**](https://manual.nexusformat.org/classes/applications/NXmpes.html){:target="_blank" rel="noopener"}: a general-purpose photoemission standard
- [**`NXxps`**](https://manual.nexusformat.org/classes/applications/NXxps.html){:target="_blank" rel="noopener"}: its XPS-specific specialization.

Both are accepted as standards by the NeXus International Advisory Committee (NIAC) and are available in the
[official NeXus definitions repository](https://github.com/nexusformat/definitions){:target="_blank" rel="noopener"}.

## Motivation and inspiration

Photoemission spectroscopy (PES) suffers from a fragmented data landscape: dozens of proprietary and open file formats, each encoding the same physical quantities with different names, units, and conventions.
The existing [VAMAS / ISO 14976](https://www.iso.org/standard/24269.html){:target="_blank" rel="noopener"} standard covers basic data transfer but does not capture the full richness of modern photoemission experiments — instrument geometry, calibration provenance, peak-fitting results, and
multi-dimensional momentum-resolved datasets.

`NXmpes` and `NXxps` were designed to close this gap. They draw their vocabulary from [ISO 18115-1:2023](https://www.iso.org/standard/82860.html){:target="_blank" rel="noopener"} (Surface Chemical Analysis — Vocabulary), IUPAC recommendations for PES nomenclature, and the established VAMAS data model, elevating these community standards into a machine-readable,
self-describing HDF5 format.

The definitions were developed within the [FAIRmat project](https://www.fairmat-nfdi.eu/){:target="_blank" rel="noopener"} of the German National Research Data Infrastructure (NFDI) in close collaboration with PES
researchers and instrument vendors. The process was deliberately bottom-up: existing data from real instruments drove the schema, not an abstract top-down design. In 2025, the application definitions as well as a group of related base classes was proposed to NIAC for standardization and, following a successful vote, accepted as officially supported standards.

## Instrument-agnostic by design

Both application definitions are intentionally vendor-neutral. There are no vendor- or format-specific concepts anywhere in the schema. All vendor-specific concepts are handled in the parser layer of `pynxtools-xps`
(see [Explanation > Parser architecture](parser_architecture.md)); the application definitions express only the physical experiment.

This separation means that a dataset from a Scienta PEAK system and one from a SPECS PHOIBOS instrument can be stored in identically structured NeXus files, differing only in their metadata values — not in their layout.

The technique scope is equally broad. `NXmpes` covers the full photoemission family. This includes, among others:

- X-ray Photoelectron Spectroscopy (XPS/ESCA)
- Angle-Resolved XPS (ARXPS)
- Ultraviolet Photoelectron Spectroscopy (UPS)
- Hard X-ray Photoelectron Spectroscopy (HAXPES)
- Near-Ambient Pressure XPS (NAPXPS)
- Angle-Resolved Photoemission Spectroscopy (ARPES)
- Photoemission Electron Microscopy (PEEM)
- Two-Photon Photoemission (2PPE)
- Time-resolved variants (tr-ARPES, tr-XPS)

`NXxps` focuses on the XPS subset of this family.

## `NXmpes` — the general photoemission base

`NXmpes` provides the minimal but complete metadata framework for any photoemission experiment. Its philosophy is to capture what every PES measurement shares: a photon source, an electron analyzer, a sample, and a measured spectrum — while leaving room for multi-dimensional data axes
and post-acquisition calibration records.

### Structure

```text
NXmpes
└── ENTRY (NXentry)
    ├── definition = "NXmpes"
    ├── title, start_time, end_time
    ├── method (XPS · UPS · ARPES · …)
    ├── transitions (ISO 18115-1 notation: "C 1s", "Fe 2p3/2", …)
    ├── INSTRUMENT (NXinstrument)
    │   ├── beam_probe (NXbeam) — incident photon beam
    │   └── ELECTRONANALYZER (NXelectronanalyzer)
    │       ├── COLLECTIONCOLUMN (NXcollectioncolumn)
    │       ├── ENERGYDISPERSION (NXenergydispersion)
    │       └── DETECTOR (NXdetector)
    ├── SAMPLE (NXsample)
    ├── DATA (NXdata) — photoemission counts vs. energy/momentum/…
    ├── energy_referencing (NXcalibration)
    └── transmission_correction (NXcalibration)
```

### Key design choices

**`transitions` field.**
Peak labels follow strict ISO 18115-1 notation (`"C 1s"`, `"Fe 2p3/2"`, `"O KLL"`), making automated peak identification unambiguous across datasets.

**Multi-dimensional `DATA`.**
The `DATA` group axes are not fixed to energy alone. Axes can be kinetic or binding energy, momentum components ($k_x$, $k_y$, $k_z$), angular,
spatial ($x$, $y$), or pump-probe time delay. A basic XPS survey spectrum and a full 4D ARPES dataset both fit within the same group structure.

**Calibration provenance.**
The optional `energy_referencing` and `transmission_correction` groups (as well as further optional `NXprocess` groups) record the calibration and processing parameters alongside the data rather than embedding corrections silently into the intensity values. This ensures reproducibility and allows downstream reanalysis.

## `NXxps` — the XPS specialization

`NXxps` extends `NXmpes` with fields that are specific to X-ray photoelectron spectroscopy: explicit laboratory geometry, analyzer characterization, and a structured peak-fitting record. Every valid `NXxps` file is also a valid `NXmpes` file; `NXxps` only adds constraints and groups.

### Additions over `NXmpes`

| Group / field | Purpose |
| ------------- | ------- |
| `xps_coordinate_system` | Explicit lab-frame coordinate system anchored at the sample stage with $z$ pointing along the sample normal, following existing conventions in the VAMAS/ISO standard |
| `source_probe.power` | X-ray source power in watts |
| `beam_probe` angles | `beam_polar_angle_of_incidence` and `beam_azimuth_angle` relative to the sample frame |
| `ELECTRONANALYZER.work_function` | Analyzer work function (required) |
| Analyzer take-off angles | `analyzer_take_off_polar_angle` and `analyzer_take_off_azimuth_angle` |
| Sample orientation | `sample_rotation_angle`, `sample_normal_polar_angle_of_tilt`, `sample_normal_tilt_azimuth_angle` |
| `FIT` group | A complete, self-contained record of a peak-fitting analysis for one spectral region |

### The FIT group

The `FIT` group is unique to `NXxps` and encodes quantitative XPS analysis alongside the raw data.

```text
NXmpes
└── FIT (NXfit)
    ├── label — region identifier
    ├── figure_of_merit
    ├── data (NXdata) — measured spectrum + fit + residual
    ├── peakPEAK (NXpeak) — one per fitted component
    │   ├── label, total_area, relative_atomic_concentration
    │   └── function (NXfit_function) — function_type · formula · parameters
    ├── backgroundBACKGROUND (NXpeak) — background model
    │   └── function (NXfit_function) — Shirley · Tougaard · Linear · Step
    ├── global_fit_function (NXfit_function)
    └── error_function (NXfit_function)
```

Peak function types supported by `NXfit_function`: Gaussian, Lorentzian, Voigt, Gaussian–Lorentzian Sum/Product, Asymmetric Lorentzian, Doniach–Šunjić, Asymmetric Finite Lorentzian.

Storing the fit alongside the measured data keeps quantification results FAIR: parameters, residuals, and quality metrics are preserved in the same file and can be re-evaluated or compared across experiments without relying on external analysis software state.

### The XPS coordinate system

`NXxps` defines a coordinate system based on the sample stage, which is the same coordinate system as in the [ISO standard](https://www.iso.org/standard/24269.html) for surface chemical analysis. You can learn more about this coordinate system in [Explanation > The XPS coordinate system](coordinate_system.md).

## Relationship between `NXmpes` and `NXxps`

`NXxps` is a strict specialization of `NXmpes`: it inherits the full `NXmpes` structure and adds XPS-specific constraints and groups.

```text
NXmpes (general photoemission)  ── extended by ──►  NXxps (XPS specialization)
```

In practice:

- **Use `NXmpes`** when the conversion targets a multi-technique context (e.g. ARPES + XPS),
  when only the core photoemission metadata is available, or when compatibility with MPES-family
  tooling is the priority.
- **Use `NXxps`** as the default for XPS workflows: it adds the coordinate system, analyzer work
  function, x-ray source power, and the peak-fitting record that characterize quantitative XPS
  analysis.

`pynxtools-xps` supports both targets.
The choice is made at conversion time via the `--nxdl` flag:

```shell
dataconverter example.vms eln_data.yaml --reader xps --nxdl NXxps  --output example.nxs
dataconverter example.vms eln_data.yaml --reader xps --nxdl NXmpes --output example.nxs
```

`NXxps` is recommended for most XPS workflows; `NXmpes` is the appropriate choice when the richer XPS-specific structure is not needed or not supported by the source data.

!!! note
    There exist further specializations of `NXmpes` for specific photoemission spectroscopy sub-techniques. An example is [`NXmpes_arpes`](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXmpes_arpes.html){:target="_blank" rel="noopener"} for angle-resolved PES (ARPES). Learn more about these in the [NeXus > MPES application definitions overview](https://manual.nexusformat.org/classes/applications/mpes-structure.html#appdef-mpes-structure){:target="_blank" rel="noopener"}.

## Further reading

- [NeXus > NXmpes](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXmpes.html){:target="_blank" rel="noopener"} — full field reference
- [NeXus > NXxps](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXxps.html){:target="_blank" rel="noopener"} — full field reference
- [NeXus > MPES application definitions overview](https://manual.nexusformat.org/classes/applications/mpes-structure.html#appdef-mpes-structure){:target="_blank" rel="noopener"} — the broader NXmpes family (NXmpes_arpes, NXxps, legacy NXarpes)
- [NeXus > Base classes used in `NXmpes`](https://manual.nexusformat.org/classes/base_classes/mpes-structure.html){:target="_blank" rel="noopener"} — NXelectronanalyzer, NXcollectioncolumn, NXenergydispersion, and others
- [Explanation > Motivation](motivation.md) — why `pynxtools-xps` exists and the XPS software landscape it addresses
- [Explanation > Parser architecture](parser_architecture.md) — how raw vendor data is converted to NXmpes / NXxps
