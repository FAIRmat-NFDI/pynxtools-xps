---
hide: toc
---

# Documentation for pynxtools-xps

pynxtools-xps is a free, and open-source data software for harmonizing X-ray photolectron spectroscopy data and metadata for research data management using [NeXus](https://www.nexusformat.org/), implemented with the goal to make scientific research data FAIR (findable, accessible, interoperable and reusable).

pynxtools-xps, which is a plugin for [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools), provides a tool for reading data from various propietary and open data formats from technology partners and the wider XPS community and standardizing it such that it is compliant with the NeXus application definition [`NXmpes`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html) and <!--[`NXxps`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html), which is an extension of `NXmpes` -->. pynxtools-xps is developed both as a standalone reader and as a tool within [NOMAD](https://nomad-lab.eu/), which is the open-source data management platform for materials science we are developing with [FAIRmat](https://www.fairmat-nfdi.eu/fairmat/).

pynxtools-xps solves the challenge of using heterogeneous and unfindable data formats which is common in X-ray Photoelectron Spectroscopy. In addition, it provides an interface for writing readers for different file formats to be mapped to NeXus.

pynxtools-xps is useful for scientists from the XPS community that deal with heterogeneous data, for technology partners and data providers looking for ways to make their data FAIRer, and for research groups that want to organize thier data using NeXus and NOMAD.

<div markdown="block" class="home-grid">
<div markdown="block"> 

### Tutorial

A series of tutorials giving you an overview on how to store or convert your XPS data to NeXus compliant files.

- [Installation guide](tutorial/installation.md)
- [Standalone usage](tutorial/standalone.md)
- [How to use a NeXus/HDF5 file](tutorial/nexusio.md)
- [Usage in NOMAD](tutorial/nomad.md)

</div>
<div markdown="block">

### How-to guides

How-to guides provide step-by-step instructions for a wide range of tasks, with the overarching topics:

- [How to create your own reader for your XPS data](how-tos/build_a_reader.md)

</div>

<div markdown="block">

### Learn

The explanation section provides background knowledge on the implementation design, how the data is structured, how data processing can be incorporated, how the integration works in NOMAD, and more.

- [Design principles and implementation](explanation/implementation.md)
- [NXmpes](explanation/appdefs.md) <!-- - [NXmpes and NXxps](explanation/appdefs.md) -->
- [The XPS coordinate system](explanation/coordinate_system.md)
- [How to map pieces of information to NeXus](explanation/contextualization.md)
<!-- - [Data processing](explanation/data_processing.md) -->
<!-- - - [NOMAD integration](explanation/nomad_integration.md) -->

</div>
<div markdown="block">

### Reference

Here you can learn which specific measurement setups and file formats from technology partners pynxtools-xps currently supports.

The reader decides which data parser to use based on the file extension of the files provided. For the main XPS files, the following file extensions are supported:

- .ibw: [Igor Binary Wave Format](https://www.wavemetrics.com/) files, exported from [Scienta Omicron](https://scientaomicron.com/en)
- .npl: VAMAS files, ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)), both in regular and irregular format
- .spe, .pro: [Phi MultiPak](https://www.phi.com/surface-analysis-equipment/genesis.html#software:multi-pak-data-reduction-software/) files, propietary format of PHI Electronics
- .sle: [SpecsLabProdigy](https://www.specs-group.com/nc/specs/products/detail/prodigy/) files, propietary format of SPECS GmbH (1 and v4)
- .xml: SpecsLab 2files, XML format from SPECS GmbH (v1.6)
- .vms: VAMAS files, ISO standard data transfer format ([ISO 14976](https://www.iso.org/standard/24269.html)), both in regular and irregular format
- .xy: SpecsLabProdigy export format in XY format (including all export settings)
- .txt:
    - exported by [Scienta Omicron](https://scientaomicron.com/en) instruments
    - exported by [CasaXPS](https://www.casaxps.com/) analysis software

You can find more information regarding the readers for data from different technology partners here:

- [Data exported by SPECS spectrometers](reference/specs.md)
- [Data exported by Scienta Omicron spectrometers](reference/scienta.md)
- [Data exported by Phi spectrometers](reference/phi.md)
- [VAMAS ISO Standard format](reference/vms.md)

</div>
</div>

<h2>Project and community</h2>

- [NOMAD code guidelines](https://nomad-lab.eu/prod/v1/staging/docs/reference/code_guidelines.html) 

Any questions or suggestions? [Get in touch!](https://www.fair-di.eu/fairmat/about-fairmat/team-fairmat)

[The work is funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 460197019 (FAIRmat).](https://gepris.dfg.de/gepris/projekt/460197019?language=en)
