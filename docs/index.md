---
hide: toc
---

# Documentation for pynxtools-xps: 

pynxtools-xps is a free, and open-source data software for harmonizing X-ray photolectron spectroscopy data and metadata for research data management using [NeXus](https://www.nexusformat.org/), implemented with the goal to make scientific research data FAIR (findable, accessible, interoperable and reusable).

pynxtools-xps, which is a plugin for [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools), provides a tool for reading data from various propietary and open data formats from technology partners and the wider XPS community and standardizing it such that it is compliant with the NeXus application definitions [`NXmpes`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html) and [`NXxps`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html), which is an extension of `NXmpes`. pynxtools-xps is devloped both as a standalone reader and as a tool within [NOMAD](https://nomad-lab.eu/), which is the open-source data management platform for materials science we are developing with [FAIRmat](https://www.fairmat-nfdi.eu/fairmat/).

pynxtools-xps solves the challenge of using heterogeneous and unfindable data formats which is common in X-ray Photoelectron Spectroscopy. In addition, it provides an interface for writing readers for different file formats to be mapped to NeXus.

pynxtools-xps is useful for scientists from the XPS community that deal with heterogeneous data, for technology partners and data providers looking for ways to make their data FAIRer, and for research groups that want to organize thier data using NeXus and NOMAD.

<div markdown="block" class="home-grid">
<div markdown="block">

### Tutorial

A series of tutorials giving you an overview on how to store or convert your XPS data to NeXus compliant files.

- [Installation guide](tutorial/installation.md)
- [Standalone usage](tutorial/standalone.md)
- [How to use a NeXus/HDF5 file](tutorial/nexusio.md)
- [Usage on NOMAD](tutorial/nomad.md)

</div>
<div markdown="block">

### How-to guides

How-to guides provide step-by-step instructions for a wide range of tasks, with the overarching topics:

- [How to create your own reader for your XPS data](how-tos/build_a_reader.md)

</div>

<div markdown="block">

### Learn

The explanation section provides background knowledge on the implementation design, how the data is structured, how data processing can be incorporated, how the integration works in NOMAD, and more.

- [Implementation design](explanation/implementation.md)
- [NXmpes and NXxps](explanation/appdefs.md)
- [How to map pieces of information to NeXus](explanation/contextualization.md)
<!-- - [Data processing](explanation/data_processing.md) -->
<!-- - - [NOMAD integration](explanation/nomad_integration.md) -->

</div>
<div markdown="block">

### Reference

Here you can learn which specific measurement setups and file formats from technology partners pynxtools-xps currently supports.

- [Data exported by SPECS spectrometers](reference/specs.md)
- [Data exported by Scienta Omicron spectrometers](reference/scienta.md)
- [Data exported by Phi spectrometers](reference/phi.md)
- [VAMAS ISO Standard format](reference/vms.md)

</div>
</div>

<h2>Project and community</h2>
- [Code guidelines](reference/code_guidelines.md)

Thinking about using NOMAD for your next project? Get in touch!
