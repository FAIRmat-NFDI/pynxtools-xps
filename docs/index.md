---
hide: toc
---

# Documentation for pynxtools-xps

`pynxtools-xps` is a free, open-source data converter for X-ray photoelectron spectroscopy (XPS) using [NeXus](https://www.nexusformat.org/). It reads heterogeneous vendor formats and produces standardized
[`NXmpes`](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXmpes.html){:target="_blank" rel="noopener"}
and [`NXxps`](https://fairmat-nfdi.github.io/nexus_definitions/classes/applications/NXxps.html){:target="_blank" rel="noopener"}
NeXus files, making XPS data FAIR (findable, accessible, interoperable, and reusable).
It works both as a standalone reader tool and as a plugin for
[NOMAD](https://nomad-lab.eu/){:target="_blank" rel="noopener"}, the open-source data management platform for materials science we are developing with [FAIRmat](https://www.fairmat-nfdi.eu/fairmat/).

pynxtools-xps solves the challenge of using heterogeneous and undocumented data formats which is common in X-ray Photoelectron Spectroscopy. In addition, it provides an interface for writing readers for different file formats to be mapped to NeXus.

pynxtools-xps is useful for scientists from the XPS community that deal with heterogeneous data, for technology partners and data providers looking for ways to make their data FAIRer, and for research groups that want to organize their data using NeXus and NOMAD.

<div markdown="block" class="home-grid">
<div markdown="block">

### Tutorials

A series of tutorials giving you an overview on how to store or convert your XPS data to NeXus compliant files.

- [Installation guide](tutorial/installation.md)
- [Standalone data conversion](tutorial/standalone.md)
- [Usage in NOMAD](tutorial/nomad_usage.md)
- [Development guide](tutorial/contributing.md)

</div>
<div markdown="block">

### How-To Guides

How-to guides provide step-by-step instructions for a wide range of tasks, with the overarching topics:

- [Build a new parser for your XPS data](how-tos/build_a_parser.md)

</div>

<div markdown="block">

### Explanation

The explanation section provides background knowledge on the implementation design, how the data is structured, how data processing can be incorporated, how the integration works in NOMAD, and more.

- [Motivation and software design](explanation/motivation.md)
- [NXmpes and NXxps](explanation/appdefs.md)
- [Parser architecture](explanation/parser_architecture.md)
- [The XPS coordinate system](explanation/coordinate_system.md)
- [Mapping of data processing performed in CasaXPS](explanation/data_processing.md)
- [NOMAD integration](explanation/nomad_integration.md)
- [Note on versioning](explanation/versioning.md)

</div>
<div markdown="block">

### Reference

The following file formats are supported. The reader selects the correct parser automatically based on the file extension.

{{ supported_formats_table() }}

</div>
</div>

<h2> Contact </h2>

For questions or suggestions:

- Open an issue on the [`pynxtools-xps` GitHub](https://github.com/FAIRmat-NFDI/pynxtools-xps/issues){:target="_blank" rel="noopener"}
- Join our [Discord channel](https://discord.gg/Gyzx3ukUw8){:target="_blank" rel="noopener"}
- Get in contact with our [lead developers](contact.md).

<h2>Project and community</h2>

- [NOMAD code guidelines](https://nomad-lab.eu/prod/v1/staging/docs/reference/code_guidelines.html){:target="_blank" rel="noopener"}

[The work is funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 460197019 (FAIRmat).](https://gepris.dfg.de/gepris/projekt/460197019?language=en){:target="_blank" rel="noopener"}
