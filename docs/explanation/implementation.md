# Purpose and aim of pynxtools-xps
pynxtools-xps aims for the implementation of [FAIR principles of data stewardship](https://doi.org/10.1162/dint_r_00024) in photoelectron spectroscopy (PES). In many experimental fields, there has been a push towards such standardization and interoperability in recent yeards; however, there has been a distinct lack of such efforts in PES.

While there exists a widely adopted [ISO standard](https://www.iso.org/standard/24269.html) for data transfer in surface chemical analysis, it does not fully cover all of the information that is obtained in modern photoemission experiments. Within the FAIRmat project of the German National Research Data Infrastructure Germany (NFDI), we have spent considerable effort towards building developing an extensive and elaborated standard ([NXmpes](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html) with its specialization [NXxps](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html)) for harmonizing PES data using [NeXus](https://www.nexusformat.org/), a community-driven data-modeling framework for experiments.

The goal of pynxtools-xps is to provide a mapping from the diverse proprietary and open-source software solutions used in the XPS community towards this NeXus standard. The software implements a suggestion how diverse (meta)data from the research field of photoelectron spectroscopy can be parsed and normalized to enable users to compare data. The software can parse data provided by many different data providers and frequently used serialization and formatting. This data is then mapped onto the NeXus standard.

As part of [pynxtools and its plugin infrastructure](https://github.com/FAIRmat-NFDI/pynxtools), pynxtools-xps is fully integrated into the NOMAD research data management systems (RDMS), with the aim of facilitating harmonization of XPS data and enabling development of data-centric software tools and services.

# Software landscape in photoelectron spectroscopy - a mixture of proprietary and open-source solutions
As in many other experimental fields, the software landscape in photoelectron spectroscopy is extremely diverse, ranging from fully integrated software solution from technology partners, that integrate the measurement, post-processing, and data analysis, to custom-written software for specific uses cases.
While propietary software is often easy to use for end users, such software often writes to proprietary serialization formats (file or database entries). Not only are these formats not openly readable, but they are often not well-documented and the content and meaning of the semantic concepts is very often not documented publicly. While open-source software typically writes to more openly (and sometimes better documented) formats, they tend to loose much of the metadata that commercial vendors can provide with their data. pynxtools-xps aims at both endpoints of this spectrum and everything in-between: it provides an easy-to-use framework for writing a standardization parser for small, encapsulated solutions, while also providing the possibility of mapping the full richness of data and metadata acquired in a high-end XPS laboratory onto NeXus.

# Implementation design

pynxtools-xps is a community-based tool that provides a bottom-up approach for mapping XPS data onto the NXmpes and NXxps standards. Specificallly, the software contains example parsers for data that was measured by PES researchers in a wide array of experimental setups. The goal is not neccessarily to implement a fully comprehensive mapping of all possible existing file formats, but rather help the individual researcher or technology partner to start reading their data into the NeXus standard.

Therefore, the following design patterns guide our implementation:

- We do not consider that our work is complete (from the perspective of the idea in mind that a user can expect to drag-and-drop arbitrary content).
- We consider ontology matching a team effort that can only be achieved with technology partners and scientists working together.
- Our work is open to suggestions by the PES community, always realizing that just being able to read from a specific file alone is not solving the challenge that pynxtools-xps addresses.
- We provide specific tangible examples of (meta)data semantic mapping for specific file formats that are frequently used in XPS. These include the main formats of the leading vendors of PES spectrometers.
- The tool itself is build such that is easily extendable.
- The goal is to continously grow the number of parsers available for different communities. We therefore encourage researchers and technology partners to get in contact in order to get started with standardization in NeXus and NOMAD.
