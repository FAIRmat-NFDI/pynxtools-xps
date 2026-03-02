# Standalone data conversion

## Who is this tutorial for?

This tutorial is for researchers who have a raw XPS measurement file and want to
produce a standardized NeXus output without using NOMAD.

## What should you know before this tutorial?

- `pynxtools-xps` must be installed — see [Tutorial > Installation guide](installation.md).
- Basic command-line familiarity is sufficient. No Python coding is required.

## What you will know at the end of this tutorial?

You will understand:

- What inputs the dataconverter needs and why
- How to run the conversion for each supported file format
- How to inspect the resulting NeXus file

## Steps

### 1. Understand the three inputs

Every conversion takes three inputs and produces one output:

```text
Raw measurement file  ─┐
ELN yaml              ─┼─►  dataconverter (reader: xps, nxdl: NXxps)  ─►  output.nxs
Optional params/config ─┘
```

**Raw measurement file** — the vendor file from your instrument (`.sle`, `.vms`, `.xml`, …).

**ELN yaml** — a YAML file you fill in with instrument, sample, and user metadata that
the raw file does not contain. `NXxps` requires fields such as the analyzer work function,
X-ray source details, and sample description; vendor files rarely provide all of them.
You do not need to write this file from scratch — ready-to-use templates are in the
[`examples/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/){:target="_blank" rel="noopener"}
directory of the repository, one per format.

**Optional params file or config JSON** — some formats accept additional conversion
options (e.g., filtering alignment spectra in SLE files). These are passed with
`--params-file` or `-c`.

### 2. Get an ELN template for your format

Copy the ELN template closest to your format from the table below:

| Format | ELN template |
| ------ | ------------ |
| VAMAS `.vms` / `.npl` | [`examples/vms/eln_data_vms.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/vms/eln_data_vms.yaml){:target="_blank" rel="noopener"} |
| SPECS SLE `.sle` | [`examples/specs/sle/eln_data_sle.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/specs/sle/eln_data_sle.yaml){:target="_blank" rel="noopener"} |
| SPECS XML `.xml` | [`examples/specs/xml/eln_data_xml.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/specs/xml/eln_data_xml.yaml){:target="_blank" rel="noopener"} |
| SPECS XY `.xy` | [`examples/specs/xy/eln_data_xy.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/specs/xy/eln_data_xy.yaml){:target="_blank" rel="noopener"} |
| Scienta TXT `.txt` | [`examples/scienta/txt/eln_data_scienta_txt.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/scienta/txt/eln_data_scienta_txt.yaml){:target="_blank" rel="noopener"} |
| Scienta IBW `.ibw` | [`examples/scienta/ibw/eln_data_scienta_ibw.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/scienta/ibw/eln_data_scienta_ibw.yaml){:target="_blank" rel="noopener"} |
| PHI `.spe` / `.pro` | [`examples/phi/eln_data_phi.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/phi/eln_data_phi.yaml){:target="_blank" rel="noopener"} |

Fill in the relevant fields (instrument details, sample description, user information) before running the conversion.

### 3. Run the conversion

The general command is:

```console
user@box:~$ dataconverter <measurement-file> <eln_data.yaml> \
    --reader xps \
    --nxdl NXxps \
    --output output.nxs
```

Replace `<measurement-file>` and `<eln_data.yaml>` with your actual file paths.
Use `--nxdl NXmpes` instead of `NXxps` if you need the general photoemission schema
(see [Explanation > NXmpes and NXxps](../explanation/appdefs.md)).

The format-specific examples below use the files from the `examples/` directory.

#### VAMAS (`.vms`, `.npl`)

```console
user@box:~$ dataconverter regular.vms eln_data_vms.yaml \
    --reader xps \
    --nxdl NXxps \
    --output regular.nxs
```

Both regular (equally spaced energy axis) and irregular VAMAS files are supported.

#### SPECS SLE (`.sle`)

The SLE example uses a `params.yaml` file that bundles all arguments:

```console
user@box:~$ dataconverter --params-file params.yaml
```

The `params.yaml` in [`examples/specs/sle/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/examples/specs/sle){:target="_blank" rel="noopener"}
supports a `remove_align` keyword. Setting it to `true` removes alignment spectra from the
output, which can significantly speed up conversion for large files:

```yaml
# params.yaml excerpt
remove_align: true
```

#### SPECS XML (`.xml`)

```console
user@box:~$ dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_xml.yaml \
    --reader xps \
    --nxdl NXxps \
    --output output.nxs
```

#### SPECS XY (`.xy`)

```console
user@box:~$ dataconverter MgFe2O4.xy eln_data_xy.yaml \
    --reader xps \
    --nxdl NXxps \
    --output MgFe2O4.nxs
```

#### Scienta TXT and IBW (`.txt`, `.ibw`)

Scienta instruments export one file per scan region. Pass all scan files as positional
arguments before the ELN yaml:

```console
user@box:~$ dataconverter \
    Cu-HHTP_0001.txt Cu-HHTP_0002.txt Cu-HHTP_0003.txt \
    eln_data_scienta_txt.yaml \
    --reader xps \
    --nxdl NXxps \
    --output Cu-HHTP.nxs
```

The same multi-file pattern applies to IBW files:

```console
user@box:~$ dataconverter \
    Cu-HHTP_0001Cu-HHTP__001.ibw Cu-HHTP_0002Cu-HHTP__002.ibw \
    eln_data_scienta_ibw.yaml \
    --reader xps \
    --nxdl NXxps \
    --output Cu-HHTP.nxs
```

#### PHI (`.spe`, `.pro`)

```console
user@box:~$ dataconverter SnO2_10nm.spe eln_data_phi.yaml \
    --reader xps \
    --nxdl NXxps \
    --output SnO2_10nm.nxs
```

### 4. Inspect the output

The output is a standard HDF5 file with the `.nxs` extension. You can inspect it in several ways:

**Python (`h5py`):**

```python
import h5py

with h5py.File("output.nxs", "r") as f:
    f.visit(print)  # print the full group/dataset tree
```

**Browser (H5Web):**

Open [h5web.panosc.eu](https://h5web.panosc.eu/){:target="_blank" rel="noopener"} and
upload your `.nxs` file for interactive exploration of datasets, attributes, and groups.

**NOMAD:**

Upload the `.nxs` file to a NOMAD instance to get full parsing into NOMAD Metainfo,
Elasticsearch indexing, and the H5Web viewer embedded in the FILES tab.
See [Tutorial > Usage in NOMAD](nomad_usage.md) for instructions.

## Further reading

- [Tutorial > Installation guide](installation.md) — setting up `pynxtools-xps`
- [Tutorial > Usage in NOMAD](nomad_usage.md) — using the converter inside NOMAD
- [Explanation > Parser architecture](../explanation/parser_architecture.md) — how raw files are mapped to NeXus
- [Explanation > NXmpes and NXxps](../explanation/appdefs.md) — the output schemas
- [Reference > SPECS formats](../reference/specs.md) — SLE, XML, XY format details
- [Reference > Scienta Omicron formats](../reference/scienta.md) — TXT, IBW, HDF5 format details
- [Reference > PHI formats](../reference/phi.md) — SPE, PRO format details
- [Reference > VAMAS format](../reference/vms.md) — VMS, NPL format details
