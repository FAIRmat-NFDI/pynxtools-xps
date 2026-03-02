# Data from Kratos instruments

[Kratos Analytical](https://www.kratos.com/) XPS instruments (AXIS Ultra, AXIS Nova, AXIS Supra, and related models) export data in VAMAS ISO 14976 format (`.vms` files). The reader handles Kratos files through the standard
[VAMAS parser](vms.md), with automatic detection and parsing of the Kratos-specific metadata stored in the VAMAS comment blocks.

## Supported versions

The Kratos reader supports any VAMAS file that contains a recognized Kratos metadata block. Version detection is not enforced for this format.

## Supported metadata

The Kratos comment parser (see
[`parsers/kratos/`](https://github.com/FAIRmat-NFDI/pynxtools-xps/tree/main/src/pynxtools_xps/parsers/kratos))
extracts the following instrument-specific fields from the VAMAS comment blocks when present:

- Lens mode (`lens`)
- X-ray source (anode material, power, deflection)
- Charge neutralizer state
- Stage position (`tilt`)
- Energy range (`start`, `end`, `centre`, `width`)
- Step size and dwell time
- Sweep count and acquisition date

## Example conversion

Kratos `.vms` files are converted using the same command as any other VAMAS file:

```console
user@box:~$ dataconverter my_kratos_data.vms eln_data.yaml --reader xps --nxdl NXxps --output my_kratos_data.nxs
```

An ELN YAML file is required to supply any required NeXus fields not present in the
VAMAS file. Use
[`examples/vms/eln_data_vms.yaml`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/examples/vms/eln_data_vms.yaml)
as a starting template.

!!! note
    Kratos-specific metadata is read automatically whenever the reader detects Kratos comment blocks in a VAMAS file. No additional flags or config options are required.

## Further reading

- [Reference > VAMAS ISO standard](vms.md) — full description of the VAMAS format
- [Explanation > Parser architecture](../explanation/parser_architecture.md)