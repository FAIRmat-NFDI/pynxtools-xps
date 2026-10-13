# SPECS parameter-history logs (.slh and .csv)

SpecsLab Prodigy can log instrument parameters — pressure, temperature, voltages, currents — over the course of a measurement session. `pynxtools-xps` reads these logs from two export formats and merges the readings into the matching entry of an `NXxps` conversion.

Parsers:
[`SPECSMetadataSLHParser`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/specs/slh/parser.py) and [`SPECSMetadataCSVParser`](https://github.com/FAIRmat-NFDI/pynxtools-xps/blob/main/src/pynxtools_xps/parsers/specs/slh/parser.py). Both are metadata-only parsers: they add data to entries produced by a primary parser, not entries of their own. The only compatible primary parser is [`SPECSSLEParser`](../reference/specs.md).

## File formats

### .slh

An `.slh` file is a SQLite database. One file logs one group of parameters (for example`xray`, `nap_parameters`, `phoibos_voltages`) for a whole session — not any single spectrum's own acquisition window.

Three tables matter:

| Table | Holds |
| --- | --- |
| `ParameterHistory` | Which parameters the file logs: name, device, base timestamp |
| `ParameterInfo` | Unit and scaling factor for each parameter |
| `NumericalHistoryData` | The observations: `(Offset_s, Observation)` pairs, relative to the base timestamp |

Two points to know if you read the raw database directly:

- `Observation` is stored unscaled. Multiply by `ParameterInfo.Scaling` to get the real value.
- `NumericalHistoryData` rows are not guaranteed to be stored in chronological `Offset_s`
  order.

### .csv

A flat export of the same underlying log data: one row per timestamp, one column per parameter. Despite the `.csv` extension, the delimiter is `;`, not `,`:

```csv
"Universal Time";"Local Date";"Local Time of Day";"Relative Time [s]";"CoilCurrent Current [µA] (Phoibos)";"Pressure [mbar] (Ceravac_AC)"
2025-Apr-03 06:57:15.625973;2025-04-03;08:57:15.625973;0.000135;0;1.2e-06
2025-Apr-03 06:57:16.625973;2025-04-03;08:57:16.625973;1.000135;;1.1e-06
```

Each parameter column header encodes name, unit, and device: `"<param> [<unit>] (<device>)"`. A blank cell means "unchanged since the last row", not zero — each parameter's series is built only from its own non-blank rows.

The timestamp comes from `Universal Time` if present, otherwise from `Local Date` + `Local Time of Day` combined. If neither column is present, the file is skipped and a warning is logged.

## Format detection

`SPECSMetadataSLHParser.matches_file` checks the SQLite magic bytes, then confirms all three tables above are present. `.sle` files are SQLite too, but never have a `ParameterHistory` table, so the table check is what tells the two formats apart.

`SPECSMetadataCSVParser.matches_file` reads only the header line: it confirms the `;` delimiter and that at least one of the three timestamp columns is present. This currently rejects a comma-delimited CSV.

## Resampling

Both parsers resample each parameter's observations to a fixed interval, forward-filling gaps, so a value that hasn't changed since the last sample still has a reading at every step. Default interval: `1s`. Set per file type through reader kwargs:

- `slh_resample_interval` for `.slh` files
- `csv_resample_interval` for `.csv` files

Pass `None` to keep the raw, irregularly spaced timestamps instead.

## Merging into entries

`update_main_file_data` attaches each parameter's readings to the matching `NXentry`, sliced to that entry's own acquisition window:

- **Window**: `[time_stamp, time_stamp + dwell_time × n_values × total_scans]`, taken from the entry's own metadata (already parsed from the `.sle` file).
- Only observations that fall inside the window are attached.
- If nothing falls inside the window, that parameter is left out of the entry — there is no fallback to a "last known value" scalar.

Attached values are stored on the entry as `<param>_<device>`, with `<param>_<device>/@units` and `<param>_<device>/time` alongside. A device that's only
sometimes logged is handled naturally: its data is simply absent from entries outside its `.slh` file's own session.

## Usage

Give the `.sle` file and any number of `.slh`/`.csv` files to the reader together:

```console
pynx convert my_experiment.sle xray.slh nap_parameters.slh phoibos_voltages.slh \
    eln_data.yaml --reader xps --nxdl NXxps --output EX1559_S1710.nxs
```

To set the resample interval, use a `params.yaml` file instead (same convention as the SLE parser's `remove_align`; see [Reference > SPECS](../reference/specs.md)):

```yaml
dataconverter:
  reader: xps
  nxdl: NXxps
  input-file:
    - EX1559_S1710.sle
    - xray.slh
    - eln_data.yaml
  slh_resample_interval: "5s"
  output: EX1559_S1710.nxs
```

```console
pynx convert --params-file params.yaml
```

## Version support

{{ parser_version_table("specs.slh.parser", "SPECSMetadataSLHParser") }}

`SPECSMetadataCSVParser` declares no version constraint: nothing in a `.csv` export carries a version marker to check.

## Further reading

- [Reference > SPECS](../reference/specs.md)
- [Explanation > Parser architecture](parser_architecture.md)
