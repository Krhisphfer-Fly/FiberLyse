# Importing data

FiberLyse separates **how a file is read** from **how the photometry signal is analyzed**. The import setup translates different CSV organizations into the same internal study/reference representation.

## Neurophotometrics / Bonsai

The built-in profile expects the familiar structure containing:

- `SystemTimestamp` for time;
- `LedState` for row/state identity;
- fluorescence columns named `G`, `G0`, `G1`, etc.

FiberLyse uses the LED-state information to separate reference and study measurements and handles the configured recording-start/exclusion logic before downstream analysis.

## Separate signal and reference columns

Choose this layout when every row contains:

- one time value;
- a study-signal column;
- a separate comparison/reference column.

The wizard can define multiple study/reference pairs from one CSV.

Typical example:

```text
time (sec);405_channel;465_channel;dF/F
0,0009;66,5;407,6;4,2
...
```

A sensible raw-data mapping for this kind of file could be:

- Time → `time (sec)`
- Study → `465_channel`
- Reference → `405_channel`

A pre-existing processed `dF/F` column does not need to be selected; FiberLyse can calculate its own ΔF/F from the mapped raw study/reference channels.

## Interleaved/alternating rows

Choose this layout when study and reference measurements take turns in the rows and another column identifies what each row represents.

The setup asks for:

- time column;
- row/state/type column;
- study row value;
- reference row value;
- fluorescence-value column(s).

For LedState-like files, values such as `1`, `2`, and `7` can be used when they are the actual coding in the recording. Always verify the mapping from the acquisition system rather than assuming wavelength meaning from the number alone.

## CSV separator and decimal mark

Before column mapping, FiberLyse can determine or let the user choose how the CSV itself is written.

Supported separator choices include common formats such as:

- comma;
- semicolon;
- tab;
- pipe.

Decimal numbers can use a period or comma. The same selected CSV settings are used for the actual analysis, not only the preview.

## Time units

Custom profiles can define the units used by the time column. FiberLyse converts the imported times to seconds internally.

## Recording start and initial exclusion

The custom setup can define where the recording begins and how much initial time should be excluded before downstream analysis. This is useful when the first part of a recording is known to contain setup/non-biological acquisition rather than experimental signal.

## Saved profiles

Custom CSV mappings can be remembered for later sessions. Saved mappings are based on the configured column names and layout, so reopen **CSV setup…** if another exporter changes its column naming or CSV format.

## Validation recommendation

After creating a new import profile:

1. inspect the raw study/reference traces;
2. confirm the detected sampling rate is plausible;
3. confirm the recording begins where expected;
4. inspect the fitted-reference graph;
5. only then interpret normalized/frequency outputs.
