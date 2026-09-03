# Troubleshooting

## The GUI does not open / Tkinter is missing

FiberLyse uses Tkinter. On Windows, install Python from python.org using the standard installer. Some company-managed or minimal Python distributions omit Tk/Tcl.

A quick test is:

```bash
python -m tkinter
```

If a small Tk window does not open, fix the Python/Tk installation before troubleshooting FiberLyse itself.

## `ModuleNotFoundError` for NumPy, pandas, Matplotlib, or SciPy

Run:

```bash
python -m pip install -r requirements.txt
```

Use the same Python interpreter for both installation and launching FiberLyse.

## My CSV appears as one giant column

Open **CSV setup…** and check the CSV-format step.

Common combinations include:

- comma separator + period decimal;
- semicolon separator + comma decimal;
- tab separator + period decimal.

Adjust the separator and decimal choices until the preview shows the expected columns.

## FiberLyse cannot find my columns

A saved custom profile refers to the column names used when it was created. If another acquisition/export program changed the header names, reopen **CSV setup…** and save a new mapping.

## The detected sampling rate looks wrong

Inspect the selected time column and its units. The effective rate is calculated from the study-signal timestamps after import/separation.

For interleaved recordings, remember that the per-study-signal sampling rate may differ from the acquisition frame/update rate because study and reference rows take turns.

For downsampled files, the detected rate should describe the **current file**, not the original acquisition system.

## A frequency band says unavailable

Check the current detected signal rate. A requested band cannot extend beyond the current Nyquist frequency (`Fs / 2`).

A band can also be unavailable if there is not enough continuous finite data to support several cycles at the lower edge.

## The fit says no valid fit

The selected fit window may contain too few finite paired study/reference points. Inspect:

- the fit interval;
- artifact holes;
- whether study/reference mapping is correct;
- whether the reference is aligned over that time range.

Select a scientifically appropriate interval with enough usable data rather than relying on an automatic fallback.

## The beginning of my recording is missing

Check the import profile's recording-start and initial-exclusion settings. The built-in/custom workflow can intentionally remove a configured initial period before downstream analysis.

## A graph contains gaps

Gaps can be intentional. Artifact removal and non-interpolated analysis preserve missing sections rather than drawing or filtering directly through them.

## Exported values do not look identical to a displayed batch axis

Batch ΔF/F views can use display transformations such as percent units and zeroing the first finite displayed point. Keep the export metadata with batch results so the transformation is explicit.

## Reporting a bug

Open a GitHub issue and include:

- FiberLyse version;
- OS;
- import layout/profile;
- relevant analysis settings;
- exact error text;
- screenshot if useful;
- the smallest shareable example CSV that reproduces the issue.
