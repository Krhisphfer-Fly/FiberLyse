# Getting started with FiberLyse

This guide covers a first installation and analysis using the public `FiberLyse.py` script.

## 1. Install Python

Install a current Python 3 release from python.org. On Windows, use the standard installer so Tkinter/Tk is included.

To check the installation, open PowerShell or a terminal and run:

```powershell
py --version
```

or:

```bash
python --version
```

## 2. Download FiberLyse

Download the GitHub repository ZIP or clone the repository. Extract the ZIP before running the application.

## 3. Install dependencies

Open a terminal in the extracted FiberLyse folder:

```powershell
py -m pip install -r requirements.txt
```

The runtime dependencies are NumPy, pandas, Matplotlib, and SciPy.

## 4. Start the application

```powershell
py FiberLyse.py
```

The main FiberLyse window should open.

## 5. Choose how FiberLyse should read the CSV

For Neurophotometrics/Bonsai-style CSV files, use the built-in import mapping.

For other CSV layouts, open **CSV setup…** and choose the layout that matches the recording:

- signal and reference are separate columns; or
- signal and reference alternate/take turns in rows.

Use an example file to map the time column, fluorescence columns, row/state values, recording start rule, and initial exclusion.

## 6. Add files and run analysis

Use **Add CSV(s)…**, confirm the analysis controls, then choose **Run analysis**.

FiberLyse creates separate views for the analysis stages. Inspect the raw data, artifact-removal result, and fitted reference before interpreting normalized traces.

## 7. Check the sampling rate

FiberLyse calculates the effective study-signal sampling rate from timestamps. The displayed detected signal rate is the rate used for frequency analysis.

If the file was downsampled before export, the original acquisition rate may be stored as metadata in Frequency settings, but it does not replace the sampling rate of the current file.

## 8. Export

Use the plot export controls or the centralized **Export…** workflow to save figures and related data.

For scientific reproducibility, keep the exported metadata alongside figures used in analysis or publication.

## Next guides

- [Importing data](importing-data.md)
- [Analysis workflow](analysis-workflow.md)
- [Frequency analysis](frequency-analysis.md)
- [Troubleshooting](troubleshooting.md)
