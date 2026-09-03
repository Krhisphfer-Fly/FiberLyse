# FiberLyse visual assets

This folder is for the images used by the GitHub repository and documentation.

## Recommended files

| Filename | Suggested content |
|---|---|
| `main-interface.png` | V25 main window with a representative analyzed recording loaded. |
| `csv-wizard.png` | CSV setup wizard showing a clean example mapping. |
| `artifact-removal.png` | Artifact-removal tab with raw/cleaned traces visible. |
| `control-fit.png` | Fit tab showing study signal and fitted reference. |
| `normalization.png` | ΔF/F or zF normalization view. |
| `frequency-analysis.png` | V24+/V25 frequency overview with PSD and valid/unavailable bands visible. |
| `batch-analysis.png` | Batch compare or Group A/Group B average view. |
| `fiberlyse-social-preview.png` | Optional GitHub social-preview graphic. |

## Screenshot style

For a cohesive repository page:

1. Use the same FiberLyse theme for the main screenshot set (light mode is usually clearest for documentation).
2. Use the same approximate window size for each capture, ideally around 1400–1600 px wide.
3. Load non-sensitive example data with clear traces.
4. Avoid showing personal folder paths, animal IDs, usernames, network drives, or unpublished metadata.
5. Crop away desktop clutter while keeping enough window chrome to make it obvious this is a desktop application.
6. Do not heavily annotate the screenshot itself; let the README caption explain what is being shown.

## Adding the images to the README

After uploading the image files, replace the commented screenshot block in the main `README.md` with something like:

```html
<p align="center">
  <img src="assets/main-interface.png" width="49%" alt="FiberLyse main interface">
  <img src="assets/frequency-analysis.png" width="49%" alt="FiberLyse frequency analysis">
</p>
```

Then add smaller feature-specific images further down only where they help explain the workflow.

## GitHub social preview

A separate wide social-preview image can be uploaded in the repository settings so links to FiberLyse look polished when shared. A 2:1 image around 1280 × 640 px works well for that purpose.
