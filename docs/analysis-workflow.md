# Analysis workflow

This page describes the main scientific stages exposed by FiberLyse V25.

## 1. Import and recording-time handling

The selected import profile maps the source CSV to study and reference channels with timestamps. Recording-start handling and configured initial exclusion occur before the downstream artifact, fitting, normalization, frequency, AUC, and batch steps.

## 2. Artifact detection

FiberLyse detects abrupt changes from the derivative of each fluorescence trace.

The normal robust threshold uses the median absolute deviation (MAD) of the derivative. If the MAD is degenerate/zero, the implementation contains a standard-deviation fallback rather than silently disabling detection in that special case.

Options include:

- enable/disable artifact removal;
- threshold factor;
- additional sample padding around detected regions;
- require artifacts to be shared in time between study and reference channels.

Shared-artifact matching is performed before the extra padding is applied.

## 3. Artifact holes and interpolation

Detected artifact samples become missing values in the hole-preserving representation.

When linear interpolation is enabled, FiberLyse fills interior holes using time-aware interpolation. Edge holes are preserved rather than extrapolated beyond the first/last finite sample.

The frequency-analysis workflow uses the non-interpolated ΔF/F representation so removed sections are not silently reconstructed for band filtering.

## 4. Reference alignment

Reference samples are aligned to the study timestamps using the configured alignment mode. The normal mode is nearest-neighbor mapping without extrapolating beyond the available reference time support.

## 5. Linear reference fit

Within the selected fit interval(s), FiberLyse fits:

```text
study ≈ a × aligned_reference + b
```

where `a` is the regression coefficient and `b` is the intercept.

The Fit tab displays the study trace, fitted reference, fit window, coefficient, intercept, and R² so the fit can be inspected visually.

If the selected fit interval does not contain enough usable paired points, FiberLyse returns an unavailable fit rather than silently substituting a different interval.

## 6. ΔF/F

FiberLyse calculates:

```text
ΔF = study − fitted_reference
ΔF/F = ΔF / fitted_reference
```

The implementation requires a finite, positive fitted denominator before calculating ΔF/F.

Individual ΔF/F views use fractional units. Some batch displays transform ΔF/F to percent units for presentation; exported metadata should be retained with batch figures.

## 7. zF

Two zF options are available:

- global zF using the complete available ΔF/F distribution;
- interval-based zF using mean and standard deviation estimated from a selected time interval.

An invalid interval does not silently fall back to global zF.

## 8. Smoothing

Smoothed normalization traces use a centered moving average. Smoothing is applied independently to contiguous finite segments so NaN/artifact holes are not bridged.

The current V25 main smoothing control remains sample-based.

## 9. Frequency analysis

Frequency analysis works from non-interpolated ΔF/F and uses the sampling rate measured from timestamps. See [Frequency analysis](frequency-analysis.md).

## 10. AUC

The active graph can calculate interval AUC for visible labeled line traces. Available results include:

- signed area;
- absolute area;
- positive area;
- negative area;
- usable coverage;
- mean relative to a chosen baseline.

Missing sections are not integrated straight across as though data were present.

## 11. Batch compare and average

Batch comparison overlays selected traces. Batch averaging aligns recordings over shared time support without extrapolating outside each recording's timestamps and can display Group A / Group B means, SEM shading, and optional individual traces.

## 12. Export

FiberLyse can export plot data and analysis metadata. Keep metadata with figures when exact settings, import mappings, fit windows, normalization choices, and frequency parameters need to be reproduced.
