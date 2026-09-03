# Frequency analysis

FiberLyse's frequency tools are designed to work with recordings produced at different sampling rates.

## Current sampling rate

The analysis sampling rate is calculated from the timestamps of the current study signal:

```text
Fs ≈ 1 / median(diff(time))
```

This current-file rate is the value used for digital filtering and frequency limits.

## Nyquist frequency

The highest representable frequency is approximately:

```text
Nyquist = Fs / 2
```

If a requested band's upper edge is outside the current file's measurable range, FiberLyse marks the band unavailable instead of silently changing the requested band.

This is especially important for data that were downsampled before being exported to CSV.

## Band modes

### FiberLyse standard bands

The standard ranges remain constant between recordings so experiments can be compared using the same definitions. A range that is impossible for a low-sampling-rate file is shown as unavailable.

### Automatic bands for this recording

Automatic mode creates logarithmically spaced bands inside the measurable range of the current recording.

This is useful for exploration across hardware with different sampling rates, but automatically changing band definitions should be considered carefully when comparing experimental groups.

### Custom bands

Custom mode lets the user enter explicit low/high frequency ranges in Hz. Each band is validated against the current sampling rate and available continuous data.

## Band-pass filtering

With SciPy available, FiberLyse uses Butterworth band-pass filtering in contiguous usable sections. The filter does not intentionally bridge NaN holes.

The implementation also considers large timestamp discontinuities when creating continuous sections for frequency processing so an acquisition gap is not treated as one normal sample interval.

## Continuous-data requirement

Low-frequency bands require enough continuous data to contain multiple cycles. FiberLyse applies a minimum-length rule based on the lower frequency edge rather than returning a filtered trace from an obviously too-short segment.

## Power spectral density

When SciPy is available, the frequency overview uses Welch PSD on the longest continuous finite segment.

The PSD window is specified in **seconds**, not samples, so the same setting has the same time meaning for recordings acquired at different rates.

## Original acquisition rate and downsampling metadata

Frequency settings can store information such as:

- original acquisition sampling rate;
- downsampling factor;
- whether anti-alias/low-pass filtering before downsampling is known.

These values are metadata only. They do **not** restore samples or frequencies absent from the current file.

Example:

```text
Current CSV Fs:          ~10.14 Hz
Current Nyquist:         ~5.07 Hz
Original acquisition Fs: ~1014 Hz
Downsampling factor:     ~100
```

The current CSV can still only support frequencies below its ~5.07 Hz Nyquist limit.

## Interpolation

The frequency tab intentionally uses the non-interpolated ΔF/F representation. This keeps artifact holes visible to the filtering workflow rather than treating reconstructed points as originally measured frequency content.

## Reporting frequency analyses

For reproducibility, report or retain:

- current effective Fs;
- requested band(s);
- filter type/order;
- PSD method and window duration when used;
- downsampling history when known;
- whether the requested band was available for every recording in a group comparison.
