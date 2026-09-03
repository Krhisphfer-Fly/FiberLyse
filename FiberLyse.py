# FiberLyse V24.1 - universal CSV import + adaptive frequency analysis; timestamp-only sampling rate
from __future__ import annotations;import sys;import os;import re;import argparse;import threading;import zipfile;from xml.sax.saxutils import escape as _xml_escape;from dataclasses import dataclass;from typing import Dict, List, Tuple, Optional, Callable, Any;import numpy as np;import pandas as pd
try:import tkinter as tk;from tkinter import ttk, filedialog, messagebox, simpledialog
except Exception as e:raise SystemExit(f'Tkinter is not available in this Python environment.\n\nOn Windows, Tkinter usually comes with the official python.org installer.\nIf your company Python build excludes Tkinter, you’ll need IT to install a Python distribution that includes Tk/Tcl.\n\nOriginal error: {e}')
import matplotlib;matplotlib.use('TkAgg');from matplotlib.figure import Figure;from matplotlib.widgets import SpanSelector;from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk;from matplotlib import colors as mcolors;from matplotlib.font_manager import FontProperties
import importlib
try:
    _scipy_signal = importlib.import_module('scipy.signal')
    butter = getattr(_scipy_signal, 'butter', None)
    sosfiltfilt = getattr(_scipy_signal, 'sosfiltfilt', None)
    filtfilt = getattr(_scipy_signal, 'filtfilt', None)
    _HAVE_SCIPY_SIGNAL = True if (_scipy_signal is not None) else False
except Exception:
    butter = None
    sosfiltfilt = None
    filtfilt = None
    _HAVE_SCIPY_SIGNAL = False
ISO_STATE = 1;EXC_STATE = 2;ANALYSIS_EXCLUDE_INITIAL_SECONDS = 10.0 * 60.0;DEFAULT_FIT_WINDOWS = [(ANALYSIS_EXCLUDE_INITIAL_SECONDS, 6500.0)];DEFAULT_ALIGN_MODE = 'nearest';DEFAULT_ACQ_FPS_HZ = 40.0;DEFAULT_BUTTER_ORDER = 2;DEFAULT_SMOOTH_WINDOW = 20;DEFAULT_ZF_INTERVAL_START_S = ANALYSIS_EXCLUDE_INITIAL_SECONDS;DEFAULT_ZF_INTERVAL_END_S = 6500.0;DEFAULT_USE_LINEAR_INTERP = True
# Compatibility alias retained for late V21 cache code. V22 excludes this period upstream for ALL analyses.
BATCH_HIDE_INITIAL_SECONDS = ANALYSIS_EXCLUDE_INITIAL_SECONDS
FIBERLYSE_VERSION = '22'
def _fontsize_points(v: Any, fallback: float) -> float:
    try:return float(FontProperties(size=v).get_size_in_points())
    except Exception:
        try:return float(v)
        except Exception:return float(fallback)
DEFAULT_AXIS_LABEL_FONTSIZE = _fontsize_points(matplotlib.rcParams.get('axes.labelsize', 10), fallback=10.0);DEFAULT_GRAPH_TITLE_FONTSIZE = _fontsize_points(matplotlib.rcParams.get('axes.titlesize', 12), fallback=12.0);DEFAULT_TICK_LABEL_FONTSIZE = _fontsize_points(matplotlib.rcParams.get('xtick.labelsize', matplotlib.rcParams.get('ytick.labelsize', 8)), fallback=8.0);NORMALIZATION_OPTION_DF_OVER_F = 'ΔF/F';NORMALIZATION_OPTION_ZF_GLOBAL = 'zF (global, GUI)';NORMALIZATION_OPTION_ZF_INTERVAL = 'zF - interval based';NORMALIZATION_ALL_OPTIONS = [NORMALIZATION_OPTION_DF_OVER_F, NORMALIZATION_OPTION_ZF_GLOBAL, NORMALIZATION_OPTION_ZF_INTERVAL];NORM_DFF = NORMALIZATION_OPTION_DF_OVER_F;NORM_ZF_GLOBAL = NORMALIZATION_OPTION_ZF_GLOBAL;NORM_ZF_INTERVAL = NORMALIZATION_OPTION_ZF_INTERVAL;NORM_CHOICES = NORMALIZATION_ALL_OPTIONS;FREQ_BANDS: List[Tuple[float, float]] = [(5.0, 10.0), (2.5, 5.0), (1.25, 2.5), (0.6, 1.25), (0.3, 0.6), (0.15, 0.3)]
def find_g_columns(columns) -> List[str]:
    """Return fluorescence columns named G, G0, G1, ... in source order.

    V21 accepted every column beginning with ``G`` which could accidentally
    ingest metadata such as ``Gain`` or ``Group``. V22 intentionally accepts
    only the Bonsai/Neurophotometrics-style G channel names.
    """
    return [c for c in columns if re.fullmatch(r'G\d*', str(c))]
def detect_artifacts_by_derivative(t: np.ndarray, y: np.ndarray, factor: float, method: str='mad', pad: int=1) -> np.ndarray:
    t = np.asarray(t, dtype=float);y = np.asarray(y, dtype=float);valid = np.isfinite(t) & np.isfinite(y)
    if valid.sum() < 3:return np.zeros_like(y, dtype=bool)
    tv = t[valid];yv = y[valid];dt = np.diff(tv);dy = np.diff(yv);dt_safe = np.where(dt == 0, np.nan, dt);slopes = dy / dt_safe;slopes_valid = slopes[np.isfinite(slopes)]
    if slopes_valid.size < 3:return np.zeros_like(y, dtype=bool)
    method = (method or 'mad').strip().lower()
    if method != 'mad':raise ValueError("Only 'mad' artifact method is supported (sd removed).")
    center = float(np.median(slopes_valid));mad = float(np.median(np.abs(slopes_valid - center)));scale = float(1.4826 * mad)
    if not np.isfinite(scale) or scale <= 0:
        # Degenerate robust case (commonly a flat/quantized trace). A zero MAD
        # previously disabled artifact detection completely. Use SD only as a
        # fallback scale when it is finite/non-zero; MAD remains the normal path.
        try:scale = float(np.std(slopes_valid, ddof=1))
        except Exception:scale = np.nan
        if not np.isfinite(scale) or scale <= 0:return np.zeros_like(y, dtype=bool)
    bad_seg = np.abs(slopes - center) > factor * scale;art_valid = np.zeros_like(yv, dtype=bool);bad_idx = np.where(bad_seg)[0]
    for k in bad_idx:i0 = max(0, k - pad);i1 = min(len(art_valid) - 1, k + 1 + pad);art_valid[i0:i1 + 1] = True
    art = np.zeros_like(y, dtype=bool);art[np.where(valid)[0]] = art_valid;return art
def shared_artifacts_by_time(t_iso: np.ndarray, art_iso: np.ndarray, t_exc: np.ndarray, art_exc: np.ndarray, tol: Optional[float]=None) -> Tuple[np.ndarray, np.ndarray]:
    t_iso = np.asarray(t_iso, dtype=float);t_exc = np.asarray(t_exc, dtype=float);art_iso = np.asarray(art_iso, dtype=bool);art_exc = np.asarray(art_exc, dtype=bool)
    if not art_iso.any() or not art_exc.any():return (np.zeros_like(art_iso, dtype=bool), np.zeros_like(art_exc, dtype=bool))
    if tol is None:dt_iso = np.median(np.diff(t_iso)) if t_iso.size > 1 else np.inf;dt_exc = np.median(np.diff(t_exc)) if t_exc.size > 1 else np.inf;base_dt = min(dt_iso, dt_exc);tol = 0.5 * base_dt if np.isfinite(base_dt) else 0.0
    iso_idx = np.where(art_iso)[0];exc_idx = np.where(art_exc)[0];iso_times = t_iso[iso_idx];exc_times = t_exc[exc_idx];iso_shared = np.zeros_like(art_iso, dtype=bool);exc_shared = np.zeros_like(art_exc, dtype=bool)
    for ei, te in zip(exc_idx, exc_times):
        if np.any(np.abs(iso_times - te) <= tol):exc_shared[ei] = True
    for ii, ti in zip(iso_idx, iso_times):
        if np.any(np.abs(exc_times - ti) <= tol):iso_shared[ii] = True
    return (iso_shared, exc_shared)
def expand_artifact_mask(mask: np.ndarray, pad: int) -> np.ndarray:
    """Expand an already-detected sample mask by ``pad`` samples on each side.

    A point-to-point derivative event naturally implicates the two samples
    bounding that transition. ``pad`` therefore means *additional* neighboring
    samples per side. Shared-artifact matching is performed before this
    expansion so padding cannot create a false cross-channel match.
    """
    mask = np.asarray(mask, dtype=bool)
    try:
        p = max(0, int(pad))
    except Exception:
        p = 0
    if p == 0 or mask.size == 0 or not mask.any():
        return mask.copy()
    out = mask.copy()
    idx = np.where(mask)[0]
    for i in idx:
        out[max(0, i - p):min(mask.size, i + p + 1)] = True
    return out

def remove_with_holes(y: np.ndarray, artifact_mask: np.ndarray) -> np.ndarray:y = np.asarray(y, dtype=float).copy();y[np.asarray(artifact_mask, dtype=bool)] = np.nan;return y
def linear_interpolate_by_time(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Linearly fill interior NaN holes using time; preserve edge holes."""
    t = np.asarray(t, dtype=float);y = np.asarray(y, dtype=float);out = y.copy()
    valid = np.isfinite(t) & np.isfinite(out)
    if valid.sum() < 2:return out
    idx = np.where(valid)[0];first = int(idx[0]);last = int(idx[-1])
    interior = np.arange(first, last + 1, dtype=int)
    finite_t = t[valid]
    if np.any(np.diff(finite_t) <= 0):return out
    out[interior] = np.interp(t[interior], finite_t, out[valid]).astype(float)
    return out
def align_iso_to_exc_no_interp(t_iso: np.ndarray, y_iso: np.ndarray, t_exc: np.ndarray, mode: str='nearest') -> np.ndarray:
    """Map isosbestic samples to excitatory timestamps without extrapolation."""
    t_iso = np.asarray(t_iso, dtype=float);y_iso = np.asarray(y_iso, dtype=float);t_exc = np.asarray(t_exc, dtype=float)
    out = np.full_like(t_exc, np.nan, dtype=float)
    if t_iso.size == 0:return out
    finite_iso_t = np.isfinite(t_iso)
    if not finite_iso_t.all() or np.any(np.diff(t_iso) <= 0):return out
    in_support = np.isfinite(t_exc) & (t_exc >= t_iso[0]) & (t_exc <= t_iso[-1])
    if not np.any(in_support):return out
    tq = t_exc[in_support];idx = np.searchsorted(t_iso, tq, side='left');mode = (mode or 'nearest').strip().lower()
    if mode == 'prev':j = idx - 1
    elif mode == 'next':j = idx
    elif mode == 'nearest':
        j_prev = np.clip(idx - 1, 0, len(t_iso) - 1);j_next = np.clip(idx, 0, len(t_iso) - 1)
        d_prev = np.abs(tq - t_iso[j_prev]);d_next = np.abs(tq - t_iso[j_next]);j = np.where(d_next < d_prev, j_next, j_prev)
    else:raise ValueError("align mode must be 'prev', 'next', or 'nearest'")
    j = np.clip(j, 0, len(t_iso) - 1);out[in_support] = y_iso[j]
    return out
def fit_linear(y: np.ndarray, x: np.ndarray) -> Tuple[float, float]:X = np.vstack([x, np.ones_like(x)]).T;(a, b), *_ = np.linalg.lstsq(X, y, rcond=None);return (float(a), float(b))
def r2_score(y: np.ndarray, yhat: np.ndarray) -> float:ss_res = np.nansum((y - yhat) ** 2);ss_tot = np.nansum((y - np.nanmean(y)) ** 2);return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
def _fill_nans_linear_1d(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 0:return x.copy()
    out = x.copy();finite = np.isfinite(out)
    if not np.any(finite):return out
    idx = np.where(finite)[0]
    if idx.size == 1:out[~finite] = out[idx[0]];return out
    xp = idx.astype(float);fp = out[idx].astype(float);x_all = np.arange(out.size, dtype=float);out = np.interp(x_all, xp, fp).astype(float);return out
def _moving_average_centered_fast(x: np.ndarray, window_size: int) -> np.ndarray:
    x = np.asarray(x, dtype=float);n = x.size
    if n == 0:return x.copy()
    w = int(max(1, window_size))
    if w == 1:return x.copy()
    left = (w - 1) // 2;right = w // 2;idx = np.arange(n, dtype=int);start = idx - left;end = idx + right;start = np.clip(start, 0, n - 1);end = np.clip(end, 0, n - 1);prefix = np.concatenate([[0.0], np.cumsum(x, dtype=float)]);sums = prefix[end + 1] - prefix[start];counts = (end - start + 1).astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):y = sums / counts
    return np.asarray(y, dtype=float)
def smooth_like_batch(x: np.ndarray, window_size: int) -> np.ndarray:
    """Centered moving average that never bridges NaN holes.

    V21 switched between filtfilt and a centered moving average depending on
    SciPy/data length. V22 deliberately uses one method everywhere so the same
    input and window always produce the same result.
    """
    x = np.asarray(x, dtype=float);n = x.size
    if n == 0:return x.copy()
    try:w = max(1, int(window_size))
    except Exception:w = DEFAULT_SMOOTH_WINDOW
    if w == 1:return x.copy()
    out = np.full_like(x, np.nan, dtype=float)
    finite = np.isfinite(x)
    for i0, i1 in _contiguous_true_runs(finite):
        seg = x[i0:i1 + 1]
        out[i0:i1 + 1] = _moving_average_centered_fast(seg, w)
    return out
def _contiguous_true_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0:return []
    idx = np.where(mask)[0]
    if idx.size == 0:return []
    runs: List[Tuple[int, int]] = [];s = idx[0];p = idx[0]
    for k in idx[1:]:
        if k == p + 1:p = k
        else:runs.append((s, p));s = p = k
    runs.append((s, p));return runs
def bandpass_butterworth_segmentwise_no_interp(x: np.ndarray, low_hz: float, high_hz: float, fs: float, order: int=DEFAULT_BUTTER_ORDER) -> np.ndarray:
    x = np.asarray(x, dtype=float);out = np.full_like(x, np.nan, dtype=float)
    if not _HAVE_SCIPY_SIGNAL or butter is None or sosfiltfilt is None:return out
    if not np.isfinite(fs) or fs <= 0:return out
    nyq = fs / 2.0;high_hz = min(float(high_hz), nyq - 1e-06)
    if low_hz <= 0 or high_hz <= 0 or low_hz >= high_hz:return out
    low = float(low_hz) / nyq;high = float(high_hz) / nyq
    try:sos = butter(int(order), [low, high], btype='band', output='sos')
    except Exception:return out
    finite = np.isfinite(x)
    for i0, i1 in _contiguous_true_runs(finite):
        seg = x[i0:i1 + 1]
        min_samples = max(8, int(np.ceil((3.0 * float(fs)) / float(low_hz))) + 1)
        if seg.size < min_samples:continue
        try:yseg = sosfiltfilt(sos, seg)
        except Exception:
            try:yseg = sosfiltfilt(sos, seg, padlen=0)
            except Exception:continue
        out[i0:i1 + 1] = yseg
    return out
def bandpass_fft_no_interp(x: np.ndarray, low_hz: float, high_hz: float, fs: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if not np.isfinite(fs) or fs <= 0:return np.full_like(x, np.nan, dtype=float)
    if not np.all(np.isfinite(x)):return np.full_like(x, np.nan, dtype=float)
    nyq = fs / 2.0;high_hz = min(float(high_hz), nyq - 1e-06)
    if low_hz <= 0 or high_hz <= 0 or low_hz >= high_hz:return np.full_like(x, np.nan, dtype=float)
    n = x.size
    if n < 4:return np.full_like(x, np.nan, dtype=float)
    X = np.fft.rfft(x);freqs = np.fft.rfftfreq(n, d=1.0 / fs);mask = (freqs >= float(low_hz)) & (freqs <= float(high_hz));Xf = np.where(mask, X, 0.0);y = np.fft.irfft(Xf, n=n);return np.asarray(y, dtype=float)
def estimate_fs_from_t(t: np.ndarray) -> float:
    t = np.asarray(t, dtype=float);dt = np.diff(t);dt = dt[np.isfinite(dt) & (dt > 0)]
    if dt.size == 0:return np.nan
    med = float(np.median(dt))
    if med <= 0 or not np.isfinite(med):return np.nan
    return float(1.0 / med)
def _nanstd_safe(x: np.ndarray, ddof: int) -> float:
    x = np.asarray(x, dtype=float);x = x[np.isfinite(x)]
    if x.size < 2:return np.nan
    try:return float(np.std(x, ddof=int(ddof)))
    except Exception:return np.nan
def zscore_global_gui(dff: np.ndarray, ddof: int=0) -> np.ndarray:
    dff = np.asarray(dff, dtype=float);mu = float(np.nanmean(dff));sigma = _nanstd_safe(dff, ddof=ddof)
    if np.isfinite(sigma) and sigma > 1e-12:return (dff - mu) / sigma
    return np.full_like(dff, np.nan, dtype=float)
def zscore_interval_based(dff: np.ndarray, t: np.ndarray, start_s: float, end_s: float, ddof: int=1) -> np.ndarray:
    dff = np.asarray(dff, dtype=float);t = np.asarray(t, dtype=float);a = float(min(start_s, end_s));b = float(max(start_s, end_s))
    mask = (t >= a) & (t <= b) & np.isfinite(dff)
    if mask.sum() < 2:return np.full_like(dff, np.nan, dtype=float)
    mu = float(np.nanmean(dff[mask]));sigma = _nanstd_safe(dff[mask], ddof=ddof)
    if np.isfinite(sigma) and sigma > 1e-12:return (dff - mu) / sigma
    return np.full_like(dff, np.nan, dtype=float)
def get_norm_array(res: 'ChannelResult', norm_mode: str) -> np.ndarray:
    if norm_mode == NORM_DFF:return np.asarray(res.dFF, dtype=float)
    if norm_mode == NORM_ZF_GLOBAL:return np.asarray(res.zF_global, dtype=float)
    if norm_mode == NORM_ZF_INTERVAL:return np.asarray(res.zF_interval, dtype=float)
    return np.asarray(res.dFF, dtype=float)
def get_smoothed_norm_array(res: 'ChannelResult', norm_mode: str, window_size: int) -> np.ndarray:
    try:w = int(window_size)
    except Exception:w = DEFAULT_SMOOTH_WINDOW
    w = max(1, w);cache = getattr(res, '_smooth_cache', None)
    if cache is None or not isinstance(cache, dict):cache = {};setattr(res, '_smooth_cache', cache)
    key = (str(norm_mode), int(w))
    if key in cache:return cache[key]
    y = get_norm_array(res, norm_mode);y_s = smooth_like_batch(y, window_size=w)
    if len(cache) > 12:cache.clear()
    cache[key] = y_s;return y_s
@dataclass
class ChannelResult:gcol: str;source_path: str;analysis_exclude_initial_seconds: float;artifact_enabled: bool;artifact_factor: float;artifact_pad: int;require_shared: bool;align_mode: str;t_iso: np.ndarray;t_exc: np.ndarray;iso_raw: np.ndarray;exc_raw: np.ndarray;art_iso: np.ndarray;art_exc: np.ndarray;iso_clean_holes: np.ndarray;exc_clean_holes: np.ndarray;iso_clean_interp: np.ndarray;exc_clean_interp: np.ndarray;iso_on_exc_holes: np.ndarray;iso_on_exc_interp: np.ndarray;use_interpolation: bool;iso_clean: np.ndarray;exc_clean: np.ndarray;iso_on_exc: np.ndarray;windows: List[Tuple[float, float]];acq_fps_hz: float;eff_fs_hz: float;slope: float;intercept: float;r2: float;fitted_iso_on_exc: np.ndarray;residual: np.ndarray;dF: np.ndarray;dFF: np.ndarray;slope_nointerp: float;intercept_nointerp: float;r2_nointerp: float;fitted_iso_on_exc_nointerp: np.ndarray;residual_nointerp: np.ndarray;dF_nointerp: np.ndarray;dFF_nointerp: np.ndarray;smooth_window: int;zf_interval_start_s: float;zf_interval_end_s: float;zF_global: np.ndarray;zF_interval: np.ndarray
def recompute_normalizations(res: ChannelResult) -> None:
    res.zF_global = zscore_global_gui(res.dFF, ddof=0);res.zF_interval = zscore_interval_based(res.dFF, res.t_exc, start_s=res.zf_interval_start_s, end_s=res.zf_interval_end_s, ddof=1);cache = getattr(res, '_smooth_cache', None)
    if isinstance(cache, dict):cache.clear()
    else:setattr(res, '_smooth_cache', {})
def _compute_fit_and_downstream(t_exc: np.ndarray, exc: np.ndarray, iso_on_exc: np.ndarray, windows: List[Tuple[float, float]]) -> Tuple[float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    t_exc = np.asarray(t_exc, dtype=float);exc = np.asarray(exc, dtype=float);iso_on_exc = np.asarray(iso_on_exc, dtype=float)
    fit_mask = np.zeros_like(t_exc, dtype=bool)
    for a, b in windows:fit_mask |= (t_exc >= float(min(a, b))) & (t_exc <= float(max(a, b)))
    finite = np.isfinite(exc) & np.isfinite(iso_on_exc);fit_idx = fit_mask & finite
    if fit_idx.sum() >= 2:
        a_fit, b_fit = fit_linear(exc[fit_idx], iso_on_exc[fit_idx]);fitted = a_fit * iso_on_exc + b_fit;r2 = r2_score(exc[fit_idx], fitted[fit_idx])
    else:
        a_fit, b_fit, r2 = (np.nan, np.nan, np.nan);fitted = np.full_like(iso_on_exc, np.nan, dtype=float)
    residual = exc - fitted;dF = residual
    with np.errstate(divide='ignore', invalid='ignore'):
        dFF = np.divide(dF, fitted, out=np.full_like(dF, np.nan), where=np.isfinite(fitted) & (fitted > 1e-12))
    return (float(a_fit), float(b_fit), float(r2), fitted, residual, dF, dFF)
def recompute_fit_and_downstream(res: ChannelResult) -> None:
    a, b, r2, fitted, resid, dF, dFF = _compute_fit_and_downstream(res.t_exc, res.exc_clean, res.iso_on_exc, res.windows);res.slope = a;res.intercept = b;res.r2 = r2;res.fitted_iso_on_exc = fitted;res.residual = resid;res.dF = dF;res.dFF = dFF;a2, b2, r22, fitted2, resid2, dF2, dFF2 = _compute_fit_and_downstream(res.t_exc, res.exc_clean_holes, res.iso_on_exc_holes, res.windows);res.slope_nointerp = a2;res.intercept_nointerp = b2;res.r2_nointerp = r22;res.fitted_iso_on_exc_nointerp = fitted2;res.residual_nointerp = resid2;res.dF_nointerp = dF2;res.dFF_nointerp = dFF2;recompute_normalizations(res);res._data_version = int(getattr(res, '_data_version', 0)) + 1;fcache = getattr(res, '_freq_cache', None)
    if isinstance(fcache, dict):fcache.clear()
    else:setattr(res, '_freq_cache', {})
def set_interpolation_mode(res: ChannelResult, use_interpolation: bool) -> None:
    res.use_interpolation = bool(use_interpolation)
    if res.use_interpolation:res.iso_clean = res.iso_clean_interp;res.exc_clean = res.exc_clean_interp;res.iso_on_exc = res.iso_on_exc_interp
    else:res.iso_clean = res.iso_clean_holes;res.exc_clean = res.exc_clean_holes;res.iso_on_exc = res.iso_on_exc_holes
    recompute_fit_and_downstream(res)
def recompute_artifact_pipeline_inplace(res: ChannelResult, artifact_enabled: bool, artifact_factor: float, artifact_pad: int, require_shared: bool, align_mode: str, use_linear_interp: bool) -> None:
    res.artifact_enabled = bool(artifact_enabled);res.artifact_factor = float(artifact_factor);res.artifact_pad = int(artifact_pad);res.require_shared = bool(require_shared);res.align_mode = str(align_mode)
    if artifact_enabled:
        base_iso = detect_artifacts_by_derivative(res.t_iso, res.iso_raw, factor=float(artifact_factor), method='mad', pad=0)
        base_exc = detect_artifacts_by_derivative(res.t_exc, res.exc_raw, factor=float(artifact_factor), method='mad', pad=0)
        if require_shared:base_iso, base_exc = shared_artifacts_by_time(res.t_iso, base_iso, res.t_exc, base_exc)
        art_iso = expand_artifact_mask(base_iso, int(artifact_pad));art_exc = expand_artifact_mask(base_exc, int(artifact_pad))
    else:
        art_iso = np.zeros_like(res.iso_raw, dtype=bool);art_exc = np.zeros_like(res.exc_raw, dtype=bool)
    res.art_iso = np.asarray(art_iso, dtype=bool);res.art_exc = np.asarray(art_exc, dtype=bool)
    res.iso_clean_holes = remove_with_holes(res.iso_raw, res.art_iso);res.exc_clean_holes = remove_with_holes(res.exc_raw, res.art_exc)
    res.iso_clean_interp = linear_interpolate_by_time(res.t_iso, res.iso_clean_holes);res.exc_clean_interp = linear_interpolate_by_time(res.t_exc, res.exc_clean_holes)
    res.iso_on_exc_holes = align_iso_to_exc_no_interp(res.t_iso, res.iso_clean_holes, res.t_exc, mode=align_mode);res.iso_on_exc_interp = align_iso_to_exc_no_interp(res.t_iso, res.iso_clean_interp, res.t_exc, mode=align_mode)
    res.iso_on_exc_holes = np.asarray(res.iso_on_exc_holes, dtype=float);res.iso_on_exc_interp = np.asarray(res.iso_on_exc_interp, dtype=float);res.iso_on_exc_holes[res.art_exc] = np.nan
    set_interpolation_mode(res, use_interpolation=bool(use_linear_interp))
def analyze_csv(csv_path: str, artifact_enabled: bool, artifact_factor: float, artifact_method: str, artifact_pad: int, require_shared: bool, align_mode: str, fit_windows: List[Tuple[float, float]], acq_fps_hz: Optional[float], smooth_window: int, zf_interval_start_s: float, zf_interval_end_s: float, use_linear_interp: bool) -> Dict[str, ChannelResult]:
    header = pd.read_csv(csv_path, nrows=0)
    for col in ['SystemTimestamp', 'LedState']:
        if col not in header.columns:raise ValueError(f"Missing required column '{col}'.")
    g_cols = find_g_columns(header.columns)
    if not g_cols:raise ValueError('No G* fluorescence columns found (expected names such as G0, G1...).')
    usecols = ['SystemTimestamp', 'LedState'] + list(g_cols);df = pd.read_csv(csv_path, usecols=usecols)
    if df.empty:raise ValueError('CSV contains no data rows.')
    try:timestamps = pd.to_numeric(df['SystemTimestamp'], errors='raise').to_numpy(dtype=float)
    except Exception as e:raise ValueError(f'SystemTimestamp must be numeric: {e}')
    if not np.all(np.isfinite(timestamps)):raise ValueError('SystemTimestamp contains NaN or infinite values.')
    led = pd.to_numeric(df['LedState'], errors='coerce').to_numpy()
    start_candidates = np.flatnonzero(led == 7)
    start_idx = int(start_candidates[0]) if start_candidates.size else 0
    # Remove everything before the recording-start marker rather than retaining negative-time samples.
    df = df.iloc[start_idx:].reset_index(drop=True);timestamps = timestamps[start_idx:];led = led[start_idx:]
    if timestamps.size < 2:raise ValueError('Not enough samples after the recording-start marker.')
    dt = np.diff(timestamps)
    if np.any(~np.isfinite(dt)) or np.any(dt <= 0):
        n_bad = int(np.sum(~np.isfinite(dt) | (dt <= 0)))
        raise ValueError(f'SystemTimestamp must be strictly increasing after recording start; found {n_bad} duplicate/decreasing interval(s).')
    t0 = float(timestamps[0]);t_full = timestamps - t0
    # V22 scientific rule: the initial non-biological period is excluded BEFORE
    # LED separation, artifact statistics, fitting, normalization, filtering,
    # AUC, and batch analysis. Surviving timestamps retain recording-time values.
    keep = np.isfinite(t_full) & (t_full >= float(ANALYSIS_EXCLUDE_INITIAL_SECONDS))
    if int(np.sum(keep)) < 4:
        raise ValueError(f'No usable recording remains after excluding the first {ANALYSIS_EXCLUDE_INITIAL_SECONDS:g} s.')
    df = df.loc[keep].reset_index(drop=True);t_full = t_full[keep];led = led[keep]
    mask_iso = led == ISO_STATE;mask_exc = led == EXC_STATE;t_iso = t_full[mask_iso];t_exc = t_full[mask_exc]
    if t_iso.size == 0 or t_exc.size == 0:raise ValueError('No iso/exc samples remain after the initial exclusion. Check LedState coding and recording length.')
    measured_eff_fs = estimate_fs_from_t(t_exc)
    # V24.1: the current signal sampling rate is defined only by the timestamps
    # present in the file. Acquisition/camera FPS is provenance metadata, not an
    # analysis input, and is therefore never used as a sampling-rate fallback.
    supplied_acq = np.nan
    eff_fs = measured_eff_fs if np.isfinite(measured_eff_fs) and measured_eff_fs > 0 else np.nan
    acq_fps_val = np.nan
    results: Dict[str, ChannelResult] = {}
    for gcol in g_cols:
        iso_raw = pd.to_numeric(df.loc[mask_iso, gcol], errors='coerce').to_numpy(dtype=float);exc_raw = pd.to_numeric(df.loc[mask_exc, gcol], errors='coerce').to_numpy(dtype=float)
        if iso_raw.size == 0 or exc_raw.size == 0:continue
        if artifact_enabled:
            base_iso = detect_artifacts_by_derivative(t_iso, iso_raw, factor=artifact_factor, method='mad', pad=0);base_exc = detect_artifacts_by_derivative(t_exc, exc_raw, factor=artifact_factor, method='mad', pad=0)
            if require_shared:base_iso, base_exc = shared_artifacts_by_time(t_iso, base_iso, t_exc, base_exc)
            art_iso = expand_artifact_mask(base_iso, artifact_pad);art_exc = expand_artifact_mask(base_exc, artifact_pad)
        else:
            art_iso = np.zeros_like(iso_raw, dtype=bool);art_exc = np.zeros_like(exc_raw, dtype=bool)
        iso_clean_holes = remove_with_holes(iso_raw, art_iso);exc_clean_holes = remove_with_holes(exc_raw, art_exc)
        iso_clean_interp = linear_interpolate_by_time(t_iso, iso_clean_holes);exc_clean_interp = linear_interpolate_by_time(t_exc, exc_clean_holes)
        iso_on_exc_holes = align_iso_to_exc_no_interp(t_iso, iso_clean_holes, t_exc, mode=align_mode);iso_on_exc_interp = align_iso_to_exc_no_interp(t_iso, iso_clean_interp, t_exc, mode=align_mode)
        iso_on_exc_holes = np.asarray(iso_on_exc_holes, dtype=float);iso_on_exc_interp = np.asarray(iso_on_exc_interp, dtype=float);iso_on_exc_holes[np.asarray(art_exc, dtype=bool)] = np.nan
        if use_linear_interp:iso_clean = iso_clean_interp;exc_clean = exc_clean_interp;iso_on_exc = iso_on_exc_interp
        else:iso_clean = iso_clean_holes;exc_clean = exc_clean_holes;iso_on_exc = iso_on_exc_holes
        res = ChannelResult(gcol=gcol, source_path=os.path.abspath(str(csv_path)), analysis_exclude_initial_seconds=float(ANALYSIS_EXCLUDE_INITIAL_SECONDS), artifact_enabled=bool(artifact_enabled), artifact_factor=float(artifact_factor), artifact_pad=int(artifact_pad), require_shared=bool(require_shared), align_mode=str(align_mode), t_iso=t_iso, t_exc=t_exc, iso_raw=iso_raw, exc_raw=exc_raw, art_iso=np.asarray(art_iso, dtype=bool), art_exc=np.asarray(art_exc, dtype=bool), iso_clean_holes=iso_clean_holes, exc_clean_holes=exc_clean_holes, iso_clean_interp=iso_clean_interp, exc_clean_interp=exc_clean_interp, iso_on_exc_holes=iso_on_exc_holes, iso_on_exc_interp=iso_on_exc_interp, use_interpolation=bool(use_linear_interp), iso_clean=iso_clean, exc_clean=exc_clean, iso_on_exc=iso_on_exc, windows=list(fit_windows), acq_fps_hz=float(acq_fps_val), eff_fs_hz=float(eff_fs), slope=np.nan, intercept=np.nan, r2=np.nan, fitted_iso_on_exc=np.full_like(t_exc, np.nan, dtype=float), residual=np.full_like(t_exc, np.nan, dtype=float), dF=np.full_like(t_exc, np.nan, dtype=float), dFF=np.full_like(t_exc, np.nan, dtype=float), slope_nointerp=np.nan, intercept_nointerp=np.nan, r2_nointerp=np.nan, fitted_iso_on_exc_nointerp=np.full_like(t_exc, np.nan, dtype=float), residual_nointerp=np.full_like(t_exc, np.nan, dtype=float), dF_nointerp=np.full_like(t_exc, np.nan, dtype=float), dFF_nointerp=np.full_like(t_exc, np.nan, dtype=float), smooth_window=max(1, int(smooth_window)), zf_interval_start_s=float(zf_interval_start_s), zf_interval_end_s=float(zf_interval_end_s), zF_global=np.full_like(t_exc, np.nan, dtype=float), zF_interval=np.full_like(t_exc, np.nan, dtype=float))
        recompute_fit_and_downstream(res);results[gcol] = res
    if not results:raise ValueError('No usable G* channels found to plot.')
    return results
class PlotTabTk(ttk.Frame):
    def __init__(self, master, tab_name: str, default_filename_prefix: str='', figsize: Tuple[float, float]=(7.2, 4.6), dpi: int=110):super().__init__(master);self.tab_name = tab_name;self.default_filename_prefix = default_filename_prefix;self.export_provider: Optional[Callable[[], Dict[str, pd.DataFrame]]] = None;self.fig = Figure(figsize=figsize, dpi=dpi);self.ax = self.fig.add_subplot(111);self.canvas = FigureCanvasTkAgg(self.fig, master=self);self.canvas_widget = self.canvas.get_tk_widget();self.toolbar = NavigationToolbar2Tk(self.canvas, self);self.toolbar.update();btn_row = ttk.Frame(self);btn_row.pack(side=tk.BOTTOM, fill=tk.X, padx=6, pady=6);self.save_btn = ttk.Button(btn_row, text='Save this graph…', command=self.save_plot);self.save_btn.pack(side=tk.RIGHT);self.export_btn = ttk.Button(btn_row, text='Export data (Excel)…', command=self.export_excel);self.export_btn.pack(side=tk.RIGHT, padx=(0, 8));self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True);self.axis_label_fontsize: Optional[float] = None;self.graph_title_fontsize: Optional[float] = None;self.tick_label_fontsize: Optional[float] = None;self._suptitle_override: Optional[str] = None;self._title_overrides: Dict[int, str] = {};self._xlabel_overrides: Dict[int, str] = {};self._ylabel_overrides: Dict[int, str] = {};self._legend_title_overrides: Dict[int, str] = {};self._legend_label_overrides: Dict[int, Dict[str, str]] = {};self._color_overrides: Dict[int, Dict[str, str]] = {};self._artist_color_overrides: Dict[int, str] = {};self._time_markers: List[Dict[str, Any]] = [];self._time_marker_artists: List[Any] = [];self._cid_button_press = self.canvas.mpl_connect('button_press_event', self._on_mpl_button_press)
    def _ax_key(self, ax) -> int:
        try:return int(self.fig.axes.index(ax))
        except Exception:return id(ax)
    def save_plot(self):
        suggested = f'{self.default_filename_prefix}_{self.tab_name}'.strip('_');path = filedialog.asksaveasfilename(title='Save graph', defaultextension='.png', initialfile=f'{suggested}.png', filetypes=[('PNG', '*.png'), ('SVG', '*.svg'), ('PDF', '*.pdf'), ('All files', '*.*')])
        if not path:return
        self.fig.savefig(path, bbox_inches='tight')
    @staticmethod
    def _safe_sheet_name(name: str, used: set) -> str:
        name = re.sub('[:\\\\/?*\\[\\]]+', '_', str(name)).strip()
        if not name:name = 'Sheet'
        name = name[:31];base = name;k = 2
        while name in used:suffix = f'_{k}';name = (base[:max(1, 31 - len(suffix))] + suffix)[:31];k += 1
        used.add(name);return name
    def _payload_from_artists_fallback(self) -> Dict[str, pd.DataFrame]:
        payload: Dict[str, pd.DataFrame] = {};used: set = set()
        for ai, ax in enumerate(self.fig.axes, start=1):
            for li, line in enumerate(ax.get_lines(), start=1):
                x = np.asarray(line.get_xdata(), dtype=float);y = np.asarray(line.get_ydata(), dtype=float);label = line.get_label()
                if not label or label.startswith('_'):label = f'line{li}'
                sheet = self._safe_sheet_name(f'ax{ai}_{label}', used);payload[sheet] = pd.DataFrame({'x': x, 'y': y})
            for ci, coll in enumerate(getattr(ax, 'collections', []), start=1):
                try:offs = coll.get_offsets()
                except Exception:continue
                if offs is None or len(offs) == 0:continue
                offs = np.asarray(offs, dtype=float)
                if offs.ndim != 2 or offs.shape[1] < 2:continue
                sheet = self._safe_sheet_name(f'ax{ai}_scatter{ci}', used);payload[sheet] = pd.DataFrame({'x': offs[:, 0], 'y': offs[:, 1]})
        if not payload:payload[self._safe_sheet_name('empty', used)] = pd.DataFrame({'note': ['No plottable artists found.']})
        return payload
    @staticmethod
    def _excel_col_name(n_1based: int) -> str:
        n = int(n_1based);s = ''
        while n > 0:n, r = divmod(n - 1, 26);s = chr(65 + r) + s
        return s
    @staticmethod
    def _is_nan(v: Any) -> bool:
        try:return v is None or (isinstance(v, float) and (not np.isfinite(v))) or (isinstance(v, np.floating) and (not np.isfinite(v)))
        except Exception:return v is None
    def _write_worksheet_xml(self, zf: zipfile.ZipFile, sheet_path: str, df: pd.DataFrame) -> None:
        cols = list(df.columns);ncols = len(cols)
        def write_cell(fh, row_idx_1based: int, col_idx_1based: int, value: Any, force_str: bool=False):
            if self._is_nan(value):return
            col_letter = self._excel_col_name(col_idx_1based);cell_ref = f'{col_letter}{row_idx_1based}'
            if isinstance(value, (bool, np.bool_)) and (not force_str):v = '1' if bool(value) else '0';fh.write(f'<c r="{cell_ref}" t="b"><v>{v}</v></c>'.encode('utf-8'));return
            if not force_str and isinstance(value, (int, np.integer)):fh.write(f'<c r="{cell_ref}"><v>{int(value)}</v></c>'.encode('utf-8'));return
            if not force_str and isinstance(value, (float, np.floating)):
                if np.isfinite(value):fh.write(f'<c r="{cell_ref}"><v>{float(value)}</v></c>'.encode('utf-8'))
                return
            s = str(value)
            if len(s) > 32767:s = s[:32767]
            s = _xml_escape(s);fh.write(f'<c r="{cell_ref}" t="inlineStr"><is><t xml:space="preserve">{s}</t></is></c>'.encode('utf-8'))
        with zf.open(sheet_path, 'w') as fh:
            fh.write(b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>');fh.write(b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">');fh.write(b'<sheetData>');fh.write(b'<row r="1">')
            for j, c in enumerate(cols, start=1):write_cell(fh, 1, j, c, force_str=True)
            fh.write(b'</row>')
            for i, row in enumerate(df.itertuples(index=False, name=None), start=2):
                fh.write(f'<row r="{i}">'.encode('utf-8'))
                for j in range(ncols):write_cell(fh, i, j + 1, row[j], force_str=False)
                fh.write(b'</row>')
            fh.write(b'</sheetData></worksheet>')
    def _write_xlsx_minimal(self, path: str, payload: Dict[str, pd.DataFrame]) -> None:
        used: set = set();sheets: List[Tuple[str, pd.DataFrame]] = []
        for sheet_name, df in payload.items():
            safe = self._safe_sheet_name(sheet_name, used)
            if not isinstance(df, pd.DataFrame):df = pd.DataFrame(df)
            sheets.append((safe, df))
        ct_lines = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>', '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>', '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>', '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>', '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>']
        for i in range(1, len(sheets) + 1):ct_lines.append(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
        ct_lines.append('</Types>');content_types_xml = '\n'.join(ct_lines);rels_xml = '\n'.join(['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">', '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>', '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>', '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>', '</Relationships>']);core_xml = '\n'.join(['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">', '<dc:title>Fiberlyse Export</dc:title>', '<dc:creator>Fiberlyse</dc:creator>', '<cp:lastModifiedBy>Fiberlyse</cp:lastModifiedBy>', '</cp:coreProperties>']);app_xml = '\n'.join(['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">', '<Application>Fiberlyse</Application>', '</Properties>']);wb_lines = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">', '<sheets>']
        for i, (sheet_name, _df) in enumerate(sheets, start=1):name_esc = _xml_escape(sheet_name);wb_lines.append(f'<sheet name="{name_esc}" sheetId="{i}" r:id="rId{i}"/>')
        wb_lines += ['</sheets>', '</workbook>'];workbook_xml = '\n'.join(wb_lines);wb_rels_lines = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
        for i in range(1, len(sheets) + 1):wb_rels_lines.append(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>')
        wb_rels_lines.append(f'<Relationship Id="rId{len(sheets) + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>');wb_rels_lines.append('</Relationships>');workbook_rels_xml = '\n'.join(wb_rels_lines);styles_xml = '\n'.join(['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">', '<fonts count="1"><font><sz val="11"/><color theme="1"/><name val="Calibri"/><family val="2"/></font></fonts>', '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>', '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>', '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>', '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>', '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>', '</styleSheet>'])
        with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr('[Content_Types].xml', content_types_xml.encode('utf-8'));zf.writestr('_rels/.rels', rels_xml.encode('utf-8'));zf.writestr('docProps/core.xml', core_xml.encode('utf-8'));zf.writestr('docProps/app.xml', app_xml.encode('utf-8'));zf.writestr('xl/workbook.xml', workbook_xml.encode('utf-8'));zf.writestr('xl/_rels/workbook.xml.rels', workbook_rels_xml.encode('utf-8'));zf.writestr('xl/styles.xml', styles_xml.encode('utf-8'))
            for i, (_sheet_name, df) in enumerate(sheets, start=1):self._write_worksheet_xml(zf, f'xl/worksheets/sheet{i}.xml', df)
    def export_excel(self):
        suggested = f'{self.default_filename_prefix}_{self.tab_name}'.strip('_');path = filedialog.asksaveasfilename(title='Export data (Excel)', defaultextension='.xlsx', initialfile=f'{suggested}.xlsx', filetypes=[('Excel workbook', '*.xlsx'), ('All files', '*.*')])
        if not path:return
        try:
            if callable(self.export_provider):payload = self.export_provider()
            else:payload = self._payload_from_artists_fallback()
            if not isinstance(payload, dict) or not payload:raise ValueError('Export provider returned no data.')
            self._write_xlsx_minimal(path, payload);messagebox.showinfo('Export complete', f'Saved Excel file:\n{path}')
        except Exception as e:messagebox.showerror('Export failed', f'Could not export Excel data:\n\n{e}')
    def set_axis_label_fontsize(self, fontsize: Optional[float]) -> None:
        if fontsize is None:self.axis_label_fontsize = None;self.redraw();return
        try:fs = float(fontsize)
        except Exception:return
        if not np.isfinite(fs) or fs <= 0:return
        self.axis_label_fontsize = fs;self.redraw()
    def set_graph_title_fontsize(self, fontsize: Optional[float]) -> None:
        if fontsize is None:self.graph_title_fontsize = None;self.redraw();return
        try:fs = float(fontsize)
        except Exception:return
        if not np.isfinite(fs) or fs <= 0:return
        self.graph_title_fontsize = fs;self.redraw()
    def set_tick_label_fontsize(self, fontsize: Optional[float]) -> None:
        if fontsize is None:self.tick_label_fontsize = None;self.redraw();return
        try:fs = float(fontsize)
        except Exception:return
        if not np.isfinite(fs) or fs <= 0:return
        self.tick_label_fontsize = fs;self.redraw()
    def _toolbar_is_active(self) -> bool:
        try:return bool(getattr(self.toolbar, 'mode', ''))
        except Exception:return False
    @staticmethod
    def _normalize_hex_color(s: str) -> str:
        s = (s or '').strip()
        if not s:raise ValueError('Empty color')
        if not s.startswith('#'):s = '#' + s
        if re.fullmatch('#[0-9a-fA-F]{3}', s):r, g, b = (s[1], s[2], s[3]);return f'#{r}{r}{g}{g}{b}{b}'
        if re.fullmatch('#[0-9a-fA-F]{4}', s):r, g, b, a = (s[1], s[2], s[3], s[4]);return f'#{r}{r}{g}{g}{b}{b}{a}{a}'
        if re.fullmatch('#[0-9a-fA-F]{6}', s) or re.fullmatch('#[0-9a-fA-F]{8}', s):return s
        raise ValueError('Expected hex like #RRGGBB (or #RRGGBBAA)')
    @staticmethod
    def _artist_contains(artist, event) -> bool:
        try:hit, _ = artist.contains(event);return bool(hit)
        except Exception:return False
    @staticmethod
    def _resolve_label(label: str, mapping: Dict[str, str], max_hops: int=12) -> str:
        cur = str(label)
        for _ in range(max_hops):
            nxt = mapping.get(cur)
            if not nxt or nxt == cur:return cur
            cur = nxt
        return cur
    @staticmethod
    def _iter_axes_color_targets(ax):
        for art in list(getattr(ax, 'lines', [])):yield art
        for art in list(getattr(ax, 'collections', [])):yield art
        for art in list(getattr(ax, 'patches', [])):
            if art is getattr(ax, 'patch', None):continue
            yield art
    @staticmethod
    def _set_artist_color(artist, hex_color: str) -> None:
        if hasattr(artist, 'set_color'):
            try:artist.set_color(hex_color);return
            except Exception:pass
        if hasattr(artist, 'set_facecolor'):
            try:artist.set_facecolor(hex_color)
            except Exception:pass
        if hasattr(artist, 'set_edgecolor'):
            try:artist.set_edgecolor(hex_color)
            except Exception:pass
    @staticmethod
    def _artist_current_hex(artist) -> Optional[str]:
        try:
            if hasattr(artist, 'get_color'):return mcolors.to_hex(artist.get_color(), keep_alpha=False)
        except Exception:pass
        try:
            if hasattr(artist, 'get_facecolor'):
                fc = artist.get_facecolor()
                if fc is not None and len(fc):return mcolors.to_hex(fc[0], keep_alpha=False)
        except Exception:pass
        return None
    def _legend_handles(self, leg):
        handles = getattr(leg, 'legend_handles', None)
        if handles is None:handles = getattr(leg, 'legendHandles', [])
        return list(handles) if handles is not None else []
    def _hit_test_editable_text(self, event):
        fig = self.fig
        for ax in fig.axes:
            leg = ax.get_legend()
            if leg:
                for txt in leg.get_texts():
                    if self._artist_contains(txt, event):return ('legend_text', ax, txt)
                lt = leg.get_title()
                if lt and lt.get_text() is not None and self._artist_contains(lt, event):return ('legend_title', ax, lt)
        for ax in fig.axes:
            t = ax.title
            if t and t.get_text() is not None and self._artist_contains(t, event):return ('graph_title', ax, t)
            xl = ax.xaxis.label
            if xl and xl.get_text() is not None and self._artist_contains(xl, event):return ('xlabel', ax, xl)
            yl = ax.yaxis.label
            if yl and yl.get_text() is not None and self._artist_contains(yl, event):return ('ylabel', ax, yl)
        st = getattr(fig, '_suptitle', None)
        if st is not None and self._artist_contains(st, event):return ('suptitle', None, st)
        return None
    def _find_legend_entry_at_event(self, ax, event):
        leg = ax.get_legend()
        if not leg:return None
        texts = list(leg.get_texts());handles = self._legend_handles(leg)
        for i, txt in enumerate(texts):
            if self._artist_contains(txt, event):return ('legend', txt.get_text(), i)
        for i, h in enumerate(handles):
            if self._artist_contains(h, event):label = texts[i].get_text() if i < len(texts) else getattr(h, 'get_label', lambda: '')();return ('legend', str(label), i)
        return None
    def _find_axes_artist_at_event(self, ax, event):
        cands = list(self._iter_axes_color_targets(ax));cands.sort(key=lambda a: float(getattr(a, 'get_zorder', lambda: 0.0)()), reverse=True)
        for art in cands:
            if not getattr(art, 'get_visible', lambda: True)():continue
            if self._artist_contains(art, event):return art
        return None
    def _update_label_override(self, ax_key: int, old: str, new: str) -> None:
        m = self._legend_label_overrides.setdefault(ax_key, {})
        for k in list(m.keys()):
            if m[k] == old:m[k] = new
        m[old] = new;cmap = self._color_overrides.get(ax_key)
        if cmap and old in cmap and (new != old):cmap[new] = cmap.pop(old)
    def _rename_axes_artists_label(self, ax, old: str, new: str) -> None:
        for art in self._iter_axes_color_targets(ax):
            try:
                if hasattr(art, 'get_label') and hasattr(art, 'set_label'):
                    if str(art.get_label()) == str(old):art.set_label(str(new))
            except Exception:pass
    def _apply_user_overrides(self) -> None:
        fig = self.fig;st = getattr(fig, '_suptitle', None)
        if self._suptitle_override is not None:
            if st is None:st = fig.suptitle(self._suptitle_override)
            else:st.set_text(self._suptitle_override)
        if self.graph_title_fontsize is not None and st is not None:
            try:st.set_fontsize(self.graph_title_fontsize)
            except Exception:pass
        for i, ax in enumerate(fig.axes):
            ax_key = int(i)
            if self.axis_label_fontsize is not None:
                try:ax.xaxis.label.set_fontsize(self.axis_label_fontsize)
                except Exception:pass
                try:ax.yaxis.label.set_fontsize(self.axis_label_fontsize)
                except Exception:pass
            if self.tick_label_fontsize is not None:
                try:ax.tick_params(axis='both', labelsize=self.tick_label_fontsize)
                except Exception:pass
            if self.graph_title_fontsize is not None:
                try:ax.title.set_fontsize(self.graph_title_fontsize)
                except Exception:pass
            if ax_key in self._title_overrides:
                try:ax.title.set_text(self._title_overrides[ax_key])
                except Exception:pass
            if ax_key in self._xlabel_overrides:
                try:ax.set_xlabel(self._xlabel_overrides[ax_key])
                except Exception:pass
            if ax_key in self._ylabel_overrides:
                try:ax.set_ylabel(self._ylabel_overrides[ax_key])
                except Exception:pass
            leg = ax.get_legend()
            if leg and ax_key in self._legend_title_overrides:
                try:leg.set_title(self._legend_title_overrides[ax_key])
                except Exception:
                    try:leg.get_title().set_text(self._legend_title_overrides[ax_key])
                    except Exception:pass
            label_map = self._legend_label_overrides.get(ax_key, {})
            if label_map:
                for art in self._iter_axes_color_targets(ax):
                    try:
                        if hasattr(art, 'get_label') and hasattr(art, 'set_label'):
                            lab = str(art.get_label())
                            if lab in label_map:art.set_label(self._resolve_label(lab, label_map))
                    except Exception:pass
                if leg:
                    for txt in leg.get_texts():
                        try:txt.set_text(self._resolve_label(txt.get_text(), label_map))
                        except Exception:pass
            cmap = self._color_overrides.get(ax_key, {})
            if cmap:
                for art in self._iter_axes_color_targets(ax):
                    try:
                        if hasattr(art, 'get_label'):
                            lab = str(art.get_label())
                            if lab in cmap:self._set_artist_color(art, cmap[lab])
                    except Exception:pass
                if leg:
                    texts = list(leg.get_texts());handles = self._legend_handles(leg)
                    for j, txt in enumerate(texts):
                        lab = str(txt.get_text())
                        if lab in cmap and j < len(handles):self._set_artist_color(handles[j], cmap[lab])
        for ax in fig.axes:
            for art in ax.get_children():
                c = self._artist_color_overrides.get(id(art))
                if c:self._set_artist_color(art, c)
    def _on_mpl_button_press(self, event):
        if event is None:return
        if self._toolbar_is_active():return
        if getattr(event, 'dblclick', False) and getattr(event, 'button', None) == 1:
            hit = self._hit_test_editable_text(event)
            if not hit:return
            kind, ax, txt = hit;old_text = str(txt.get_text());new_text = simpledialog.askstring('Edit text', 'Enter new text:', initialvalue=old_text, parent=self.winfo_toplevel())
            if new_text is None:return
            new_text = str(new_text)
            if kind == 'suptitle':self._suptitle_override = new_text;txt.set_text(new_text);self.redraw();return
            if ax is None:return
            ax_key = self._ax_key(ax)
            if kind == 'graph_title':self._title_overrides[ax_key] = new_text;txt.set_text(new_text);self.redraw();return
            if kind == 'xlabel':
                self._xlabel_overrides[ax_key] = new_text
                try:ax.set_xlabel(new_text)
                except Exception:txt.set_text(new_text)
                self.redraw();return
            if kind == 'ylabel':
                self._ylabel_overrides[ax_key] = new_text
                try:ax.set_ylabel(new_text)
                except Exception:txt.set_text(new_text)
                self.redraw();return
            if kind == 'legend_title':self._legend_title_overrides[ax_key] = new_text;txt.set_text(new_text);self.redraw();return
            if kind == 'legend_text':self._update_label_override(ax_key, old_text, new_text);self._rename_axes_artists_label(ax, old_text, new_text);txt.set_text(new_text);self.redraw();return
        if getattr(event, 'button', None) == 3:
            ax = getattr(event, 'inaxes', None)
            if ax is not None:leg_hit = self._find_legend_entry_at_event(ax, event)
            else:leg_hit = None
            if leg_hit and ax is not None:
                _kind, label, _i = leg_hit;ax_key = self._ax_key(ax);current = None
                for art in self._iter_axes_color_targets(ax):
                    try:
                        if hasattr(art, 'get_label') and str(art.get_label()) == str(label):current = self._artist_current_hex(art);break
                    except Exception:pass
                initial = current or self._color_overrides.get(ax_key, {}).get(label, '#000000');s = simpledialog.askstring('Set color', f"Hex color for '{label}' (#RRGGBB or #RRGGBBAA):", initialvalue=initial, parent=self.winfo_toplevel())
                if s is None:return
                try:hex_color = self._normalize_hex_color(s)
                except Exception as e:messagebox.showerror('Invalid color', f'{e}');return
                self._color_overrides.setdefault(ax_key, {})[str(label)] = hex_color;self.redraw();return
            if ax is None:return
            art = self._find_axes_artist_at_event(ax, event)
            if art is None:return
            label = None
            try:
                if hasattr(art, 'get_label'):
                    label = str(art.get_label())
                    if label.startswith('_'):label = None
            except Exception:label = None
            initial = self._artist_current_hex(art) or '#000000';s = simpledialog.askstring('Set color', 'Hex color (#RRGGBB or #RRGGBBAA):', initialvalue=initial, parent=self.winfo_toplevel())
            if s is None:return
            try:hex_color = self._normalize_hex_color(s)
            except Exception as e:messagebox.showerror('Invalid color', f'{e}');return
            if label is not None:self._color_overrides.setdefault(self._ax_key(ax), {})[label] = hex_color
            else:self._artist_color_overrides[id(art)] = hex_color
            self.redraw()
    @staticmethod
    def _axis_looks_like_time(ax) -> bool:
        try:xl = str(ax.get_xlabel() or '').lower()
        except Exception:xl = ''
        if 'time' in xl or 'sec' in xl or '(s' in xl:return True
        try:
            for line in list(getattr(ax, 'lines', [])):
                try:x = np.asarray(line.get_xdata(), dtype=float)
                except Exception:continue
                if x.size < 10:continue
                xf = x[np.isfinite(x)]
                if xf.size < 10:continue
                if float(np.nanmax(xf) - np.nanmin(xf)) < 1.0:continue
                dx = np.diff(xf)
                if dx.size and float(np.nanmean(dx >= 0)) > 0.8:return True
        except Exception:pass
        return False
    def _clear_time_marker_artists(self) -> None:
        for art in list(getattr(self, '_time_marker_artists', [])):
            try:art.remove()
            except Exception:pass
        self._time_marker_artists = []
    def _apply_time_markers(self) -> None:
        self._clear_time_marker_artists();marks = list(getattr(self, '_time_markers', []) or [])
        if not marks:
            for ax in list(getattr(self.fig, 'axes', [])):
                try:
                    leg = ax.get_legend()
                    if leg is not None:leg.remove()
                    handles, labels = ax.get_legend_handles_labels();keep = [(h, lab) for h, lab in zip(handles, labels) if lab and not str(lab).startswith('_')]
                    if keep:ax.legend(*zip(*keep), loc='best')
                except Exception:pass
            return
        styles = [dict(linestyle=':', linewidth=1.3), dict(linestyle=(0, (1, 1)), linewidth=1.3), dict(linestyle=(0, (1, 3)), linewidth=1.3), dict(linestyle=(0, (3, 1, 1, 1)), linewidth=1.3)]
        for ax in list(getattr(self.fig, 'axes', [])):
            if not self._axis_looks_like_time(ax):continue
            added_any = False
            for m in marks:
                try:x = float(m.get('x'))
                except Exception:continue
                if not np.isfinite(x):continue
                label = str(m.get('label') or f't={x:g}s')
                try:si = int(m.get('style_idx', 0))
                except Exception:si = 0
                si = max(0, min(3, si))
                try:art = ax.axvline(x, label=label, color='0.2', alpha=0.85, zorder=9, **styles[si]);self._time_marker_artists.append(art);added_any = True
                except Exception:pass
            if added_any or ax.get_legend() is not None:
                try:
                    if ax.get_legend() is not None:ax.get_legend().remove()
                    ax.legend(loc='best')
                except Exception:pass
    def add_time_marker(self, x_s: float, label: Optional[str]=None) -> None:
        try:x = float(x_s)
        except Exception:return
        if not np.isfinite(x):return
        cur = list(getattr(self, '_time_markers', []) or [])
        if len(cur) >= 4:
            try:messagebox.showwarning('Time markers', 'You can add up to 4 time markers per plot.\n\nRemove one with Ctrl+I, then Ctrl+Backspace (in the dialog).')
            except Exception:pass
            return
        si = len(cur);lab = str(label).strip() if label is not None else ''
        if not lab:lab = f't={x:g}s'
        cur.append({'x': x, 'label': lab, 'style_idx': int(si)});self._time_markers = cur;self.redraw()
    def remove_last_time_marker(self) -> bool:
        cur = list(getattr(self, '_time_markers', []) or [])
        if not cur:return False
        cur.pop(-1);self._time_markers = cur;self.redraw();return True
    def clear_time_markers(self) -> None:self._time_markers = [];self.redraw()
    def redraw(self):self._apply_user_overrides();self._apply_time_markers();self._apply_user_overrides();self.canvas.draw_idle()
class ChannelTabsTk(ttk.Frame):
    def __init__(self, master, res: ChannelResult, norm_mode: str, parent_app=None):super().__init__(master);self.res = res;self.norm_mode = norm_mode;self.parent_app = parent_app;self.tabs = ttk.Notebook(self);self.tabs.pack(fill=tk.BOTH, expand=True);source_base = _basename_no_ext(getattr(res, 'source_path', '')) if getattr(res, 'source_path', '') else '';prefix = f'{source_base}_{res.gcol}'.strip('_');self.tab_raw = PlotTabTk(self.tabs, 'Raw', prefix);self.tab_art = PlotTabTk(self.tabs, 'ArtifactRemover', prefix);self.tab_fit = PlotTabTk(self.tabs, 'Fit', prefix);self.tab_norm = PlotTabTk(self.tabs, 'Normalization', prefix);self.tab_norm_smooth = PlotTabTk(self.tabs, 'Normalization_smoothed', prefix);self.tab_freq = PlotTabTk(self.tabs, 'Freq analysis', prefix, figsize=(8.6, 6.4));self._configure_tabs();self._freq_drawn_version = None;self._fit_initialized = False;self.tabs.bind('<<NotebookTabChanged>>', self._on_inner_tab_changed);self._attach_exporters();self._draw_current_tab()
    def set_axis_label_fontsize(self, fontsize: Optional[float]) -> None:
        for tab in [self.tab_raw, self.tab_art, self.tab_fit, self.tab_norm, self.tab_norm_smooth, self.tab_freq]:
            try:tab.set_axis_label_fontsize(fontsize)
            except Exception:pass
    def set_tick_label_fontsize(self, fontsize: Optional[float]) -> None:
        for tab in [self.tab_raw, self.tab_art, self.tab_fit, self.tab_norm, self.tab_norm_smooth, self.tab_freq]:
            try:tab.set_tick_label_fontsize(fontsize)
            except Exception:pass
    def set_graph_title_fontsize(self, fontsize: Optional[float]) -> None:
        for tab in [self.tab_raw, self.tab_art, self.tab_fit, self.tab_norm, self.tab_norm_smooth, self.tab_freq]:
            try:tab.set_graph_title_fontsize(fontsize)
            except Exception:pass
    def _on_inner_tab_changed(self, _event=None) -> None:self._draw_current_tab()
    def _draw_current_tab(self) -> None:
        sel = self.tabs.select()
        if not sel:return
        if sel == str(self.tab_raw):self._draw_raw();return
        if sel == str(self.tab_art):self._draw_artifact();return
        if sel == str(self.tab_fit):
            if not self._fit_initialized:self._draw_fit_and_attach_selector();self._fit_initialized = True
            else:self.refresh_after_pipeline_change()
            return
        if sel == str(self.tab_norm):self._draw_norm();return
        if sel == str(self.tab_norm_smooth):self._draw_norm_smooth();return
        if self.tab_freq is not None and sel == str(self.tab_freq):
            current_version = int(getattr(self.res, '_data_version', 0))
            if self._freq_drawn_version != current_version:self._draw_frequency();self._freq_drawn_version = current_version
            else:self.tab_freq.redraw()
            return
    def get_active_plot_tab(self) -> Optional[PlotTabTk]:
        sel = self.tabs.select()
        if not sel:return None
        if sel == str(self.tab_raw):return self.tab_raw
        if sel == str(self.tab_art):return self.tab_art
        if sel == str(self.tab_fit):return self.tab_fit
        if sel == str(self.tab_norm):return self.tab_norm
        if sel == str(self.tab_norm_smooth):return self.tab_norm_smooth
        if self.tab_freq is not None and sel == str(self.tab_freq):return self.tab_freq
        return None
    def _norm_tab_texts(self) -> Tuple[str, str]:
        if self.norm_mode == NORM_DFF:return ('ΔF/F', 'ΔF/F smoothed')
        if self.norm_mode == NORM_ZF_GLOBAL:return ('zF (global)', 'zF smoothed')
        if self.norm_mode == NORM_ZF_INTERVAL:return ('zF - interval based', 'zF smoothed')
        return ('Normalization', 'Smoothed')
    def _configure_tabs(self):
        norm_txt, smooth_txt = self._norm_tab_texts()
        if len(self.tabs.tabs()) == 0:
            self.tabs.add(self.tab_raw, text='Raw');self.tabs.add(self.tab_art, text='Artifact remover');self.tabs.add(self.tab_fit, text='Fit');self.tabs.add(self.tab_norm, text=norm_txt);self.tabs.add(self.tab_norm_smooth, text=smooth_txt)
            if self.norm_mode == NORM_DFF:self.tabs.add(self.tab_freq, text='Freq analysis')
            return
        self.tabs.tab(self.tab_norm, text=norm_txt);self.tabs.tab(self.tab_norm_smooth, text=smooth_txt);freq_present = str(self.tab_freq) in self.tabs.tabs()
        if self.norm_mode == NORM_DFF:
            if not freq_present:self.tabs.add(self.tab_freq, text='Freq analysis')
            else:self.tabs.tab(self.tab_freq, text='Freq analysis')
        elif freq_present:
            if self.tabs.select() == str(self.tab_freq):self.tabs.select(self.tab_norm)
            self.tabs.forget(self.tab_freq)
    def set_norm_mode(self, norm_mode: str):self.norm_mode = norm_mode;self._configure_tabs();self._draw_current_tab()
    def _draw_static_tabs(self):self._draw_raw();self._draw_artifact();self._draw_fit_and_attach_selector()
    def refresh_after_pipeline_change(self):
        sel = self.tabs.select()
        if not sel:return
        if sel == str(self.tab_art):self._draw_artifact();return
        if sel == str(self.tab_fit):
            if not self._fit_initialized:self._draw_fit_and_attach_selector();self._fit_initialized = True
            else:
                try:self._line_exc.set_ydata(self.res.exc_clean);self._line_fit.set_ydata(self.res.fitted_iso_on_exc);self._fit_info.set_text(self._fit_info_text());self.tab_fit.redraw()
                except Exception:self._draw_fit_and_attach_selector();self._fit_initialized = True
            return
        if sel == str(self.tab_norm):self._draw_norm();return
        if sel == str(self.tab_norm_smooth):self._draw_norm_smooth();return
        if self.norm_mode == NORM_DFF and self.tab_freq is not None and (sel == str(self.tab_freq)):self._draw_frequency();self._freq_drawn_version = int(getattr(self.res, '_data_version', 0));return
    def _draw_norm_dependent_tabs(self):
        self._draw_norm();self._draw_norm_smooth()
        if self.norm_mode == NORM_DFF:self._draw_frequency()
    def _draw_raw(self):ax = self.tab_raw.ax;ax.clear();ax.plot(self.res.t_exc, self.res.exc_raw, label='Excitatory raw');ax.plot(self.res.t_iso, self.res.iso_raw, label='Isosbestic raw');ax.set_title(f'{self.res.gcol} - Raw');ax.set_xlabel('Time (s)');ax.set_ylabel('Signal');ax.legend(loc='best');self.tab_raw.redraw()
    def _draw_artifact(self):
        ax = self.tab_art.ax;ax.clear();ax.plot(self.res.t_exc, self.res.exc_raw, color='0.75', linewidth=1.0, label='Exc raw');ax.plot(self.res.t_iso, self.res.iso_raw, color='0.75', linewidth=1.0, linestyle='--', label='Iso raw');ax.plot(self.res.t_exc, self.res.exc_clean_holes, linewidth=1.2, label='Exc cleaned (holes)');ax.plot(self.res.t_iso, self.res.iso_clean_holes, linewidth=1.2, label='Iso cleaned (holes)')
        if self.res.use_interpolation:
            def _interp_overlay(y_holes: np.ndarray, y_interp: np.ndarray) -> np.ndarray:
                y_holes = np.asarray(y_holes, dtype=float);y_interp = np.asarray(y_interp, dtype=float);overlay = np.full_like(y_interp, np.nan, dtype=float);nan_mask = ~np.isfinite(y_holes)
                for i0, i1 in _contiguous_true_runs(nan_mask):j0 = max(0, i0 - 1);j1 = min(len(overlay) - 1, i1 + 1);overlay[j0:j1 + 1] = y_interp[j0:j1 + 1]
                return overlay
            exc_fill = _interp_overlay(self.res.exc_clean_holes, self.res.exc_clean_interp);iso_fill = _interp_overlay(self.res.iso_clean_holes, self.res.iso_clean_interp);ax.plot(self.res.t_exc, exc_fill, linewidth=2.0, color='tab:orange', label='Exc interpolated (filled parts)');ax.plot(self.res.t_iso, iso_fill, linewidth=2.0, color='tab:orange', linestyle='--', label='Iso interpolated (filled parts)')
        if np.any(self.res.art_exc):ax.scatter(self.res.t_exc[self.res.art_exc], self.res.exc_raw[self.res.art_exc], s=12, color='red', label='Shared artifacts', zorder=5)
        if np.any(self.res.art_iso):ax.scatter(self.res.t_iso[self.res.art_iso], self.res.iso_raw[self.res.art_iso], s=12, color='red', label='_nolegend_', zorder=5)
        ax.set_title(f'{self.res.gcol} - Artifact remover (red=shared, orange=filled)');ax.set_xlabel('Time (s)');ax.set_ylabel('Signal');ax.legend(loc='best');self.tab_art.redraw()
    def _draw_fit_and_attach_selector(self):
        ax = self.tab_fit.ax;ax.clear();self._fit_window_artists = []
        for a, b in self.res.windows:span = ax.axvspan(a, b, alpha=0.15);v1 = ax.axvline(a, linestyle=':', linewidth=1.0, alpha=0.7);v2 = ax.axvline(b, linestyle=':', linewidth=1.0, alpha=0.7);self._fit_window_artists.extend([span, v1, v2])
        self._line_exc, = ax.plot(self.res.t_exc, self.res.exc_clean, label='Exc (active)');self._line_fit, = ax.plot(self.res.t_exc, self.res.fitted_iso_on_exc, label='Fitted iso → exc (active)');self._fit_info = ax.text(0.99, 0.02, self._fit_info_text(), transform=ax.transAxes, ha='right', va='bottom', fontsize=11, bbox=dict(boxstyle='round,pad=0.25', alpha=0.08, lw=0.6));ax.set_title(f'{self.res.gcol} - Fit (drag to select new window)');ax.set_xlabel('Time (s)');ax.set_ylabel('Signal');ax.legend(loc='best')
        def onselect(xmin, xmax):
            if xmin == xmax:return
            left, right = sorted([float(xmin), float(xmax)]);self.update_fit_window(left, right)
        self._span_selector = SpanSelector(ax, onselect, 'horizontal', useblit=True);self.tab_fit.redraw()
    def _fit_info_text(self) -> str:
        a = self.res.slope;b = self.res.intercept;r2 = self.res.r2;tag = 'interp' if self.res.use_interpolation else 'holes'
        if not (np.isfinite(a) and np.isfinite(b)):return f'{tag}\nNO VALID FIT\n(selected window has <2 usable points)'
        return f'{tag}\ncoefficient={a:.4f}\nintercept={b:.4f}\nR²={r2:.4f}'
    def _get_norm_series(self) -> Tuple[np.ndarray, str, str, str]:
        if self.norm_mode == NORM_DFF:return (self.res.dFF, 'ΔF/F', 'ΔF/F', f'{self.res.gcol} - ΔF/F')
        if self.norm_mode == NORM_ZF_GLOBAL:return (self.res.zF_global, 'zF', 'zF (global)', f'{self.res.gcol} - zF (global, GUI)')
        if self.norm_mode == NORM_ZF_INTERVAL:a = self.res.zf_interval_start_s;b = self.res.zf_interval_end_s;return (self.res.zF_interval, 'zF', f'zF (interval stats: {min(a, b):g}–{max(a, b):g} s)', f'{self.res.gcol} - zF - interval based')
        return (self.res.dFF, 'Signal', 'Signal', f'{self.res.gcol} - Normalization')
    def _get_norm_series_for_result(self, res: ChannelResult) -> Tuple[np.ndarray, str, str, str]:
        if self.norm_mode == NORM_DFF:return (res.dFF, 'ΔF/F', 'ΔF/F', f'{res.gcol} - ΔF/F')
        if self.norm_mode == NORM_ZF_GLOBAL:return (res.zF_global, 'zF', 'zF (global)', f'{res.gcol} - zF (global, GUI)')
        if self.norm_mode == NORM_ZF_INTERVAL:a = res.zf_interval_start_s;b = res.zf_interval_end_s;return (res.zF_interval, 'zF', f'zF (interval stats: {min(a, b):g}–{max(a, b):g} s)', f'{res.gcol} - zF - interval based')
        return (res.dFF, 'Signal', 'Signal', f'{res.gcol} - Normalization')
    def _draw_norm(self):ax = self.tab_norm.ax;ax.clear();y, ylabel, label, title = self._get_norm_series();ax.plot(self.res.t_exc, y, label=label);ax.set_title(title);ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel);ax.legend(loc='best');self.tab_norm.redraw()
    def _draw_norm_smooth(self):
        ax = self.tab_norm_smooth.ax;ax.clear();y, ylabel, label, title = self._get_norm_series();w = int(getattr(self.res, 'smooth_window', DEFAULT_SMOOTH_WINDOW));y_s = get_smoothed_norm_array(self.res, self.norm_mode, window_size=w);line_raw, = ax.plot(self.res.t_exc, y, alpha=0.25, linewidth=1.0, zorder=1, label=f'{label} (raw)');ax.plot(self.res.t_exc, y_s, linewidth=1.4, zorder=3, color=line_raw.get_color(), label=f'{ylabel} smoothed (win={w})');ax.set_title(f'{title} - smoothed');ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel);ax.legend(loc='best')
        if self.parent_app is not None and hasattr(self.parent_app, '_results'):
            try:
                all_results = self.parent_app._results;global_y_min = np.inf;global_y_max = -np.inf
                for mid, res_other in all_results.items():
                    y_other, _, _, _ = self._get_norm_series_for_result(res_other);y_s_other = get_smoothed_norm_array(res_other, self.norm_mode, window_size=w);finite_vals = y_s_other[np.isfinite(y_s_other)]
                    if len(finite_vals) > 0:global_y_min = min(global_y_min, np.nanmin(finite_vals));global_y_max = max(global_y_max, np.nanmax(finite_vals))
                if np.isfinite(global_y_min) and np.isfinite(global_y_max):padding = (global_y_max - global_y_min) * 0.05;ax.set_ylim(global_y_min - padding, global_y_max + padding)
            except Exception:pass
        self.tab_norm_smooth.redraw()
    def _draw_frequency(self):
        fig = self.tab_freq.fig;fig.clear();axs = fig.subplots(3, 2);axes = np.asarray(axs).flatten();t = np.asarray(self.res.t_exc, dtype=float);dff = np.asarray(self.res.dFF_nointerp, dtype=float);fs = float(getattr(self.res, 'eff_fs_hz', np.nan));acq = float(getattr(self.res, 'acq_fps_hz', 0.0))
        if not np.isfinite(fs) or fs <= 0:fs = estimate_fs_from_t(t)
        have_holes = np.any(~np.isfinite(dff))
        if _HAVE_SCIPY_SIGNAL:method = 'Butterworth (segment-wise, no interpolation)'
        else:method = 'FFT (only if no NaNs)' if not have_holes else 'Unavailable (SciPy signal required for NaN holes)'
        supt = f'{self.res.gcol} – Band-limited ΔF/F (NO interpolation)\nUses non-interpolated ΔF/F; holes preserved. Acq FPS={acq:.3g} Hz → eff fs={fs:.3g} Hz | {method}';fig.suptitle(supt, y=0.995, fontsize=11)
        for i, (low_hz, high_hz) in enumerate(FREQ_BANDS):
            ax = axes[i];ax.clear();label = f'{low_hz:g}–{high_hz:g} Hz'
            if not np.any(np.isfinite(dff)):ax.text(0.5, 0.5, 'ΔF/F is all NaN', ha='center', va='center', transform=ax.transAxes);ax.set_title(label);ax.set_xlabel('Time (s)');ax.set_ylabel('ΔF/F');continue
            if _HAVE_SCIPY_SIGNAL:y_band = bandpass_butterworth_segmentwise_no_interp(dff, low_hz, high_hz, fs, order=DEFAULT_BUTTER_ORDER)
            else:y_band = bandpass_fft_no_interp(dff, low_hz, high_hz, fs)
            if not np.any(np.isfinite(y_band)):ax.text(0.5, 0.5, 'Band unavailable\n(check SciPy / NaNs / Nyquist)', ha='center', va='center', transform=ax.transAxes)
            else:ax.plot(t, y_band, label=label, linewidth=1.0);ax.legend(loc='best')
            ax.set_title(label);ax.set_xlabel('Time (s)');ax.set_ylabel('ΔF/F')
        for j in range(len(FREQ_BANDS), len(axes)):axes[j].axis('off')
        fig.tight_layout(rect=(0, 0.0, 1, 0.93));self.tab_freq.redraw()
    def update_fit_window(self, start_s: float, end_s: float) -> None:
        if not isinstance(self.res.windows, list) or not self.res.windows:self.res.windows = [(float(start_s), float(end_s))]
        else:self.res.windows[0] = (float(start_s), float(end_s))
        recompute_fit_and_downstream(self.res)
        if not self._fit_initialized:self._draw_fit_and_attach_selector();self._fit_initialized = True;return
        try:self._line_exc.set_ydata(self.res.exc_clean);self._line_fit.set_ydata(self.res.fitted_iso_on_exc);self._fit_info.set_text(self._fit_info_text());self.tab_fit.redraw()
        except Exception:self._draw_fit_and_attach_selector();self._fit_initialized = True
    def _attach_exporters(self):self.tab_raw.export_provider = self._export_raw;self.tab_art.export_provider = self._export_artifact;self.tab_fit.export_provider = self._export_fit;self.tab_norm.export_provider = self._export_norm;self.tab_norm_smooth.export_provider = self._export_norm_smoothed;self.tab_freq.export_provider = self._export_frequency
    def _meta_df(self) -> pd.DataFrame:
        wtxt = '; '.join([f'[{min(a, b):g},{max(a, b):g}]' for a, b in self.res.windows]) if self.res.windows else ''
        rows = [('fiberlyse_version', FIBERLYSE_VERSION), ('source_path', getattr(self.res, 'source_path', '')), ('gcol', self.res.gcol), ('analysis_exclude_initial_seconds', float(getattr(self.res, 'analysis_exclude_initial_seconds', ANALYSIS_EXCLUDE_INITIAL_SECONDS))), ('artifact_enabled', bool(getattr(self.res, 'artifact_enabled', False))), ('artifact_factor', float(getattr(self.res, 'artifact_factor', np.nan))), ('artifact_pad_extra_samples_per_side', int(getattr(self.res, 'artifact_pad', 0))), ('require_shared_artifacts', bool(getattr(self.res, 'require_shared', False))), ('align_mode', str(getattr(self.res, 'align_mode', DEFAULT_ALIGN_MODE))), ('use_interpolation', bool(self.res.use_interpolation)), ('smooth_window', int(getattr(self.res, 'smooth_window', DEFAULT_SMOOTH_WINDOW))), ('acq_fps_user_hz', float(getattr(self.res, 'acq_fps_hz', np.nan))), ('eff_fs_from_timestamps_hz', float(getattr(self.res, 'eff_fs_hz', np.nan))), ('fit_windows_s', wtxt), ('fit_coefficient_active', float(self.res.slope)), ('intercept_active', float(self.res.intercept)), ('r2_active', float(self.res.r2)), ('zf_interval_start_s', float(getattr(self.res, 'zf_interval_start_s', DEFAULT_ZF_INTERVAL_START_S))), ('zf_interval_end_s', float(getattr(self.res, 'zf_interval_end_s', DEFAULT_ZF_INTERVAL_END_S))), ('norm_mode', str(self.norm_mode))]
        return pd.DataFrame(rows, columns=['key', 'value'])
    def _export_raw(self) -> Dict[str, pd.DataFrame]:exc = pd.DataFrame({'t_s': self.res.t_exc, 'exc_raw': self.res.exc_raw});iso = pd.DataFrame({'t_s': self.res.t_iso, 'iso_raw': self.res.iso_raw});return {'exc_raw': exc, 'iso_raw': iso, 'meta': self._meta_df()}
    def _export_artifact(self) -> Dict[str, pd.DataFrame]:exc = pd.DataFrame({'t_s': self.res.t_exc, 'exc_raw': self.res.exc_raw, 'exc_clean_holes': self.res.exc_clean_holes, 'exc_clean_interp': self.res.exc_clean_interp, 'artifact_exc': self.res.art_exc.astype(bool)});iso = pd.DataFrame({'t_s': self.res.t_iso, 'iso_raw': self.res.iso_raw, 'iso_clean_holes': self.res.iso_clean_holes, 'iso_clean_interp': self.res.iso_clean_interp, 'artifact_iso': self.res.art_iso.astype(bool)});return {'exc_artifact': exc, 'iso_artifact': iso, 'meta': self._meta_df()}
    def _export_fit(self) -> Dict[str, pd.DataFrame]:df = pd.DataFrame({'t_s': self.res.t_exc, 'exc_active': self.res.exc_clean, 'iso_on_exc_active': self.res.iso_on_exc, 'fitted_iso_on_exc_active': self.res.fitted_iso_on_exc, 'residual_active': self.res.residual, 'dF_active': self.res.dF, 'dFF_active': self.res.dFF, 'zF_global': self.res.zF_global, 'zF_interval': self.res.zF_interval});return {'fit_active': df, 'meta': self._meta_df()}
    def _export_norm(self) -> Dict[str, pd.DataFrame]:y, ylabel, label, title = self._get_norm_series();df = pd.DataFrame({'t_s': self.res.t_exc, 'value': y});return {'normalization': df, 'meta': self._meta_df()}
    def _export_norm_smoothed(self) -> Dict[str, pd.DataFrame]:y, ylabel, label, title = self._get_norm_series();w = int(getattr(self.res, 'smooth_window', DEFAULT_SMOOTH_WINDOW));y_s = get_smoothed_norm_array(self.res, self.norm_mode, window_size=w);df = pd.DataFrame({'t_s': self.res.t_exc, 'value_raw': y, f'value_smoothed_win{w}': y_s});return {'normalization_smoothed': df, 'meta': self._meta_df()}
    def _export_frequency(self) -> Dict[str, pd.DataFrame]:
        t = np.asarray(self.res.t_exc, dtype=float);dff = np.asarray(self.res.dFF_nointerp, dtype=float);fs = float(getattr(self.res, 'eff_fs_hz', np.nan));acq = float(getattr(self.res, 'acq_fps_hz', 0.0))
        if not np.isfinite(fs) or fs <= 0:fs = estimate_fs_from_t(t)
        out = {'t_s': t, 'dFF_nointerp': dff}
        for low_hz, high_hz in FREQ_BANDS:
            col = f'band_{low_hz:g}_{high_hz:g}_Hz'
            if _HAVE_SCIPY_SIGNAL:y_band = bandpass_butterworth_segmentwise_no_interp(dff, low_hz, high_hz, fs, order=DEFAULT_BUTTER_ORDER)
            else:y_band = bandpass_fft_no_interp(dff, low_hz, high_hz, fs)
            out[col] = y_band
        df = pd.DataFrame(out);meta = self._meta_df();extra = pd.DataFrame([('freq_acq_fps_hz', acq), ('freq_eff_fs_hz', fs), ('method', 'butter_sos_segmentwise' if _HAVE_SCIPY_SIGNAL else 'fft')], columns=['key', 'value']);meta2 = pd.concat([meta, extra], ignore_index=True);return {'freq_bands': df, 'meta': meta2}
def _fiberlyse_first_finite_value(y: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(y)
    if not finite.any():return np.nan
    idx = int(np.flatnonzero(finite)[0])
    return float(y[idx])
def _fiberlyse_zero_start_percent_dff(y: np.ndarray) -> np.ndarray:
    y_pct = np.asarray(y, dtype=float).copy() * 100.0
    baseline = _fiberlyse_first_finite_value(y_pct)
    if np.isfinite(baseline):y_pct = y_pct - baseline
    return y_pct
def _fiberlyse_batch_trim_minutes_label() -> str:
    try:
        minutes = float(BATCH_HIDE_INITIAL_SECONDS) / 60.0
        return f'{minutes:g} min'
    except Exception:
        return '10 min'
def _fiberlyse_trim_batch_initial_period(t: np.ndarray, y: np.ndarray, start_s: float=BATCH_HIDE_INITIAL_SECONDS) -> Tuple[np.ndarray, np.ndarray]:
    """Compatibility helper: V22 already excludes the initial period upstream."""
    t = np.asarray(t, dtype=float).reshape(-1);y = np.asarray(y, dtype=float).reshape(-1);n = min(t.size, y.size)
    return (t[:n], y[:n])
def get_batch_ylabel_for_norm_mode(mode: str) -> str:
    if mode == NORM_DFF:return '%ΔF/F'
    if mode == NORM_ZF_GLOBAL:return 'zF'
    if mode == NORM_ZF_INTERVAL:return 'zF'
    return 'Signal'
def get_batch_series_for_norm_mode(res: ChannelResult, mode: str, smooth_win: int, do_smooth: bool=False) -> Tuple[np.ndarray, np.ndarray, str, str]:
    """Return a batch-only display/export series.

    V22 receives traces whose first 10 minutes were already excluded upstream,
    before any scientific processing. For these batch views
    only, ΔF/F is converted from fractional units to percent ΔF/F and shifted
    so the first finite displayed sample is 0%. Individual channel tabs keep
    the original full-length ΔF/F and zF traces.
    """
    t = np.asarray(res.t_exc, dtype=float)
    trim_txt = _fiberlyse_batch_trim_minutes_label()
    if mode == NORM_DFF:
        if do_smooth:y = get_smoothed_norm_array(res, mode, window_size=smooth_win)
        else:y = np.asarray(res.dFF, dtype=float)
        t, y = _fiberlyse_trim_batch_initial_period(t, y)
        return (t, _fiberlyse_zero_start_percent_dff(y), '%ΔF/F', f'%ΔF/F (first {trim_txt} excluded upstream; zeroed to first finite sample)')
    if mode == NORM_ZF_GLOBAL:
        if do_smooth:y = get_smoothed_norm_array(res, mode, window_size=smooth_win)
        else:y = np.asarray(res.zF_global, dtype=float)
        t, y = _fiberlyse_trim_batch_initial_period(t, y)
        return (t, y, 'zF', f'zF (global; first {trim_txt} excluded upstream)')
    if mode == NORM_ZF_INTERVAL:
        if do_smooth:y = get_smoothed_norm_array(res, mode, window_size=smooth_win)
        else:y = np.asarray(res.zF_interval, dtype=float)
        t, y = _fiberlyse_trim_batch_initial_period(t, y)
        return (t, y, 'zF', f'zF (interval based; first {trim_txt} excluded upstream)')
    if do_smooth:y = get_smoothed_norm_array(res, NORM_DFF, window_size=smooth_win)
    else:y = np.asarray(res.dFF, dtype=float)
    t, y = _fiberlyse_trim_batch_initial_period(t, y)
    return (t, _fiberlyse_zero_start_percent_dff(y), '%ΔF/F', f'%ΔF/F (first {trim_txt} excluded upstream; zeroed to first finite sample)')
class BatchCompareTk(ttk.Frame):
    def __init__(self, master, app: 'MainAppTk'):super().__init__(master);self.app = app;self.selected_ids: List[str] = [];self._available_ids: List[str] = [];self._display: Dict[str, str] = {};self.columnconfigure(1, weight=1);self.rowconfigure(0, weight=1);controls = ttk.Frame(self, padding=8);controls.grid(row=0, column=0, sticky='nsw');plot_holder = ttk.Frame(self);plot_holder.grid(row=0, column=1, sticky='nsew', padx=(0, 8), pady=(8, 8));plot_holder.rowconfigure(0, weight=1);plot_holder.columnconfigure(0, weight=1);ttk.Label(controls, text='Available recordings/channels:').grid(row=0, column=0, sticky='w');self.lst_available = tk.Listbox(controls, selectmode=tk.EXTENDED, height=10, width=28);sb_av = ttk.Scrollbar(controls, orient='vertical', command=self.lst_available.yview);self.lst_available.config(yscrollcommand=sb_av.set);self.lst_available.grid(row=1, column=0, sticky='nsew');sb_av.grid(row=1, column=1, sticky='ns');btns = ttk.Frame(controls);btns.grid(row=2, column=0, columnspan=2, pady=(6, 10), sticky='ew');self.btn_add = ttk.Button(btns, text='Add →', command=self.add_selected);self.btn_add.pack(side=tk.LEFT, fill=tk.X, expand=True);self.btn_remove = ttk.Button(btns, text='← Remove', command=self.remove_selected);self.btn_remove.pack(side=tk.LEFT, padx=(6, 0), fill=tk.X, expand=True);ttk.Label(controls, text='Selected for compare:').grid(row=3, column=0, sticky='w');self.lst_selected = tk.Listbox(controls, selectmode=tk.EXTENDED, height=10, width=28);sb_sel = ttk.Scrollbar(controls, orient='vertical', command=self.lst_selected.yview);self.lst_selected.config(yscrollcommand=sb_sel.set);self.lst_selected.grid(row=4, column=0, sticky='nsew');sb_sel.grid(row=4, column=1, sticky='ns');opts = ttk.Frame(controls);opts.grid(row=5, column=0, columnspan=2, sticky='ew', pady=(10, 0));self.var_smoothed = tk.BooleanVar(value=False);ttk.Checkbutton(opts, text='Plot smoothed trace', variable=self.var_smoothed, command=self.refresh_plot).pack(side=tk.TOP, anchor='w');self.btn_plot = ttk.Button(opts, text='Plot selected', command=self.refresh_plot);self.btn_plot.pack(side=tk.TOP, fill=tk.X, pady=(8, 0));self.btn_clear = ttk.Button(opts, text='Clear selection', command=self.clear_selection);self.btn_clear.pack(side=tk.TOP, fill=tk.X, pady=(6, 0));hint = ttk.Label(controls, text='Tip: Ctrl+I adds a line or interval event annotation\n(remove last with Ctrl+I then Ctrl+Backspace).', foreground='gray35', justify='left');hint.grid(row=6, column=0, columnspan=2, sticky='w', pady=(12, 0));self.plot = PlotTabTk(plot_holder, 'Compare', default_filename_prefix='batch_compare', figsize=(8.2, 5.2));self.plot.grid(row=0, column=0, sticky='nsew');self.plot.export_provider = self._export_compare;self.lst_available.bind('<Double-Button-1>', lambda _e: self.add_selected());self.lst_selected.bind('<Double-Button-1>', lambda _e: self.remove_selected());self.refresh_available();self.refresh_plot()
    def refresh_available(self) -> None:
        items = self.app.get_mouse_items();self._display = {mid: label for mid, label in items};self.selected_ids = [mid for mid in self.selected_ids if mid in self._display];self._available_ids = [mid for mid, _lab in items if mid not in set(self.selected_ids)];self.lst_available.delete(0, tk.END)
        for mid in self._available_ids:self.lst_available.insert(tk.END, self._display.get(mid, mid))
        self.lst_selected.delete(0, tk.END)
        for mid in self.selected_ids:self.lst_selected.insert(tk.END, self._display.get(mid, mid))
    def add_selected(self) -> None:
        sel = list(self.lst_available.curselection())
        if not sel:return
        for i in sel:
            if 0 <= int(i) < len(self._available_ids):
                mid = self._available_ids[int(i)]
                if mid not in self.selected_ids:self.selected_ids.append(mid)
        self.refresh_available();self.refresh_plot()
    def remove_selected(self) -> None:
        sel = sorted(list(self.lst_selected.curselection()), reverse=True)
        if not sel:return
        for i in sel:
            if 0 <= int(i) < len(self.selected_ids):self.selected_ids.pop(int(i))
        self.refresh_available();self.refresh_plot()
    def clear_selection(self) -> None:self.selected_ids = [];self.refresh_available();self.refresh_plot()
    def set_axis_label_fontsize(self, fontsize: Optional[float]) -> None:
        try:self.plot.set_axis_label_fontsize(fontsize)
        except Exception:pass
    def set_graph_title_fontsize(self, fontsize: Optional[float]) -> None:
        try:self.plot.set_graph_title_fontsize(fontsize)
        except Exception:pass
    def set_tick_label_fontsize(self, fontsize: Optional[float]) -> None:
        try:self.plot.set_tick_label_fontsize(fontsize)
        except Exception:pass
    def refresh_plot(self) -> None:
        ax = self.plot.ax;ax.clear();mode = self.app.get_norm_mode();smooth_win = self.app.get_smooth_window();do_smooth = bool(self.var_smoothed.get());ylabel = get_batch_ylabel_for_norm_mode(mode);plotted = 0;skipped = 0
        for mid in self.selected_ids:
            res = self.app.get_mouse_result(mid)
            if res is None:continue
            t, y, ylabel, _series_label = get_batch_series_for_norm_mode(res, mode, smooth_win=smooth_win, do_smooth=do_smooth)
            if t.size == 0 or y.size == 0 or not np.any(np.isfinite(y)):
                skipped += 1;continue
            ax.plot(t, y, label=self._display.get(mid, mid));plotted += 1
        trim_txt = _fiberlyse_batch_trim_minutes_label();title_extra = f' – first {trim_txt} excluded upstream'
        if mode == NORM_DFF:title_extra += '; traces zeroed to 0%'
        ax.set_title(f"Batch compare ({plotted} trace{('s' if plotted != 1 else '')}){title_extra}");ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel)
        if plotted:ax.legend(loc='best')
        elif skipped:ax.text(0.5, 0.5, f'No finite data after upstream exclusion of first {trim_txt}.', ha='center', va='center', transform=ax.transAxes)
        self.plot.redraw()
    def _export_compare(self) -> Dict[str, pd.DataFrame]:
        mode = self.app.get_norm_mode();smooth_win = self.app.get_smooth_window();do_smooth = bool(self.var_smoothed.get());payload: Dict[str, pd.DataFrame] = {}
        for mid in self.selected_ids:
            res = self.app.get_mouse_result(mid)
            if res is None:continue
            t, y, _ylabel, series_label = get_batch_series_for_norm_mode(res, mode, smooth_win=smooth_win, do_smooth=do_smooth)
            payload[self._display.get(mid, mid)] = pd.DataFrame({'t_s': t, 'value': y, 'series': series_label})
        payload['meta'] = pd.DataFrame([('mode', mode), ('analysis_exclude_initial_seconds', float(ANALYSIS_EXCLUDE_INITIAL_SECONDS)), ('analysis_exclusion_note', 'Initial non-biological samples are removed upstream before artifact detection, fitting, normalization, filtering, AUC, and batch analysis.'), ('batch_dff_transform', 'ΔF/F only: percent units, first finite displayed sample subtracted so trace starts at 0%'), ('smoothed', do_smooth), ('smooth_window', smooth_win), ('n_selected_traces', len(self.selected_ids))], columns=['key', 'value']);return payload
class BatchAverageTk(ttk.Frame):
    def __init__(self, master, app: 'MainAppTk'):
        super().__init__(master)
        self.app = app
        self.group_a: List[str] = []
        self.group_b: List[str] = []
        self._available_ids: List[str] = []
        self._display: Dict[str, str] = {}
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        # The group-selection controls can be taller than a 13-inch laptop
        # viewport. Keep the plot fixed on the right, and make only the
        # left control panel vertically scrollable so Group B remains reachable.
        controls_outer = ttk.Frame(self, padding=0)
        controls_outer.grid(row=0, column=0, sticky='nsw')
        controls_outer.rowconfigure(0, weight=1)
        controls_outer.columnconfigure(0, weight=1)
        controls_canvas = tk.Canvas(controls_outer, borderwidth=0, highlightthickness=0, width=320)
        controls_scroll = ttk.Scrollbar(controls_outer, orient='vertical', command=controls_canvas.yview)
        controls_canvas.configure(yscrollcommand=controls_scroll.set)
        controls_canvas.grid(row=0, column=0, sticky='nsew')
        controls_scroll.grid(row=0, column=1, sticky='ns')
        controls = ttk.Frame(controls_canvas, padding=8)
        controls_window = controls_canvas.create_window((0, 0), window=controls, anchor='nw')
        self._batch_avg_controls_outer = controls_outer
        self._batch_avg_controls_canvas = controls_canvas
        self._batch_avg_controls_scrollbar = controls_scroll
        self._batch_avg_controls_inner = controls

        def _update_scrollregion(_event=None):
            try:
                controls_canvas.configure(scrollregion=controls_canvas.bbox('all'))
            except Exception:
                pass

        def _sync_inner_width(event=None):
            try:
                width = int(getattr(event, 'width', 0) or controls_canvas.winfo_width())
                controls_canvas.itemconfigure(controls_window, width=max(1, width))
            except Exception:
                pass

        def _scroll_controls(event):
            try:
                bbox = controls_canvas.bbox('all')
                if not bbox:
                    return None
                content_h = int(bbox[3] - bbox[1])
                canvas_h = int(controls_canvas.winfo_height())
                if content_h <= canvas_h:
                    return None
                num = getattr(event, 'num', None)
                if num == 4:
                    units = -3
                elif num == 5:
                    units = 3
                else:
                    delta = int(getattr(event, 'delta', 0) or 0)
                    if delta == 0:
                        return None
                    units = -3 if delta > 0 else 3
                controls_canvas.yview_scroll(units, 'units')
                return 'break'
            except Exception:
                return None

        def _bind_panel_wheel(widget):
            try:
                widget.bind('<MouseWheel>', _scroll_controls)
                widget.bind('<Button-4>', _scroll_controls)
                widget.bind('<Button-5>', _scroll_controls)
            except Exception:
                pass
            try:
                for child in widget.winfo_children():
                    _bind_panel_wheel(child)
            except Exception:
                pass

        controls.bind('<Configure>', _update_scrollregion, add='+')
        controls_canvas.bind('<Configure>', _sync_inner_width, add='+')

        plot_holder = ttk.Frame(self)
        plot_holder.grid(row=0, column=1, sticky='nsew', padx=(0, 8), pady=(8, 8))
        plot_holder.rowconfigure(0, weight=1)
        plot_holder.columnconfigure(0, weight=1)

        ttk.Label(controls, text='Available recordings/channels:').grid(row=0, column=0, sticky='w')
        self.lst_available = tk.Listbox(controls, selectmode=tk.EXTENDED, height=9, width=28)
        sb_av = ttk.Scrollbar(controls, orient='vertical', command=self.lst_available.yview)
        self.lst_available.config(yscrollcommand=sb_av.set)
        self.lst_available.grid(row=1, column=0, sticky='nsew')
        sb_av.grid(row=1, column=1, sticky='ns')

        btns = ttk.Frame(controls)
        btns.grid(row=2, column=0, columnspan=2, pady=(6, 10), sticky='ew')
        ttk.Button(btns, text='Add to A \u2192', command=self.add_to_a).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(btns, text='Add to B \u2192', command=self.add_to_b).pack(side=tk.LEFT, padx=(6, 0), fill=tk.X, expand=True)

        ttk.Label(controls, text='Group A:').grid(row=3, column=0, sticky='w')
        self.lst_a = tk.Listbox(controls, selectmode=tk.EXTENDED, height=6, width=28)
        sb_a = ttk.Scrollbar(controls, orient='vertical', command=self.lst_a.yview)
        self.lst_a.config(yscrollcommand=sb_a.set)
        self.lst_a.grid(row=4, column=0, sticky='nsew')
        sb_a.grid(row=4, column=1, sticky='ns')
        ttk.Button(controls, text='Remove from A', command=self.remove_from_a).grid(row=5, column=0, columnspan=2, sticky='ew', pady=(4, 10))

        ttk.Label(controls, text='Group B:').grid(row=6, column=0, sticky='w')
        self.lst_b = tk.Listbox(controls, selectmode=tk.EXTENDED, height=6, width=28)
        sb_b = ttk.Scrollbar(controls, orient='vertical', command=self.lst_b.yview)
        self.lst_b.config(yscrollcommand=sb_b.set)
        self.lst_b.grid(row=7, column=0, sticky='nsew')
        sb_b.grid(row=7, column=1, sticky='ns')
        ttk.Button(controls, text='Remove from B', command=self.remove_from_b).grid(row=8, column=0, columnspan=2, sticky='ew', pady=(4, 10))

        opts = ttk.Frame(controls)
        opts.grid(row=9, column=0, columnspan=2, sticky='ew')
        self.var_show_individual = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts, text='Show individual traces', variable=self.var_show_individual, command=self.refresh_plot).pack(side=tk.TOP, anchor='w')
        self.var_show_sem = tk.BooleanVar(value=True)
        ttk.Checkbutton(opts, text='Show SEM shading', variable=self.var_show_sem, command=self.refresh_plot).pack(side=tk.TOP, anchor='w')
        self.var_smoothed = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts, text='Average smoothed traces', variable=self.var_smoothed, command=self.refresh_plot).pack(side=tk.TOP, anchor='w')
        self.btn_plot = ttk.Button(opts, text='Plot averages', command=self.refresh_plot)
        self.btn_plot.pack(side=tk.TOP, fill=tk.X, pady=(8, 0))
        self.btn_clear = ttk.Button(opts, text='Clear groups', command=self.clear_groups)
        self.btn_clear.pack(side=tk.TOP, fill=tk.X, pady=(6, 0))
        hint = ttk.Label(controls, text='Tip: Ctrl+I adds a line or interval event annotation.', foreground='gray35', justify='left')
        hint.grid(row=10, column=0, columnspan=2, sticky='w', pady=(12, 0))

        self.plot = PlotTabTk(plot_holder, 'Average', default_filename_prefix='batch_average', figsize=(8.2, 5.2))
        self.plot.grid(row=0, column=0, sticky='nsew')
        self.plot.export_provider = self._export_average

        _bind_panel_wheel(controls_outer)
        _update_scrollregion()
        self.refresh_available()
        self.refresh_plot()
    def set_axis_label_fontsize(self, fontsize: Optional[float]) -> None:
        try:self.plot.set_axis_label_fontsize(fontsize)
        except Exception:pass
    def set_graph_title_fontsize(self, fontsize: Optional[float]) -> None:
        try:self.plot.set_graph_title_fontsize(fontsize)
        except Exception:pass
    def set_tick_label_fontsize(self, fontsize: Optional[float]) -> None:
        try:self.plot.set_tick_label_fontsize(fontsize)
        except Exception:pass
    def refresh_available(self) -> None:
        items = self.app.get_mouse_items();self._display = {mid: label for mid, label in items};all_ids = set(self._display.keys());self.group_a = [mid for mid in self.group_a if mid in all_ids];self.group_b = [mid for mid in self.group_b if mid in all_ids];used = set(self.group_a) | set(self.group_b);self._available_ids = [mid for mid, _lab in items if mid not in used];self.lst_available.delete(0, tk.END)
        for mid in self._available_ids:self.lst_available.insert(tk.END, self._display.get(mid, mid))
        self.lst_a.delete(0, tk.END)
        for mid in self.group_a:self.lst_a.insert(tk.END, self._display.get(mid, mid))
        self.lst_b.delete(0, tk.END)
        for mid in self.group_b:self.lst_b.insert(tk.END, self._display.get(mid, mid))
    def _add_from_available(self, target: List[str]) -> None:
        sel = list(self.lst_available.curselection())
        if not sel:return
        for i in sel:
            if 0 <= int(i) < len(self._available_ids):
                mid = self._available_ids[int(i)]
                if mid not in target:target.append(mid)
        self.refresh_available();self.refresh_plot()
    def add_to_a(self) -> None:self._add_from_available(self.group_a)
    def add_to_b(self) -> None:self._add_from_available(self.group_b)
    def remove_from_a(self) -> None:
        sel = sorted(list(self.lst_a.curselection()), reverse=True)
        if not sel:return
        for i in sel:
            if 0 <= int(i) < len(self.group_a):self.group_a.pop(int(i))
        self.refresh_available();self.refresh_plot()
    def remove_from_b(self) -> None:
        sel = sorted(list(self.lst_b.curselection()), reverse=True)
        if not sel:return
        for i in sel:
            if 0 <= int(i) < len(self.group_b):self.group_b.pop(int(i))
        self.refresh_available();self.refresh_plot()
    def clear_groups(self) -> None:self.group_a = [];self.group_b = [];self.refresh_available();self.refresh_plot()
    @staticmethod
    def _interp_to_common_grid(t_list: List[np.ndarray], y_list: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
        """Align recordings only across their shared time support.

        Interpolation is used only to place finite contiguous segments onto the
        common grid. NaN/artifact holes are never bridged and values are never
        extrapolated beyond a recording's first/last timestamp.
        """
        pairs: List[Tuple[np.ndarray, np.ndarray]] = []
        for t, y in zip(t_list, y_list):
            t = np.asarray(t, dtype=float);y = np.asarray(y, dtype=float)
            if t.size < 2 or y.size != t.size:continue
            tf = t[np.isfinite(t)]
            if tf.size < 2 or np.any(np.diff(tf) <= 0):continue
            pairs.append((t, y))
        if not pairs:return (np.array([], dtype=float), np.empty((0, 0), dtype=float))
        overlap_start = max(float(np.nanmin(t)) for t, _ in pairs);overlap_end = min(float(np.nanmax(t)) for t, _ in pairs)
        if not np.isfinite(overlap_start) or not np.isfinite(overlap_end) or overlap_end <= overlap_start:return (np.array([], dtype=float), np.empty((0, 0), dtype=float))
        candidates = []
        for t, _y in pairs:
            mask = np.isfinite(t) & (t >= overlap_start) & (t <= overlap_end);tg = t[mask]
            if tg.size >= 2:candidates.append(tg)
        if not candidates:return (np.array([], dtype=float), np.empty((0, 0), dtype=float))
        # Use the sparsest real timestamp grid to avoid inventing unnecessary samples.
        t0 = np.asarray(min(candidates, key=len), dtype=float)
        Y = []
        for t, y in pairs:
            yi = np.full_like(t0, np.nan, dtype=float);finite = np.isfinite(t) & np.isfinite(y)
            for i0, i1 in _contiguous_true_runs(finite):
                ts = t[i0:i1 + 1];ys = y[i0:i1 + 1]
                if ts.size < 2 or np.any(np.diff(ts) <= 0):continue
                on_seg = (t0 >= ts[0]) & (t0 <= ts[-1])
                if np.any(on_seg):yi[on_seg] = np.interp(t0[on_seg], ts, ys)
            if np.any(np.isfinite(yi)):Y.append(yi)
        if not Y:return (t0, np.empty((0, t0.size), dtype=float))
        return (t0, np.vstack(Y))
    def _group_stats(self, ids: List[str], mode: str, smooth_win: int, do_smooth: bool) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        t_list: List[np.ndarray] = [];y_list: List[np.ndarray] = []
        for mid in ids:
            res = self.app.get_mouse_result(mid)
            if res is None:continue
            t, y, _ylabel, _series_label = get_batch_series_for_norm_mode(res, mode, smooth_win=smooth_win, do_smooth=do_smooth)
            if t.size < 2 or y.size != t.size:continue
            t_list.append(t);y_list.append(y)
        t0, Y = self._interp_to_common_grid(t_list, y_list)
        if mode == NORM_DFF and Y.size:
            rows = []
            for row in Y:
                baseline = _fiberlyse_first_finite_value(row)
                rows.append(row - baseline if np.isfinite(baseline) else row)
            Y = np.vstack(rows)
        if Y.size == 0:return (t0, np.full_like(t0, np.nan, dtype=float), np.full_like(t0, np.nan, dtype=float), 0)
        mean = np.nanmean(Y, axis=0);n_eff = np.sum(np.isfinite(Y), axis=0);std = np.nanstd(Y, axis=0, ddof=1)
        with np.errstate(divide='ignore', invalid='ignore'):sem = std / np.sqrt(np.where(n_eff > 0, n_eff, np.nan))
        return (t0, mean, sem, int(Y.shape[0]))
    def refresh_plot(self) -> None:
        ax = self.plot.ax
        ax.clear()
        mode = self.app.get_norm_mode()
        smooth_win = self.app.get_smooth_window()
        do_smooth = bool(self.var_smoothed.get())
        ylabel = get_batch_ylabel_for_norm_mode(mode)
        show_individual = bool(self.var_show_individual.get())
        show_sem = bool(self.var_show_sem.get())

        tA, mA, sA, nA = self._group_stats(self.group_a, mode, smooth_win, do_smooth)
        tB, mB, sB, nB = self._group_stats(self.group_b, mode, smooth_win, do_smooth)

        mean_label = 'smoothed mean' if do_smooth else 'mean'
        individual_label = 'smoothed individual traces' if do_smooth else 'individual traces'
        sem_alpha = 0.28
        lineA = None
        lineB = None

        if nA > 0 and tA.size:
            lineA, = ax.plot(tA, mA, linewidth=2.0, zorder=3, label=f'Group A {mean_label} (traces={nA})')
            if show_sem and np.any(np.isfinite(sA)):
                c = lineA.get_color()
                ax.fill_between(
                    tA,
                    mA - sA,
                    mA + sA,
                    facecolor=c,
                    edgecolor='none',
                    alpha=sem_alpha,
                    zorder=1,
                    label='Group A SEM',
                )

        if nB > 0 and tB.size:
            lineB, = ax.plot(tB, mB, linewidth=2.0, zorder=3, label=f'Group B {mean_label} (traces={nB})')
            if show_sem and np.any(np.isfinite(sB)):
                c = lineB.get_color()
                ax.fill_between(
                    tB,
                    mB - sB,
                    mB + sB,
                    facecolor=c,
                    edgecolor='none',
                    alpha=sem_alpha,
                    zorder=1,
                    label='Group B SEM',
                )

        if show_individual:
            labeled_a = False
            color_a = lineA.get_color() if lineA is not None else None
            for mid in self.group_a:
                res = self.app.get_mouse_result(mid)
                if res is None:
                    continue
                t, y, ylabel, _series_label = get_batch_series_for_norm_mode(res, mode, smooth_win=smooth_win, do_smooth=do_smooth)
                if t.size and y.size and np.any(np.isfinite(y)):
                    label = f'Group A {individual_label}' if not labeled_a else '_nolegend_'
                    kwargs = dict(alpha=0.25, linewidth=1.0, zorder=2, label=label)
                    if color_a is not None:
                        kwargs['color'] = color_a
                    ax.plot(t, y, **kwargs)
                    labeled_a = True

            labeled_b = False
            color_b = lineB.get_color() if lineB is not None else None
            for mid in self.group_b:
                res = self.app.get_mouse_result(mid)
                if res is None:
                    continue
                t, y, ylabel, _series_label = get_batch_series_for_norm_mode(res, mode, smooth_win=smooth_win, do_smooth=do_smooth)
                if t.size and y.size and np.any(np.isfinite(y)):
                    label = f'Group B {individual_label}' if not labeled_b else '_nolegend_'
                    kwargs = dict(alpha=0.25, linewidth=1.0, zorder=2, label=label)
                    if color_b is not None:
                        kwargs['color'] = color_b
                    ax.plot(t, y, **kwargs)
                    labeled_b = True

        trim_txt = _fiberlyse_batch_trim_minutes_label()
        title_extra = f' \u2013 first {trim_txt} excluded upstream'
        if mode == NORM_DFF:
            title_extra += '; traces zeroed to 0%'
        ax.set_title(f'Batch average (Group A vs Group B){title_extra}')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(ylabel)

        if nA > 0 or nB > 0:
            handles, labels = ax.get_legend_handles_labels()
            filtered = [(h, lab) for h, lab in zip(handles, labels) if lab and not str(lab).startswith('_')]
            if filtered:
                handles, labels = zip(*filtered)
                ax.legend(handles, labels, loc='best')
        elif self.group_a or self.group_b:
            ax.text(0.5, 0.5, f'No finite data after upstream exclusion of first {trim_txt}.', ha='center', va='center', transform=ax.transAxes)

        self.plot.redraw()
    def _export_average(self) -> Dict[str, pd.DataFrame]:mode = self.app.get_norm_mode();smooth_win = self.app.get_smooth_window();do_smooth = bool(self.var_smoothed.get());tA, mA, sA, nA = self._group_stats(self.group_a, mode, smooth_win, do_smooth);tB, mB, sB, nB = self._group_stats(self.group_b, mode, smooth_win, do_smooth);payload: Dict[str, pd.DataFrame] = {};payload['groupA_mean_sem'] = pd.DataFrame({'t_s': tA, 'mean': mA, 'sem': sA});payload['groupB_mean_sem'] = pd.DataFrame({'t_s': tB, 'mean': mB, 'sem': sB});payload['groupA_members'] = pd.DataFrame({'recording_channel': [self._display.get(mid, mid) for mid in self.group_a]});payload['groupB_members'] = pd.DataFrame({'recording_channel': [self._display.get(mid, mid) for mid in self.group_b]});payload['meta'] = pd.DataFrame([('mode', mode), ('analysis_exclude_initial_seconds', float(ANALYSIS_EXCLUDE_INITIAL_SECONDS)), ('analysis_exclusion_note', 'Initial non-biological samples are removed upstream before artifact detection, fitting, normalization, filtering, AUC, and batch analysis.'), ('batch_dff_transform', 'ΔF/F only: percent units, first finite displayed sample subtracted so each trace starts at 0% before averaging'), ('avg_smoothed', do_smooth), ('smooth_window', smooth_win), ('n_groupA_traces', nA), ('n_groupB_traces', nB)], columns=['key', 'value']);return payload
def _basename_no_ext(path: str) -> str:base = os.path.basename(str(path));name, _ext = os.path.splitext(base);name = name.strip();return name if name else base if base else 'file'
def _unique_aliases(paths: List[str]) -> List[str]:
    used: set = set();aliases: List[str] = []
    for p in paths:
        base = _basename_no_ext(p);alias = base;k = 2
        while alias in used:alias = f'{base}_{k}';k += 1
        used.add(alias);aliases.append(alias)
    return aliases
def get_series_for_norm_mode(res: ChannelResult, mode: str) -> Tuple[np.ndarray, np.ndarray, str]:
    if mode == NORM_DFF:return (np.asarray(res.t_exc, dtype=float), np.asarray(res.dFF, dtype=float), 'ΔF/F')
    if mode == NORM_ZF_GLOBAL:return (np.asarray(res.t_exc, dtype=float), np.asarray(res.zF_global, dtype=float), 'zF')
    if mode == NORM_ZF_INTERVAL:return (np.asarray(res.t_exc, dtype=float), np.asarray(res.zF_interval, dtype=float), 'zF')
    return (np.asarray(res.t_exc, dtype=float), np.asarray(res.dFF, dtype=float), 'Signal')
class MainAppTk:
    def __init__(self, initial_csvs: Optional[List[str]]=None, autorun: bool=False):
        self.root = tk.Tk();self.root.title('Fiberlyse');self.root.resizable(True, True);self.root.geometry('1200x700');self.csv_paths: List[str] = []
        if initial_csvs:
            for p in list(initial_csvs):
                if p and p not in self.csv_paths:self.csv_paths.append(p)
        self._results: Optional[Dict[str, ChannelResult]] = None;self._mouse_display: Dict[str, str] = {};self._mouse_order: List[str] = [];self._analysis_thread: Optional[threading.Thread] = None;self._channel_widgets: Dict[str, ChannelTabsTk] = {};self.compare_widget: Optional[BatchCompareTk] = None;self.average_widget: Optional[BatchAverageTk] = None;self._result_file_order: List[str] = [];self._file_alias_by_key: Dict[str, str] = {};self._file_path_by_key: Dict[str, str] = {};self._axis_label_fs_override: Optional[float] = None;self._graph_title_fs_override: Optional[float] = None;self._tick_label_fs_override: Optional[float] = None;top = ttk.Frame(self.root, padding=8);top.pack(side=tk.TOP, fill=tk.X);self.lbl_file = ttk.Label(top, text='');self.lbl_file.grid(row=0, column=0, sticky='w', padx=(0, 8));self.btn_add = ttk.Button(top, text='Add CSV(s)…', command=self.add_csvs);self.btn_add.grid(row=0, column=1, sticky='e', padx=(0, 6));self.btn_clear = ttk.Button(top, text='Clear files', command=self.clear_csvs);self.btn_clear.grid(row=0, column=2, sticky='e', padx=(0, 6));self.btn_run = ttk.Button(top, text='Run analysis', command=self.run_analysis);self.btn_run.grid(row=0, column=3, sticky='e');self._update_files_label();row2 = ttk.Frame(self.root, padding=(8, 0, 8, 4));row2.pack(side=tk.TOP, fill=tk.X);self.var_artifact_enabled = tk.BooleanVar(value=True);self.chk_artifact_enabled = ttk.Checkbutton(row2, text='Enable artifact remover (MAD)', variable=self.var_artifact_enabled, command=self.on_artifact_enabled_toggled);self.chk_artifact_enabled.grid(row=0, column=0, columnspan=2, sticky='w', padx=(0, 12));ttk.Label(row2, text='Factor:').grid(row=0, column=2, sticky='w');self.var_factor = tk.StringVar(value='11.9');self.spin_factor = ttk.Spinbox(row2, from_=0.1, to=1000.0, increment=0.1, textvariable=self.var_factor, width=8);self.spin_factor.grid(row=0, column=3, padx=(6, 12), sticky='w');ttk.Label(row2, text='Extra pad (samples/side):').grid(row=0, column=4, sticky='w');self.var_pad = tk.StringVar(value='1');self.spin_pad = ttk.Spinbox(row2, from_=0, to=50, increment=1, textvariable=self.var_pad, width=6);self.spin_pad.grid(row=0, column=5, padx=(6, 12), sticky='w');self.var_shared = tk.BooleanVar(value=True);self.chk_shared = ttk.Checkbutton(row2, text='Require shared artifacts', variable=self.var_shared);self.chk_shared.grid(row=0, column=6, padx=(0, 12), sticky='w');ttk.Label(row2, text='Acq FPS (Hz):').grid(row=0, column=7, sticky='w');self.var_acq_fps = tk.StringVar(value=str(DEFAULT_ACQ_FPS_HZ));self.spin_acq_fps = ttk.Spinbox(row2, from_=0.0, to=10000.0, increment=0.1, textvariable=self.var_acq_fps, width=8);self.spin_acq_fps.grid(row=0, column=8, padx=(6, 6), sticky='w');ttk.Label(row2, text='(filter fs measured from timestamps)').grid(row=0, column=9, sticky='w');ttk.Label(row2, text='Smooth win (samples):').grid(row=0, column=10, sticky='w');self.var_smooth_win = tk.StringVar(value=str(DEFAULT_SMOOTH_WINDOW));self.spin_smooth_win = ttk.Spinbox(row2, from_=1, to=100000, increment=1, textvariable=self.var_smooth_win, width=8);self.spin_smooth_win.grid(row=0, column=11, padx=(6, 12), sticky='w');self.var_interp = tk.BooleanVar(value=DEFAULT_USE_LINEAR_INTERP);self.chk_interp = ttk.Checkbutton(row2, text='Linear interpolate holes (after artifact removal)', variable=self.var_interp, command=self.on_interp_toggled);self.chk_interp.grid(row=1, column=0, columnspan=6, sticky='w', pady=(4, 0));row3 = ttk.Frame(self.root, padding=(8, 0, 8, 8));row3.pack(side=tk.TOP, fill=tk.X);ttk.Label(row3, text='Normalization view:').grid(row=0, column=0, sticky='w');self.var_norm = tk.StringVar(value=NORM_DFF);self.cmb_norm = ttk.Combobox(row3, values=NORM_CHOICES, textvariable=self.var_norm, state='readonly', width=24);self.cmb_norm.grid(row=0, column=1, padx=(6, 12), sticky='w');self.cmb_norm.bind('<<ComboboxSelected>>', lambda _e: self.on_norm_mode_changed());self.lbl_interval = ttk.Label(row3, text='Interval (s):');self.var_interval_start = tk.StringVar(value=str(DEFAULT_ZF_INTERVAL_START_S));self.var_interval_end = tk.StringVar(value=str(DEFAULT_ZF_INTERVAL_END_S));self.spin_interval_start = ttk.Spinbox(row3, from_=-1000000000.0, to=1000000000.0, increment=1.0, textvariable=self.var_interval_start, width=10);self.lbl_interval_to = ttk.Label(row3, text='to');self.spin_interval_end = ttk.Spinbox(row3, from_=-1000000000.0, to=1000000000.0, increment=1.0, textvariable=self.var_interval_end, width=10);self.btn_apply_norm = ttk.Button(row3, text='Apply interval', command=self.apply_normalization)
        for w in [self.spin_interval_start, self.spin_interval_end]:w.bind('<Return>', lambda _e: self.apply_normalization())
        self.lbl_interval.grid(row=0, column=2, sticky='w');self.spin_interval_start.grid(row=0, column=3, padx=(6, 4), sticky='w');self.lbl_interval_to.grid(row=0, column=4, sticky='w');self.spin_interval_end.grid(row=0, column=5, padx=(4, 12), sticky='w');self.btn_apply_norm.grid(row=0, column=6, sticky='w');self.var_axis_label_fs = tk.StringVar(value=f'{DEFAULT_AXIS_LABEL_FONTSIZE:g}');self.var_graph_title_fs = tk.StringVar(value=f'{DEFAULT_GRAPH_TITLE_FONTSIZE:g}');self.var_tick_label_fs = tk.StringVar(value=f'{DEFAULT_TICK_LABEL_FONTSIZE:g}');self.update_norm_controls_visibility();view_row = ttk.Frame(self.root, padding=(8, 0, 8, 4));view_row.pack(side=tk.TOP, fill=tk.X);view_row.columnconfigure(2, weight=1);ttk.Label(view_row, text='File:').grid(row=0, column=0, sticky='w');self.var_view_file = tk.StringVar(value='');self.cmb_view_file = ttk.Combobox(view_row, textvariable=self.var_view_file, values=[], width=16, state='disabled');self.cmb_view_file.grid(row=0, column=1, padx=(6, 12), sticky='w');self.cmb_view_file.bind('<<ComboboxSelected>>', self.on_view_file_changed);self.lbl_view_hint = ttk.Label(view_row, text='Run analysis to load the file selector.');self.lbl_view_hint.grid(row=0, column=2, sticky='w');self.outer_tabs = ttk.Notebook(self.root);self.outer_tabs.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=(0, 6));self._mouse_frames: Dict[str, ttk.Frame] = {};self._frame_to_mid: Dict[str, str] = {};self.outer_tabs.bind('<<NotebookTabChanged>>', self._on_outer_tab_changed);self.status = tk.StringVar(value='Ready.');status_bar = ttk.Label(self.root, textvariable=self.status, relief=tk.SUNKEN, anchor='w');status_bar.pack(side=tk.BOTTOM, fill=tk.X);self._sync_artifact_controls_state();self._install_time_marker_hotkeys();self._install_file_map_hotkeys()
        if self.csv_paths and autorun:self.root.after(0, self.run_analysis)
    def _install_file_map_hotkeys(self) -> None:self.root.bind_all('<Control-j>', self._on_ctrl_j, add='+');self.root.bind_all('<Control-J>', self._on_ctrl_j, add='+')
    def _on_ctrl_j(self, _event=None) -> None:
        try:self._show_file_number_map_dialog()
        except Exception as e:
            try:print(f'_on_ctrl_j error: {e}', file=sys.stderr)
            except Exception:pass
            try:self.status.set(f'Error showing file map: {e}')
            except Exception:pass
    def _show_file_number_map_dialog(self) -> None:
        paths = list(self.csv_paths)
        if not paths:messagebox.showinfo('File number mapping (Ctrl+J)', "No CSV files selected.\n\nUse 'Add CSV(s)…' first.");return
        top = tk.Toplevel(self.root);top.title('File number mapping (Ctrl+J)');top.transient(self.root);top.grab_set();frm = ttk.Frame(top, padding=10);frm.pack(fill=tk.BOTH, expand=True);ttk.Label(frm, text='These are the names used in the file drop-down:', justify='left').pack(anchor='w');txt_frame = ttk.Frame(frm);txt_frame.pack(fill=tk.BOTH, expand=True, pady=(8, 8));txt = tk.Text(txt_frame, width=90, height=min(18, 3 + len(paths) * 2), wrap='none');sb = ttk.Scrollbar(txt_frame, orient='vertical', command=txt.yview);txt.configure(yscrollcommand=sb.set);txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=True);sb.pack(side=tk.RIGHT, fill=tk.Y);lines = []
        for i, p in enumerate(paths, start=1):display_name = f'File {i}';basename = os.path.basename(p);lines.append(f'{display_name}: {basename}\n    {p}\n')
        txt.insert('1.0', '\n'.join(lines).rstrip() + '\n');txt.configure(state='disabled');btn_row = ttk.Frame(frm);btn_row.pack(fill=tk.X)
        def copy_to_clipboard():
            try:self.root.clipboard_clear();self.root.clipboard_append('\n'.join(lines).rstrip() + '\n')
            except Exception:pass
        ttk.Button(btn_row, text='Copy', command=copy_to_clipboard).pack(side=tk.LEFT);ttk.Button(btn_row, text='Close', command=top.destroy).pack(side=tk.RIGHT);self.root.wait_window(top)
    @staticmethod
    def _mouse_file_key(mouse_id: str) -> str:return str(mouse_id).split(':', 1)[0]
    @staticmethod
    def _mouse_channel_label(mouse_id: str) -> str:parts = str(mouse_id).split(':', 1);return parts[1] if len(parts) == 2 and parts[1] else str(mouse_id)
    def _set_result_file_choices(self, file_meta: List[Tuple[str, str, str]]) -> None:
        self._result_file_order = [file_key for file_key, _alias, _path in file_meta];self._file_alias_by_key = {file_key: alias for file_key, alias, _path in file_meta};self._file_path_by_key = {file_key: path for file_key, _alias, path in file_meta};values = [self._file_alias_by_key[file_key] for file_key in self._result_file_order]
        if not values:self.var_view_file.set('');self.cmb_view_file.config(values=(), state='disabled');self.lbl_view_hint.config(text='Run analysis to load the file selector.');return
        current = (self.var_view_file.get() or '').strip()
        if current not in values:current = values[0]
        self.cmb_view_file.config(values=tuple(values), state='readonly');self.var_view_file.set(current)
        if len(values) == 1:self.lbl_view_hint.config(text="Showing the selected file's G channels below.")
        else:self.lbl_view_hint.config(text='Pick a file to show its G channels below.')
    def _get_selected_view_file_key(self) -> Optional[str]:
        if not self._result_file_order:return None
        selected_alias = (self.var_view_file.get() or '').strip()
        for file_key in self._result_file_order:
            if self._file_alias_by_key.get(file_key) == selected_alias:return file_key
        fallback = self._result_file_order[0];self.var_view_file.set(self._file_alias_by_key.get(fallback, ''));return fallback
    def on_view_file_changed(self, _event=None) -> None:
        if not self._results:return
        self._refresh_visible_tabs();file_key = self._get_selected_view_file_key();alias = self._file_alias_by_key.get(file_key, '') if file_key else ''
        if alias:self.status.set(f'Viewing file: {alias}.')
    def _update_files_label(self):
        if not self.csv_paths:self.lbl_file.config(text='No CSV files selected.');return
        if len(self.csv_paths) == 1:self.lbl_file.config(text=self.csv_paths[0]);return
        first = os.path.basename(self.csv_paths[0]);self.lbl_file.config(text=f'{len(self.csv_paths)} CSV files selected (first: {first})')
    def add_csvs(self):
        paths = filedialog.askopenfilenames(title='Select CSV file(s)', filetypes=[('CSV', '*.csv'), ('All files', '*.*')])
        if not paths:return
        for p in paths:
            if p and p not in self.csv_paths:self.csv_paths.append(p)
        self._update_files_label();self.status.set(f'Selected {len(self.csv_paths)} file(s).')
    def clear_csvs(self):self.csv_paths = [];self._update_files_label();self.status.set('Cleared file list.')
    def get_mouse_items(self) -> List[Tuple[str, str]]:
        items: List[Tuple[str, str]] = []
        for mid in list(self._mouse_order):items.append((mid, self._mouse_display.get(mid, mid)))
        return items
    def get_mouse_result(self, mouse_id: str) -> Optional[ChannelResult]:
        if not self._results:return None
        return self._results.get(mouse_id)
    def get_mouse_display(self, mouse_id: str) -> str:return self._mouse_display.get(mouse_id, mouse_id)
    def get_norm_mode(self) -> str:return self._read_norm_mode()
    def get_smooth_window(self) -> int:
        try:return max(1, int(float(self.var_smooth_win.get())))
        except Exception:return DEFAULT_SMOOTH_WINDOW
    def _install_time_marker_hotkeys(self) -> None:self.root.bind_all('<Control-i>', self._on_ctrl_i, add='+');self.root.bind_all('<Control-I>', self._on_ctrl_i, add='+')
    def _on_ctrl_i(self, _event=None) -> None:
        tab = self.find_active_plot_tab()
        if tab is None:return
        self._show_time_marker_dialog(tab)
    def _show_time_marker_dialog(self, tab: 'PlotTabTk') -> None:
        top = tk.Toplevel(self.root);top.title('Time marker (Ctrl+I)');top.transient(self.root);top.grab_set();time_var = tk.StringVar(value=getattr(self, '_last_marker_time_str', ''));info = 'Add a vertical time marker to the *currently active* plot.\n\n• To ADD: type a time in seconds and press Enter / OK.\n• To REMOVE the last marker: press Ctrl + Backspace.\n\nTip: marker lines appear in the legend, so you can rename them by double-clicking.';frm = ttk.Frame(top, padding=10);frm.pack(fill=tk.BOTH, expand=True);ttk.Label(frm, text=info, justify='left').pack(anchor='w');ttk.Label(frm, text='Time (seconds):').pack(anchor='w', pady=(10, 0));entry = ttk.Entry(frm, textvariable=time_var, width=28);entry.pack(anchor='w', pady=(2, 8));entry.focus_set();btn_row = ttk.Frame(frm);btn_row.pack(fill=tk.X, pady=(4, 0))
        def close_dialog() -> None:
            try:top.grab_release()
            except Exception:pass
            top.destroy()
        def do_add(_e=None) -> str:
            s = (time_var.get() or '').strip()
            if not s:messagebox.showerror('Missing time', 'Please type a time in seconds (for example: 12.5).');return 'break'
            try:t = float(s)
            except Exception:messagebox.showerror('Invalid time', f"Could not read '{s}' as a number of seconds.");return 'break'
            if not np.isfinite(t):messagebox.showerror('Invalid time', 'Time must be a finite number.');return 'break'
            self._last_marker_time_str = s;tab.add_time_marker(t);close_dialog();return 'break'
        def do_remove(_e=None) -> str:tab.remove_last_time_marker();close_dialog();return 'break'
        def do_cancel(_e=None) -> str:close_dialog();return 'break'
        ttk.Button(btn_row, text='OK (add)', command=do_add).pack(side=tk.RIGHT);ttk.Button(btn_row, text='Cancel', command=do_cancel).pack(side=tk.RIGHT, padx=(0, 8));top.bind('<Return>', do_add);top.bind('<Escape>', do_cancel);top.bind('<Control-BackSpace>', do_remove);entry.bind('<Control-BackSpace>', do_remove);self.root.wait_window(top)
    def _sync_artifact_controls_state(self):
        enabled = bool(self.var_artifact_enabled.get());self.spin_factor.config(state='normal' if enabled else 'disabled');self.spin_pad.config(state='normal' if enabled else 'disabled')
        if enabled:self.chk_shared.state(['!disabled'])
        else:self.chk_shared.state(['disabled'])
    def _read_norm_mode(self) -> str:mode = (self.var_norm.get() or NORM_DFF).strip();return mode if mode in NORM_CHOICES else NORM_DFF
    def update_norm_controls_visibility(self):
        mode = self._read_norm_mode();show_interval = mode == NORM_ZF_INTERVAL
        if show_interval:self.lbl_interval.grid();self.spin_interval_start.grid();self.lbl_interval_to.grid();self.spin_interval_end.grid();self.btn_apply_norm.grid()
        else:self.lbl_interval.grid_remove();self.spin_interval_start.grid_remove();self.lbl_interval_to.grid_remove();self.spin_interval_end.grid_remove();self.btn_apply_norm.grid_remove()
    def on_norm_mode_changed(self):self.update_norm_controls_visibility();self.apply_normalization()
    def apply_axis_label_fontsize(self):
        try:fs = float(self.var_axis_label_fs.get())
        except Exception as e:messagebox.showerror('Invalid font size', f'Could not parse axis label font size:\n\n{e}');return
        if not np.isfinite(fs) or fs <= 0:messagebox.showerror('Invalid font size', 'Axis label font size must be a positive number.');return
        self._axis_label_fs_override = fs
        for widget in self._channel_widgets.values():widget.set_axis_label_fontsize(fs)
        if self.compare_widget is not None:self.compare_widget.set_axis_label_fontsize(fs)
        if self.average_widget is not None:self.average_widget.set_axis_label_fontsize(fs)
        self.status.set(f'Applied axis label font size = {fs:g} (x/y labels).')
    def apply_graph_title_fontsize(self):
        try:fs = float(self.var_graph_title_fs.get())
        except Exception as e:messagebox.showerror('Invalid font size', f'Could not parse graph title font size:\n\n{e}');return
        if not np.isfinite(fs) or fs <= 0:messagebox.showerror('Invalid font size', 'Graph title font size must be a positive number.');return
        self._graph_title_fs_override = fs
        for widget in self._channel_widgets.values():widget.set_graph_title_fontsize(fs)
        if self.compare_widget is not None:self.compare_widget.set_graph_title_fontsize(fs)
        if self.average_widget is not None:self.average_widget.set_graph_title_fontsize(fs)
        self.status.set(f'Applied graph title font size = {fs:g} (plot titles).')
    def apply_tick_label_fontsize(self):
        try:fs = float(self.var_tick_label_fs.get())
        except Exception as e:messagebox.showerror('Invalid font size', f'Could not parse tick label font size:\n\n{e}');return
        if not np.isfinite(fs) or fs <= 0:messagebox.showerror('Invalid font size', 'Tick label font size must be a positive number.');return
        self._tick_label_fs_override = fs
        for widget in self._channel_widgets.values():widget.set_tick_label_fontsize(fs)
        if self.compare_widget is not None:self.compare_widget.set_tick_label_fontsize(fs)
        if self.average_widget is not None:self.average_widget.set_tick_label_fontsize(fs)
        self.status.set(f'Applied tick label font size = {fs:g} (axis numbers).')
    def on_artifact_enabled_toggled(self):
        enabled = bool(self.var_artifact_enabled.get());self._sync_artifact_controls_state()
        if not self._results:self.status.set(f"Artifact remover {('ENABLED' if enabled else 'DISABLED')} (will apply on next Run).");return
        try:artifact_factor = float(self.var_factor.get());artifact_pad = int(float(self.var_pad.get()));require_shared = bool(self.var_shared.get());align_mode = DEFAULT_ALIGN_MODE;use_interp = bool(self.var_interp.get())
        except Exception as e:messagebox.showerror('Invalid settings', f'Could not parse artifact settings:\n\n{e}');return
        self.status.set('Updating artifact pipeline…')
        def worker():
            try:
                for res in self._results.values():recompute_artifact_pipeline_inplace(res, artifact_enabled=enabled, artifact_factor=artifact_factor, artifact_pad=artifact_pad, require_shared=require_shared, align_mode=align_mode, use_linear_interp=use_interp)
                self.root.after(0, self._after_artifact_pipeline_updated)
            except Exception as e:self.root.after(0, lambda: self.on_analysis_failed(str(e)))
        threading.Thread(target=worker, daemon=True).start()
    def _after_artifact_pipeline_updated(self):
        for widget in self._channel_widgets.values():widget.refresh_after_pipeline_change()
        if self.compare_widget is not None:self.compare_widget.refresh_plot()
        if self.average_widget is not None:self.average_widget.refresh_plot()
        self.status.set(f"Artifact remover {('ENABLED' if self.var_artifact_enabled.get() else 'DISABLED')} (MAD). Interp={('ON' if self.var_interp.get() else 'OFF')}.")
    def on_interp_toggled(self):
        if not self._results:return
        use_interp = bool(self.var_interp.get())
        for res in self._results.values():set_interpolation_mode(res, use_interp)
        for widget in self._channel_widgets.values():widget.refresh_after_pipeline_change()
        if self.compare_widget is not None:self.compare_widget.refresh_plot()
        if self.average_widget is not None:self.average_widget.refresh_plot()
        self.status.set(f"Interpolation {('ON' if use_interp else 'OFF')} (frequency analysis still uses NO interpolation).")
    def apply_normalization(self):
        mode = self._read_norm_mode()
        try:smooth_win = max(1, int(float(self.var_smooth_win.get())))
        except Exception:smooth_win = DEFAULT_SMOOTH_WINDOW
        interval_start: Optional[float] = None;interval_end: Optional[float] = None
        if mode == NORM_ZF_INTERVAL:
            try:interval_start = float(self.var_interval_start.get());interval_end = float(self.var_interval_end.get())
            except Exception as e:messagebox.showerror('Invalid interval', f'Could not parse interval start/end:\n\n{e}');return
        if not self._results:
            if mode == NORM_DFF:self.status.set('Normalization view: ΔF/F')
            elif mode == NORM_ZF_GLOBAL:self.status.set('Normalization view: zF (global, GUI)')
            elif mode == NORM_ZF_INTERVAL:self.status.set('Normalization view: zF - interval based (set interval and apply)')
            return
        for res in self._results.values():
            res.smooth_window = smooth_win
            if mode == NORM_ZF_INTERVAL and interval_start is not None and (interval_end is not None):res.zf_interval_start_s = interval_start;res.zf_interval_end_s = interval_end
            recompute_normalizations(res)
        for widget in self._channel_widgets.values():widget.set_norm_mode(mode)
        if self.compare_widget is not None:self.compare_widget.refresh_plot()
        if self.average_widget is not None:self.average_widget.refresh_plot()
        if mode == NORM_DFF:self.status.set('Normalization view: ΔF/F.')
        elif mode == NORM_ZF_GLOBAL:self.status.set('Normalization view: zF (global, GUI).')
        elif mode == NORM_ZF_INTERVAL:a = float(interval_start);b = float(interval_end);self.status.set(f'Normalization view: zF - interval based. Interval = {min(a, b):g}–{max(a, b):g} s')
    def run_analysis(self):
        if not self.csv_paths:messagebox.showwarning('No files', 'Please add one or more CSV files first.');return
        try:artifact_enabled = bool(self.var_artifact_enabled.get());artifact_factor = float(self.var_factor.get());artifact_pad = int(float(self.var_pad.get()));require_shared = bool(self.var_shared.get());align_mode = DEFAULT_ALIGN_MODE;acq_fps = float(self.var_acq_fps.get());smooth_win = max(1, int(float(self.var_smooth_win.get())));mode = self._read_norm_mode();interval_start = float(self.var_interval_start.get());interval_end = float(self.var_interval_end.get());use_interp = bool(self.var_interp.get())
        except Exception as e:messagebox.showerror('Invalid settings', f'Could not parse settings:\n\n{e}');return
        fit_windows = list(DEFAULT_FIT_WINDOWS);paths = list(self.csv_paths);file_labels = [f'File {i}' for i in range(1, len(paths) + 1)];file_meta = [(f'F{i}', label, path) for i, (path, label) in enumerate(zip(paths, file_labels), start=1)];self.btn_run.config(state=tk.DISABLED);self.btn_add.config(state=tk.DISABLED);self.btn_clear.config(state=tk.DISABLED);self.status.set('Analyzing…')
        def worker():
            try:
                aggregated: Dict[str, ChannelResult] = {};display: Dict[str, str] = {};order: List[str] = []
                for i, (file_key, alias, path) in enumerate(file_meta, start=1):
                    self.root.after(0, lambda i=i, alias=alias: self.status.set(f'Analyzing {i}/{len(paths)}: {alias}…'));per_file = analyze_csv(path, artifact_enabled=artifact_enabled, artifact_factor=artifact_factor, artifact_method='mad', artifact_pad=artifact_pad, require_shared=require_shared, align_mode=align_mode, fit_windows=fit_windows, acq_fps_hz=acq_fps if np.isfinite(acq_fps) and acq_fps > 0 else None, smooth_window=smooth_win, zf_interval_start_s=interval_start, zf_interval_end_s=interval_end, use_linear_interp=use_interp)
                    for gcol in sorted(per_file.keys()):mid = f'{file_key}:{gcol}';aggregated[mid] = per_file[gcol];display[mid] = f'{alias}:{gcol}';order.append(mid)
                self.root.after(0, lambda: self.on_analysis_finished(aggregated, display, order, mode, file_meta))
            except Exception as e:self.root.after(0, lambda: self.on_analysis_failed(str(e)))
        self._analysis_thread = threading.Thread(target=worker, daemon=True);self._analysis_thread.start()
    def on_analysis_finished(self, results: Dict[str, ChannelResult], display: Dict[str, str], order: List[str], mode: str, file_meta: List[Tuple[str, str, str]]):
        self._results = results;self._mouse_display = dict(display);self._mouse_order = list(order);self._set_result_file_choices(file_meta);self.build_tabs(norm_mode=mode)
        try:any_key = next(iter(results.keys()));r0 = results[any_key];self.status.set(f"Done. Files={len(file_meta)} | Mice={len(results)} | Mode={mode}. Smooth win={r0.smooth_window}. Artifacts={('ON' if self.var_artifact_enabled.get() else 'OFF')}. Interp={('ON' if r0.use_interpolation else 'OFF')}. Freq eff fs={r0.eff_fs_hz:.3g} Hz (Acq FPS={r0.acq_fps_hz:.3g}).")
        except Exception:self.status.set('Done.')
        self.btn_run.config(state=tk.NORMAL);self.btn_add.config(state=tk.NORMAL);self.btn_clear.config(state=tk.NORMAL);self.on_norm_mode_changed()
    def on_analysis_failed(self, msg: str):messagebox.showerror('Analysis failed', msg);self.status.set('Analysis failed.');self.btn_run.config(state=tk.NORMAL);self.btn_add.config(state=tk.NORMAL);self.btn_clear.config(state=tk.NORMAL)
    def build_tabs(self, norm_mode: str):
        outer_tabs = self.outer_tabs
        for tab_id in outer_tabs.tabs():outer_tabs.forget(tab_id)
        for frame in list(getattr(self, '_mouse_frames', {}).values()):
            try:frame.destroy()
            except Exception:pass
        if self.compare_widget is not None:
            try:self.compare_widget.destroy()
            except Exception:pass
        if self.average_widget is not None:
            try:self.average_widget.destroy()
            except Exception:pass
        self._channel_widgets.clear();self._mouse_frames = {};self._frame_to_mid = {};self.compare_widget = None;self.average_widget = None
        if not self._results:return
        mids = [mid for mid in self._mouse_order if mid in self._results]
        for mid in mids:frame = ttk.Frame(outer_tabs);display_label = self._mouse_display.get(mid, mid);ttk.Label(frame, text=f'Click this tab to load plots for {display_label} (lazy-loaded for speed).', foreground='#666666').pack(anchor='w', padx=10, pady=10);self._mouse_frames[mid] = frame
        self._refresh_visible_tabs()
    def _refresh_visible_tabs(self) -> None:
        outer_tabs = self.outer_tabs
        for tab_id in outer_tabs.tabs():outer_tabs.forget(tab_id)
        self._frame_to_mid = {};file_key = self._get_selected_view_file_key();mids: List[str] = []
        for mid in self._mouse_order:
            if mid not in self._results:continue
            if file_key is not None and self._mouse_file_key(mid) != file_key:continue
            mids.append(mid)
        for mid in mids:
            frame = self._mouse_frames.get(mid)
            if frame is None:continue
            tab_label = self._mouse_channel_label(mid) if file_key is not None else self._mouse_display.get(mid, mid);outer_tabs.add(frame, text=tab_label);self._frame_to_mid[str(frame)] = mid
        if mids:outer_tabs.select(self._mouse_frames[mids[0]]);self._ensure_mouse_widget(mids[0])
    def _on_outer_tab_changed(self, _event=None) -> None:
        sel = self.outer_tabs.select()
        if not sel:return
        mid = self._frame_to_mid.get(sel)
        if mid:self._ensure_mouse_widget(mid)
    def _ensure_mouse_widget(self, mid: str) -> None:
        if mid in self._channel_widgets:return
        frame = self._mouse_frames.get(mid)
        if frame is None:return
        for child in frame.winfo_children():child.destroy()
        res = self._results[mid];widget = ChannelTabsTk(frame, res, norm_mode=self._read_norm_mode(), parent_app=self);widget.pack(fill=tk.BOTH, expand=True)
        if self._axis_label_fs_override is not None:widget.set_axis_label_fontsize(self._axis_label_fs_override)
        if self._graph_title_fs_override is not None:widget.set_graph_title_fontsize(self._graph_title_fs_override)
        if self._tick_label_fs_override is not None:widget.set_tick_label_fontsize(self._tick_label_fs_override)
        self._channel_widgets[mid] = widget
    def find_active_plot_tab(self) -> Optional[PlotTabTk]:
        sel = self.outer_tabs.select()
        if not sel:return None
        mid = self._frame_to_mid.get(sel)
        if not mid:return None
        self._ensure_mouse_widget(mid);cw = self._channel_widgets.get(mid)
        if cw is None:return None
        return cw.get_active_plot_tab()
def parse_cli():parser = argparse.ArgumentParser(description='Fiberlyse GUI (batch + markers)');parser.add_argument('--csv', nargs='*', default=[], help='Optional CSV file(s) to load automatically (space-separated).');return parser.parse_args()
def main():
    args = parse_cli();app = MainAppTk()
    if getattr(args, 'csv', None):
        app.csv_paths = list(args.csv)
        try:app._update_files_label()
        except Exception:pass
        try:app.run_analysis()
        except Exception:pass
    app.root.mainloop()
def _fiberlyse_auc_insert_zero_crossings(x, y):
    x = np.asarray(x, dtype=float);y = np.asarray(y, dtype=float)
    if x.size < 2 or y.size < 2 or x.size != y.size:return (x, y)
    xs = [float(x[0])];ys = [float(y[0])]
    for k in range(1, x.size):
        x0 = float(x[k - 1]);x1 = float(x[k]);y0 = float(y[k - 1]);y1 = float(y[k])
        if np.isfinite(y0) and np.isfinite(y1) and (y0 != 0.0) and (y1 != 0.0):
            if y0 < 0.0 < y1 or y1 < 0.0 < y0:
                denom = y1 - y0
                if abs(denom) > 1e-15:
                    frac = -y0 / denom;xc = x0 + frac * (x1 - x0)
                    if min(x0, x1) < xc < max(x0, x1):xs.append(float(xc));ys.append(0.0)
        xs.append(x1);ys.append(y1)
    return (np.asarray(xs, dtype=float), np.asarray(ys, dtype=float))
def _fiberlyse_auc_stats_for_xy(x, y, start_x, end_x, baseline_y=0.0):
    try:lo = float(min(start_x, end_x));hi = float(max(start_x, end_x));baseline = float(baseline_y)
    except Exception:return None
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:return None
    if not np.isfinite(baseline):baseline = 0.0
    x = np.asarray(x, dtype=float).reshape(-1);y = np.asarray(y, dtype=float).reshape(-1);n = min(x.size, y.size)
    if n < 2:return None
    x = x[:n];y = y[:n];finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 2:return None
    trapz = getattr(np, 'trapezoid', None)
    if trapz is None:trapz = np.trapz
    signed_auc = 0.0;abs_auc = 0.0;positive_auc = 0.0;negative_auc = 0.0;coverage_x = 0.0;n_points_used = 0
    for i0, i1 in _contiguous_true_runs(finite):
        xs = x[i0:i1 + 1].astype(float);ys = y[i0:i1 + 1].astype(float)
        if xs.size < 2:continue
        if np.any(np.diff(xs) < 0):order = np.argsort(xs);xs = xs[order];ys = ys[order]
        unique_x, unique_idx = np.unique(xs, return_index=True);xs = unique_x.astype(float);ys = ys[unique_idx].astype(float)
        if xs.size < 2:continue
        left = max(lo, float(xs[0]));right = min(hi, float(xs[-1]))
        if right <= left:continue
        inside = (xs >= left) & (xs <= right);xi = xs[inside].astype(float);yi = ys[inside].astype(float)
        if xi.size == 0:xi = np.array([left, right], dtype=float);yi = np.interp(xi, xs, ys).astype(float)
        else:
            if not np.isclose(float(xi[0]), left, rtol=0.0, atol=1e-12):xi = np.insert(xi, 0, left);yi = np.insert(yi, 0, float(np.interp(left, xs, ys)))
            if not np.isclose(float(xi[-1]), right, rtol=0.0, atol=1e-12):xi = np.append(xi, right);yi = np.append(yi, float(np.interp(right, xs, ys)))
        if xi.size < 2:continue
        yb = yi - baseline;xi_z, yb_z = _fiberlyse_auc_insert_zero_crossings(xi, yb);signed_auc += float(trapz(yb_z, xi_z));abs_auc += float(trapz(np.abs(yb_z), xi_z));positive_auc += float(trapz(np.maximum(yb_z, 0.0), xi_z));negative_auc += float(trapz(np.minimum(yb_z, 0.0), xi_z));coverage_x += float(right - left);n_points_used += int(xi_z.size)
    if coverage_x <= 0.0 or n_points_used < 2:return None
    requested_x = hi - lo;coverage_percent = 100.0 * coverage_x / requested_x if requested_x > 0 else np.nan;mean_minus_baseline = signed_auc / coverage_x if coverage_x > 0 else np.nan;return {'signed_auc': float(signed_auc), 'abs_auc': float(abs_auc), 'positive_auc': float(positive_auc), 'negative_auc': float(negative_auc), 'coverage_x': float(coverage_x), 'coverage_percent': float(coverage_percent), 'mean_minus_baseline': float(mean_minus_baseline), 'n_points_used': float(n_points_used)}
def _fiberlyse_plot_calculate_auc_for_interval(self, start_x, end_x, baseline_y=0.0):
    rows = []
    try:self._apply_user_overrides()
    except Exception:pass
    for ai, ax in enumerate(list(getattr(self.fig, 'axes', [])), start=1):
        axis_title = str(ax.get_title() or '').strip()
        if not axis_title:axis_title = f'Axis {ai}'
        for line in list(getattr(ax, 'lines', [])):
            try:
                if not line.get_visible():continue
            except Exception:continue
            label = str(line.get_label() or '').strip()
            if not label or label.startswith('_'):continue
            try:x = np.asarray(line.get_xdata(), dtype=float);y = np.asarray(line.get_ydata(), dtype=float)
            except Exception:continue
            stats = _fiberlyse_auc_stats_for_xy(x=x, y=y, start_x=start_x, end_x=end_x, baseline_y=baseline_y)
            if stats is None:continue
            row = {'axis': axis_title, 'trace': label};row.update(stats);rows.append(row)
    return rows
def _fiberlyse_auc_fmt(v):
    try:vf = float(v)
    except Exception:return 'n/a'
    if not np.isfinite(vf):return 'n/a'
    return f'{vf:.8g}'
def _fiberlyse_plot_format_auc_results_text(self, rows, start_x, end_x, baseline_y):
    lo = float(min(start_x, end_x));hi = float(max(start_x, end_x));lines = [f'AUC results for plot: {self.tab_name}', f'Interval on x-axis: {lo:g} to {hi:g}', f'Baseline: B = {baseline_y:g}', '', 'PLAIN-LANGUAGE EXPLANATION', '', 'Signed AUC:', 'Net area after baseline subtraction. Area above baseline is positive; area below baseline is negative. Positive and negative parts can cancel each other out.', '', 'Absolute AUC:', 'Total area away from baseline, regardless of direction. Both increases and decreases count as positive area.', '', 'Positive AUC:', 'Only the area above baseline. Any part of the signal below baseline is treated as zero.', '', 'Negative AUC:', 'Only the area below baseline. This value is usually negative.', '', 'Coverage:', 'The amount of your selected x/time interval that had usable finite data. NaN holes or missing regions are not integrated across.', '', 'Coverage %:', 'The percentage of the requested interval that actually had usable data.', '', 'Mean - baseline:', 'The average baseline-corrected signal over the usable part of the selected interval.', '', 'Points used:', 'The number of x/y points used for integration, including inserted interval-boundary points and zero-crossing points.', '', 'MATHEMATICAL FORMULAS', '', 'Notation:', 'B = baseline', 'yᵦ(x) = y(x) - B', 'U = usable parts of [x_start, x_end] after removing NaN/missing segments', 'Δxᵢ = xᵢ₊₁ - xᵢ', '', 'Trapezoidal rule used by the calculator:', '∫ f(x) dx ≈ Σᵢ ½ · [f(xᵢ) + f(xᵢ₊₁)] · Δxᵢ', '', 'Signed AUC:', 'A_signed ≈ Σᵢ ½ · [yᵦᵢ + yᵦᵢ₊₁] · Δxᵢ', '', 'Absolute AUC:', 'A_abs ≈ Σᵢ ½ · [|yᵦᵢ| + |yᵦᵢ₊₁|] · Δxᵢ', '', 'Positive AUC:', 'A_pos ≈ Σᵢ ½ · [max(yᵦᵢ, 0) + max(yᵦᵢ₊₁, 0)] · Δxᵢ', '', 'Negative AUC:', 'A_neg ≈ Σᵢ ½ · [min(yᵦᵢ, 0) + min(yᵦᵢ₊₁, 0)] · Δxᵢ', '', 'Coverage:', 'C = Σₖ (bₖ - aₖ)', 'where [aₖ, bₖ] are the usable finite data segments inside the selected interval.', '', 'Coverage %:', 'C_% = 100 · C / (x_end - x_start)', '', 'Mean - baseline:', 'μᵦ = A_signed / C', '', 'On time plots, x is time in seconds, so AUC units are y-units · seconds.', '', '\t'.join(['Axis', 'Trace', 'Signed AUC', 'Absolute AUC', 'Positive AUC', 'Negative AUC', 'Coverage', 'Coverage %', 'Mean - baseline', 'Points used'])]
    for r in rows:lines.append('\t'.join([str(r.get('axis', '')), str(r.get('trace', '')), _fiberlyse_auc_fmt(r.get('signed_auc')), _fiberlyse_auc_fmt(r.get('abs_auc')), _fiberlyse_auc_fmt(r.get('positive_auc')), _fiberlyse_auc_fmt(r.get('negative_auc')), _fiberlyse_auc_fmt(r.get('coverage_x')), _fiberlyse_auc_fmt(r.get('coverage_percent')), _fiberlyse_auc_fmt(r.get('mean_minus_baseline')), _fiberlyse_auc_fmt(r.get('n_points_used'))]))
    return '\n'.join(lines)
def _fiberlyse_plot_default_auc_interval_strings(self):
    ax_default = None
    try:
        for ax in list(getattr(self.fig, 'axes', [])):
            if self._axis_looks_like_time(ax):ax_default = ax;break
    except Exception:ax_default = None
    if ax_default is None:
        try:ax_default = self.fig.axes[0] if self.fig.axes else None
        except Exception:ax_default = None
    if ax_default is None:return ('0', '1')
    try:
        x0, x1 = ax_default.get_xlim();lo = min(float(x0), float(x1));hi = max(float(x0), float(x1))
        if np.isfinite(lo) and np.isfinite(hi) and (hi > lo):return (f'{lo:g}', f'{hi:g}')
    except Exception:pass
    return ('0', '1')
def _fiberlyse_tree_sort(tree, col, reverse=False):
    data = []
    for item in tree.get_children(''):
        val = tree.set(item, col)
        try:key = float(val)
        except Exception:key = str(val).lower()
        data.append((key, item))
    data.sort(reverse=bool(reverse))
    for index, (_key, item) in enumerate(data):tree.move(item, '', index)
    tree.heading(col, command=lambda: _fiberlyse_tree_sort(tree, col, not reverse))
def _fiberlyse_plot_show_auc_results_window(self, rows, start_x, end_x, baseline_y):
    parent = self.winfo_toplevel();text_out = self._format_auc_results_text(rows, start_x, end_x, baseline_y);lo = float(min(start_x, end_x));hi = float(max(start_x, end_x));requested = hi - lo;top = tk.Toplevel(parent);top.title('AUC results');top.transient(parent);_fiberlyse_set_toplevel_pixel_geometry(top, parent, 1160, 840);main = ttk.Frame(top, padding=10);main.pack(fill=tk.BOTH, expand=True);main.rowconfigure(1, weight=1);main.rowconfigure(2, weight=0);main.columnconfigure(0, weight=1);summary = f'Plot: {self.tab_name} | Interval: {lo:g} to {hi:g} (requested width {requested:g}) | Baseline: B = {baseline_y:g}';ttk.Label(main, text=summary, justify='left').grid(row=0, column=0, sticky='w');table_frame = ttk.LabelFrame(main, text='AUC table');table_frame.grid(row=1, column=0, sticky='nsew', pady=(8, 8));table_frame.rowconfigure(0, weight=1);table_frame.columnconfigure(0, weight=1);columns = ['axis', 'trace', 'signed_auc', 'abs_auc', 'positive_auc', 'negative_auc', 'coverage_x', 'coverage_percent', 'mean_minus_baseline', 'n_points_used'];headings = {'axis': 'Axis', 'trace': 'Trace', 'signed_auc': 'Signed AUC', 'abs_auc': 'Absolute AUC', 'positive_auc': 'Positive AUC', 'negative_auc': 'Negative AUC', 'coverage_x': 'Coverage', 'coverage_percent': 'Coverage %', 'mean_minus_baseline': 'Mean - baseline', 'n_points_used': 'Points used'};widths = {'axis': 220, 'trace': 190, 'signed_auc': 105, 'abs_auc': 105, 'positive_auc': 105, 'negative_auc': 105, 'coverage_x': 90, 'coverage_percent': 90, 'mean_minus_baseline': 120, 'n_points_used': 90};tree = ttk.Treeview(table_frame, columns=columns, show='headings', height=12);ysb = ttk.Scrollbar(table_frame, orient='vertical', command=tree.yview);xsb = ttk.Scrollbar(table_frame, orient='horizontal', command=tree.xview);tree.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set);tree.grid(row=0, column=0, sticky='nsew');ysb.grid(row=0, column=1, sticky='ns');xsb.grid(row=1, column=0, sticky='ew')
    for col in columns:tree.heading(col, text=headings[col], command=lambda c=col: _fiberlyse_tree_sort(tree, c, False));tree.column(col, width=widths.get(col, 100), minwidth=60, stretch=col in ['axis', 'trace'])
    for r in rows:
        vals = []
        for col in columns:
            if col in ['axis', 'trace']:vals.append(str(r.get(col, '')))
            else:vals.append(_fiberlyse_auc_fmt(r.get(col)))
        tree.insert('', tk.END, values=vals)
    explanation_frame = ttk.LabelFrame(main, text='Explanation and formulas');explanation_frame.grid(row=2, column=0, sticky='ew', pady=(0, 8));explanation_frame.rowconfigure(0, weight=1);explanation_frame.columnconfigure(0, weight=1);explanation = 'PLAIN-LANGUAGE EXPLANATION\n\nSigned AUC:\nNet area after baseline subtraction. Area above baseline is positive; area below baseline is negative. Positive and negative parts can cancel each other out.\n\nAbsolute AUC:\nTotal area away from baseline, regardless of direction. Both increases and decreases count as positive area.\n\nPositive AUC:\nOnly the area above baseline. Any part of the signal below baseline is treated as zero.\n\nNegative AUC:\nOnly the area below baseline. This value is usually negative.\n\nCoverage:\nThe amount of your selected x/time interval that had usable finite data. NaN holes or missing regions are not integrated across.\n\nCoverage %:\nThe percentage of the requested interval that actually had usable data.\n\nMean - baseline:\nThe average baseline-corrected signal over the usable part of the selected interval.\n\nPoints used:\nThe number of x/y points used for integration, including inserted interval-boundary points and zero-crossing points.\n\nMATHEMATICAL FORMULAS\n\nNotation:\nB = baseline\nyᵦ(x) = y(x) - B\nU = usable parts of [x_start, x_end] after removing NaN/missing segments\nΔxᵢ = xᵢ₊₁ - xᵢ\n\nTrapezoidal rule:\n∫ f(x) dx ≈ Σᵢ ½ · [f(xᵢ) + f(xᵢ₊₁)] · Δxᵢ\n\nSigned AUC:\nA_signed ≈ Σᵢ ½ · [yᵦᵢ + yᵦᵢ₊₁] · Δxᵢ\n\nAbsolute AUC:\nA_abs ≈ Σᵢ ½ · [|yᵦᵢ| + |yᵦᵢ₊₁|] · Δxᵢ\n\nPositive AUC:\nA_pos ≈ Σᵢ ½ · [max(yᵦᵢ, 0) + max(yᵦᵢ₊₁, 0)] · Δxᵢ\n\nNegative AUC:\nA_neg ≈ Σᵢ ½ · [min(yᵦᵢ, 0) + min(yᵦᵢ₊₁, 0)] · Δxᵢ\n\nCoverage:\nC = Σₖ (bₖ - aₖ)\nwhere [aₖ, bₖ] are the usable finite data segments inside the selected interval.\n\nCoverage %:\nC_% = 100 · C / (x_end - x_start)\n\nMean - baseline:\nμᵦ = A_signed / C\n\nOn time plots, x is time in seconds, so AUC units are y-units · seconds.';txt = tk.Text(explanation_frame, height=18, wrap='word');exp_scroll = ttk.Scrollbar(explanation_frame, orient='vertical', command=txt.yview);txt.configure(yscrollcommand=exp_scroll.set);txt.grid(row=0, column=0, sticky='nsew', padx=(8, 0), pady=6);exp_scroll.grid(row=0, column=1, sticky='ns', padx=(4, 8), pady=6);txt.insert('1.0', explanation);txt.configure(state='disabled');btn_row = ttk.Frame(main);btn_row.grid(row=3, column=0, sticky='ew')
    def export_auc_table_csv():
        try:
            path = filedialog.asksaveasfilename(title='Export AUC table (CSV)', defaultextension='.csv', initialfile=f'{self.default_filename_prefix}_{self.tab_name}_AUC.csv'.strip('_'), filetypes=[('CSV', '*.csv'), ('All files', '*.*')], parent=top)
            if not path:return
            csv_columns = [headings[col] for col in columns];csv_rows = []
            for item in tree.get_children(''):
                vals = list(tree.item(item, 'values'))
                csv_rows.append({csv_columns[i]: (vals[i] if i < len(vals) else '') for i in range(len(csv_columns))})
            pd.DataFrame(csv_rows, columns=csv_columns).to_csv(path, index=False)
            messagebox.showinfo('Export complete', f'Saved CSV file:\n{path}', parent=top)
        except Exception as e:messagebox.showerror('Export failed', f'Could not export AUC table CSV:\n\n{e}', parent=top)
    ttk.Button(btn_row, text='Export AUC table (CSV)...', command=export_auc_table_csv).pack(side=tk.LEFT);ttk.Button(btn_row, text='Close', command=top.destroy).pack(side=tk.RIGHT)
def _fiberlyse_plot_show_auc_dialog(self):
    parent = self.winfo_toplevel()
    if not getattr(self.fig, 'axes', []):messagebox.showinfo('AUC', 'No axes/plot found in this graph.');return
    default_start, default_end = self._default_auc_interval_strings();start_var = tk.StringVar(value=str(getattr(self, '_last_auc_start_str', default_start) or default_start));end_var = tk.StringVar(value=str(getattr(self, '_last_auc_end_str', default_end) or default_end));baseline_var = tk.StringVar(value=str(getattr(self, '_last_auc_baseline_str', '0') or '0'));top = tk.Toplevel(parent);top.title('AUC interval');top.transient(parent);top.grab_set();frm = ttk.Frame(top, padding=10);frm.pack(fill=tk.BOTH, expand=True);ttk.Label(frm, text='Calculate AUC for every visible line trace in the current graph.\n\nFor time plots, enter times in seconds.\nAUC is computed relative to the baseline B you enter below.\nUse B = 0 if you want ordinary area relative to zero.', justify='left').grid(row=0, column=0, columnspan=2, sticky='w');ttk.Label(frm, text='Start x/time:').grid(row=1, column=0, sticky='w', pady=(10, 2));entry_start = ttk.Entry(frm, textvariable=start_var, width=18);entry_start.grid(row=1, column=1, sticky='w', pady=(10, 2));ttk.Label(frm, text='End x/time:').grid(row=2, column=0, sticky='w', pady=2);entry_end = ttk.Entry(frm, textvariable=end_var, width=18);entry_end.grid(row=2, column=1, sticky='w', pady=2);ttk.Label(frm, text='Baseline B:').grid(row=3, column=0, sticky='w', pady=2);ttk.Entry(frm, textvariable=baseline_var, width=18).grid(row=3, column=1, sticky='w', pady=2);small_help = 'Formula setup: yᵦ(x) = y(x) - B. Use B = 0 if you do not want to subtract another baseline.';ttk.Label(frm, text=small_help, foreground='gray35', justify='left', wraplength=420).grid(row=4, column=0, columnspan=2, sticky='w', pady=(6, 0));btn_row = ttk.Frame(frm);btn_row.grid(row=5, column=0, columnspan=2, sticky='ew', pady=(10, 0))
    def close_dialog():
        try:top.grab_release()
        except Exception:pass
        try:top.destroy()
        except Exception:pass
    def do_calculate(_e=None):
        try:start_x = float((start_var.get() or '').strip());end_x = float((end_var.get() or '').strip());baseline_y = float((baseline_var.get() or '0').strip())
        except Exception as e:messagebox.showerror('Invalid AUC interval', f'Could not parse the AUC inputs:\n\n{e}');return 'break'
        if not np.isfinite(start_x) or not np.isfinite(end_x):messagebox.showerror('Invalid AUC interval', 'Start and end must be finite numbers.');return 'break'
        if start_x == end_x:messagebox.showerror('Invalid AUC interval', 'Start and end cannot be the same value.');return 'break'
        if not np.isfinite(baseline_y):baseline_y = 0.0
        rows = self.calculate_auc_for_interval(start_x=start_x, end_x=end_x, baseline_y=baseline_y);self._last_auc_start_str = (start_var.get() or '').strip();self._last_auc_end_str = (end_var.get() or '').strip();self._last_auc_baseline_str = (baseline_var.get() or '0').strip();close_dialog()
        if not rows:messagebox.showinfo('AUC', 'No visible labeled line traces overlapped that interval.\n\nScatter-only plots, histogram patches, hidden lines, and vertical marker lines are ignored.');return 'break'
        self._show_auc_results_window(rows, start_x, end_x, baseline_y);return 'break'
    def do_cancel(_e=None):close_dialog();return 'break'
    ttk.Button(btn_row, text='Calculate AUC', command=do_calculate).pack(side=tk.RIGHT);ttk.Button(btn_row, text='Cancel', command=do_cancel).pack(side=tk.RIGHT, padx=(0, 8));top.bind('<Return>', do_calculate);top.bind('<Escape>', do_cancel);entry_start.focus_set();parent.wait_window(top)
PlotTabTk.calculate_auc_for_interval = _fiberlyse_plot_calculate_auc_for_interval;PlotTabTk._format_auc_results_text = _fiberlyse_plot_format_auc_results_text;PlotTabTk._default_auc_interval_strings = _fiberlyse_plot_default_auc_interval_strings;PlotTabTk._show_auc_results_window = _fiberlyse_plot_show_auc_results_window;PlotTabTk.show_auc_dialog = _fiberlyse_plot_show_auc_dialog;_FIBERLYSE_ORIGINAL_PLOTTABTK_INIT = PlotTabTk.__init__
def _fiberlyse_plottabtk_init_with_auc(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_INIT(self, *args, **kwargs)
    try:
        existing = getattr(self, 'auc_btn', None)
        if existing is not None:
            try:
                if bool(existing.winfo_exists()):return
            except Exception:pass
        btn_parent = self.save_btn.master;self.auc_btn = ttk.Button(btn_parent, text='AUC interval...', command=self.show_auc_dialog);self.auc_btn.pack(side=tk.RIGHT, padx=(0, 8))
    except Exception:pass
PlotTabTk.__init__ = _fiberlyse_plottabtk_init_with_auc
def _fiberlyse_install_auc_hotkeys(self):self.root.bind_all('<Control-u>', self._on_ctrl_u, add='+');self.root.bind_all('<Control-U>', self._on_ctrl_u, add='+')
def _fiberlyse_on_ctrl_u(self, _event=None):
    tab = self.find_active_plot_tab()
    if tab is None:return
    tab.show_auc_dialog()
MainAppTk._install_auc_hotkeys = _fiberlyse_install_auc_hotkeys;MainAppTk._on_ctrl_u = _fiberlyse_on_ctrl_u;_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT = MainAppTk.__init__
def _fiberlyse_mainapptk_init_with_auc(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT(self, *args, **kwargs)
    try:self._install_auc_hotkeys()
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_mainapptk_init_with_auc;_FIBERLYSE_ORIGINAL_FIND_ACTIVE_PLOT_TAB = MainAppTk.find_active_plot_tab
def _fiberlyse_find_active_plot_tab_with_auc(self):
    sel = self.outer_tabs.select()
    if not sel:return None
    try:selected_widget = self.root.nametowidget(sel)
    except Exception:selected_widget = None
    try:
        if self.compare_widget is not None and selected_widget is self.compare_widget:return self.compare_widget.plot
    except Exception:pass
    try:
        if self.average_widget is not None and selected_widget is self.average_widget:return self.average_widget.plot
    except Exception:pass
    return _FIBERLYSE_ORIGINAL_FIND_ACTIVE_PLOT_TAB(self)
MainAppTk.find_active_plot_tab = _fiberlyse_find_active_plot_tab_with_auc

# ---- Fiberlyse axis tick interval slider extension ----
# Adds a per-plot dialog for adjusting x/y major tick spacing with sliders.
try:
    from matplotlib.ticker import MultipleLocator as _FiberlyseMultipleLocator, AutoLocator as _FiberlyseAutoLocator
except Exception:
    _FiberlyseMultipleLocator = None
    _FiberlyseAutoLocator = None

def _fiberlyse_axis_float_or_none(value):
    try:
        v = float(value)
    except Exception:
        return None
    if not np.isfinite(v) or v <= 0:
        return None
    return float(v)

def _fiberlyse_nice_interval_from_axis(ax, axis_name: str) -> float:
    try:
        lim = ax.get_xlim() if axis_name == 'x' else ax.get_ylim()
        span = abs(float(lim[1]) - float(lim[0]))
    except Exception:
        span = np.nan
    if not np.isfinite(span) or span <= 0:
        span = 1.0
    try:
        axis = ax.xaxis if axis_name == 'x' else ax.yaxis
        ticks = np.asarray(axis.get_majorticklocs(), dtype=float)
        ticks = ticks[np.isfinite(ticks)]
        diffs = np.diff(np.sort(np.unique(ticks)))
        diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
        if diffs.size:
            val = float(np.median(diffs))
            if np.isfinite(val) and val > 0:
                return val
    except Exception:
        pass
    return max(span / 8.0, 1e-12)

def _fiberlyse_plot_axis_interval_bounds(self, axis_name: str):
    spans = []
    for ax in list(getattr(self.fig, 'axes', [])):
        try:
            lim = ax.get_xlim() if axis_name == 'x' else ax.get_ylim()
            span = abs(float(lim[1]) - float(lim[0]))
            if np.isfinite(span) and span > 0:
                spans.append(span)
        except Exception:
            pass
    span = max(spans) if spans else 1.0
    cur = getattr(self, f'_axis_tick_interval_{axis_name}', None)
    cur = _fiberlyse_axis_float_or_none(cur)
    if cur is None:
        ax0 = None
        try:
            ax0 = self.fig.axes[0] if self.fig.axes else None
        except Exception:
            ax0 = None
        cur = _fiberlyse_nice_interval_from_axis(ax0, axis_name) if ax0 is not None else max(span / 8.0, 1e-12)
    lo = max(span / 1000.0, cur / 100.0, 1e-12)
    hi = max(span, cur * 100.0, lo * 10.0)
    cur = min(max(cur, lo), hi)
    return lo, hi, cur

def _fiberlyse_slider_value_to_pos(value: float, lo: float, hi: float) -> float:
    value = max(float(value), float(lo))
    lo = max(float(lo), 1e-12)
    hi = max(float(hi), lo * 10.0)
    try:
        llo = np.log10(lo); lhi = np.log10(hi); lv = np.log10(value)
        if not np.isfinite(llo) or not np.isfinite(lhi) or lhi <= llo:
            return 50.0
        return float(100.0 * (lv - llo) / (lhi - llo))
    except Exception:
        return 50.0

def _fiberlyse_slider_pos_to_value(pos: float, lo: float, hi: float) -> float:
    lo = max(float(lo), 1e-12)
    hi = max(float(hi), lo * 10.0)
    try:
        p = min(max(float(pos), 0.0), 100.0) / 100.0
        return float(10.0 ** (np.log10(lo) + p * (np.log10(hi) - np.log10(lo))))
    except Exception:
        return float(lo)

def _fiberlyse_plot_apply_axis_tick_intervals(self, reset_auto: bool=False) -> None:
    if _FiberlyseMultipleLocator is None:
        return
    x_interval = _fiberlyse_axis_float_or_none(getattr(self, '_axis_tick_interval_x', None))
    y_interval = _fiberlyse_axis_float_or_none(getattr(self, '_axis_tick_interval_y', None))
    for ax in list(getattr(self.fig, 'axes', [])):
        try:
            if x_interval is not None:
                ax.xaxis.set_major_locator(_FiberlyseMultipleLocator(x_interval))
            elif reset_auto and _FiberlyseAutoLocator is not None:
                ax.xaxis.set_major_locator(_FiberlyseAutoLocator())
        except Exception:
            pass
        try:
            if y_interval is not None:
                ax.yaxis.set_major_locator(_FiberlyseMultipleLocator(y_interval))
            elif reset_auto and _FiberlyseAutoLocator is not None:
                ax.yaxis.set_major_locator(_FiberlyseAutoLocator())
        except Exception:
            pass

def _fiberlyse_plot_set_axis_tick_intervals(self, x_interval=None, y_interval=None) -> None:
    self._axis_tick_interval_x = _fiberlyse_axis_float_or_none(x_interval)
    self._axis_tick_interval_y = _fiberlyse_axis_float_or_none(y_interval)
    try:
        self._apply_axis_tick_intervals(reset_auto=True)
    except Exception:
        pass
    try:
        self.canvas.draw_idle()
    except Exception:
        pass

def _fiberlyse_plot_show_axis_interval_dialog(self):
    parent = self.winfo_toplevel()
    if not getattr(self.fig, 'axes', []):
        try:messagebox.showinfo('Graph customization', 'No axes/plot found in this graph.')
        except Exception:pass
        return
    x_lo, x_hi, x_cur = self._axis_interval_bounds('x')
    y_lo, y_hi, y_cur = self._axis_interval_bounds('y')
    orig_x = getattr(self, '_axis_tick_interval_x', None)
    orig_y = getattr(self, '_axis_tick_interval_y', None)
    x_auto = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_x) is None)
    y_auto = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_y) is None)
    x_entry = tk.StringVar(value=f'{(_fiberlyse_axis_float_or_none(orig_x) or x_cur):.8g}')
    y_entry = tk.StringVar(value=f'{(_fiberlyse_axis_float_or_none(orig_y) or y_cur):.8g}')
    x_scale = tk.DoubleVar(value=_fiberlyse_slider_value_to_pos(_fiberlyse_axis_float_or_none(orig_x) or x_cur, x_lo, x_hi))
    y_scale = tk.DoubleVar(value=_fiberlyse_slider_value_to_pos(_fiberlyse_axis_float_or_none(orig_y) or y_cur, y_lo, y_hi))
    top = tk.Toplevel(parent);top.title('Graph customization');top.transient(parent);top.grab_set();top.resizable(True, False)
    frm = ttk.Frame(top, padding=10);frm.pack(fill=tk.BOTH, expand=True)
    ttk.Label(frm, text='Customize the current graph. Use Auto/Default to return an item to the plot default.', justify='left').grid(row=0, column=0, columnspan=5, sticky='w', pady=(0, 10))
    busy = {'value': False}
    def parse_entry(s, label):
        try:v = float((s or '').strip())
        except Exception:raise ValueError(f'{label} interval must be a number.')
        if not np.isfinite(v) or v <= 0:raise ValueError(f'{label} interval must be a positive finite number.')
        return float(v)
    def apply_from_controls(show_errors=False):
        try:
            xi = None if bool(x_auto.get()) else parse_entry(x_entry.get(), 'X')
            yi = None if bool(y_auto.get()) else parse_entry(y_entry.get(), 'Y')
        except Exception as e:
            if show_errors:messagebox.showerror('Invalid graph customization', str(e), parent=top)
            return False
        self.set_axis_tick_intervals(xi, yi)
        return True
    def on_scale(axis_name):
        if busy['value']:
            return
        busy['value'] = True
        try:
            if axis_name == 'x':
                val = _fiberlyse_slider_pos_to_value(x_scale.get(), x_lo, x_hi);x_entry.set(f'{val:.8g}');x_auto.set(False)
            else:
                val = _fiberlyse_slider_pos_to_value(y_scale.get(), y_lo, y_hi);y_entry.set(f'{val:.8g}');y_auto.set(False)
        finally:
            busy['value'] = False
        apply_from_controls(show_errors=False)
    def on_entry_return(_event=None):
        apply_from_controls(show_errors=True)
        try:
            xv = parse_entry(x_entry.get(), 'X')
            x_scale.set(_fiberlyse_slider_value_to_pos(xv, x_lo, x_hi))
        except Exception:
            pass
        try:
            yv = parse_entry(y_entry.get(), 'Y')
            y_scale.set(_fiberlyse_slider_value_to_pos(yv, y_lo, y_hi))
        except Exception:
            pass
        return 'break'
    def on_auto_changed():
        apply_from_controls(show_errors=False)
    def row(r, label, lo, hi, var_scale, var_entry, var_auto, axis_name):
        ttk.Label(frm, text=label).grid(row=r, column=0, sticky='w', padx=(0, 8), pady=4)
        scale = ttk.Scale(frm, from_=0.0, to=100.0, orient=tk.HORIZONTAL, variable=var_scale, command=lambda _v, n=axis_name: on_scale(n))
        scale.grid(row=r, column=1, sticky='ew', padx=(0, 8), pady=4)
        ent = ttk.Entry(frm, textvariable=var_entry, width=14)
        ent.grid(row=r, column=2, sticky='w', padx=(0, 8), pady=4)
        ent.bind('<Return>', on_entry_return)
        chk = ttk.Checkbutton(frm, text='Auto', variable=var_auto, command=on_auto_changed)
        chk.grid(row=r, column=3, sticky='w', padx=(0, 8), pady=4)
        ttk.Label(frm, text=f'range: {lo:.3g} to {hi:.3g}', foreground='gray35').grid(row=r, column=4, sticky='w', pady=4)
    frm.columnconfigure(1, weight=1)
    row(1, 'X tick interval:', x_lo, x_hi, x_scale, x_entry, x_auto, 'x')
    row(2, 'Y tick interval:', y_lo, y_hi, y_scale, y_entry, y_auto, 'y')
    hint = ttk.Label(frm, text='Tip: press Enter in a number box to apply an exact interval. The sliders use a logarithmic scale for fine control.', foreground='gray35', justify='left')
    hint.grid(row=3, column=0, columnspan=5, sticky='w', pady=(8, 0))
    btn_row = ttk.Frame(frm);btn_row.grid(row=4, column=0, columnspan=5, sticky='ew', pady=(12, 0))
    def close_dialog():
        try:top.grab_release()
        except Exception:pass
        try:top.destroy()
        except Exception:pass
    def do_reset():
        x_auto.set(True);y_auto.set(True);self.set_axis_tick_intervals(None, None)
    def do_ok(_event=None):
        if apply_from_controls(show_errors=True):close_dialog()
        return 'break'
    def do_cancel(_event=None):
        self.set_axis_tick_intervals(orig_x, orig_y);close_dialog();return 'break'
    ttk.Button(btn_row, text='Reset both to Auto', command=do_reset).pack(side=tk.LEFT)
    ttk.Button(btn_row, text='OK', command=do_ok).pack(side=tk.RIGHT)
    ttk.Button(btn_row, text='Cancel', command=do_cancel).pack(side=tk.RIGHT, padx=(0, 8))
    top.bind('<Return>', do_ok);top.bind('<Escape>', do_cancel)
    parent.wait_window(top)

def _fiberlyse_install_axis_interval_hotkeys(self):
    self.root.bind_all('<Control-k>', self._on_ctrl_k_axis_intervals, add='+')
    self.root.bind_all('<Control-K>', self._on_ctrl_k_axis_intervals, add='+')

def _fiberlyse_on_ctrl_k_axis_intervals(self, _event=None):
    tab = self.find_active_plot_tab()
    if tab is None:
        return
    tab.show_axis_interval_dialog()

PlotTabTk._axis_interval_bounds = _fiberlyse_plot_axis_interval_bounds
PlotTabTk._apply_axis_tick_intervals = _fiberlyse_plot_apply_axis_tick_intervals
PlotTabTk.set_axis_tick_intervals = _fiberlyse_plot_set_axis_tick_intervals
PlotTabTk.show_axis_interval_dialog = _fiberlyse_plot_show_axis_interval_dialog
_FIBERLYSE_ORIGINAL_PLOTTABTK_REDRAW_FOR_AXIS_INTERVALS = PlotTabTk.redraw
def _fiberlyse_plottabtk_redraw_with_axis_intervals(self):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_REDRAW_FOR_AXIS_INTERVALS(self)
    try:self._apply_axis_tick_intervals(reset_auto=False)
    except Exception:pass
    try:self.canvas.draw_idle()
    except Exception:pass
PlotTabTk.redraw = _fiberlyse_plottabtk_redraw_with_axis_intervals
_FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_AXIS_INTERVALS = PlotTabTk.__init__
def _fiberlyse_plottabtk_init_with_axis_intervals(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_AXIS_INTERVALS(self, *args, **kwargs)
    try:
        existing = getattr(self, 'axis_interval_btn', None)
        if existing is not None:
            try:
                if bool(existing.winfo_exists()):return
            except Exception:pass
        btn_parent = self.save_btn.master
        self.axis_interval_btn = ttk.Button(btn_parent, text='Graph customization...', command=self.show_axis_interval_dialog)
        self.axis_interval_btn.pack(side=tk.RIGHT, padx=(0, 8))
    except Exception:pass
PlotTabTk.__init__ = _fiberlyse_plottabtk_init_with_axis_intervals
MainAppTk._install_axis_interval_hotkeys = _fiberlyse_install_axis_interval_hotkeys
MainAppTk._on_ctrl_k_axis_intervals = _fiberlyse_on_ctrl_k_axis_intervals
_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_AXIS_INTERVALS = MainAppTk.__init__
def _fiberlyse_mainapptk_init_with_axis_intervals(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_AXIS_INTERVALS(self, *args, **kwargs)
    try:self._install_axis_interval_hotkeys()
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_mainapptk_init_with_axis_intervals
# ---- End Fiberlyse axis tick interval slider extension ----

# ---- Fiberlyse axis visible range extension ----
# Extends the axis interval dialog so each plot can also set visible x/y ranges.

def _fiberlyse_axis_range_tuple_or_none(value):
    if value is None:
        return None
    try:
        a, b = value
        a = float(a); b = float(b)
    except Exception:
        return None
    if not np.isfinite(a) or not np.isfinite(b) or a == b:
        return None
    if b < a:
        a, b = b, a
    return (float(a), float(b))

def _fiberlyse_plot_current_axis_range(self, axis_name: str):
    try:
        ax = self.fig.axes[0] if self.fig.axes else None
    except Exception:
        ax = None
    if ax is None:
        return (0.0, 1.0)
    try:
        lim = ax.get_xlim() if axis_name == 'x' else ax.get_ylim()
        a = float(lim[0]); b = float(lim[1])
        if np.isfinite(a) and np.isfinite(b) and a != b:
            if b < a:
                a, b = b, a
            return (float(a), float(b))
    except Exception:
        pass
    return (0.0, 1.0)

def _fiberlyse_plot_axis_range_slider_bounds(self, axis_name: str):
    vals = []
    for ax in list(getattr(self.fig, 'axes', [])):
        try:
            cur = ax.get_xlim() if axis_name == 'x' else ax.get_ylim()
            for v in cur:
                vf = float(v)
                if np.isfinite(vf):
                    vals.append(vf)
        except Exception:
            pass
        try:
            dl = getattr(ax, 'dataLim', None)
            if dl is not None:
                interval = dl.intervalx if axis_name == 'x' else dl.intervaly
                for v in interval:
                    vf = float(v)
                    if np.isfinite(vf):
                        vals.append(vf)
        except Exception:
            pass
        try:
            for line in list(getattr(ax, 'lines', [])):
                arr = np.asarray(line.get_xdata() if axis_name == 'x' else line.get_ydata(), dtype=float)
                arr = arr[np.isfinite(arr)]
                if arr.size:
                    vals.append(float(np.nanmin(arr))); vals.append(float(np.nanmax(arr)))
        except Exception:
            pass
        try:
            for coll in list(getattr(ax, 'collections', [])):
                offs = coll.get_offsets()
                if offs is None or len(offs) == 0:
                    continue
                offs = np.asarray(offs, dtype=float)
                if offs.ndim == 2 and offs.shape[1] >= 2:
                    arr = offs[:, 0] if axis_name == 'x' else offs[:, 1]
                    arr = arr[np.isfinite(arr)]
                    if arr.size:
                        vals.append(float(np.nanmin(arr))); vals.append(float(np.nanmax(arr)))
        except Exception:
            pass
    stored = _fiberlyse_axis_range_tuple_or_none(getattr(self, f'_axis_visible_range_{axis_name}', None))
    if stored is not None:
        vals.extend([stored[0], stored[1]])
    vals = [float(v) for v in vals if np.isfinite(v)]
    if not vals:
        vals = [0.0, 1.0]
    lo = float(min(vals)); hi = float(max(vals))
    if not np.isfinite(lo) or not np.isfinite(hi):
        lo, hi = 0.0, 1.0
    if hi <= lo:
        pad = max(abs(lo) * 0.05, 1.0)
        lo -= pad; hi += pad
    span = hi - lo
    pad = max(span * 0.05, 1e-12)
    lo -= pad; hi += pad
    if hi <= lo:
        hi = lo + 1.0
    return (float(lo), float(hi))

def _fiberlyse_plot_apply_axis_ranges(self, reset_auto: bool=False) -> None:
    x_range = _fiberlyse_axis_range_tuple_or_none(getattr(self, '_axis_visible_range_x', None))
    y_range = _fiberlyse_axis_range_tuple_or_none(getattr(self, '_axis_visible_range_y', None))
    for ax in list(getattr(self.fig, 'axes', [])):
        try:
            if x_range is not None:
                ax.set_xlim(x_range[0], x_range[1], auto=False)
            elif reset_auto:
                ax.autoscale(enable=True, axis='x')
                try:ax.autoscale_view(scalex=True, scaley=False)
                except Exception:pass
        except Exception:
            pass
        try:
            if y_range is not None:
                ax.set_ylim(y_range[0], y_range[1], auto=False)
            elif reset_auto:
                ax.autoscale(enable=True, axis='y')
                try:ax.autoscale_view(scalex=False, scaley=True)
                except Exception:pass
        except Exception:
            pass

def _fiberlyse_plot_set_axis_ranges(self, x_range=None, y_range=None) -> None:
    self._axis_visible_range_x = _fiberlyse_axis_range_tuple_or_none(x_range)
    self._axis_visible_range_y = _fiberlyse_axis_range_tuple_or_none(y_range)
    try:
        self._apply_axis_ranges(reset_auto=True)
    except Exception:
        pass
    try:
        self._apply_axis_tick_intervals(reset_auto=False)
    except Exception:
        pass
    try:
        self.canvas.draw_idle()
    except Exception:
        pass

def _fiberlyse_clamp(v, lo, hi):
    try:
        v = float(v); lo = float(lo); hi = float(hi)
        if hi < lo:
            lo, hi = hi, lo
        return float(min(max(v, lo), hi))
    except Exception:
        return v

def _fiberlyse_plot_show_axis_interval_range_dialog(self, initial_tab: str='axes', publish_legend_on_ok: bool=True):
    parent = self.winfo_toplevel()
    if not getattr(self.fig, 'axes', []):
        try:messagebox.showinfo('Graph customization', 'No axes/plot found in this graph.')
        except Exception:pass
        return

    x_lo, x_hi, x_cur = self._axis_interval_bounds('x')
    y_lo, y_hi, y_cur = self._axis_interval_bounds('y')
    xb_lo, xb_hi = self._axis_range_slider_bounds('x')
    yb_lo, yb_hi = self._axis_range_slider_bounds('y')

    orig_tick_x = getattr(self, '_axis_tick_interval_x', None)
    orig_tick_y = getattr(self, '_axis_tick_interval_y', None)
    orig_range_x = getattr(self, '_axis_visible_range_x', None)
    orig_range_y = getattr(self, '_axis_visible_range_y', None)
    orig_axis_label_fs = getattr(self, 'axis_label_fontsize', None)
    orig_graph_title_fs = getattr(self, 'graph_title_fontsize', None)
    orig_tick_label_fs = getattr(self, 'tick_label_fontsize', None)
    try:
        orig_legend_snapshot = _fiberlyse_legend_sync_snapshot(self)
    except Exception:
        orig_legend_snapshot = None

    cur_x_range = _fiberlyse_axis_range_tuple_or_none(orig_range_x) or self._current_axis_range('x')
    cur_y_range = _fiberlyse_axis_range_tuple_or_none(orig_range_y) or self._current_axis_range('y')

    x_auto = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_tick_x) is None)
    y_auto = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_tick_y) is None)
    xr_auto = tk.BooleanVar(value=_fiberlyse_axis_range_tuple_or_none(orig_range_x) is None)
    yr_auto = tk.BooleanVar(value=_fiberlyse_axis_range_tuple_or_none(orig_range_y) is None)
    axis_label_font_default = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_axis_label_fs) is None)
    graph_title_font_default = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_graph_title_fs) is None)
    tick_label_font_default = tk.BooleanVar(value=_fiberlyse_axis_float_or_none(orig_tick_label_fs) is None)

    x_entry = tk.StringVar(value=f'{(_fiberlyse_axis_float_or_none(orig_tick_x) or x_cur):.8g}')
    y_entry = tk.StringVar(value=f'{(_fiberlyse_axis_float_or_none(orig_tick_y) or y_cur):.8g}')
    x_scale = tk.DoubleVar(value=_fiberlyse_slider_value_to_pos(_fiberlyse_axis_float_or_none(orig_tick_x) or x_cur, x_lo, x_hi))
    y_scale = tk.DoubleVar(value=_fiberlyse_slider_value_to_pos(_fiberlyse_axis_float_or_none(orig_tick_y) or y_cur, y_lo, y_hi))

    xr_min_entry = tk.StringVar(value=f'{cur_x_range[0]:.8g}')
    xr_max_entry = tk.StringVar(value=f'{cur_x_range[1]:.8g}')
    yr_min_entry = tk.StringVar(value=f'{cur_y_range[0]:.8g}')
    yr_max_entry = tk.StringVar(value=f'{cur_y_range[1]:.8g}')
    xr_min_scale = tk.DoubleVar(value=_fiberlyse_clamp(cur_x_range[0], xb_lo, xb_hi))
    xr_max_scale = tk.DoubleVar(value=_fiberlyse_clamp(cur_x_range[1], xb_lo, xb_hi))
    yr_min_scale = tk.DoubleVar(value=_fiberlyse_clamp(cur_y_range[0], yb_lo, yb_hi))
    yr_max_scale = tk.DoubleVar(value=_fiberlyse_clamp(cur_y_range[1], yb_lo, yb_hi))

    def _font_entry_value(value, fallback):
        v = _fiberlyse_axis_float_or_none(value)
        if v is None:
            v = float(fallback)
        return f'{v:.8g}'

    axis_label_font_entry = tk.StringVar(value=_font_entry_value(orig_axis_label_fs, DEFAULT_AXIS_LABEL_FONTSIZE))
    graph_title_font_entry = tk.StringVar(value=_font_entry_value(orig_graph_title_fs, DEFAULT_GRAPH_TITLE_FONTSIZE))
    tick_label_font_entry = tk.StringVar(value=_font_entry_value(orig_tick_label_fs, DEFAULT_TICK_LABEL_FONTSIZE))

    top = tk.Toplevel(parent)
    top.title('Graph customization')
    top.transient(parent)
    top.grab_set()
    top.resizable(True, False)
    frm = ttk.Frame(top, padding=10)
    frm.pack(fill=tk.BOTH, expand=True)
    ttk.Label(frm, text='Customize the current graph: axes, fonts, and legend box settings. Use Auto/Default to return an item to the plot default.', justify='left').pack(anchor='w', pady=(0, 10))

    notebook = ttk.Notebook(frm)
    notebook.pack(fill=tk.BOTH, expand=True)
    axes_page = ttk.Frame(notebook, padding=8)
    legend_page = ttk.Frame(notebook, padding=8)
    notebook.add(axes_page, text='Axes / fonts')
    notebook.add(legend_page, text='Legend box')

    tick_box = ttk.LabelFrame(axes_page, text='Tick intervals')
    tick_box.pack(fill=tk.X, expand=True)
    range_box = ttk.LabelFrame(axes_page, text='Visible range')
    range_box.pack(fill=tk.X, expand=True, pady=(10, 0))
    font_box = ttk.LabelFrame(axes_page, text='Fonts')
    font_box.pack(fill=tk.X, expand=True, pady=(10, 0))
    tick_box.columnconfigure(1, weight=1)
    range_box.columnconfigure(2, weight=1)
    font_box.columnconfigure(1, weight=1)
    legend_page.columnconfigure(1, weight=1)

    busy = {'value': False}
    legend_changed = {'value': False}

    def parse_interval_entry(s, label):
        try:v = float((s or '').strip())
        except Exception:raise ValueError(f'{label} tick interval must be a number.')
        if not np.isfinite(v) or v <= 0:raise ValueError(f'{label} tick interval must be a positive finite number.')
        return float(v)

    def parse_range_entries(min_s, max_s, label):
        try:a = float((min_s or '').strip()); b = float((max_s or '').strip())
        except Exception:raise ValueError(f'{label} visible range values must be numbers.')
        if not np.isfinite(a) or not np.isfinite(b):raise ValueError(f'{label} visible range values must be finite numbers.')
        if b <= a:raise ValueError(f'{label} visible range maximum must be greater than the minimum.')
        return (float(a), float(b))

    def parse_font_entry(s, label):
        try:v = float((s or '').strip())
        except Exception:raise ValueError(f'{label} must be a number.')
        if not np.isfinite(v) or v <= 0:raise ValueError(f'{label} must be a positive finite number.')
        return float(v)

    def apply_from_controls(show_errors=False):
        try:
            xi = None if bool(x_auto.get()) else parse_interval_entry(x_entry.get(), 'X')
            yi = None if bool(y_auto.get()) else parse_interval_entry(y_entry.get(), 'Y')
            xr = None if bool(xr_auto.get()) else parse_range_entries(xr_min_entry.get(), xr_max_entry.get(), 'X')
            yr = None if bool(yr_auto.get()) else parse_range_entries(yr_min_entry.get(), yr_max_entry.get(), 'Y')
            axis_fs = None if bool(axis_label_font_default.get()) else parse_font_entry(axis_label_font_entry.get(), 'Axis label font size')
            title_fs = None if bool(graph_title_font_default.get()) else parse_font_entry(graph_title_font_entry.get(), 'Graph title font size')
            tick_fs = None if bool(tick_label_font_default.get()) else parse_font_entry(tick_label_font_entry.get(), 'Tick label font size')
        except Exception as e:
            if show_errors:messagebox.showerror('Invalid graph customization', str(e), parent=top)
            return False
        self.set_axis_tick_intervals(xi, yi)
        self.set_axis_ranges(xr, yr)
        self.set_axis_label_fontsize(axis_fs)
        self.set_graph_title_fontsize(title_fs)
        self.set_tick_label_fontsize(tick_fs)
        return True

    def update_interval_slider_from_entry(axis_name):
        try:
            if axis_name == 'x':
                xv = parse_interval_entry(x_entry.get(), 'X')
                x_scale.set(_fiberlyse_slider_value_to_pos(xv, x_lo, x_hi))
            else:
                yv = parse_interval_entry(y_entry.get(), 'Y')
                y_scale.set(_fiberlyse_slider_value_to_pos(yv, y_lo, y_hi))
        except Exception:
            pass

    def update_range_sliders_from_entries(axis_name):
        try:
            if axis_name == 'x':
                a, b = parse_range_entries(xr_min_entry.get(), xr_max_entry.get(), 'X')
                xr_min_scale.set(_fiberlyse_clamp(a, xb_lo, xb_hi)); xr_max_scale.set(_fiberlyse_clamp(b, xb_lo, xb_hi))
            else:
                a, b = parse_range_entries(yr_min_entry.get(), yr_max_entry.get(), 'Y')
                yr_min_scale.set(_fiberlyse_clamp(a, yb_lo, yb_hi)); yr_max_scale.set(_fiberlyse_clamp(b, yb_lo, yb_hi))
        except Exception:
            pass

    def on_interval_scale(axis_name):
        if busy['value']:
            return
        busy['value'] = True
        try:
            if axis_name == 'x':
                val = _fiberlyse_slider_pos_to_value(x_scale.get(), x_lo, x_hi); x_entry.set(f'{val:.8g}'); x_auto.set(False)
            else:
                val = _fiberlyse_slider_pos_to_value(y_scale.get(), y_lo, y_hi); y_entry.set(f'{val:.8g}'); y_auto.set(False)
        finally:
            busy['value'] = False
        apply_from_controls(show_errors=False)

    def on_range_scale(axis_name, which):
        if busy['value']:
            return
        busy['value'] = True
        try:
            if axis_name == 'x':
                b_lo, b_hi = xb_lo, xb_hi; eps = max((b_hi - b_lo) * 0.001, 1e-12)
                a = float(xr_min_scale.get()); b = float(xr_max_scale.get())
                if which == 'min' and a >= b:
                    b = min(b_hi, a + eps); xr_max_scale.set(b)
                    if b <= a:
                        a = max(b_lo, b - eps); xr_min_scale.set(a)
                elif which == 'max' and b <= a:
                    a = max(b_lo, b - eps); xr_min_scale.set(a)
                    if b <= a:
                        b = min(b_hi, a + eps); xr_max_scale.set(b)
                xr_min_entry.set(f'{a:.8g}'); xr_max_entry.set(f'{b:.8g}'); xr_auto.set(False)
            else:
                b_lo, b_hi = yb_lo, yb_hi; eps = max((b_hi - b_lo) * 0.001, 1e-12)
                a = float(yr_min_scale.get()); b = float(yr_max_scale.get())
                if which == 'min' and a >= b:
                    b = min(b_hi, a + eps); yr_max_scale.set(b)
                    if b <= a:
                        a = max(b_lo, b - eps); yr_min_scale.set(a)
                elif which == 'max' and b <= a:
                    a = max(b_lo, b - eps); yr_min_scale.set(a)
                    if b <= a:
                        b = min(b_hi, a + eps); yr_max_scale.set(b)
                yr_min_entry.set(f'{a:.8g}'); yr_max_entry.set(f'{b:.8g}'); yr_auto.set(False)
        finally:
            busy['value'] = False
        apply_from_controls(show_errors=False)

    def on_font_spin(default_var):
        try:default_var.set(False)
        except Exception:pass
        apply_from_controls(show_errors=False)

    def on_entry_return(_event=None):
        ok = apply_from_controls(show_errors=True)
        if ok:
            update_interval_slider_from_entry('x'); update_interval_slider_from_entry('y')
            update_range_sliders_from_entries('x'); update_range_sliders_from_entries('y')
        return 'break'

    def on_auto_changed():
        apply_from_controls(show_errors=False)

    def tick_row(r, label, lo, hi, var_scale, var_entry, var_auto, axis_name):
        ttk.Label(tick_box, text=label).grid(row=r, column=0, sticky='w', padx=(8, 8), pady=4)
        scale = ttk.Scale(tick_box, from_=0.0, to=100.0, orient=tk.HORIZONTAL, variable=var_scale, command=lambda _v, n=axis_name: on_interval_scale(n))
        scale.grid(row=r, column=1, sticky='ew', padx=(0, 8), pady=4)
        ent = ttk.Entry(tick_box, textvariable=var_entry, width=14)
        ent.grid(row=r, column=2, sticky='w', padx=(0, 8), pady=4)
        ent.bind('<Return>', on_entry_return)
        ttk.Checkbutton(tick_box, text='Auto', variable=var_auto, command=on_auto_changed).grid(row=r, column=3, sticky='w', padx=(0, 8), pady=4)
        ttk.Label(tick_box, text=f'slider: {lo:.3g} to {hi:.3g}', foreground='gray35').grid(row=r, column=4, sticky='w', padx=(0, 8), pady=4)

    def range_row(r, label, b_lo, b_hi, min_scale, max_scale, min_entry, max_entry, auto_var, axis_name):
        ttk.Label(range_box, text=label).grid(row=r, column=0, rowspan=2, sticky='w', padx=(8, 8), pady=4)
        ttk.Label(range_box, text='Min:').grid(row=r, column=1, sticky='e', padx=(0, 6), pady=2)
        smin = ttk.Scale(range_box, from_=b_lo, to=b_hi, orient=tk.HORIZONTAL, variable=min_scale, command=lambda _v, n=axis_name: on_range_scale(n, 'min'))
        smin.grid(row=r, column=2, sticky='ew', padx=(0, 8), pady=2)
        e_min = ttk.Entry(range_box, textvariable=min_entry, width=14)
        e_min.grid(row=r, column=3, sticky='w', padx=(0, 8), pady=2)
        e_min.bind('<Return>', on_entry_return)
        ttk.Checkbutton(range_box, text='Auto', variable=auto_var, command=on_auto_changed).grid(row=r, column=4, rowspan=2, sticky='w', padx=(0, 8), pady=2)
        ttk.Label(range_box, text=f'slider: {b_lo:.3g} to {b_hi:.3g}', foreground='gray35').grid(row=r, column=5, rowspan=2, sticky='w', padx=(0, 8), pady=2)
        ttk.Label(range_box, text='Max:').grid(row=r + 1, column=1, sticky='e', padx=(0, 6), pady=2)
        smax = ttk.Scale(range_box, from_=b_lo, to=b_hi, orient=tk.HORIZONTAL, variable=max_scale, command=lambda _v, n=axis_name: on_range_scale(n, 'max'))
        smax.grid(row=r + 1, column=2, sticky='ew', padx=(0, 8), pady=2)
        e_max = ttk.Entry(range_box, textvariable=max_entry, width=14)
        e_max.grid(row=r + 1, column=3, sticky='w', padx=(0, 8), pady=2)
        e_max.bind('<Return>', on_entry_return)

    def font_row(r, label, var_entry, default_var, min_size, max_size):
        ttk.Label(font_box, text=label).grid(row=r, column=0, sticky='w', padx=(8, 8), pady=4)
        spn = ttk.Spinbox(font_box, from_=min_size, to=max_size, increment=1, textvariable=var_entry, width=10, command=lambda dv=default_var: on_font_spin(dv))
        spn.grid(row=r, column=1, sticky='w', padx=(0, 8), pady=4)
        def on_font_return(_event=None, dv=default_var):
            try:dv.set(False)
            except Exception:pass
            return on_entry_return(_event)
        spn.bind('<Return>', on_font_return)
        spn.bind('<KeyRelease>', lambda _event, dv=default_var: dv.set(False))
        ttk.Checkbutton(font_box, text='Default', variable=default_var, command=on_auto_changed).grid(row=r, column=2, sticky='w', padx=(0, 8), pady=4)

    tick_row(0, 'X tick interval:', x_lo, x_hi, x_scale, x_entry, x_auto, 'x')
    tick_row(1, 'Y tick interval:', y_lo, y_hi, y_scale, y_entry, y_auto, 'y')
    range_row(0, 'X visible range:', xb_lo, xb_hi, xr_min_scale, xr_max_scale, xr_min_entry, xr_max_entry, xr_auto, 'x')
    range_row(2, 'Y visible range:', yb_lo, yb_hi, yr_min_scale, yr_max_scale, yr_min_entry, yr_max_entry, yr_auto, 'y')
    font_row(0, 'Axis label font:', axis_label_font_entry, axis_label_font_default, 6, 60)
    font_row(1, 'Graph title font:', graph_title_font_entry, graph_title_font_default, 6, 80)
    font_row(2, 'Tick label font:', tick_label_font_entry, tick_label_font_default, 4, 40)

    ttk.Label(axes_page, text='Tip: for exact limits or font sizes, type numbers into the boxes and press Enter. Ctrl+K opens Graph customization for the active plot.', foreground='gray35', justify='left').pack(anchor='w', pady=(8, 0))

    legend_choices = []
    try:
        legend_choices = self._legend_axes_choices()
    except Exception:
        legend_choices = []

    if not legend_choices:
        ttk.Label(legend_page, text='No legend is visible on the current graph. Create or show a legend on the plot, then open Graph customization again.', justify='left', foreground='gray35', wraplength=560).grid(row=0, column=0, sticky='w')
    else:
        ttk.Label(legend_page, text='Move and resize the legend for the current plot. You can also drag legends directly on the graph.', justify='left').grid(row=0, column=0, columnspan=4, sticky='w', pady=(0, 8))
        legend_labels = [c[0] for c in legend_choices]
        legend_axis_var = tk.StringVar(value=legend_labels[0])
        ttk.Label(legend_page, text='Legend:').grid(row=1, column=0, sticky='w', padx=(0, 8), pady=4)
        legend_cmb = ttk.Combobox(legend_page, values=legend_labels, textvariable=legend_axis_var, state='readonly', width=42)
        legend_cmb.grid(row=1, column=1, columnspan=3, sticky='ew', pady=4)
        legend_text_size = tk.DoubleVar(value=10.0)
        legend_title_size = tk.DoubleVar(value=10.0)
        legend_alpha_var = tk.DoubleVar(value=0.8)
        legend_auto_pos = tk.BooleanVar(value=False)
        legend_x_var = tk.DoubleVar(value=0.5)
        legend_y_var = tk.DoubleVar(value=0.5)
        legend_busy = {'value': False}

        def legend_current_choice():
            lab = legend_axis_var.get()
            for item in legend_choices:
                if item[0] == lab:
                    return item
            return legend_choices[0]

        def legend_load_current():
            legend_busy['value'] = True
            try:
                _label, ax, leg = legend_current_choice()
                state = _fiberlyse_plot_get_legend_state(self, ax)
                fs = _fiberlyse_legend_float(state.get('fontsize', None), None)
                if fs is None:
                    fs = _fiberlyse_legend_text_fontsize(leg, 10.0)
                tfs = _fiberlyse_legend_float(state.get('title_fontsize', None), None)
                if tfs is None:
                    tfs = _fiberlyse_legend_title_fontsize(leg, fs)
                alpha = _fiberlyse_legend_float(state.get('frame_alpha', None), None)
                if alpha is None:
                    try:alpha = float(leg.get_frame().get_alpha())
                    except Exception:alpha = 0.8
                    if alpha is None or not np.isfinite(alpha):alpha = 0.8
                anchor = state.get('anchor_axes', None)
                if anchor is None:
                    anchor = _fiberlyse_legend_current_anchor_axes(ax, leg)
                legend_text_size.set(_fiberlyse_legend_clamp(fs, 5.0, 48.0))
                legend_title_size.set(_fiberlyse_legend_clamp(tfs, 5.0, 56.0))
                legend_alpha_var.set(_fiberlyse_legend_clamp(alpha, 0.0, 1.0))
                legend_auto_pos.set(not bool(state.get('manual_position', False)))
                legend_x_var.set(_fiberlyse_legend_clamp(anchor[0], -0.75, 1.75))
                legend_y_var.set(_fiberlyse_legend_clamp(anchor[1], -0.75, 1.75))
            finally:
                legend_busy['value'] = False

        def legend_apply_live(show_errors=False):
            if legend_busy['value']:
                return True
            try:
                _label, ax, _leg = legend_current_choice()
                state = _fiberlyse_plot_get_legend_state(self, ax)
                state['fontsize'] = _fiberlyse_legend_clamp(legend_text_size.get(), 5.0, 48.0)
                state['title_fontsize'] = _fiberlyse_legend_clamp(legend_title_size.get(), 5.0, 56.0)
                state['frame_alpha'] = _fiberlyse_legend_clamp(legend_alpha_var.get(), 0.0, 1.0)
                if bool(legend_auto_pos.get()):
                    state['manual_position'] = False
                else:
                    state['manual_position'] = True
                    state['anchor_axes'] = (_fiberlyse_legend_clamp(legend_x_var.get(), -0.75, 1.75), _fiberlyse_legend_clamp(legend_y_var.get(), -0.75, 1.75))
                self._apply_legend_box_overrides()
                self.canvas.draw_idle()
                legend_changed['value'] = True
                return True
            except Exception as e:
                if show_errors:
                    messagebox.showerror('Legend box', f'Could not apply legend settings:\n\n{e}', parent=top)
                return False

        def legend_slider_row(row, label, var, lo, hi, fmt, callback):
            ttk.Label(legend_page, text=label).grid(row=row, column=0, sticky='w', padx=(0, 8), pady=4)
            scale = ttk.Scale(legend_page, from_=lo, to=hi, orient=tk.HORIZONTAL, variable=var, command=lambda _v: callback())
            scale.grid(row=row, column=1, sticky='ew', padx=(0, 8), pady=4)
            entry_var = tk.StringVar(value=fmt(var.get()))
            entry = ttk.Entry(legend_page, textvariable=entry_var, width=10)
            entry.grid(row=row, column=2, sticky='w', padx=(0, 8), pady=4)
            def sync_entry(*_args):
                try:entry_var.set(fmt(var.get()))
                except Exception:pass
            def entry_return(_e=None):
                try:var.set(float((entry_var.get() or '').strip()))
                except Exception:
                    messagebox.showerror('Legend box', f'Could not read {label} as a number.', parent=top)
                    sync_entry();return 'break'
                callback();sync_entry();return 'break'
            try:var.trace_add('write', sync_entry)
            except Exception:pass
            entry.bind('<Return>', entry_return)
            return scale, entry

        legend_slider_row(2, 'Text / box size:', legend_text_size, 5.0, 48.0, lambda v: f'{float(v):.1f}', lambda: legend_apply_live(False))
        legend_slider_row(3, 'Title size:', legend_title_size, 5.0, 56.0, lambda v: f'{float(v):.1f}', lambda: legend_apply_live(False))
        legend_slider_row(4, 'Frame opacity:', legend_alpha_var, 0.0, 1.0, lambda v: f'{float(v):.2f}', lambda: legend_apply_live(False))
        ttk.Checkbutton(legend_page, text='Automatic position', variable=legend_auto_pos, command=lambda: legend_apply_live(False)).grid(row=5, column=0, columnspan=4, sticky='w', pady=(6, 2))
        legend_slider_row(6, 'Position X:', legend_x_var, -0.75, 1.75, lambda v: f'{float(v):.3f}', lambda: (legend_auto_pos.set(False), legend_apply_live(False)))
        legend_slider_row(7, 'Position Y:', legend_y_var, -0.75, 1.75, lambda v: f'{float(v):.3f}', lambda: (legend_auto_pos.set(False), legend_apply_live(False)))
        ttk.Label(legend_page, text='Mouse: left-drag the legend to move. Drag a legend corner, or Shift/Ctrl-drag inside it, to resize.', foreground='gray35', justify='left').grid(row=8, column=0, columnspan=4, sticky='w', pady=(8, 0))
        legend_btn_row = ttk.Frame(legend_page)
        legend_btn_row.grid(row=9, column=0, columnspan=4, sticky='ew', pady=(12, 0))

        def legend_on_axis_changed(_event=None):
            legend_load_current()
        legend_cmb.bind('<<ComboboxSelected>>', legend_on_axis_changed)

        def legend_reset_selected():
            _label, ax, _leg = legend_current_choice()
            try:
                key = self._ax_key(ax)
            except Exception:
                key = id(ax)
            states = getattr(self, '_legend_box_overrides', {})
            if isinstance(states, dict) and key in states:
                states.pop(key, None)
            try:
                leg = ax.get_legend()
                if leg is not None:
                    try:
                        if hasattr(leg, 'set_loc'):
                            leg.set_loc('best')
                        else:
                            leg._loc = 0
                    except Exception:
                        try:leg._loc = 0
                        except Exception:pass
                    try:leg.set_bbox_to_anchor(None)
                    except Exception:pass
            except Exception:
                pass
            legend_changed['value'] = True
            self.redraw()
            legend_load_current()

        ttk.Button(legend_btn_row, text='Reset selected legend', command=legend_reset_selected).pack(side=tk.LEFT)
        legend_load_current()

    try:
        if str(initial_tab or '').strip().lower().startswith('legend'):
            notebook.select(legend_page)
    except Exception:
        pass

    btn_row = ttk.Frame(frm)
    btn_row.pack(fill=tk.X, pady=(12, 0))

    def close_dialog():
        try:top.grab_release()
        except Exception:pass
        try:top.destroy()
        except Exception:pass

    def do_reset_axes_fonts():
        x_auto.set(True); y_auto.set(True); xr_auto.set(True); yr_auto.set(True)
        axis_label_font_default.set(True); graph_title_font_default.set(True); tick_label_font_default.set(True)
        axis_label_font_entry.set(f'{DEFAULT_AXIS_LABEL_FONTSIZE:g}')
        graph_title_font_entry.set(f'{DEFAULT_GRAPH_TITLE_FONTSIZE:g}')
        tick_label_font_entry.set(f'{DEFAULT_TICK_LABEL_FONTSIZE:g}')
        self.set_axis_tick_intervals(None, None)
        self.set_axis_ranges(None, None)
        self.set_axis_label_fontsize(None)
        self.set_graph_title_fontsize(None)
        self.set_tick_label_fontsize(None)

    def do_ok(_event=None):
        if apply_from_controls(show_errors=True):
            if bool(publish_legend_on_ok) and bool(legend_changed.get('value')):
                try:_fiberlyse_legend_sync_publish(self, reason='graph_customization')
                except Exception:pass
            close_dialog()
        return 'break'

    def do_cancel(_event=None):
        self.set_axis_tick_intervals(orig_tick_x, orig_tick_y)
        self.set_axis_ranges(orig_range_x, orig_range_y)
        self.set_axis_label_fontsize(orig_axis_label_fs)
        self.set_graph_title_fontsize(orig_graph_title_fs)
        self.set_tick_label_fontsize(orig_tick_label_fs)
        if orig_legend_snapshot is not None:
            try:_fiberlyse_legend_sync_apply_to_plot(self, orig_legend_snapshot, redraw=True)
            except Exception:pass
        close_dialog()
        return 'break'

    ttk.Button(btn_row, text='Reset axes/fonts to Auto/default', command=do_reset_axes_fonts).pack(side=tk.LEFT)
    ttk.Button(btn_row, text='OK', command=do_ok).pack(side=tk.RIGHT)
    ttk.Button(btn_row, text='Cancel', command=do_cancel).pack(side=tk.RIGHT, padx=(0, 8))
    top.bind('<Return>', do_ok)
    top.bind('<Escape>', do_cancel)
    try:top.protocol('WM_DELETE_WINDOW', do_cancel)
    except Exception:pass
    parent.wait_window(top)

PlotTabTk._current_axis_range = _fiberlyse_plot_current_axis_range
PlotTabTk._axis_range_slider_bounds = _fiberlyse_plot_axis_range_slider_bounds
PlotTabTk._apply_axis_ranges = _fiberlyse_plot_apply_axis_ranges
PlotTabTk.set_axis_ranges = _fiberlyse_plot_set_axis_ranges
PlotTabTk.show_axis_interval_dialog = _fiberlyse_plot_show_axis_interval_range_dialog

_FIBERLYSE_ORIGINAL_PLOTTABTK_REDRAW_FOR_AXIS_RANGES = PlotTabTk.redraw
def _fiberlyse_plottabtk_redraw_with_axis_ranges(self):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_REDRAW_FOR_AXIS_RANGES(self)
    try:self._apply_axis_ranges(reset_auto=False)
    except Exception:pass
    try:self._apply_axis_tick_intervals(reset_auto=False)
    except Exception:pass
    try:self.canvas.draw_idle()
    except Exception:pass
PlotTabTk.redraw = _fiberlyse_plottabtk_redraw_with_axis_ranges

_FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_AXIS_RANGES = PlotTabTk.__init__
def _fiberlyse_plottabtk_init_with_axis_ranges(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_AXIS_RANGES(self, *args, **kwargs)
    try:
        btn = getattr(self, 'axis_interval_btn', None)
        if btn is not None and bool(btn.winfo_exists()):
            btn.config(text='Graph customization...', command=self.show_axis_interval_dialog)
    except Exception:
        pass
PlotTabTk.__init__ = _fiberlyse_plottabtk_init_with_axis_ranges
# ---- End Fiberlyse axis visible range extension ----



# ---- Fiberlyse draggable and resizable legend extension ----
# Adds per-plot legend movement and readability controls.
# Mouse controls:
#   * Left-drag inside a legend to move it.
#   * Drag a legend corner, or Shift/Ctrl-left-drag inside it, to resize the legend text/box.
#   * Use the "Graph customization..." button or Ctrl+L for exact legend settings.

def _fiberlyse_legend_float(value, default=None):
    try:
        v = float(value)
    except Exception:
        return default
    if not np.isfinite(v):
        return default
    return float(v)

def _fiberlyse_legend_clamp(value, lo, hi):
    v = _fiberlyse_legend_float(value, lo)
    try:
        return float(min(max(v, float(lo)), float(hi)))
    except Exception:
        return float(lo)

def _fiberlyse_legend_event_has_resize_modifier(event) -> bool:
    try:
        key = str(getattr(event, 'key', '') or '').lower()
    except Exception:
        key = ''
    return ('shift' in key) or ('control' in key) or ('ctrl' in key)

def _fiberlyse_legend_text_fontsize(leg, fallback: float=10.0) -> float:
    try:
        texts = list(leg.get_texts())
        if texts:
            return float(texts[0].get_fontsize())
    except Exception:
        pass
    try:
        return float(matplotlib.rcParams.get('legend.fontsize', fallback))
    except Exception:
        return float(fallback)

def _fiberlyse_legend_title_fontsize(leg, fallback: Optional[float]=None) -> float:
    if fallback is None:
        fallback = _fiberlyse_legend_text_fontsize(leg, 10.0)
    try:
        title = leg.get_title()
        if title is not None:
            fs = float(title.get_fontsize())
            if np.isfinite(fs) and fs > 0:
                return fs
    except Exception:
        pass
    return float(fallback)

def _fiberlyse_legend_current_anchor_axes(ax, leg, renderer=None):
    try:
        if renderer is None:
            renderer = leg.figure.canvas.get_renderer()
    except Exception:
        renderer = None
    try:
        bbox = leg.get_window_extent(renderer)
        cx = float((bbox.x0 + bbox.x1) / 2.0)
        cy = float((bbox.y0 + bbox.y1) / 2.0)
        xy = ax.transAxes.inverted().transform((cx, cy))
        return (float(xy[0]), float(xy[1]))
    except Exception:
        return (0.5, 0.5)

def _fiberlyse_legend_contains_event(self, event):
    if event is None or getattr(event, 'x', None) is None or getattr(event, 'y', None) is None:
        return None
    renderer = None
    try:
        renderer = self.canvas.get_renderer()
    except Exception:
        try:
            self.canvas.draw()
            renderer = self.canvas.get_renderer()
        except Exception:
            renderer = None
    axes = list(getattr(self.fig, 'axes', []))
    for ax in reversed(axes):
        try:
            leg = ax.get_legend()
        except Exception:
            leg = None
        if leg is None:
            continue
        try:
            if not leg.get_visible():
                continue
        except Exception:
            pass
        try:
            bbox = leg.get_window_extent(renderer)
            if bbox is not None and bbox.contains(float(event.x), float(event.y)):
                return (ax, leg, bbox, renderer)
        except Exception:
            continue
    return None

def _fiberlyse_legend_near_resize_area(bbox, event, margin: float=16.0) -> bool:
    try:
        x = float(event.x); y = float(event.y)
        mx = float(margin)
        near_right = abs(x - float(bbox.x1)) <= mx
        near_left = abs(x - float(bbox.x0)) <= mx
        near_top = abs(y - float(bbox.y1)) <= mx
        near_bottom = abs(y - float(bbox.y0)) <= mx
        return (near_right or near_left) and (near_top or near_bottom)
    except Exception:
        return False

def _fiberlyse_plot_get_legend_state(self, ax=None):
    states = getattr(self, '_legend_box_overrides', None)
    if not isinstance(states, dict):
        states = {}
        self._legend_box_overrides = states
    if ax is None:
        return states
    try:
        key = self._ax_key(ax)
    except Exception:
        key = id(ax)
    state = states.setdefault(key, {})
    if not isinstance(state, dict):
        state = {}
        states[key] = state
    return state

def _fiberlyse_plot_set_legend_cursor(self, cursor: str) -> None:
    try:
        self.canvas_widget.configure(cursor=cursor)
    except Exception:
        try:
            self.canvas_widget.configure(cursor='')
        except Exception:
            pass

def _fiberlyse_plot_apply_legend_box_overrides(self) -> None:
    states = getattr(self, '_legend_box_overrides', None)
    if not isinstance(states, dict):
        states = {}
        self._legend_box_overrides = states
    for ax in list(getattr(self.fig, 'axes', [])):
        try:
            leg = ax.get_legend()
        except Exception:
            leg = None
        if leg is None:
            continue
        try:
            key = self._ax_key(ax)
        except Exception:
            key = id(ax)
        state = states.get(key, {}) if isinstance(states.get(key, {}), dict) else {}
        if bool(state.get('manual_position')):
            try:
                anchor = state.get('anchor_axes', None)
                if anchor is not None:
                    x, y = float(anchor[0]), float(anchor[1])
                    try:
                        if hasattr(leg, 'set_loc'):
                            leg.set_loc('center')
                        else:
                            leg._loc = 10
                    except Exception:
                        try:leg._loc = 10
                        except Exception:pass
                    leg.set_bbox_to_anchor((x, y), transform=ax.transAxes)
            except Exception:
                pass
        elif 'manual_position' in state:
            try:
                if hasattr(leg, 'set_loc'):
                    leg.set_loc('best')
                else:
                    leg._loc = 0
            except Exception:
                try:leg._loc = 0
                except Exception:pass
            try:
                leg.set_bbox_to_anchor(None)
            except Exception:
                pass
        fs = _fiberlyse_legend_float(state.get('fontsize', None), None)
        if fs is not None and fs > 0:
            try:
                for txt in leg.get_texts():
                    txt.set_fontsize(fs)
            except Exception:
                pass
        tfs = _fiberlyse_legend_float(state.get('title_fontsize', None), None)
        if tfs is not None and tfs > 0:
            try:
                title = leg.get_title()
                if title is not None:
                    title.set_fontsize(tfs)
            except Exception:
                pass
        alpha = _fiberlyse_legend_float(state.get('frame_alpha', None), None)
        if alpha is not None:
            try:
                leg.get_frame().set_alpha(_fiberlyse_legend_clamp(alpha, 0.0, 1.0))
            except Exception:
                pass
        face = state.get('facecolor', None)
        if face:
            try:
                leg.get_frame().set_facecolor(face)
            except Exception:
                pass
        try:
            leg.set_draggable(False)
        except Exception:
            pass

def _fiberlyse_plot_on_legend_button_press(self, event):
    if event is None:
        return
    try:
        if self._toolbar_is_active():
            return
    except Exception:
        pass
    try:
        if getattr(event, 'button', None) != 1:
            return
        if bool(getattr(event, 'dblclick', False)):
            return
    except Exception:
        return
    hit = _fiberlyse_legend_contains_event(self, event)
    if hit is None:
        return
    ax, leg, bbox, renderer = hit
    action = 'resize' if (_fiberlyse_legend_event_has_resize_modifier(event) or _fiberlyse_legend_near_resize_area(bbox, event)) else 'move'
    try:
        anchor = _fiberlyse_legend_current_anchor_axes(ax, leg, renderer)
        pointer_axes = ax.transAxes.inverted().transform((float(event.x), float(event.y)))
        offset = (float(anchor[0] - pointer_axes[0]), float(anchor[1] - pointer_axes[1]))
    except Exception:
        anchor = (0.5, 0.5)
        offset = (0.0, 0.0)
    fs0 = _fiberlyse_legend_text_fontsize(leg, 10.0)
    tfs0 = _fiberlyse_legend_title_fontsize(leg, fs0)
    try:
        cx = float((bbox.x0 + bbox.x1) / 2.0); cy = float((bbox.y0 + bbox.y1) / 2.0)
        d0 = float(max(8.0, ((float(event.x) - cx) ** 2.0 + (float(event.y) - cy) ** 2.0) ** 0.5))
    except Exception:
        cx = cy = 0.0; d0 = 80.0
    self._legend_drag_state = {
        'action': action,
        'ax': ax,
        'legend': leg,
        'start_x': float(event.x),
        'start_y': float(event.y),
        'anchor0': anchor,
        'offset': offset,
        'fontsize0': fs0,
        'title_fontsize0': tfs0,
        'center_px': (cx, cy),
        'dist0': d0,
    }
    _fiberlyse_plot_get_legend_state(self, ax)['manual_position'] = True
    _fiberlyse_plot_get_legend_state(self, ax)['anchor_axes'] = tuple(anchor)
    _fiberlyse_plot_set_legend_cursor(self, 'sizing' if action == 'resize' else 'fleur')

def _fiberlyse_plot_on_legend_motion(self, event):
    state = getattr(self, '_legend_drag_state', None)
    if not state:
        hit = _fiberlyse_legend_contains_event(self, event)
        if hit is None:
            _fiberlyse_plot_set_legend_cursor(self, '')
        else:
            _ax, _leg, bbox, _renderer = hit
            cur = 'sizing' if (_fiberlyse_legend_event_has_resize_modifier(event) or _fiberlyse_legend_near_resize_area(bbox, event)) else 'fleur'
            _fiberlyse_plot_set_legend_cursor(self, cur)
        return
    if event is None or getattr(event, 'x', None) is None or getattr(event, 'y', None) is None:
        return
    ax = state.get('ax')
    if ax is None:
        return
    action = str(state.get('action', 'move'))
    st = _fiberlyse_plot_get_legend_state(self, ax)
    if action == 'resize':
        try:
            cx, cy = state.get('center_px', (float(event.x), float(event.y)))
            d0 = float(state.get('dist0', 80.0))
            d = float(max(2.0, ((float(event.x) - float(cx)) ** 2.0 + (float(event.y) - float(cy)) ** 2.0) ** 0.5))
            scale = _fiberlyse_legend_clamp(d / max(d0, 1e-9), 0.35, 3.25)
        except Exception:
            dx = float(event.x) - float(state.get('start_x', event.x))
            scale = _fiberlyse_legend_clamp(1.0 + dx / 220.0, 0.35, 3.25)
        fs = _fiberlyse_legend_clamp(float(state.get('fontsize0', 10.0)) * scale, 5.0, 48.0)
        tfs = _fiberlyse_legend_clamp(float(state.get('title_fontsize0', fs)) * scale, 5.0, 56.0)
        st['fontsize'] = fs
        st['title_fontsize'] = tfs
        st['manual_position'] = True
        st['anchor_axes'] = tuple(state.get('anchor0', (0.5, 0.5)))
    else:
        try:
            pointer_axes = ax.transAxes.inverted().transform((float(event.x), float(event.y)))
            off = state.get('offset', (0.0, 0.0))
            x = float(pointer_axes[0]) + float(off[0])
            y = float(pointer_axes[1]) + float(off[1])
        except Exception:
            x, y = state.get('anchor0', (0.5, 0.5))
        x = _fiberlyse_legend_clamp(x, -0.75, 1.75)
        y = _fiberlyse_legend_clamp(y, -0.75, 1.75)
        st['manual_position'] = True
        st['anchor_axes'] = (x, y)
    try:
        self._apply_legend_box_overrides()
    except Exception:
        pass
    try:
        self.canvas.draw_idle()
    except Exception:
        pass

def _fiberlyse_plot_on_legend_button_release(self, event):
    if getattr(self, '_legend_drag_state', None):
        self._legend_drag_state = None
        _fiberlyse_plot_set_legend_cursor(self, '')
        try:
            self._apply_legend_box_overrides()
            self.canvas.draw_idle()
        except Exception:
            pass

def _fiberlyse_plot_legend_axes_choices(self):
    choices = []
    for i, ax in enumerate(list(getattr(self.fig, 'axes', [])), start=1):
        try:
            leg = ax.get_legend()
        except Exception:
            leg = None
        if leg is None:
            continue
        title = ''
        try:title = str(ax.get_title() or '').strip()
        except Exception:title = ''
        label = f'Axis {i}' + (f': {title}' if title else '')
        choices.append((label, ax, leg))
    return choices

def _fiberlyse_plot_show_legend_box_dialog(self):
    return self.show_axis_interval_dialog(initial_tab='legend', publish_legend_on_ok=False)

PlotTabTk._legend_axes_choices = _fiberlyse_plot_legend_axes_choices
PlotTabTk._apply_legend_box_overrides = _fiberlyse_plot_apply_legend_box_overrides
PlotTabTk.show_legend_box_dialog = _fiberlyse_plot_show_legend_box_dialog

_FIBERLYSE_ORIGINAL_APPLY_USER_OVERRIDES_FOR_LEGENDS = PlotTabTk._apply_user_overrides
def _fiberlyse_apply_user_overrides_with_legends(self):
    _FIBERLYSE_ORIGINAL_APPLY_USER_OVERRIDES_FOR_LEGENDS(self)
    try:self._apply_legend_box_overrides()
    except Exception:pass
PlotTabTk._apply_user_overrides = _fiberlyse_apply_user_overrides_with_legends

_FIBERLYSE_ORIGINAL_PLOTTABTK_REDRAW_FOR_LEGENDS = PlotTabTk.redraw
def _fiberlyse_plottabtk_redraw_with_legends(self):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_REDRAW_FOR_LEGENDS(self)
    try:self._apply_legend_box_overrides()
    except Exception:pass
    try:self.canvas.draw_idle()
    except Exception:pass
PlotTabTk.redraw = _fiberlyse_plottabtk_redraw_with_legends

_FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_LEGENDS = PlotTabTk.__init__
def _fiberlyse_plottabtk_init_with_legends(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_LEGENDS(self, *args, **kwargs)
    try:self._legend_box_overrides = {}
    except Exception:pass
    try:self._legend_drag_state = None
    except Exception:pass
    try:self.legend_box_btn = None
    except Exception:pass
    try:self._cid_legend_press = self.canvas.mpl_connect('button_press_event', self._on_legend_button_press)
    except Exception:pass
    try:self._cid_legend_motion = self.canvas.mpl_connect('motion_notify_event', self._on_legend_motion)
    except Exception:pass
    try:self._cid_legend_release = self.canvas.mpl_connect('button_release_event', self._on_legend_button_release)
    except Exception:pass
PlotTabTk.__init__ = _fiberlyse_plottabtk_init_with_legends
PlotTabTk._on_legend_button_press = _fiberlyse_plot_on_legend_button_press
PlotTabTk._on_legend_motion = _fiberlyse_plot_on_legend_motion
PlotTabTk._on_legend_button_release = _fiberlyse_plot_on_legend_button_release


def _fiberlyse_install_legend_hotkeys(self):
    self.root.bind_all('<Control-l>', self._on_ctrl_l_legend_box, add='+')
    self.root.bind_all('<Control-L>', self._on_ctrl_l_legend_box, add='+')

def _fiberlyse_on_ctrl_l_legend_box(self, _event=None):
    tab = self.find_active_plot_tab()
    if tab is None:return
    tab.show_legend_box_dialog()

MainAppTk._install_legend_hotkeys = _fiberlyse_install_legend_hotkeys
MainAppTk._on_ctrl_l_legend_box = _fiberlyse_on_ctrl_l_legend_box
_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_LEGENDS = MainAppTk.__init__
def _fiberlyse_mainapptk_init_with_legends(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_LEGENDS(self, *args, **kwargs)
    try:self._install_legend_hotkeys()
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_mainapptk_init_with_legends
# ---- End Fiberlyse draggable and resizable legend extension ----


# ---- Begin Fiberlyse synced legend settings across analysis files ----
# When a legend box is moved/resized or edited on one plot, copy the same
# legend settings to the same plot type for every loaded file/channel in the
# current analysis. The stored setting is also applied to lazy-loaded tabs.
try:
    import copy as _fiberlyse_legend_sync_copy
except Exception:
    _fiberlyse_legend_sync_copy = None

def _fiberlyse_legend_sync_deepcopy(obj, fallback=None):
    try:
        if _fiberlyse_legend_sync_copy is not None:
            return _fiberlyse_legend_sync_copy.deepcopy(obj)
    except Exception:
        pass
    try:
        if isinstance(obj, dict):
            return {k: _fiberlyse_legend_sync_deepcopy(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_fiberlyse_legend_sync_deepcopy(v) for v in obj]
        if isinstance(obj, tuple):
            return tuple(_fiberlyse_legend_sync_deepcopy(v) for v in obj)
    except Exception:
        pass
    if obj is None:
        return fallback
    return obj

def _fiberlyse_legend_sync_find_app(widget):
    try:
        root = widget.winfo_toplevel()
        app = getattr(root, '_fiberlyse_app', None)
        if app is not None:
            return app
    except Exception:
        pass
    return None

def _fiberlyse_legend_sync_has_legend(plot_tab) -> bool:
    try:
        for ax in list(getattr(plot_tab.fig, 'axes', [])):
            try:
                leg = ax.get_legend()
            except Exception:
                leg = None
            if leg is not None:
                return True
    except Exception:
        pass
    return False

def _fiberlyse_legend_sync_snapshot(plot_tab):
    try:
        plot_tab._apply_legend_box_overrides()
    except Exception:
        pass
    return {
        'legend_box_overrides': _fiberlyse_legend_sync_deepcopy(getattr(plot_tab, '_legend_box_overrides', {}), {}),
        'legend_title_overrides': _fiberlyse_legend_sync_deepcopy(getattr(plot_tab, '_legend_title_overrides', {}), {}),
        'legend_label_overrides': _fiberlyse_legend_sync_deepcopy(getattr(plot_tab, '_legend_label_overrides', {}), {}),
        'color_overrides': _fiberlyse_legend_sync_deepcopy(getattr(plot_tab, '_color_overrides', {}), {}),
    }

def _fiberlyse_legend_sync_apply_to_plot(plot_tab, snapshot, redraw: bool=True) -> None:
    if not isinstance(snapshot, dict):
        return
    try:
        setattr(plot_tab, '_legend_box_overrides', _fiberlyse_legend_sync_deepcopy(snapshot.get('legend_box_overrides', {}), {}))
    except Exception:
        pass
    for attr, key in [
        ('_legend_title_overrides', 'legend_title_overrides'),
        ('_legend_label_overrides', 'legend_label_overrides'),
        ('_color_overrides', 'color_overrides'),
    ]:
        try:
            setattr(plot_tab, attr, _fiberlyse_legend_sync_deepcopy(snapshot.get(key, {}), {}))
        except Exception:
            pass
    if redraw:
        try:
            plot_tab.redraw()
            return
        except Exception:
            pass
        try:
            plot_tab._apply_user_overrides()
        except Exception:
            pass
        try:
            plot_tab._apply_legend_box_overrides()
        except Exception:
            pass
        try:
            plot_tab.canvas.draw_idle()
        except Exception:
            pass

def _fiberlyse_legend_sync_iter_plot_tabs(app):
    seen = set()
    try:
        widgets = list(getattr(app, '_channel_widgets', {}).values())
    except Exception:
        widgets = []
    attr_names = ['tab_raw', 'tab_art', 'tab_fit', 'tab_norm', 'tab_norm_smooth', 'tab_freq']
    for cw in widgets:
        for name in attr_names:
            try:
                tab = getattr(cw, name, None)
            except Exception:
                tab = None
            if tab is None:
                continue
            if id(tab) in seen:
                continue
            seen.add(id(tab))
            yield tab
    try:
        cmpw = getattr(app, 'compare_widget', None)
        if cmpw is not None and getattr(cmpw, 'plot', None) is not None and id(cmpw.plot) not in seen:
            seen.add(id(cmpw.plot));yield cmpw.plot
    except Exception:
        pass
    try:
        avgw = getattr(app, 'average_widget', None)
        if avgw is not None and getattr(avgw, 'plot', None) is not None and id(avgw.plot) not in seen:
            seen.add(id(avgw.plot));yield avgw.plot
    except Exception:
        pass

def _fiberlyse_legend_sync_publish(plot_tab, reason: str='legend') -> None:
    try:
        if bool(getattr(plot_tab, '_fiberlyse_legend_sync_suppress', False)):
            return
    except Exception:
        pass
    app = _fiberlyse_legend_sync_find_app(plot_tab)
    if app is None:
        return
    tab_name = str(getattr(plot_tab, 'tab_name', '') or '')
    if not tab_name:
        return
    snapshot = _fiberlyse_legend_sync_snapshot(plot_tab)
    try:
        styles = getattr(app, '_fiberlyse_synced_legend_styles_by_tab', None)
        if not isinstance(styles, dict):
            styles = {}
            setattr(app, '_fiberlyse_synced_legend_styles_by_tab', styles)
        styles[tab_name] = _fiberlyse_legend_sync_deepcopy(snapshot, {})
    except Exception:
        pass
    changed = 0
    for other in _fiberlyse_legend_sync_iter_plot_tabs(app):
        if other is plot_tab:
            continue
        try:
            if str(getattr(other, 'tab_name', '') or '') != tab_name:
                continue
            setattr(other, '_fiberlyse_legend_sync_suppress', True)
            _fiberlyse_legend_sync_apply_to_plot(other, snapshot, redraw=True)
            changed += 1
        except Exception:
            pass
        finally:
            try:setattr(other, '_fiberlyse_legend_sync_suppress', False)
            except Exception:pass
    try:
        if getattr(app, 'status', None) is not None:
            if changed > 0:
                app.status.set(f"Legend settings synced for '{tab_name}' across loaded files/channels.")
            else:
                app.status.set(f"Legend settings saved for '{tab_name}' and will be used on other files/channels when opened.")
    except Exception:
        pass

def _fiberlyse_legend_sync_apply_saved_style(plot_tab) -> None:
    app = _fiberlyse_legend_sync_find_app(plot_tab)
    if app is None:
        return
    try:
        styles = getattr(app, '_fiberlyse_synced_legend_styles_by_tab', {})
        if not isinstance(styles, dict):
            return
        tab_name = str(getattr(plot_tab, 'tab_name', '') or '')
        snapshot = styles.get(tab_name)
        if snapshot is None:
            return
        setattr(plot_tab, '_fiberlyse_legend_sync_suppress', True)
        _fiberlyse_legend_sync_apply_to_plot(plot_tab, snapshot, redraw=False)
    except Exception:
        pass
    finally:
        try:setattr(plot_tab, '_fiberlyse_legend_sync_suppress', False)
        except Exception:pass

_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_LEGEND_SYNC = MainAppTk.__init__
def _fiberlyse_mainapptk_init_with_legend_sync(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_LEGEND_SYNC(self, *args, **kwargs)
    try:setattr(self.root, '_fiberlyse_app', self)
    except Exception:pass
    try:
        if not isinstance(getattr(self, '_fiberlyse_synced_legend_styles_by_tab', None), dict):
            self._fiberlyse_synced_legend_styles_by_tab = {}
    except Exception:
        pass
MainAppTk.__init__ = _fiberlyse_mainapptk_init_with_legend_sync

_FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_LEGEND_SYNC = PlotTabTk.__init__
def _fiberlyse_plottabtk_init_with_legend_sync(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_PLOTTABTK_INIT_FOR_LEGEND_SYNC(self, *args, **kwargs)
    try:_fiberlyse_legend_sync_apply_saved_style(self)
    except Exception:pass
PlotTabTk.__init__ = _fiberlyse_plottabtk_init_with_legend_sync

_FIBERLYSE_ORIGINAL_LEGEND_BUTTON_RELEASE_FOR_SYNC = PlotTabTk._on_legend_button_release
def _fiberlyse_plot_on_legend_button_release_with_sync(self, event):
    try:was_dragging = bool(getattr(self, '_legend_drag_state', None))
    except Exception:was_dragging = False
    result = _FIBERLYSE_ORIGINAL_LEGEND_BUTTON_RELEASE_FOR_SYNC(self, event)
    if was_dragging:
        try:_fiberlyse_legend_sync_publish(self, reason='drag')
        except Exception:pass
    return result
PlotTabTk._on_legend_button_release = _fiberlyse_plot_on_legend_button_release_with_sync

_FIBERLYSE_ORIGINAL_SHOW_LEGEND_BOX_DIALOG_FOR_SYNC = PlotTabTk.show_legend_box_dialog
def _fiberlyse_plot_show_legend_box_dialog_with_sync(self, *args, **kwargs):
    had_legend = _fiberlyse_legend_sync_has_legend(self)
    result = _FIBERLYSE_ORIGINAL_SHOW_LEGEND_BOX_DIALOG_FOR_SYNC(self, *args, **kwargs)
    try:
        if had_legend or _fiberlyse_legend_sync_has_legend(self):
            _fiberlyse_legend_sync_publish(self, reason='dialog')
    except Exception:
        pass
    return result
PlotTabTk.show_legend_box_dialog = _fiberlyse_plot_show_legend_box_dialog_with_sync

_FIBERLYSE_ORIGINAL_ON_MPL_BUTTON_PRESS_FOR_LEGEND_SYNC = PlotTabTk._on_mpl_button_press
def _fiberlyse_plot_on_mpl_button_press_with_legend_sync(self, event):
    try:
        in_legend = _fiberlyse_legend_contains_event(self, event) is not None
    except Exception:
        in_legend = False
    result = _FIBERLYSE_ORIGINAL_ON_MPL_BUTTON_PRESS_FOR_LEGEND_SYNC(self, event)
    if in_legend:
        try:_fiberlyse_legend_sync_publish(self, reason='legend_text_or_color')
        except Exception:pass
    return result
PlotTabTk._on_mpl_button_press = _fiberlyse_plot_on_mpl_button_press_with_legend_sync

# Refresh existing lazy-created channel widgets with any saved legend style after creation.
_FIBERLYSE_ORIGINAL_ENSURE_MOUSE_WIDGET_FOR_LEGEND_SYNC = MainAppTk._ensure_mouse_widget
def _fiberlyse_ensure_mouse_widget_with_legend_sync(self, mid: str) -> None:
    existed = False
    try:existed = mid in getattr(self, '_channel_widgets', {})
    except Exception:existed = False
    _FIBERLYSE_ORIGINAL_ENSURE_MOUSE_WIDGET_FOR_LEGEND_SYNC(self, mid)
    if existed:
        return
    try:
        widget = getattr(self, '_channel_widgets', {}).get(mid)
    except Exception:
        widget = None
    if widget is None:
        return
    try:
        for tab in _fiberlyse_legend_sync_iter_plot_tabs(self):
            try:_fiberlyse_legend_sync_apply_saved_style(tab)
            except Exception:pass
    except Exception:
        pass
MainAppTk._ensure_mouse_widget = _fiberlyse_ensure_mouse_widget_with_legend_sync
# ---- End Fiberlyse synced legend settings across analysis files ----


# ---- FiberLyse modernization / scalable queue extension ----
# Implements the research recommendations in-place: unlimited file queue,
# friendly aliases, centralized shortcuts, a left navigation panel, progress /
# retry / cancel controls, theme tokens, tooltips, and active-file-only plot
# materialization while preserving existing analysis and plotting commands.
FIBERLYSE_UNLIMITED_FILES = True

class FiberlyseShortcutRegistry:
    def __init__(self, root):
        self.root = root
        self.actions: Dict[str, Dict[str, Any]] = {}
        self._bound: set = set()

    def register(self, action: str, label: str, sequences, callback) -> None:
        if isinstance(sequences, str):
            sequences = [sequences]
        entry = self.actions.setdefault(str(action), {'label': str(label), 'sequences': [], 'callback': callback})
        entry['label'] = str(label)
        entry['callback'] = callback
        for seq in list(sequences):
            seq = str(seq)
            if seq not in entry['sequences']:
                entry['sequences'].append(seq)
            key = (str(action), seq)
            if key in self._bound:
                continue
            try:
                self.root.bind_all(seq, callback, add='+')
                self._bound.add(key)
            except Exception:
                pass

    def as_rows(self) -> List[Tuple[str, str, str]]:
        rows = []
        for action, info in sorted(self.actions.items()):
            rows.append((action, str(info.get('label', action)), ', '.join(info.get('sequences', []))))
        return rows

class FiberlyseToolTip:
    def __init__(self, widget, text: str, delay_ms: int=550):
        self.widget = widget
        self.text = str(text)
        self.delay_ms = int(delay_ms)
        self._after_id = None
        self._tip = None
        try:
            widget.bind('<Enter>', self._schedule, add='+')
            widget.bind('<Leave>', self._hide, add='+')
            widget.bind('<ButtonPress>', self._hide, add='+')
            widget.bind('<FocusIn>', self._schedule, add='+')
            widget.bind('<FocusOut>', self._hide, add='+')
        except Exception:
            pass

    def _schedule(self, _event=None):
        self._hide()
        try:
            self._after_id = self.widget.after(self.delay_ms, self._show)
        except Exception:
            self._after_id = None

    def _show(self):
        self._after_id = None
        if self._tip is not None:
            return
        try:
            x = self.widget.winfo_rootx() + 18
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 8
            top = tk.Toplevel(self.widget)
            top.wm_overrideredirect(True)
            top.wm_geometry(f'+{x}+{y}')
            lbl = ttk.Label(top, text=self.text, justify='left', wraplength=360, padding=(8, 5), style='Tooltip.TLabel')
            lbl.pack(fill=tk.BOTH, expand=True)
            self._tip = top
        except Exception:
            self._tip = None

    def _hide(self, _event=None):
        try:
            if self._after_id is not None:
                self.widget.after_cancel(self._after_id)
        except Exception:
            pass
        self._after_id = None
        tip = self._tip
        self._tip = None
        if tip is not None:
            try:tip.destroy()
            except Exception:pass

def _fiberlyse_shortcuts(self) -> FiberlyseShortcutRegistry:
    reg = getattr(self, 'shortcut_registry', None)
    if reg is None or not isinstance(reg, FiberlyseShortcutRegistry):
        reg = FiberlyseShortcutRegistry(self.root)
        self.shortcut_registry = reg
    return reg

def _fiberlyse_install_file_map_hotkeys_modern(self) -> None:
    _fiberlyse_shortcuts(self).register('file_map', 'Show file queue / file-number mapping', ['<Control-j>', '<Control-J>'], self._on_ctrl_j)

def _fiberlyse_install_time_marker_hotkeys_modern(self) -> None:
    _fiberlyse_shortcuts(self).register('time_marker', 'Add / remove time marker on active plot', ['<Control-i>', '<Control-I>'], self._on_ctrl_i)

def _fiberlyse_install_auc_hotkeys_modern(self) -> None:
    _fiberlyse_shortcuts(self).register('auc_interval', 'Calculate AUC interval on active plot', ['<Control-u>', '<Control-U>'], self._on_ctrl_u)

def _fiberlyse_install_axis_interval_hotkeys_modern(self) -> None:
    _fiberlyse_shortcuts(self).register('axis_controls', 'Open graph customization dialog', ['<Control-k>', '<Control-K>'], self._on_ctrl_k_axis_intervals)

def _fiberlyse_install_legend_hotkeys_modern(self) -> None:
    _fiberlyse_shortcuts(self).register('legend_box', 'Open graph customization legend controls', ['<Control-l>', '<Control-L>'], self._on_ctrl_l_legend_box)

# Replace scattered bind_all installers before MainAppTk.__init__ runs, so the
# existing wrapper chain still installs the same hotkeys through one registry.
MainAppTk._install_file_map_hotkeys = _fiberlyse_install_file_map_hotkeys_modern
MainAppTk._install_time_marker_hotkeys = _fiberlyse_install_time_marker_hotkeys_modern
MainAppTk._install_auc_hotkeys = _fiberlyse_install_auc_hotkeys_modern
MainAppTk._install_axis_interval_hotkeys = _fiberlyse_install_axis_interval_hotkeys_modern
MainAppTk._install_legend_hotkeys = _fiberlyse_install_legend_hotkeys_modern

def _fiberlyse_norm_path(path: str) -> str:
    try:
        return os.path.abspath(os.path.expanduser(str(path)))
    except Exception:
        return str(path)

def _fiberlyse_existing_records(self) -> Dict[str, Dict[str, Any]]:
    recs = getattr(self, '_file_queue_records', None)
    if not isinstance(recs, dict):
        recs = {}
        self._file_queue_records = recs
    return recs

def _fiberlyse_sync_queue_from_paths(self) -> None:
    # Normalize and de-duplicate while preserving order.
    paths: List[str] = []
    seen: set = set()
    for p in list(getattr(self, 'csv_paths', []) or []):
        npth = _fiberlyse_norm_path(p)
        if not npth or npth in seen:
            continue
        seen.add(npth)
        paths.append(npth)
    self.csv_paths = paths
    key_by_path = getattr(self, '_file_key_by_path', None)
    if not isinstance(key_by_path, dict):
        key_by_path = {}
        self._file_key_by_path = key_by_path
    # Remove keys for paths no longer in the queue.
    for p in list(key_by_path.keys()):
        if p not in seen:
            key_by_path.pop(p, None)
    next_idx = 1
    used_keys = set(key_by_path.values())
    for p in paths:
        if p in key_by_path:
            continue
        while f'F{next_idx}' in used_keys:
            next_idx += 1
        key_by_path[p] = f'F{next_idx}'
        used_keys.add(f'F{next_idx}')
    aliases = _unique_aliases(paths)
    old = _fiberlyse_existing_records(self)
    new_records: Dict[str, Dict[str, Any]] = {}
    for p, alias in zip(paths, aliases):
        prev = old.get(p, {}) if isinstance(old.get(p, {}), dict) else {}
        new_records[p] = {
            'path': p,
            'key': key_by_path.get(p),
            'alias': alias,
            'status': prev.get('status', 'queued'),
            'channels': prev.get('channels', ''),
            'note': prev.get('note', ''),
        }
    self._file_queue_records = new_records
    self._queue_iid_to_path = {rec['key']: p for p, rec in new_records.items() if rec.get('key')}
    self._file_path_by_key = {rec['key']: p for p, rec in new_records.items() if rec.get('key')}
    self._file_alias_by_key = {rec['key']: rec['alias'] for rec in new_records.values() if rec.get('key')}
    self._result_file_order = [rec['key'] for rec in new_records.values() if rec.get('key')]

def _fiberlyse_apply_modern_theme(self) -> None:
    try:
        style = ttk.Style(self.root)
        try:
            if 'clam' in style.theme_names():
                style.theme_use('clam')
        except Exception:
            pass
        self.root.option_add('*Font', 'TkDefaultFont 10')
        style.configure('TFrame', background='#f6f7f9')
        style.configure('Fiber.TFrame', background='#f6f7f9')
        style.configure('Panel.TFrame', background='#ffffff', relief='flat')
        style.configure('TLabel', background='#f6f7f9')
        style.configure('Panel.TLabel', background='#ffffff')
        style.configure('Title.TLabel', background='#ffffff', font=('TkDefaultFont', 11, 'bold'))
        style.configure('Hint.TLabel', background='#ffffff', foreground='#5f6b7a')
        style.configure('Tooltip.TLabel', background='#222222', foreground='#ffffff', relief='solid', borderwidth=1)
        style.configure('Accent.TButton', padding=(10, 4))
        style.configure('Queue.Treeview', rowheight=24, fieldbackground='#ffffff', background='#ffffff')
        style.configure('Queue.Treeview.Heading', font=('TkDefaultFont', 9, 'bold'))
        style.configure('TNotebook.Tab', padding=(12, 5))
        style.configure('TLabelframe', background='#f6f7f9')
        style.configure('TLabelframe.Label', background='#f6f7f9', font=('TkDefaultFont', 10, 'bold'))
    except Exception:
        pass

def _fiberlyse_attach_tooltip(widget, text: str) -> None:
    try:
        widget._fiberlyse_tooltip = FiberlyseToolTip(widget, text)
    except Exception:
        pass

def _fiberlyse_install_tooltips(self) -> None:
    tips = [
        (getattr(self, 'btn_add', None), 'Add CSV files to the queue. There is no fixed file limit; practical limits depend on memory and file size.'),
        (getattr(self, 'btn_clear', None), 'Clear the queue, loaded results, and visible plots.'),
        (getattr(self, 'btn_run', None), 'Analyze the full queue. The left panel tracks queued, analyzing, analyzed, failed, and canceled files.'),
        (getattr(self, 'chk_artifact_enabled', None), 'Enable derivative-based MAD artifact detection.'),
        (getattr(self, 'spin_factor', None), 'MAD multiplier for artifact detection. Larger values are less aggressive.'),
        (getattr(self, 'spin_pad', None), 'Number of neighboring samples to also mark around each detected artifact.'),
        (getattr(self, 'chk_shared', None), 'Only remove artifacts that occur in both isosbestic and excitatory streams at matching times.'),
        (getattr(self, 'spin_acq_fps', None), 'Original camera acquisition FPS. Effective per-state sampling is FPS / 2 for alternating LED states.'),
        (getattr(self, 'spin_smooth_win', None), 'Smoothing window in samples for the smoothed normalization view and batch plots.'),
        (getattr(self, 'chk_interp', None), 'Fill artifact holes for the active analysis path. Frequency analysis still uses the no-interpolation version.'),
        (getattr(self, 'cmb_norm', None), 'Choose the normalization displayed in the normalization, smoothed, compare, and average views.'),
        (getattr(self, 'spin_interval_start', None), 'Start time for interval-based zF baseline statistics.'),
        (getattr(self, 'spin_interval_end', None), 'End time for interval-based zF baseline statistics.'),
    ]
    for widget, text in tips:
        if widget is not None:
            _fiberlyse_attach_tooltip(widget, text)

def _fiberlyse_refresh_files_label(self) -> None:
    try:
        n = len(getattr(self, 'csv_paths', []) or [])
        if n <= 0:
            self.lbl_file.config(text='No CSV files in queue.')
            return
        recs = _fiberlyse_existing_records(self)
        analyzed = sum(1 for r in recs.values() if r.get('status') == 'analyzed')
        failed = sum(1 for r in recs.values() if r.get('status') == 'failed')
        self.lbl_file.config(text=f'{n} CSV file(s) in queue | analyzed {analyzed} | failed {failed}')
    except Exception:
        try:
            self.lbl_file.config(text='File queue ready.')
        except Exception:
            pass

def _fiberlyse_refresh_queue_tree(self) -> None:
    tree = getattr(self, 'queue_tree', None)
    if tree is None:
        return
    try:
        _fiberlyse_sync_queue_from_paths(self)
        recs = _fiberlyse_existing_records(self)
        filt = ''
        try:filt = (self.var_queue_filter.get() or '').strip().lower()
        except Exception:filt = ''
        tree.delete(*tree.get_children(''))
        self._queue_iid_to_path = {}
        for p in list(getattr(self, 'csv_paths', []) or []):
            rec = recs.get(p, {})
            alias = str(rec.get('alias') or _basename_no_ext(p))
            hay = f"{alias} {p} {rec.get('status','')} {rec.get('note','')}".lower()
            if filt and filt not in hay:
                continue
            key = str(rec.get('key') or f'F{len(self._queue_iid_to_path)+1}')
            self._queue_iid_to_path[key] = p
            status = str(rec.get('status') or 'queued')
            channels = rec.get('channels', '')
            note = rec.get('note', '')
            if isinstance(note, str) and len(note) > 80:
                note = note[:77] + '...'
            tree.insert('', tk.END, iid=key, text=alias, values=(status, channels, note), tags=(status,))
        for tag, fg in [('queued', '#5f6b7a'), ('analyzing', '#8a5a00'), ('analyzed', '#0f6b3f'), ('failed', '#a12020'), ('canceled', '#7a4b00')]:
            try:tree.tag_configure(tag, foreground=fg)
            except Exception:pass
        active = getattr(self, '_active_file_key', None)
        if active and active in tree.get_children(''):
            try:
                tree.selection_set(active)
                tree.see(active)
            except Exception:
                pass
        _fiberlyse_refresh_files_label(self)
    except Exception:
        pass

def _fiberlyse_mark_queue_status(self, path: str, status: str, channels: Any='', note: str='') -> None:
    try:
        p = _fiberlyse_norm_path(path)
        recs = _fiberlyse_existing_records(self)
        if p in recs:
            recs[p]['status'] = str(status)
            if channels != '':
                recs[p]['channels'] = channels
            recs[p]['note'] = str(note or '')
        _fiberlyse_refresh_queue_tree(self)
    except Exception:
        pass

def _fiberlyse_add_paths_to_queue(self, paths) -> int:
    _fiberlyse_sync_queue_from_paths(self)
    existing = set(getattr(self, 'csv_paths', []) or [])
    added = 0
    first_new = None
    for p in list(paths or []):
        npth = _fiberlyse_norm_path(p)
        if not npth or npth in existing:
            continue
        existing.add(npth)
        self.csv_paths.append(npth)
        if first_new is None:
            first_new = npth
        added += 1
    _fiberlyse_sync_queue_from_paths(self)
    recs = _fiberlyse_existing_records(self)
    for p in list(getattr(self, 'csv_paths', []) or []):
        if recs.get(p, {}).get('status') not in ('analyzed', 'failed', 'analyzing', 'canceled'):
            recs[p]['status'] = 'queued'
    if first_new is not None:
        try:self._active_file_key = recs[first_new]['key']
        except Exception:pass
    _fiberlyse_refresh_queue_tree(self)
    return added

def _fiberlyse_selected_queue_paths(self) -> List[str]:
    tree = getattr(self, 'queue_tree', None)
    out: List[str] = []
    if tree is not None:
        try:
            mapping = getattr(self, '_queue_iid_to_path', {})
            for iid in tree.selection():
                p = mapping.get(str(iid))
                if p:
                    out.append(p)
        except Exception:
            pass
    if not out:
        key = getattr(self, '_active_file_key', None)
        try:
            p = getattr(self, '_file_path_by_key', {}).get(key)
            if p:
                out.append(p)
        except Exception:
            pass
    return out

def _fiberlyse_install_modern_shell(self) -> None:
    _fiberlyse_apply_modern_theme(self)
    try:self.root.title('FiberLyse - Fiber Photometry Analyzer')
    except Exception:pass
    try:self.root.minsize(1050, 660)
    except Exception:pass
    # Hide the old compact file selector row. The left queue panel replaces it,
    # while the legacy combobox is still kept in memory for compatibility.
    try:
        self._legacy_view_row = self.cmb_view_file.master
        self._legacy_view_row.pack_forget()
    except Exception:
        pass
    # Replace the old full-width notebook with a scalable queue/workspace split.
    try:
        old_outer = self.outer_tabs
        old_outer.pack_forget()
        old_outer.destroy()
    except Exception:
        pass
    try:
        self.workspace_pane = ttk.Panedwindow(self.root, orient=tk.HORIZONTAL)
        self.workspace_pane.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=(0, 6))
        self.queue_panel = ttk.Frame(self.workspace_pane, style='Panel.TFrame', padding=(8, 8))
        self.workspace_holder = ttk.Frame(self.workspace_pane, style='Fiber.TFrame')
        self.workspace_pane.add(self.queue_panel, weight=0)
        self.workspace_pane.add(self.workspace_holder, weight=1)
        self.outer_tabs = ttk.Notebook(self.workspace_holder)
        self.outer_tabs.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._frame_to_mid = {}
        self.outer_tabs.bind('<<NotebookTabChanged>>', self._on_outer_tab_changed)
    except Exception:
        # If the shell cannot be rebuilt, fall back to a fresh notebook on root.
        self.outer_tabs = ttk.Notebook(self.root)
        self.outer_tabs.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=(0, 6))
        self.outer_tabs.bind('<<NotebookTabChanged>>', self._on_outer_tab_changed)
        return
    # Queue panel.
    panel = self.queue_panel
    panel.columnconfigure(0, weight=1)
    ttk.Label(panel, text='Files', style='Title.TLabel').grid(row=0, column=0, sticky='w')
    ttk.Label(panel, text='Queue, analyze, retry, and switch files without creating one top-level tab per file.', style='Hint.TLabel', wraplength=260, justify='left').grid(row=1, column=0, sticky='ew', pady=(2, 8))
    search_row = ttk.Frame(panel, style='Panel.TFrame')
    search_row.grid(row=2, column=0, sticky='ew', pady=(0, 6))
    search_row.columnconfigure(1, weight=1)
    ttk.Label(search_row, text='Search:', style='Panel.TLabel').grid(row=0, column=0, sticky='w', padx=(0, 6))
    self.var_queue_filter = tk.StringVar(value='')
    ent_filter = ttk.Entry(search_row, textvariable=self.var_queue_filter, width=20)
    ent_filter.grid(row=0, column=1, sticky='ew')
    try:self.var_queue_filter.trace_add('write', lambda *_: _fiberlyse_refresh_queue_tree(self))
    except Exception:pass
    tree_frame = ttk.Frame(panel, style='Panel.TFrame')
    tree_frame.grid(row=3, column=0, sticky='nsew')
    panel.rowconfigure(3, weight=1)
    tree_frame.rowconfigure(0, weight=1)
    tree_frame.columnconfigure(0, weight=1)
    self.queue_tree = ttk.Treeview(tree_frame, columns=('status', 'channels', 'note'), show='tree headings', selectmode='extended', style='Queue.Treeview', height=16)
    self.queue_tree.heading('#0', text='File')
    self.queue_tree.heading('status', text='Status')
    self.queue_tree.heading('channels', text='Ch')
    self.queue_tree.heading('note', text='Note')
    self.queue_tree.column('#0', width=150, minwidth=90, stretch=True)
    self.queue_tree.column('status', width=76, minwidth=70, stretch=False)
    self.queue_tree.column('channels', width=42, minwidth=34, stretch=False, anchor='center')
    self.queue_tree.column('note', width=105, minwidth=60, stretch=True)
    ysb = ttk.Scrollbar(tree_frame, orient='vertical', command=self.queue_tree.yview)
    self.queue_tree.configure(yscrollcommand=ysb.set)
    self.queue_tree.grid(row=0, column=0, sticky='nsew')
    ysb.grid(row=0, column=1, sticky='ns')
    self.queue_tree.bind('<<TreeviewSelect>>', self._on_queue_tree_select, add='+')
    self.queue_tree.bind('<Double-Button-1>', self._on_queue_tree_double_click, add='+')
    btn_grid = ttk.Frame(panel, style='Panel.TFrame')
    btn_grid.grid(row=4, column=0, sticky='ew', pady=(8, 0))
    for c in range(2):
        btn_grid.columnconfigure(c, weight=1)
    self.btn_add_folder = ttk.Button(btn_grid, text='Add folder...', command=self.add_folder_csvs)
    self.btn_analyze_selected = ttk.Button(btn_grid, text='Analyze selected', command=self.run_analysis_selected)
    self.btn_retry_failed = ttk.Button(btn_grid, text='Retry', command=self.retry_failed_files)
    self.btn_remove_selected = ttk.Button(btn_grid, text='Remove selected', command=self.remove_selected_files)
    self.btn_open_folder = ttk.Button(btn_grid, text='Open folder', command=self.open_selected_folder)
    self.btn_shortcuts = ttk.Button(btn_grid, text='Shortcuts', command=self.show_shortcuts_dialog)
    self.btn_add_folder.grid(row=0, column=0, sticky='ew', padx=(0, 4), pady=2)
    self.btn_analyze_selected.grid(row=0, column=1, sticky='ew', padx=(4, 0), pady=2)
    self.btn_retry_failed.grid(row=1, column=0, sticky='ew', padx=(0, 4), pady=2)
    self.btn_remove_selected.grid(row=1, column=1, sticky='ew', padx=(4, 0), pady=2)
    self.btn_open_folder.grid(row=2, column=0, sticky='ew', padx=(0, 4), pady=2)
    self.btn_shortcuts.grid(row=2, column=1, sticky='ew', padx=(4, 0), pady=2)
    ttk.Label(panel, text='Made by Christoffer Salling', style='Hint.TLabel', wraplength=260, justify='left').grid(row=5, column=0, sticky='ew', pady=(8, 0))
    # Progress and cancel in the existing top toolbar.
    try:
        top = self.btn_run.master
        self.progress_var = tk.DoubleVar(value=0.0)
        self.progress_bar = ttk.Progressbar(top, variable=self.progress_var, maximum=100.0, length=180, mode='determinate')
        self.progress_bar.grid(row=0, column=4, sticky='e', padx=(14, 6))
        self.btn_cancel = ttk.Button(top, text='Cancel', command=self.cancel_analysis, state=tk.DISABLED)
        self.btn_cancel.grid(row=0, column=5, sticky='e')
        try:top.columnconfigure(0, weight=1)
        except Exception:pass
    except Exception:
        pass
    try:self.btn_run.config(style='Accent.TButton')
    except Exception:pass
    for widget, text in [
        (self.queue_tree, 'Select a file to show its analyzed channels. Double-click a queued or failed file to analyze it.'),
        (self.btn_add_folder, 'Add every CSV found in a folder and its subfolders.'),
        (self.btn_analyze_selected, 'Analyze only the selected queue rows.'),
        (self.btn_retry_failed, 'Re-run files currently marked failed.'),
        (self.btn_remove_selected, 'Remove selected files and their loaded results from the queue.'),
        (self.btn_open_folder, 'Open the containing folder for the selected file.'),
        (self.btn_shortcuts, 'Show the preserved keyboard shortcut registry.'),
    ]:
        _fiberlyse_attach_tooltip(widget, text)
    _fiberlyse_install_tooltips(self)
    _fiberlyse_refresh_queue_tree(self)

def _fiberlyse_on_queue_tree_select(self, _event=None) -> None:
    tree = getattr(self, 'queue_tree', None)
    if tree is None:
        return
    try:
        sel = list(tree.selection())
        if not sel:
            return
        iid = str(sel[0])
        p = getattr(self, '_queue_iid_to_path', {}).get(iid)
        rec = _fiberlyse_existing_records(self).get(p, {}) if p else {}
        key = rec.get('key') or iid
        self._active_file_key = key
        alias = rec.get('alias') or self._file_alias_by_key.get(key, key)
        try:self.var_view_file.set(alias)
        except Exception:pass
        if getattr(self, '_results', None) is not None:
            self._refresh_visible_tabs()
        try:self.status.set(f'Selected {alias} ({rec.get("status", "queued")}).')
        except Exception:pass
    except Exception:
        pass

def _fiberlyse_on_queue_tree_double_click(self, _event=None) -> None:
    try:
        paths = _fiberlyse_selected_queue_paths(self)
        if not paths:
            return
        rec = _fiberlyse_existing_records(self).get(paths[0], {})
        if rec.get('status') in ('queued', 'failed', 'canceled'):
            self.run_analysis_selected()
    except Exception:
        pass

def _fiberlyse_set_buttons_for_analysis(self, running: bool) -> None:
    state_main = tk.DISABLED if running else tk.NORMAL
    try:self.btn_run.config(state=state_main)
    except Exception:pass
    try:self.btn_add.config(state=state_main)
    except Exception:pass
    try:self.btn_clear.config(state=state_main)
    except Exception:pass
    for name in ['btn_add_folder', 'btn_analyze_selected', 'btn_retry_failed', 'btn_remove_selected']:
        try:getattr(self, name).config(state=state_main)
        except Exception:pass
    try:self.btn_cancel.config(state=(tk.NORMAL if running else tk.DISABLED))
    except Exception:pass

def _fiberlyse_parse_analysis_settings(self) -> Optional[Dict[str, Any]]:
    try:
        return {
            'artifact_enabled': bool(self.var_artifact_enabled.get()),
            'artifact_factor': float(self.var_factor.get()),
            'artifact_pad': int(float(self.var_pad.get())),
            'require_shared': bool(self.var_shared.get()),
            'align_mode': DEFAULT_ALIGN_MODE,
            'acq_fps': np.nan,  # V24.1: analysis rate comes from timestamps only
            'smooth_win': max(1, int(float(self.var_smooth_win.get()))),
            'mode': self._read_norm_mode(),
            'interval_start': float(self.var_interval_start.get()),
            'interval_end': float(self.var_interval_end.get()),
            'use_interp': bool(self.var_interp.get()),
            'fit_windows': list(DEFAULT_FIT_WINDOWS),
        }
    except Exception as e:
        try:messagebox.showerror('Invalid settings', f'Could not parse settings:\n\n{e}')
        except Exception:pass
        return None

def _fiberlyse_remove_results_for_file(local_results: Dict[str, ChannelResult], local_display: Dict[str, str], file_key: str) -> None:
    prefix = f'{file_key}:'
    for mid in list(local_results.keys()):
        if str(mid).startswith(prefix):
            local_results.pop(mid, None)
            local_display.pop(mid, None)

def _fiberlyse_file_meta_all(self) -> List[Tuple[str, str, str]]:
    _fiberlyse_sync_queue_from_paths(self)
    out = []
    for p in list(getattr(self, 'csv_paths', []) or []):
        rec = _fiberlyse_existing_records(self).get(p, {})
        key = rec.get('key')
        alias = rec.get('alias') or _basename_no_ext(p)
        if key:
            out.append((key, alias, p))
    return out

def _fiberlyse_run_analysis_paths(self, selected_paths: Optional[List[str]]=None) -> None:
    _fiberlyse_sync_queue_from_paths(self)
    all_paths = list(getattr(self, 'csv_paths', []) or [])
    if not all_paths:
        try:messagebox.showwarning('No files', 'Please add one or more CSV files first.')
        except Exception:pass
        return
    paths_to_analyze = [_fiberlyse_norm_path(p) for p in (selected_paths if selected_paths is not None else all_paths)]
    paths_to_analyze = [p for p in paths_to_analyze if p in set(all_paths)]
    if not paths_to_analyze:
        try:messagebox.showinfo('Analyze selected', 'No selected files are available to analyze.')
        except Exception:pass
        return
    settings = _fiberlyse_parse_analysis_settings(self)
    if settings is None:
        return
    recs_snapshot = {p: dict(r) for p, r in _fiberlyse_existing_records(self).items()}
    file_meta_all = _fiberlyse_file_meta_all(self)
    mode = settings['mode']
    self._analysis_cancel_requested = False
    _fiberlyse_set_buttons_for_analysis(self, True)
    try:
        self.progress_bar.config(maximum=max(1, len(paths_to_analyze)))
        self.progress_var.set(0.0)
    except Exception:
        pass
    try:self.status.set(f'Analyzing {len(paths_to_analyze)} queued file(s)...')
    except Exception:pass
    for p in paths_to_analyze:
        _fiberlyse_mark_queue_status(self, p, 'queued', channels=recs_snapshot.get(p, {}).get('channels', ''), note='waiting')

    def worker():
        local_results: Dict[str, ChannelResult] = dict(getattr(self, '_results', {}) or {})
        local_display: Dict[str, str] = dict(getattr(self, '_mouse_display', {}) or {})
        successes = 0
        failures = 0
        canceled = False
        for idx, path in enumerate(paths_to_analyze, start=1):
            if bool(getattr(self, '_analysis_cancel_requested', False)):
                canceled = True
                self.root.after(0, lambda p=path: _fiberlyse_mark_queue_status(self, p, 'canceled', note='not run'))
                continue
            rec = recs_snapshot.get(path, {})
            file_key = str(rec.get('key') or f'F{idx}')
            alias = str(rec.get('alias') or _basename_no_ext(path))
            self.root.after(0, lambda idx=idx, alias=alias: self.status.set(f'Analyzing {idx}/{len(paths_to_analyze)}: {alias}...'))
            self.root.after(0, lambda p=path: _fiberlyse_mark_queue_status(self, p, 'analyzing', note='running'))
            try:
                per_file = analyze_csv(path, artifact_enabled=settings['artifact_enabled'], artifact_factor=settings['artifact_factor'], artifact_method='mad', artifact_pad=settings['artifact_pad'], require_shared=settings['require_shared'], align_mode=settings['align_mode'], fit_windows=settings['fit_windows'], acq_fps_hz=settings['acq_fps'] if np.isfinite(settings['acq_fps']) and settings['acq_fps'] > 0 else None, smooth_window=settings['smooth_win'], zf_interval_start_s=settings['interval_start'], zf_interval_end_s=settings['interval_end'], use_linear_interp=settings['use_interp'])
                _fiberlyse_remove_results_for_file(local_results, local_display, file_key)
                for gcol in sorted(per_file.keys()):
                    mid = f'{file_key}:{gcol}'
                    local_results[mid] = per_file[gcol]
                    local_display[mid] = f'{alias}:{gcol}'
                successes += 1
                self.root.after(0, lambda p=path, n=len(per_file): _fiberlyse_mark_queue_status(self, p, 'analyzed', channels=n, note=''))
            except Exception as exc:
                _fiberlyse_remove_results_for_file(local_results, local_display, file_key)
                failures += 1
                msg = str(exc)
                self.root.after(0, lambda p=path, msg=msg: _fiberlyse_mark_queue_status(self, p, 'failed', channels='', note=msg))
            finally:
                self.root.after(0, lambda v=idx: self.progress_var.set(float(v)) if hasattr(self, 'progress_var') else None)
        # Rebuild display order by queue order, then channel name.
        order: List[str] = []
        for _key, _alias, _path in file_meta_all:
            mids = sorted([mid for mid in local_results.keys() if str(mid).startswith(str(_key) + ':')], key=lambda s: str(s).split(':', 1)[-1])
            order.extend(mids)
        def finish():
            _fiberlyse_set_buttons_for_analysis(self, False)
            try:self.btn_cancel.config(state=tk.DISABLED)
            except Exception:pass
            if local_results:
                self.on_analysis_finished(local_results, local_display, order, mode, file_meta_all)
                suffix = ' Canceled.' if canceled else ''
                try:self.status.set(f'Done. Analyzed {successes} file(s); failed {failures}.{suffix}')
                except Exception:pass
            else:
                self._results = {}
                self._mouse_display = {}
                self._mouse_order = []
                try:self.build_tabs(norm_mode=mode)
                except Exception:pass
                msg = 'No files were analyzed successfully.' + (' Analysis was canceled.' if canceled else '')
                try:self.status.set(msg)
                except Exception:pass
                if failures:
                    try:messagebox.showwarning('Analysis completed with errors', msg + '\n\nSelect failed rows and use Retry failed after checking the notes column.')
                    except Exception:pass
            _fiberlyse_refresh_queue_tree(self)
        self.root.after(0, finish)
    self._analysis_thread = threading.Thread(target=worker, daemon=True)
    self._analysis_thread.start()

def _fiberlyse_run_analysis(self) -> None:
    _fiberlyse_run_analysis_paths(self, selected_paths=None)

def _fiberlyse_run_analysis_selected(self) -> None:
    paths = _fiberlyse_selected_queue_paths(self)
    if not paths:
        try:messagebox.showinfo('Analyze selected', 'Select one or more files in the queue first.')
        except Exception:pass
        return
    _fiberlyse_run_analysis_paths(self, selected_paths=paths)

def _fiberlyse_retry_failed_files(self) -> None:
    failed = [p for p, r in _fiberlyse_existing_records(self).items() if r.get('status') == 'failed']
    if not failed:
        try:messagebox.showinfo('Retry failed', 'There are no failed files to retry.')
        except Exception:pass
        return
    _fiberlyse_run_analysis_paths(self, selected_paths=failed)

def _fiberlyse_cancel_analysis(self) -> None:
    self._analysis_cancel_requested = True
    try:self.status.set('Cancel requested. Current file will finish, then the queue stops.')
    except Exception:pass
    try:self.btn_cancel.config(state=tk.DISABLED)
    except Exception:pass

def _fiberlyse_add_csvs(self) -> None:
    paths = filedialog.askopenfilenames(title='Select CSV file(s)', filetypes=[('CSV', '*.csv'), ('All files', '*.*')])
    if not paths:
        return
    added = _fiberlyse_add_paths_to_queue(self, paths)
    try:self.status.set(f'Added {added} new CSV file(s). Queue size: {len(self.csv_paths)}.')
    except Exception:pass

def _fiberlyse_add_folder_csvs(self) -> None:
    folder = filedialog.askdirectory(title='Select folder containing CSV files')
    if not folder:
        return
    found: List[str] = []
    try:
        for root_dir, _dirs, files in os.walk(folder):
            for name in files:
                if str(name).lower().endswith('.csv'):
                    found.append(os.path.join(root_dir, name))
    except Exception as e:
        try:messagebox.showerror('Add folder', f'Could not scan folder:\n\n{e}')
        except Exception:pass
        return
    if not found:
        try:messagebox.showinfo('Add folder', 'No CSV files were found in that folder.')
        except Exception:pass
        return
    added = _fiberlyse_add_paths_to_queue(self, sorted(found))
    try:self.status.set(f'Added {added} CSV file(s) from folder. Queue size: {len(self.csv_paths)}.')
    except Exception:pass

def _fiberlyse_remove_file_results_from_app(self, file_keys: set) -> None:
    if not file_keys:
        return
    for attr in ['_results', '_mouse_display']:
        d = getattr(self, attr, None)
        if isinstance(d, dict):
            for mid in list(d.keys()):
                if str(mid).split(':', 1)[0] in file_keys:
                    d.pop(mid, None)
    try:self._mouse_order = [mid for mid in self._mouse_order if str(mid).split(':', 1)[0] not in file_keys]
    except Exception:pass
    for mid in list(getattr(self, '_channel_widgets', {}).keys()):
        if str(mid).split(':', 1)[0] in file_keys:
            try:self._channel_widgets[mid].destroy()
            except Exception:pass
            self._channel_widgets.pop(mid, None)
    for mid in list(getattr(self, '_mouse_frames', {}).keys()):
        if str(mid).split(':', 1)[0] in file_keys:
            try:self._mouse_frames[mid].destroy()
            except Exception:pass
            self._mouse_frames.pop(mid, None)

def _fiberlyse_remove_selected_files(self) -> None:
    paths = _fiberlyse_selected_queue_paths(self)
    if not paths:
        try:messagebox.showinfo('Remove selected', 'Select one or more files in the queue first.')
        except Exception:pass
        return
    recs = _fiberlyse_existing_records(self)
    keys = {str(recs.get(p, {}).get('key')) for p in paths if recs.get(p, {}).get('key')}
    self.csv_paths = [p for p in list(getattr(self, 'csv_paths', []) or []) if p not in set(paths)]
    for p in paths:
        recs.pop(p, None)
    _fiberlyse_remove_file_results_from_app(self, keys)
    self._active_file_key = None
    _fiberlyse_sync_queue_from_paths(self)
    _fiberlyse_refresh_queue_tree(self)
    try:self.build_tabs(norm_mode=self._read_norm_mode())
    except Exception:pass
    try:self.status.set(f'Removed {len(paths)} file(s) from the queue.')
    except Exception:pass

def _fiberlyse_clear_csvs(self) -> None:
    if getattr(self, '_analysis_thread', None) is not None and getattr(self._analysis_thread, 'is_alive', lambda: False)():
        try:messagebox.showwarning('Analysis running', 'Cancel the current analysis before clearing the queue.')
        except Exception:pass
        return
    self.csv_paths = []
    self._results = None
    self._mouse_display = {}
    self._mouse_order = []
    self._file_queue_records = {}
    self._file_key_by_path = {}
    self._result_file_order = []
    self._file_alias_by_key = {}
    self._file_path_by_key = {}
    self._active_file_key = None
    for widget in list(getattr(self, '_channel_widgets', {}).values()):
        try:widget.destroy()
        except Exception:pass
    self._channel_widgets = {}
    for frame in list(getattr(self, '_mouse_frames', {}).values()):
        try:frame.destroy()
        except Exception:pass
    self._mouse_frames = {}
    if self.compare_widget is not None:
        try:self.compare_widget.destroy()
        except Exception:pass
    if self.average_widget is not None:
        try:self.average_widget.destroy()
        except Exception:pass
    self.compare_widget = None
    self.average_widget = None
    try:
        for tab_id in self.outer_tabs.tabs():
            self.outer_tabs.forget(tab_id)
    except Exception:pass
    _fiberlyse_refresh_queue_tree(self)
    _fiberlyse_refresh_files_label(self)
    try:self.status.set('Cleared file queue and loaded results.')
    except Exception:pass

def _fiberlyse_open_selected_folder(self) -> None:
    paths = _fiberlyse_selected_queue_paths(self)
    if not paths:
        try:messagebox.showinfo('Open folder', 'Select a file in the queue first.')
        except Exception:pass
        return
    folder = os.path.dirname(paths[0])
    try:
        import subprocess
        if sys.platform.startswith('win'):
            os.startfile(folder)  # type: ignore[attr-defined]
        elif sys.platform == 'darwin':
            subprocess.Popen(['open', folder])
        else:
            subprocess.Popen(['xdg-open', folder])
    except Exception as e:
        try:messagebox.showerror('Open folder', f'Could not open folder:\n\n{folder}\n\n{e}')
        except Exception:pass

def _fiberlyse_show_shortcuts_dialog(self) -> None:
    try:rows = _fiberlyse_shortcuts(self).as_rows()
    except Exception:rows = []
    top = tk.Toplevel(self.root)
    top.title('Keyboard shortcuts')
    top.transient(self.root)
    frm = ttk.Frame(top, padding=10)
    frm.pack(fill=tk.BOTH, expand=True)
    ttk.Label(frm, text='Preserved keyboard commands', font=('TkDefaultFont', 11, 'bold')).pack(anchor='w')
    txt = tk.Text(frm, width=72, height=max(8, min(16, len(rows) + 3)), wrap='word')
    txt.pack(fill=tk.BOTH, expand=True, pady=(8, 8))
    if rows:
        for _action, label, seqs in rows:
            txt.insert(tk.END, f'{seqs}\t{label}\n')
    else:
        txt.insert(tk.END, 'No shortcuts registered yet.\n')
    txt.configure(state='disabled')
    ttk.Button(frm, text='Close', command=top.destroy).pack(anchor='e')

def _fiberlyse_set_result_file_choices(self, file_meta: List[Tuple[str, str, str]]) -> None:
    # Keep the legacy combobox data for compatibility, but selection is driven by
    # the queue Treeview in the modern shell.
    self._result_file_order = [file_key for file_key, _alias, _path in file_meta]
    self._file_alias_by_key = {file_key: alias for file_key, alias, _path in file_meta}
    self._file_path_by_key = {file_key: path for file_key, _alias, path in file_meta}
    values = [self._file_alias_by_key[file_key] for file_key in self._result_file_order]
    try:
        self.cmb_view_file.config(values=tuple(values), state='readonly' if values else 'disabled')
    except Exception:
        pass
    if values:
        active = getattr(self, '_active_file_key', None)
        if active not in self._result_file_order:
            active = self._result_file_order[0]
            self._active_file_key = active
        try:self.var_view_file.set(self._file_alias_by_key.get(active, values[0]))
        except Exception:pass
    else:
        try:self.var_view_file.set('')
        except Exception:pass
    _fiberlyse_refresh_queue_tree(self)

def _fiberlyse_get_selected_view_file_key(self) -> Optional[str]:
    active = getattr(self, '_active_file_key', None)
    if active:
        return active
    if getattr(self, '_result_file_order', None):
        self._active_file_key = self._result_file_order[0]
        return self._active_file_key
    return None

def _fiberlyse_on_view_file_changed(self, _event=None) -> None:
    # Legacy combobox path; hidden in modern shell, but kept functional.
    selected_alias = (self.var_view_file.get() or '').strip()
    for file_key in list(getattr(self, '_result_file_order', []) or []):
        if self._file_alias_by_key.get(file_key) == selected_alias:
            self._active_file_key = file_key
            break
    if getattr(self, '_results', None):
        self._refresh_visible_tabs()

def _fiberlyse_placeholder_text(self, file_key: Optional[str]) -> str:
    if not file_key:
        return 'Select a file in the queue.'
    path = getattr(self, '_file_path_by_key', {}).get(file_key)
    rec = _fiberlyse_existing_records(self).get(path, {}) if path else {}
    alias = rec.get('alias') or self._file_alias_by_key.get(file_key, file_key)
    status = rec.get('status', 'queued')
    note = rec.get('note', '')
    if status == 'failed':
        return f'{alias}\n\nStatus: failed\n\n{note}\n\nSelect the row and click Retry failed after fixing the input or settings.'
    if status in ('queued', 'canceled'):
        return f'{alias}\n\nStatus: {status}\n\nClick Analyze selected, double-click this queue row, or run the full queue.'
    return f'{alias}\n\nNo analyzed channels are available for this file.'

def _fiberlyse_build_tabs(self, norm_mode: str):
    outer_tabs = self.outer_tabs
    try:
        for tab_id in outer_tabs.tabs():
            outer_tabs.forget(tab_id)
    except Exception:
        pass
    for widget in list(getattr(self, '_channel_widgets', {}).values()):
        try:widget.destroy()
        except Exception:pass
    for frame in list(getattr(self, '_mouse_frames', {}).values()):
        try:frame.destroy()
        except Exception:pass
    self._channel_widgets = {}
    self._mouse_frames = {}
    self._frame_to_mid = {}
    if self.compare_widget is not None:
        try:self.compare_widget.destroy()
        except Exception:pass
    if self.average_widget is not None:
        try:self.average_widget.destroy()
        except Exception:pass
    self.compare_widget = None
    self.average_widget = None
    if not self._results:
        return
    self._refresh_visible_tabs()

def _fiberlyse_refresh_visible_tabs(self) -> None:
    outer_tabs = self.outer_tabs
    try:
        for tab_id in outer_tabs.tabs():
            outer_tabs.forget(tab_id)
    except Exception:
        pass
    self._frame_to_mid = {}
    file_key = self._get_selected_view_file_key()
    mids: List[str] = []
    results = getattr(self, '_results', {}) or {}
    for mid in list(getattr(self, '_mouse_order', []) or []):
        if mid not in results:
            continue
        if file_key is not None and self._mouse_file_key(mid) != file_key:
            continue
        mids.append(mid)
    # Evict inactive plot widgets and frames so unlimited queues do not retain a
    # Matplotlib canvas for every previously visited file.
    keep = set(mids)
    for mid in list(getattr(self, '_channel_widgets', {}).keys()):
        if mid not in keep:
            try:self._channel_widgets[mid].destroy()
            except Exception:pass
            self._channel_widgets.pop(mid, None)
    for mid in list(getattr(self, '_mouse_frames', {}).keys()):
        if mid not in keep:
            try:self._mouse_frames[mid].destroy()
            except Exception:pass
            self._mouse_frames.pop(mid, None)
    if mids:
        for mid in mids:
            frame = self._mouse_frames.get(mid)
            if frame is None:
                frame = ttk.Frame(outer_tabs)
                display_label = self._mouse_display.get(mid, mid)
                ttk.Label(frame, text=f'Click this tab to load plots for {display_label} (lazy-loaded for speed).', foreground='#666666').pack(anchor='w', padx=10, pady=10)
                self._mouse_frames[mid] = frame
            tab_label = self._mouse_channel_label(mid) if file_key is not None else self._mouse_display.get(mid, mid)
            outer_tabs.add(frame, text=tab_label)
            self._frame_to_mid[str(frame)] = mid
        try:
            outer_tabs.select(self._mouse_frames[mids[0]])
            self._ensure_mouse_widget(mids[0])
        except Exception:
            pass
    else:
        frame = getattr(self, '_empty_file_frame', None)
        if frame is None or not bool(getattr(frame, 'winfo_exists', lambda: False)()):
            frame = ttk.Frame(outer_tabs, padding=18)
            self._empty_file_frame = frame
        for child in frame.winfo_children():
            try:child.destroy()
            except Exception:pass
        ttk.Label(frame, text=_fiberlyse_placeholder_text(self, file_key), justify='left', wraplength=650).pack(anchor='nw')
        outer_tabs.add(frame, text='File status')
        try:outer_tabs.select(frame)
        except Exception:pass
    # Surface the existing batch tools as first-class workspace tabs.
    if results:
        if self.compare_widget is None:
            self.compare_widget = BatchCompareTk(outer_tabs, self)
            if self._axis_label_fs_override is not None:self.compare_widget.set_axis_label_fontsize(self._axis_label_fs_override)
            if self._graph_title_fs_override is not None:self.compare_widget.set_graph_title_fontsize(self._graph_title_fs_override)
            if self._tick_label_fs_override is not None:self.compare_widget.set_tick_label_fontsize(self._tick_label_fs_override)
        if self.average_widget is None:
            self.average_widget = BatchAverageTk(outer_tabs, self)
            if self._axis_label_fs_override is not None:self.average_widget.set_axis_label_fontsize(self._axis_label_fs_override)
            if self._graph_title_fs_override is not None:self.average_widget.set_graph_title_fontsize(self._graph_title_fs_override)
            if self._tick_label_fs_override is not None:self.average_widget.set_tick_label_fontsize(self._tick_label_fs_override)
        try:self.compare_widget.refresh_available();self.compare_widget.refresh_plot()
        except Exception:pass
        try:self.average_widget.refresh_available();self.average_widget.refresh_plot()
        except Exception:pass
        try:outer_tabs.add(self.compare_widget, text='Batch compare')
        except Exception:pass
        try:outer_tabs.add(self.average_widget, text='Batch average')
        except Exception:pass

def _fiberlyse_on_analysis_finished(self, results: Dict[str, ChannelResult], display: Dict[str, str], order: List[str], mode: str, file_meta: List[Tuple[str, str, str]]):
    # Prefer the first successfully analyzed file if the selected file has no results.
    success_keys = []
    for mid in order:
        k = str(mid).split(':', 1)[0]
        if k not in success_keys:
            success_keys.append(k)
    if success_keys and getattr(self, '_active_file_key', None) not in success_keys:
        self._active_file_key = success_keys[0]
    self._results = results
    self._mouse_display = dict(display)
    self._mouse_order = list(order)
    self._set_result_file_choices(file_meta)
    self.build_tabs(norm_mode=mode)
    try:
        any_key = next(iter(results.keys()))
        r0 = results[any_key]
        self.status.set(f"Done. Files in queue={len(self.csv_paths)} | analyzed channels={len(results)} | Mode={mode}. Smooth win={r0.smooth_window}. Artifacts={('ON' if self.var_artifact_enabled.get() else 'OFF')}. Interp={('ON' if r0.use_interpolation else 'OFF')}.")
    except Exception:
        try:self.status.set('Done.')
        except Exception:pass
    _fiberlyse_set_buttons_for_analysis(self, False)
    try:self.on_norm_mode_changed()
    except Exception:pass
    _fiberlyse_refresh_queue_tree(self)

def _fiberlyse_show_file_number_map_dialog(self) -> None:
    _fiberlyse_sync_queue_from_paths(self)
    paths = list(getattr(self, 'csv_paths', []) or [])
    if not paths:
        messagebox.showinfo('File queue mapping (Ctrl+J)', "No CSV files selected.\n\nUse 'Add CSV(s)...' or 'Add folder...' first.")
        return
    top = tk.Toplevel(self.root)
    top.title('File queue mapping (Ctrl+J)')
    top.transient(self.root)
    top.grab_set()
    frm = ttk.Frame(top, padding=10)
    frm.pack(fill=tk.BOTH, expand=True)
    ttk.Label(frm, text='Friendly names, queue status, and source paths:', justify='left').pack(anchor='w')
    txt_frame = ttk.Frame(frm)
    txt_frame.pack(fill=tk.BOTH, expand=True, pady=(8, 8))
    txt = tk.Text(txt_frame, width=100, height=min(24, 4 + len(paths) * 2), wrap='none')
    sb = ttk.Scrollbar(txt_frame, orient='vertical', command=txt.yview)
    txt.configure(yscrollcommand=sb.set)
    txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    sb.pack(side=tk.RIGHT, fill=tk.Y)
    recs = _fiberlyse_existing_records(self)
    lines = []
    for i, p in enumerate(paths, start=1):
        rec = recs.get(p, {})
        lines.append(f"{i}. {rec.get('alias', _basename_no_ext(p))} [{rec.get('status', 'queued')}] channels={rec.get('channels', '')}\n    {p}\n")
    try:
        rows = _fiberlyse_shortcuts(self).as_rows()
        if rows:
            lines.append('\nShortcuts:\n')
            for _action, label, seqs in rows:
                lines.append(f'    {seqs}: {label}\n')
    except Exception:
        pass
    txt.insert('1.0', ''.join(lines).rstrip() + '\n')
    txt.configure(state='disabled')
    btn_row = ttk.Frame(frm)
    btn_row.pack(fill=tk.X)
    def copy_to_clipboard():
        try:
            self.root.clipboard_clear();self.root.clipboard_append(''.join(lines).rstrip() + '\n')
        except Exception:pass
    ttk.Button(btn_row, text='Copy', command=copy_to_clipboard).pack(side=tk.LEFT)
    ttk.Button(btn_row, text='Close', command=top.destroy).pack(side=tk.RIGHT)
    self.root.wait_window(top)

# Preserve existing extension wrappers, then add the modern shell last.
_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_MODERN_SHELL = MainAppTk.__init__
def _fiberlyse_mainapptk_init_with_modern_shell(self, initial_csvs: Optional[List[str]]=None, autorun: bool=False):
    full_initial = []
    try:
        full_initial = [_fiberlyse_norm_path(p) for p in (initial_csvs or []) if p]
    except Exception:
        full_initial = list(initial_csvs or []) if initial_csvs else []
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_MODERN_SHELL(self, initial_csvs=initial_csvs, autorun=autorun)
    try:setattr(self.root, '_fiberlyse_app', self)
    except Exception:pass
    try:
        self._active_file_key = getattr(self, '_active_file_key', None)
        self._analysis_cancel_requested = False
        if full_initial:
            self.csv_paths = []
            _fiberlyse_add_paths_to_queue(self, full_initial)
        else:
            _fiberlyse_sync_queue_from_paths(self)
        _fiberlyse_install_modern_shell(self)
        _fiberlyse_refresh_files_label(self)
    except Exception as e:
        try:
            print(f'FiberLyse modern shell failed to initialize: {e}', file=sys.stderr)
        except Exception:
            pass

# Install modern methods.
MainAppTk.__init__ = _fiberlyse_mainapptk_init_with_modern_shell
MainAppTk._update_files_label = _fiberlyse_refresh_files_label
MainAppTk.add_csvs = _fiberlyse_add_csvs
MainAppTk.add_folder_csvs = _fiberlyse_add_folder_csvs
MainAppTk.clear_csvs = _fiberlyse_clear_csvs
MainAppTk.run_analysis = _fiberlyse_run_analysis
MainAppTk.run_analysis_selected = _fiberlyse_run_analysis_selected
MainAppTk.retry_failed_files = _fiberlyse_retry_failed_files
MainAppTk.cancel_analysis = _fiberlyse_cancel_analysis
MainAppTk.remove_selected_files = _fiberlyse_remove_selected_files
MainAppTk.open_selected_folder = _fiberlyse_open_selected_folder
MainAppTk.show_shortcuts_dialog = _fiberlyse_show_shortcuts_dialog
MainAppTk._set_result_file_choices = _fiberlyse_set_result_file_choices
MainAppTk._get_selected_view_file_key = _fiberlyse_get_selected_view_file_key
MainAppTk.on_view_file_changed = _fiberlyse_on_view_file_changed
MainAppTk.build_tabs = _fiberlyse_build_tabs
MainAppTk._refresh_visible_tabs = _fiberlyse_refresh_visible_tabs
MainAppTk.on_analysis_finished = _fiberlyse_on_analysis_finished
MainAppTk._show_file_number_map_dialog = _fiberlyse_show_file_number_map_dialog
MainAppTk._on_queue_tree_select = _fiberlyse_on_queue_tree_select
MainAppTk._on_queue_tree_double_click = _fiberlyse_on_queue_tree_double_click
# ---- End FiberLyse modernization / scalable queue extension ----

# ---- FiberLyse responsive pixel layout extension ----
# Screen-aware pixel metrics are applied after the modern shell is built.  This
# keeps analysis, exports, shortcuts, and plot calculations unchanged while the
# Tk layout adapts to laptops, projectors, and small external displays.
def _fiberlyse_clamp(value: float, low: float, high: float) -> float:
    try:
        return max(float(low), min(float(high), float(value)))
    except Exception:
        return float(low)

def _fiberlyse_screen_size(root) -> Tuple[int, int]:
    try:
        root.update_idletasks()
    except Exception:
        pass
    try:
        sw = int(root.winfo_screenwidth())
        sh = int(root.winfo_screenheight())
    except Exception:
        sw, sh = 1200, 800
    if sw <= 0:sw = 1200
    if sh <= 0:sh = 800
    return sw, sh

def _fiberlyse_responsive_scale(root) -> float:
    sw, sh = _fiberlyse_screen_size(root)
    # 1200x760 is the design reference for the legacy desktop layout.  The
    # multiplier is intentionally bounded so fonts and controls remain usable.
    return _fiberlyse_clamp(min(sw / 1200.0, sh / 760.0), 0.58, 1.18)

def _fiberlyse_px(root, value: float, minimum: Optional[int]=None, maximum: Optional[int]=None) -> int:
    scale = _fiberlyse_responsive_scale(root)
    out = int(round(float(value) * scale))
    if minimum is not None:out = max(int(minimum), out)
    if maximum is not None:out = min(int(maximum), out)
    return max(1, out)

def _fiberlyse_find_label(parent, text: str):
    try:
        for child in parent.winfo_children():
            try:
                if str(child.cget('text')) == text:
                    return child
            except Exception:
                pass
    except Exception:
        pass
    return None

def _fiberlyse_set_toplevel_pixel_geometry(top, parent=None, ideal_w: int=900, ideal_h: int=650) -> None:
    try:
        root = parent.winfo_toplevel() if parent is not None and hasattr(parent, 'winfo_toplevel') else top
    except Exception:
        root = top
    sw, sh = _fiberlyse_screen_size(root)
    margin_x = 24 if sw <= 900 else 80
    margin_y = 48 if sh <= 700 else 100
    w = min(int(ideal_w), max(360, sw - margin_x))
    h = min(int(ideal_h), max(280, sh - margin_y))
    try:
        min_w = max(320, min(w, _fiberlyse_px(root, 520, minimum=320)))
        min_h = max(240, min(h, _fiberlyse_px(root, 360, minimum=240)))
        top.minsize(min_w, min_h)
    except Exception:
        pass
    try:
        x = max(0, int((sw - w) / 2))
        y = max(0, int((sh - h) / 2))
        top.geometry(f'{w}x{h}+{x}+{y}')
    except Exception:
        pass

def _fiberlyse_configure_common_listbox_pixels(root, listbox, width_px: int, height_px: int) -> None:
    try:
        listbox.configure(width=max(12, int(width_px / 8)), height=max(4, int(height_px / 21)))
    except Exception:
        pass

def _fiberlyse_make_top_toolbar_responsive(self, root_width: Optional[int]=None) -> None:
    root = self.root
    if root_width is None:
        try:root_width = int(root.winfo_width())
        except Exception:root_width = _fiberlyse_screen_size(root)[0]
    compact = root_width < 980
    very_compact = root_width < 760
    try:
        top = self.btn_run.master
    except Exception:
        return
    try:
        for col in range(8):top.columnconfigure(col, weight=0)
        top.columnconfigure(0, weight=1)
    except Exception:
        pass
    try:
        self.lbl_file.configure(wraplength=max(220, root_width - 32))
    except Exception:
        pass
    if compact:
        try:self.lbl_file.grid(row=0, column=0, columnspan=6, sticky='ew', padx=(0, 0), pady=(0, 4))
        except Exception:pass
        try:self.btn_add.grid(row=1, column=0, sticky='ew', padx=(0, 4), pady=2)
        except Exception:pass
        try:self.btn_clear.grid(row=1, column=1, sticky='ew', padx=(4, 4), pady=2)
        except Exception:pass
        try:self.btn_run.grid(row=1, column=2, sticky='ew', padx=(4, 4), pady=2)
        except Exception:pass
        try:
            self.progress_bar.configure(length=_fiberlyse_px(root, 120, minimum=80, maximum=180))
            self.progress_bar.grid(row=1, column=3, sticky='ew', padx=(4, 4), pady=2)
        except Exception:pass
        try:self.btn_cancel.grid(row=1, column=4, sticky='ew', padx=(4, 0), pady=2)
        except Exception:pass
        for c in range(5):
            try:top.columnconfigure(c, weight=1)
            except Exception:pass
    else:
        try:self.lbl_file.grid(row=0, column=0, sticky='ew', padx=(0, 8), pady=0)
        except Exception:pass
        try:self.btn_add.grid(row=0, column=1, sticky='e', padx=(0, 6), pady=0)
        except Exception:pass
        try:self.btn_clear.grid(row=0, column=2, sticky='e', padx=(0, 6), pady=0)
        except Exception:pass
        try:self.btn_run.grid(row=0, column=3, sticky='e', pady=0)
        except Exception:pass
        try:
            self.progress_bar.configure(length=_fiberlyse_px(root, 180, minimum=130, maximum=220))
            self.progress_bar.grid(row=0, column=4, sticky='e', padx=(14, 6), pady=0)
        except Exception:pass
        try:self.btn_cancel.grid(row=0, column=5, sticky='e', pady=0)
        except Exception:pass
    if very_compact:
        try:self.progress_bar.grid(row=2, column=0, columnspan=3, sticky='ew', padx=(0, 4), pady=(2, 0))
        except Exception:pass
        try:self.btn_cancel.grid(row=2, column=3, columnspan=2, sticky='ew', padx=(4, 0), pady=(2, 0))
        except Exception:pass

def _fiberlyse_make_analysis_controls_responsive(self, root_width: Optional[int]=None) -> None:
    root = self.root
    if root_width is None:
        try:root_width = int(root.winfo_width())
        except Exception:root_width = _fiberlyse_screen_size(root)[0]
    compact = root_width < 1360
    try:row2 = self.spin_factor.master
    except Exception:return
    lbl_factor = _fiberlyse_find_label(row2, 'Factor:')
    lbl_pad = _fiberlyse_find_label(row2, 'Pad (samples):')
    lbl_fps = _fiberlyse_find_label(row2, 'Acq FPS (Hz):')
    lbl_eff = _fiberlyse_find_label(row2, '(eff = FPS/2)')
    lbl_smooth = _fiberlyse_find_label(row2, 'Smooth win (samples):')
    try:
        for c in range(14):row2.columnconfigure(c, weight=0)
    except Exception:
        pass
    if compact:
        try:self.chk_artifact_enabled.grid(row=0, column=0, columnspan=6, sticky='w', padx=(0, 8), pady=(0, 2))
        except Exception:pass
        if lbl_factor is not None:
            try:lbl_factor.grid(row=1, column=0, sticky='w', padx=(0, 4), pady=2)
            except Exception:pass
        try:self.spin_factor.grid(row=1, column=1, padx=(0, 10), sticky='w', pady=2)
        except Exception:pass
        if lbl_pad is not None:
            try:lbl_pad.grid(row=1, column=2, sticky='w', padx=(0, 4), pady=2)
            except Exception:pass
        try:self.spin_pad.grid(row=1, column=3, padx=(0, 10), sticky='w', pady=2)
        except Exception:pass
        try:self.chk_shared.grid(row=1, column=4, columnspan=2, padx=(0, 0), sticky='w', pady=2)
        except Exception:pass
        if lbl_fps is not None:
            try:lbl_fps.grid(row=2, column=0, sticky='w', padx=(0, 4), pady=2)
            except Exception:pass
        try:self.spin_acq_fps.grid(row=2, column=1, padx=(0, 6), sticky='w', pady=2)
        except Exception:pass
        if lbl_eff is not None:
            try:lbl_eff.grid(row=2, column=2, sticky='w', padx=(0, 10), pady=2)
            except Exception:pass
        if lbl_smooth is not None:
            try:lbl_smooth.grid(row=2, column=3, sticky='w', padx=(0, 4), pady=2)
            except Exception:pass
        try:self.spin_smooth_win.grid(row=2, column=4, padx=(0, 10), sticky='w', pady=2)
        except Exception:pass
        try:self.chk_interp.grid(row=3, column=0, columnspan=6, sticky='w', pady=(2, 0))
        except Exception:pass
    else:
        try:self.chk_artifact_enabled.grid(row=0, column=0, columnspan=2, sticky='w', padx=(0, 12), pady=0)
        except Exception:pass
        if lbl_factor is not None:
            try:lbl_factor.grid(row=0, column=2, sticky='w', pady=0)
            except Exception:pass
        try:self.spin_factor.grid(row=0, column=3, padx=(6, 12), sticky='w', pady=0)
        except Exception:pass
        if lbl_pad is not None:
            try:lbl_pad.grid(row=0, column=4, sticky='w', pady=0)
            except Exception:pass
        try:self.spin_pad.grid(row=0, column=5, padx=(6, 12), sticky='w', pady=0)
        except Exception:pass
        try:self.chk_shared.grid(row=0, column=6, padx=(0, 12), sticky='w', pady=0)
        except Exception:pass
        if lbl_fps is not None:
            try:lbl_fps.grid(row=0, column=7, sticky='w', pady=0)
            except Exception:pass
        try:self.spin_acq_fps.grid(row=0, column=8, padx=(6, 6), sticky='w', pady=0)
        except Exception:pass
        if lbl_eff is not None:
            try:lbl_eff.grid(row=0, column=9, sticky='w', pady=0)
            except Exception:pass
        if lbl_smooth is not None:
            try:lbl_smooth.grid(row=0, column=10, sticky='w', pady=0)
            except Exception:pass
        try:self.spin_smooth_win.grid(row=0, column=11, padx=(6, 12), sticky='w', pady=0)
        except Exception:pass
        try:self.chk_interp.grid(row=1, column=0, columnspan=6, sticky='w', pady=(4, 0))
        except Exception:pass

def _fiberlyse_make_norm_controls_responsive(self, root_width: Optional[int]=None) -> None:
    root = self.root
    if root_width is None:
        try:root_width = int(root.winfo_width())
        except Exception:root_width = _fiberlyse_screen_size(root)[0]
    compact = root_width < 820
    try:row3 = self.cmb_norm.master
    except Exception:return
    lbl_norm = _fiberlyse_find_label(row3, 'Normalization view:')
    try:
        for c in range(10):row3.columnconfigure(c, weight=0)
    except Exception:
        pass
    if compact:
        if lbl_norm is not None:
            try:lbl_norm.grid(row=0, column=0, sticky='w', pady=2)
            except Exception:pass
        try:self.cmb_norm.grid(row=0, column=1, padx=(6, 0), sticky='ew', pady=2)
        except Exception:pass
        try:row3.columnconfigure(1, weight=1)
        except Exception:pass
        try:self.lbl_interval.grid(row=1, column=0, sticky='w', pady=2)
        except Exception:pass
        try:self.spin_interval_start.grid(row=1, column=1, padx=(6, 4), sticky='w', pady=2)
        except Exception:pass
        try:self.lbl_interval_to.grid(row=1, column=2, sticky='w', pady=2)
        except Exception:pass
        try:self.spin_interval_end.grid(row=1, column=3, padx=(4, 8), sticky='w', pady=2)
        except Exception:pass
        try:self.btn_apply_norm.grid(row=1, column=4, sticky='w', pady=2)
        except Exception:pass
    else:
        if lbl_norm is not None:
            try:lbl_norm.grid(row=0, column=0, sticky='w', pady=0)
            except Exception:pass
        try:self.cmb_norm.grid(row=0, column=1, padx=(6, 12), sticky='w', pady=0)
        except Exception:pass
        try:self.lbl_interval.grid(row=0, column=2, sticky='w', pady=0)
        except Exception:pass
        try:self.spin_interval_start.grid(row=0, column=3, padx=(6, 4), sticky='w', pady=0)
        except Exception:pass
        try:self.lbl_interval_to.grid(row=0, column=4, sticky='w', pady=0)
        except Exception:pass
        try:self.spin_interval_end.grid(row=0, column=5, padx=(4, 12), sticky='w', pady=0)
        except Exception:pass
        try:self.btn_apply_norm.grid(row=0, column=6, sticky='w', pady=0)
        except Exception:pass
    try:self.update_norm_controls_visibility()
    except Exception:pass

def _fiberlyse_apply_workspace_pixel_metrics(self, root_width: Optional[int]=None, root_height: Optional[int]=None) -> None:
    root = self.root
    if root_width is None:
        try:root_width = int(root.winfo_width())
        except Exception:root_width = _fiberlyse_screen_size(root)[0]
    if root_height is None:
        try:root_height = int(root.winfo_height())
        except Exception:root_height = _fiberlyse_screen_size(root)[1]
    qtree = getattr(self, 'queue_tree', None)
    try:
        ttk.Style(root).configure('Queue.Treeview', rowheight=_fiberlyse_px(root, 24, minimum=18, maximum=24))
    except Exception:
        pass
    qpanel_width = max(_fiberlyse_px(root, 220, minimum=180), min(_fiberlyse_px(root, 335, maximum=360), int(root_width * 0.34)))
    if root_width < 760:
        qpanel_width = max(170, int(root_width * 0.38))
    try:
        self.workspace_pane.sashpos(0, qpanel_width)
    except Exception:
        pass
    # On short displays, the queue panel uses a compact pixel layout: optional
    # explanatory labels are hidden first, while the actual queue and action
    # buttons remain reachable.
    try:
        for child in self.queue_panel.winfo_children():
            try:text = str(child.cget('text'))
            except Exception:text = ''
            if text.startswith('Queue, analyze') or text.startswith('Made by'):
                if root_height < 650:
                    try:child.grid_remove()
                    except Exception:pass
                else:
                    try:child.grid()
                    except Exception:pass
    except Exception:
        pass
    if qtree is not None:
        try:
            qtree.column('#0', width=max(78, int(qpanel_width * 0.42)), minwidth=66, stretch=True)
            qtree.column('status', width=max(58, int(qpanel_width * 0.20)), minwidth=54, stretch=False)
            qtree.column('channels', width=max(34, int(qpanel_width * 0.12)), minwidth=30, stretch=False, anchor='center')
            qtree.column('note', width=max(48, int(qpanel_width * 0.20)), minwidth=42, stretch=True)
            try:pane_h = int(self.workspace_pane.winfo_height())
            except Exception:pane_h = int(root_height * 0.55)
            reserved = 190 if root_height < 650 else 300
            qtree.configure(height=max(3, min(18, int((pane_h - reserved) / 24))))
        except Exception:
            pass
    try:
        self.queue_panel.configure(padding=(_fiberlyse_px(root, 8, minimum=4, maximum=10), _fiberlyse_px(root, 8, minimum=4, maximum=10)))
    except Exception:
        pass
    try:
        self.workspace_pane.pack_configure(padx=_fiberlyse_px(root, 8, minimum=2, maximum=8), pady=(0, _fiberlyse_px(root, 6, minimum=2, maximum=6)))
    except Exception:
        pass

def _fiberlyse_apply_responsive_root_geometry(self) -> None:
    root = self.root
    sw, sh = _fiberlyse_screen_size(root)
    margin_x = 24 if sw <= 900 else 72
    margin_y = 48 if sh <= 700 else 96
    target_w = min(max(520, int(sw * 0.94)), max(520, sw - margin_x))
    target_h = min(max(420, int(sh * 0.90)), max(360, sh - margin_y))
    try:
        cur = root.geometry().split('+', 1)[0]
        cur_w, cur_h = [int(v) for v in cur.split('x')[:2]]
    except Exception:
        cur_w, cur_h = (0, 0)
    min_w = max(480, min(980, sw - max(20, margin_x)))
    min_h = max(340, min(640, sh - max(40, margin_y)))
    try:root.minsize(min_w, min_h)
    except Exception:pass
    # Replace the legacy 1200x700 startup size only at startup.  Do not fight
    # the user if they later resize the window.
    if cur_w >= sw or cur_h >= sh or cur_w == 1200 or cur_h == 700:
        x = max(0, int((sw - target_w) / 2))
        y = max(0, int((sh - target_h) / 2))
        try:root.geometry(f'{target_w}x{target_h}+{x}+{y}')
        except Exception:pass

def _fiberlyse_apply_responsive_layout(self, force: bool=False) -> None:
    root = self.root
    try:
        root.update_idletasks()
    except Exception:
        pass
    try:rw = int(root.winfo_width())
    except Exception:rw = _fiberlyse_screen_size(root)[0]
    try:rh = int(root.winfo_height())
    except Exception:rh = _fiberlyse_screen_size(root)[1]
    previous = getattr(self, '_responsive_last_size', None)
    bucket = (rw // 40, rh // 40)
    if (not force) and previous == bucket:
        return
    self._responsive_last_size = bucket
    _fiberlyse_make_top_toolbar_responsive(self, rw)
    _fiberlyse_make_analysis_controls_responsive(self, rw)
    _fiberlyse_make_norm_controls_responsive(self, rw)
    _fiberlyse_apply_workspace_pixel_metrics(self, rw, rh)

def _fiberlyse_on_root_configure_responsive(self, event=None):
    try:
        if event is not None and event.widget is not self.root:
            return
    except Exception:
        pass
    try:
        if getattr(self, '_responsive_after_id', None):
            self.root.after_cancel(self._responsive_after_id)
    except Exception:
        pass
    try:
        self._responsive_after_id = self.root.after(80, lambda: _fiberlyse_apply_responsive_layout(self, force=False))
    except Exception:
        pass

def _fiberlyse_apply_responsive_to_batch_widget(widget) -> None:
    try:
        root = widget.winfo_toplevel()
    except Exception:
        return
    scale = _fiberlyse_responsive_scale(root)
    list_h_small = max(4, int(8 * scale))
    list_h_big = max(5, int(10 * scale))
    list_w = max(18, int(28 * scale))
    for name in ['lst_available', 'lst_selected', 'lst_a', 'lst_b']:
        lb = getattr(widget, name, None)
        if lb is not None:
            try:
                height = list_h_big if name in ('lst_available', 'lst_selected') else list_h_small
                lb.configure(width=list_w, height=height)
            except Exception:
                pass

def _fiberlyse_responsive_refresh_visible_tabs(self) -> None:
    _FIBERLYSE_ORIGINAL_REFRESH_VISIBLE_TABS_FOR_RESPONSIVE(self)
    try:
        if self.compare_widget is not None:_fiberlyse_apply_responsive_to_batch_widget(self.compare_widget)
        if self.average_widget is not None:_fiberlyse_apply_responsive_to_batch_widget(self.average_widget)
        _fiberlyse_apply_responsive_layout(self, force=True)
    except Exception:
        pass

def _fiberlyse_responsive_set_result_file_choices(self, file_meta: List[Tuple[str, str, str]]) -> None:
    _FIBERLYSE_ORIGINAL_SET_RESULT_FILE_CHOICES_FOR_RESPONSIVE(self, file_meta)
    try:
        root = self.root
        width_chars = max(10, min(34, int(_fiberlyse_px(root, 16, minimum=10, maximum=34))))
        self.cmb_view_file.configure(width=width_chars)
    except Exception:
        pass

def _fiberlyse_plottab_init_responsive(self, master, tab_name: str, default_filename_prefix: str='', figsize: Tuple[float, float]=(7.2, 4.6), dpi: int=110):
    try:
        root = master.winfo_toplevel()
        scale = _fiberlyse_responsive_scale(root)
        # Keep the same aspect ratio, but reduce the initial Matplotlib request
        # size on small displays so the canvas can fit inside the available pane.
        figsize = (max(4.1, float(figsize[0]) * scale), max(2.8, float(figsize[1]) * scale))
    except Exception:
        pass
    _FIBERLYSE_ORIGINAL_PLOTTAB_INIT_FOR_RESPONSIVE(self, master, tab_name, default_filename_prefix, figsize, dpi)
    try:
        root = self.winfo_toplevel()
        self.canvas_widget.configure(width=_fiberlyse_px(root, 720, minimum=360, maximum=900), height=_fiberlyse_px(root, 440, minimum=260, maximum=580))
    except Exception:
        pass

_FIBERLYSE_ORIGINAL_PLOTTAB_INIT_FOR_RESPONSIVE = PlotTabTk.__init__
PlotTabTk.__init__ = _fiberlyse_plottab_init_responsive
_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_RESPONSIVE = MainAppTk.__init__
def _fiberlyse_mainapptk_init_with_responsive_pixels(self, initial_csvs: Optional[List[str]]=None, autorun: bool=False):
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_RESPONSIVE(self, initial_csvs=initial_csvs, autorun=autorun)
    try:
        _fiberlyse_apply_responsive_root_geometry(self)
        _fiberlyse_apply_responsive_layout(self, force=True)
        self.root.bind('<Configure>', lambda e, app=self: _fiberlyse_on_root_configure_responsive(app, e), add='+')
        self.root.after(250, lambda app=self: _fiberlyse_apply_responsive_layout(app, force=True))
    except Exception as e:
        try:print(f'FiberLyse responsive layout failed to initialize: {e}', file=sys.stderr)
        except Exception:pass

MainAppTk.__init__ = _fiberlyse_mainapptk_init_with_responsive_pixels
_FIBERLYSE_ORIGINAL_REFRESH_VISIBLE_TABS_FOR_RESPONSIVE = MainAppTk._refresh_visible_tabs
MainAppTk._refresh_visible_tabs = _fiberlyse_responsive_refresh_visible_tabs
_FIBERLYSE_ORIGINAL_SET_RESULT_FILE_CHOICES_FOR_RESPONSIVE = MainAppTk._set_result_file_choices
MainAppTk._set_result_file_choices = _fiberlyse_responsive_set_result_file_choices
MainAppTk._apply_responsive_layout = _fiberlyse_apply_responsive_layout
# ---- End FiberLyse responsive pixel layout extension ----



# ---- FiberLyse deep-research performance implementation ----
# Implements the high-priority recommendations from the accompanying
# performance review without changing the scientific calculations or saved data.
# The changes below deliberately patch the existing single-file application at
# the end of the module, matching the extension style already used above.
try:
    from concurrent.futures import ThreadPoolExecutor as _FiberlyseThreadPoolExecutor
    from queue import Queue as _FiberlyseQueue, Empty as _FiberlyseEmpty
    import itertools as _fiberlyse_perf_itertools
    import traceback as _fiberlyse_perf_traceback
    from dataclasses import replace as _fiberlyse_dataclass_replace
    from matplotlib.axes import Axes as _FiberlyseAxes
except Exception:
    _FiberlyseThreadPoolExecutor = None
    _FiberlyseQueue = None
    _FiberlyseEmpty = Exception
    _fiberlyse_perf_itertools = None
    _fiberlyse_perf_traceback = None
    _fiberlyse_dataclass_replace = None
    _FiberlyseAxes = None

FIBERLYSE_PERFORMANCE_PATCH = True
FIBERLYSE_PERFORMANCE_PATCH_VERSION = '2026-06-15-deep-research'

# -- Safe result cloning -----------------------------------------------------
def _fiberlyse_perf_clone_channel_result(res: ChannelResult) -> ChannelResult:
    """Create a shallow dataclass clone whose large raw arrays are shared.

    Worker jobs assign new derived arrays to the clone, so the GUI-owned result
    is not mutated from a worker thread.  Raw arrays are intentionally shared to
    avoid doubling memory for every pending settings job.
    """
    if _fiberlyse_dataclass_replace is not None:
        try:
            clone = _fiberlyse_dataclass_replace(res)
        except Exception:
            clone = None
    else:
        clone = None
    if clone is None:
        try:
            import copy as _fiberlyse_perf_copy
            clone = _fiberlyse_perf_copy.copy(res)
        except Exception:
            clone = res
    try:
        setattr(clone, '_smooth_cache', {})
    except Exception:
        pass
    try:
        setattr(clone, '_freq_cache', {})
    except Exception:
        pass
    try:
        setattr(clone, '_data_version', int(getattr(res, '_data_version', 0)))
    except Exception:
        pass
    try:
        setattr(clone, '_pipeline_signature', getattr(res, '_pipeline_signature', None))
    except Exception:
        pass
    return clone

def _fiberlyse_perf_pipeline_signature(settings: Dict[str, Any]) -> Tuple[Any, ...]:
    return (
        bool(settings.get('artifact_enabled', True)),
        round(float(settings.get('artifact_factor', 11.9)), 8),
        int(settings.get('artifact_pad', 1)),
        bool(settings.get('require_shared', True)),
        str(settings.get('align_mode', DEFAULT_ALIGN_MODE)),
        bool(settings.get('use_interp', DEFAULT_USE_LINEAR_INTERP)),
    )

def _fiberlyse_perf_read_settings(self) -> Optional[Dict[str, Any]]:
    try:
        mode = self._read_norm_mode()
        interval_start = float(self.var_interval_start.get())
        interval_end = float(self.var_interval_end.get())
        return {
            'artifact_enabled': bool(self.var_artifact_enabled.get()),
            'artifact_factor': float(self.var_factor.get()),
            'artifact_pad': int(float(self.var_pad.get())),
            'require_shared': bool(self.var_shared.get()),
            'align_mode': DEFAULT_ALIGN_MODE,
            'acq_fps': np.nan,  # V24.1: analysis rate comes from timestamps only
            'smooth_win': max(1, int(float(self.var_smooth_win.get()))),
            'mode': mode,
            'interval_start': interval_start,
            'interval_end': interval_end,
            'use_interp': bool(self.var_interp.get()),
        }
    except Exception as e:
        try:
            messagebox.showerror('Invalid settings', f'Could not parse settings:\n\n{e}')
        except Exception:
            pass
        return None

# -- Single coordinated worker queue ---------------------------------------
def _fiberlyse_perf_init_background_jobs(self) -> None:
    if _FiberlyseThreadPoolExecutor is None or _FiberlyseQueue is None:
        return
    if getattr(self, '_fiberlyse_perf_executor', None) is not None:
        return
    self._fiberlyse_perf_executor = _FiberlyseThreadPoolExecutor(max_workers=1, thread_name_prefix='fiberlyse-settings')
    self._fiberlyse_perf_job_queue = _FiberlyseQueue()
    self._fiberlyse_perf_job_seq = _fiberlyse_perf_itertools.count(1) if _fiberlyse_perf_itertools is not None else iter(range(1, 10**12))
    self._fiberlyse_perf_latest_by_kind: Dict[str, int] = {}
    self._fiberlyse_perf_pipeline_after_id = None
    self._fiberlyse_result_source_epoch = int(getattr(self, '_fiberlyse_result_source_epoch', 0))
    try:
        self.root.after(25, self._fiberlyse_perf_poll_job_queue)
    except Exception:
        pass
    try:
        original_close = self.root.protocol('WM_DELETE_WINDOW')
    except Exception:
        original_close = None
    def _close_with_executor_shutdown():
        try:
            executor = getattr(self, '_fiberlyse_perf_executor', None)
            if executor is not None:
                executor.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
        try:
            if original_close and isinstance(original_close, str):
                # Tk returns a Tcl command name for an existing protocol.  Calling
                # it directly from Python is not reliable, so fall through to destroy.
                pass
        except Exception:
            pass
        try:
            self.root.destroy()
        except Exception:
            pass
    try:
        self.root.protocol('WM_DELETE_WINDOW', _close_with_executor_shutdown)
    except Exception:
        pass

def _fiberlyse_perf_submit_job(self, status_text: str, work_fn, apply_fn, kind: str='pipeline') -> Optional[int]:
    _fiberlyse_perf_init_background_jobs(self)
    executor = getattr(self, '_fiberlyse_perf_executor', None)
    queue_obj = getattr(self, '_fiberlyse_perf_job_queue', None)
    seq = getattr(self, '_fiberlyse_perf_job_seq', None)
    if executor is None or queue_obj is None or seq is None:
        # Fallback: run in a plain daemon thread and marshal back with after().
        def _thread_target():
            try:
                result = work_fn()
                self.root.after(0, lambda: apply_fn(result))
            except Exception:
                tb = _fiberlyse_perf_traceback.format_exc() if _fiberlyse_perf_traceback is not None else 'Background job failed.'
                self.root.after(0, lambda: self.on_analysis_failed(tb))
        threading.Thread(target=_thread_target, daemon=True).start()
        return None
    try:
        job_id = int(next(seq))
    except Exception:
        job_id = int(np.random.randint(1, 2**31 - 1))
    latest = getattr(self, '_fiberlyse_perf_latest_by_kind', None)
    if not isinstance(latest, dict):
        latest = {}
        self._fiberlyse_perf_latest_by_kind = latest
    latest[str(kind)] = job_id
    try:
        self.status.set(str(status_text))
    except Exception:
        pass
    future = executor.submit(work_fn)
    def _done(fut, job_id=job_id, kind=str(kind), apply_fn=apply_fn):
        try:
            queue_obj.put((job_id, kind, fut, apply_fn))
        except Exception:
            pass
    future.add_done_callback(_done)
    return job_id

def _fiberlyse_perf_poll_job_queue(self) -> None:
    queue_obj = getattr(self, '_fiberlyse_perf_job_queue', None)
    if queue_obj is None:
        return
    try:
        while True:
            job_id, kind, fut, apply_fn = queue_obj.get_nowait()
            latest = getattr(self, '_fiberlyse_perf_latest_by_kind', {})
            if isinstance(latest, dict) and latest.get(str(kind)) != int(job_id):
                continue
            try:
                result = fut.result()
            except Exception:
                tb = _fiberlyse_perf_traceback.format_exc() if _fiberlyse_perf_traceback is not None else 'Background job failed.'
                try:
                    messagebox.showerror('Background update failed', tb)
                except Exception:
                    pass
                try:
                    self.status.set('Background update failed.')
                except Exception:
                    pass
                continue
            latest = getattr(self, '_fiberlyse_perf_latest_by_kind', {})
            if isinstance(latest, dict) and latest.get(str(kind)) != int(job_id):
                continue
            try:
                apply_fn(result)
            except Exception:
                tb = _fiberlyse_perf_traceback.format_exc() if _fiberlyse_perf_traceback is not None else 'Background apply failed.'
                try:
                    messagebox.showerror('Background update failed', tb)
                except Exception:
                    pass
                try:
                    self.status.set('Background update failed.')
                except Exception:
                    pass
    except _FiberlyseEmpty:
        pass
    except Exception:
        pass
    finally:
        try:
            self.root.after(25, self._fiberlyse_perf_poll_job_queue)
        except Exception:
            pass

def _fiberlyse_perf_clear_average_cache(self) -> None:
    try:
        if self.average_widget is not None and isinstance(getattr(self.average_widget, '_stats_cache', None), dict):
            self.average_widget._stats_cache.clear()
    except Exception:
        pass

def _fiberlyse_perf_apply_pipeline_result(self, payload: Dict[str, Any]) -> None:
    if not isinstance(payload, dict):
        return
    expected_epoch = int(payload.get('epoch', -1))
    if expected_epoch != int(getattr(self, '_fiberlyse_result_source_epoch', 0)):
        try:self.status.set('Skipped stale settings update because analysis results changed.')
        except Exception:pass
        return
    updates = payload.get('updates', {})
    settings = payload.get('settings', {})
    reason = str(payload.get('reason', 'settings'))
    if not isinstance(updates, dict) or not updates:
        return
    if not getattr(self, '_results', None):
        return
    current_keys = set(getattr(self, '_results', {}).keys())
    for mid, res in updates.items():
        if mid in current_keys:
            self._results[mid] = res
    mode = str(settings.get('mode', self._read_norm_mode()))
    _fiberlyse_perf_clear_average_cache(self)
    for mid, widget in list(getattr(self, '_channel_widgets', {}).items()):
        try:
            if mid in self._results:
                widget.res = self._results[mid]
                widget.set_norm_mode(mode)
        except Exception:
            try:widget.refresh_after_pipeline_change()
            except Exception:pass
    try:
        if self.compare_widget is not None:
            self.compare_widget.refresh_plot()
    except Exception:
        pass
    try:
        if self.average_widget is not None:
            self.average_widget.refresh_plot()
    except Exception:
        pass
    try:
        interp_txt = 'ON' if bool(settings.get('use_interp')) else 'OFF'
        art_txt = 'ON' if bool(settings.get('artifact_enabled')) else 'OFF'
        invalid_interval = mode == NORM_ZF_INTERVAL and all(not np.any(np.isfinite(getattr(r, 'zF_interval', np.array([])))) for r in self._results.values())
        if invalid_interval:self.status.set(f'Updated {reason}, but the selected zF interval has fewer than 2 usable points; interval zF is unavailable.')
        else:self.status.set(f'Updated {reason} in background. Mode={mode}; artifacts={art_txt}; interpolation={interp_txt}.')
    except Exception:
        pass

def _fiberlyse_perf_schedule_pipeline_update(self, reason: str='settings', delay_ms: int=120) -> None:
    if not getattr(self, '_results', None):
        return
    try:
        old = getattr(self, '_fiberlyse_perf_pipeline_after_id', None)
        if old is not None:
            self.root.after_cancel(old)
    except Exception:
        pass
    def _submit_now():
        try:
            self._fiberlyse_perf_pipeline_after_id = None
        except Exception:
            pass
        settings = _fiberlyse_perf_read_settings(self)
        if settings is None:
            return
        snapshot = list((getattr(self, '_results', {}) or {}).items())
        epoch = int(getattr(self, '_fiberlyse_result_source_epoch', 0))
        pipeline_sig = _fiberlyse_perf_pipeline_signature(settings)
        def _work():
            out: Dict[str, ChannelResult] = {}
            for mid, res in snapshot:
                clone = _fiberlyse_perf_clone_channel_result(res)
                try:
                    prior_sig = getattr(clone, '_pipeline_signature', None)
                except Exception:
                    prior_sig = None
                if prior_sig != pipeline_sig:
                    recompute_artifact_pipeline_inplace(
                        clone,
                        artifact_enabled=bool(settings['artifact_enabled']),
                        artifact_factor=float(settings['artifact_factor']),
                        artifact_pad=int(settings['artifact_pad']),
                        require_shared=bool(settings['require_shared']),
                        align_mode=str(settings['align_mode']),
                        use_linear_interp=bool(settings['use_interp']),
                    )
                    try:setattr(clone, '_pipeline_signature', pipeline_sig)
                    except Exception:pass
                else:
                    # Keep the active path synchronized even if the scientific
                    # pipeline is unchanged.  This is cheap and avoids stale UI.
                    try:
                        if bool(getattr(clone, 'use_interpolation', False)) != bool(settings['use_interp']):
                            set_interpolation_mode(clone, bool(settings['use_interp']))
                    except Exception:
                        pass
                try:
                    clone.smooth_window = max(1, int(settings['smooth_win']))
                    clone.acq_fps_hz = float(settings['acq_fps']) if np.isfinite(settings['acq_fps']) and float(settings['acq_fps']) > 0 else clone.acq_fps_hz
                    clone.zf_interval_start_s = float(settings['interval_start'])
                    clone.zf_interval_end_s = float(settings['interval_end'])
                    recompute_normalizations(clone)
                except Exception:
                    pass
                out[mid] = clone
            return {'updates': out, 'settings': dict(settings), 'reason': reason, 'epoch': epoch}
        self._fiberlyse_perf_submit_job(f'Updating {reason} in background...', _work, self._fiberlyse_perf_apply_pipeline_result, kind='pipeline')
    try:
        self._fiberlyse_perf_pipeline_after_id = self.root.after(int(delay_ms), _submit_now)
    except Exception:
        _submit_now()

# -- Patch user actions that used to recompute on the Tk event thread --------
def _fiberlyse_perf_on_artifact_enabled_toggled(self):
    enabled = bool(self.var_artifact_enabled.get())
    try:self._sync_artifact_controls_state()
    except Exception:pass
    if not getattr(self, '_results', None):
        try:self.status.set(f"Artifact remover {('ENABLED' if enabled else 'DISABLED')} (will apply on next Run).")
        except Exception:pass
        return
    _fiberlyse_perf_schedule_pipeline_update(self, reason='artifact settings')

def _fiberlyse_perf_on_interp_toggled(self):
    if not getattr(self, '_results', None):
        return
    _fiberlyse_perf_schedule_pipeline_update(self, reason='interpolation setting')

def _fiberlyse_perf_apply_normalization(self):
    mode = self._read_norm_mode()
    settings = _fiberlyse_perf_read_settings(self)
    if settings is None:
        return
    if not getattr(self, '_results', None):
        if mode == NORM_DFF:
            self.status.set('Normalization view: DFF')
        elif mode == NORM_ZF_GLOBAL:
            self.status.set('Normalization view: zF (global, GUI)')
        elif mode == NORM_ZF_INTERVAL:
            self.status.set('Normalization view: zF - interval based (set interval and apply)')
        return
    _fiberlyse_perf_schedule_pipeline_update(self, reason='normalization settings')

# -- Coalesced redraw path ---------------------------------------------------
def _fiberlyse_perf_flush_redraw(self):
    try:self._fiberlyse_redraw_after_id = None
    except Exception:pass
    try:self._pending_redraw_flags = set()
    except Exception:pass
    # Apply lightweight overlays and user customizations exactly once, then ask
    # TkAgg for a single idle draw.
    try:self._apply_user_overrides()
    except Exception:pass
    try:self._apply_time_markers()
    except Exception:pass
    try:self._apply_axis_ranges(reset_auto=False)
    except Exception:pass
    try:self._apply_axis_tick_intervals(reset_auto=False)
    except Exception:pass
    try:self._apply_legend_box_overrides()
    except Exception:pass
    try:self._apply_user_overrides()
    except Exception:pass
    try:self.canvas.draw_idle()
    except Exception:pass

def _fiberlyse_perf_schedule_redraw(self, *flags):
    try:
        pending = getattr(self, '_pending_redraw_flags', None)
        if not isinstance(pending, set):
            pending = set()
            self._pending_redraw_flags = pending
        pending.update(flags or {'full'})
    except Exception:
        pass
    if getattr(self, '_fiberlyse_redraw_after_id', None) is not None:
        return
    try:
        self._fiberlyse_redraw_after_id = self.after(16, lambda tab=self: _fiberlyse_perf_flush_redraw(tab))
    except Exception:
        _fiberlyse_perf_flush_redraw(self)

def _fiberlyse_perf_legend_draw_freeze(self, event=None):
    # After the first real draw, cache the resolved location of expensive
    # legend(loc='best') placements.  Future legend() calls on the same axes use
    # the cached anchor and avoid the dense-data best-location scan.
    try:
        renderer = getattr(event, 'renderer', None)
    except Exception:
        renderer = None
    for ax in list(getattr(self.fig, 'axes', [])):
        try:
            leg = ax.get_legend()
        except Exception:
            leg = None
        if leg is None:
            continue
        try:
            if getattr(ax, '_fiberlyse_cached_legend_anchor_axes', None) is not None:
                continue
        except Exception:
            pass
        try:
            loc = getattr(leg, '_loc', None)
            if loc not in (0, 'best'):
                continue
        except Exception:
            continue
        try:
            anchor = _fiberlyse_legend_current_anchor_axes(ax, leg, renderer)
            if np.all(np.isfinite(np.asarray(anchor, dtype=float))):
                setattr(ax, '_fiberlyse_cached_legend_anchor_axes', (float(anchor[0]), float(anchor[1])))
        except Exception:
            pass

_FIBERLYSE_ORIGINAL_PLOTTAB_INIT_FOR_PERF_PATCH = PlotTabTk.__init__
def _fiberlyse_perf_plottab_init(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_PLOTTAB_INIT_FOR_PERF_PATCH(self, *args, **kwargs)
    try:self._fiberlyse_redraw_after_id = None
    except Exception:pass
    try:self._pending_redraw_flags = set()
    except Exception:pass
    try:self._cid_perf_legend_freeze = self.canvas.mpl_connect('draw_event', lambda event, tab=self: _fiberlyse_perf_legend_draw_freeze(tab, event))
    except Exception:pass

PlotTabTk.__init__ = _fiberlyse_perf_plottab_init
PlotTabTk.schedule_redraw = _fiberlyse_perf_schedule_redraw
PlotTabTk._flush_redraw = _fiberlyse_perf_flush_redraw
PlotTabTk.redraw = _fiberlyse_perf_schedule_redraw

# -- Legend placement cache --------------------------------------------------
if _FiberlyseAxes is not None and not getattr(_FiberlyseAxes.legend, '_fiberlyse_perf_cached', False):
    _FIBERLYSE_ORIGINAL_AXES_LEGEND_FOR_PERF = _FiberlyseAxes.legend
    def _fiberlyse_perf_axes_legend(self, *args, **kwargs):
        try:
            loc = kwargs.get('loc', None)
            if loc == 'best' and 'bbox_to_anchor' not in kwargs:
                anchor = getattr(self, '_fiberlyse_cached_legend_anchor_axes', None)
                if anchor is not None:
                    kwargs = dict(kwargs)
                    kwargs['loc'] = 'center'
                    kwargs['bbox_to_anchor'] = (float(anchor[0]), float(anchor[1]))
                    kwargs['bbox_transform'] = self.transAxes
        except Exception:
            pass
        return _FIBERLYSE_ORIGINAL_AXES_LEGEND_FOR_PERF(self, *args, **kwargs)
    _fiberlyse_perf_axes_legend._fiberlyse_perf_cached = True
    _FiberlyseAxes.legend = _fiberlyse_perf_axes_legend

_FIBERLYSE_ORIGINAL_APPLY_LEGEND_BOX_OVERRIDES_FOR_PERF = PlotTabTk._apply_legend_box_overrides
def _fiberlyse_perf_apply_legend_box_overrides(self):
    _FIBERLYSE_ORIGINAL_APPLY_LEGEND_BOX_OVERRIDES_FOR_PERF(self)
    states = getattr(self, '_legend_box_overrides', None)
    if not isinstance(states, dict):
        states = {}
    for ax in list(getattr(self.fig, 'axes', [])):
        try:key = self._ax_key(ax)
        except Exception:key = id(ax)
        state = states.get(key, {}) if isinstance(states.get(key, {}), dict) else {}
        if bool(state.get('manual_position')) and state.get('anchor_axes', None) is not None:
            try:
                a = state.get('anchor_axes')
                setattr(ax, '_fiberlyse_cached_legend_anchor_axes', (float(a[0]), float(a[1])))
            except Exception:
                pass
        elif state.get('manual_position') is False:
            try:delattr(ax, '_fiberlyse_cached_legend_anchor_axes')
            except Exception:pass
PlotTabTk._apply_legend_box_overrides = _fiberlyse_perf_apply_legend_box_overrides

# -- Frequency-band cache ----------------------------------------------------
def _fiberlyse_perf_get_bandpass_cached(res: ChannelResult, low_hz: float, high_hz: float, fs: float, order: int=DEFAULT_BUTTER_ORDER) -> np.ndarray:
    try:version = int(getattr(res, '_data_version', 0))
    except Exception:version = 0
    try:fs_key = round(float(fs), 12)
    except Exception:fs_key = np.nan
    method = 'butter_sos_segmentwise' if _HAVE_SCIPY_SIGNAL else 'fft'
    key = (method, version, float(low_hz), float(high_hz), fs_key, int(order))
    cache = getattr(res, '_freq_cache', None)
    if not isinstance(cache, dict):
        cache = {}
        try:setattr(res, '_freq_cache', cache)
        except Exception:pass
    if key not in cache:
        y = np.asarray(res.dFF_nointerp, dtype=float)
        if _HAVE_SCIPY_SIGNAL:
            out = bandpass_butterworth_segmentwise_no_interp(y, low_hz, high_hz, fs, order=order)
        else:
            out = bandpass_fft_no_interp(y, low_hz, high_hz, fs)
        if len(cache) > 24:
            cache.clear()
        cache[key] = out
    return cache[key]

def _fiberlyse_perf_draw_frequency(self):
    fig = self.tab_freq.fig
    fig.clear()
    axs = fig.subplots(3, 2)
    axes = np.asarray(axs).flatten()
    t = np.asarray(self.res.t_exc, dtype=float)
    dff = np.asarray(self.res.dFF_nointerp, dtype=float)
    fs = float(getattr(self.res, 'eff_fs_hz', np.nan))
    acq = float(getattr(self.res, 'acq_fps_hz', 0.0))
    if not np.isfinite(fs) or fs <= 0:
        fs = estimate_fs_from_t(t)
    have_holes = np.any(~np.isfinite(dff))
    if _HAVE_SCIPY_SIGNAL:
        method = 'Butterworth (segment-wise, no interpolation; cached)'
    else:
        method = 'FFT cached (only if no NaNs)' if not have_holes else 'Unavailable (SciPy signal required for NaN holes)'
    supt = f'{self.res.gcol} - Band-limited DFF (NO interpolation)\nUses non-interpolated DFF; holes preserved. Acq FPS={acq:.3g} Hz -> eff fs={fs:.3g} Hz | {method}'
    fig.suptitle(supt, y=0.995, fontsize=11)
    for i, (low_hz, high_hz) in enumerate(FREQ_BANDS):
        ax = axes[i]
        ax.clear()
        label = f'{low_hz:g}-{high_hz:g} Hz'
        if not np.any(np.isfinite(dff)):
            ax.text(0.5, 0.5, 'DFF is all NaN', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(label);ax.set_xlabel('Time (s)');ax.set_ylabel('DFF')
            continue
        y_band = _fiberlyse_perf_get_bandpass_cached(self.res, low_hz, high_hz, fs, order=DEFAULT_BUTTER_ORDER)
        if not np.any(np.isfinite(y_band)):
            ax.text(0.5, 0.5, 'Band unavailable\n(check SciPy / NaNs / Nyquist)', ha='center', va='center', transform=ax.transAxes)
        else:
            ax.plot(t, y_band, label=label, linewidth=1.0)
            ax.legend(loc='best')
        ax.set_title(label);ax.set_xlabel('Time (s)');ax.set_ylabel('DFF')
    for j in range(len(FREQ_BANDS), len(axes)):
        axes[j].axis('off')
    fig.tight_layout(rect=(0, 0.0, 1, 0.93))
    self.tab_freq.redraw()

def _fiberlyse_perf_export_frequency(self) -> Dict[str, pd.DataFrame]:
    t = np.asarray(self.res.t_exc, dtype=float)
    dff = np.asarray(self.res.dFF_nointerp, dtype=float)
    fs = float(getattr(self.res, 'eff_fs_hz', np.nan))
    acq = float(getattr(self.res, 'acq_fps_hz', 0.0))
    if not np.isfinite(fs) or fs <= 0:
        fs = estimate_fs_from_t(t)
    out = {'t_s': t, 'dFF_nointerp': dff}
    for low_hz, high_hz in FREQ_BANDS:
        col = f'band_{low_hz:g}_{high_hz:g}_Hz'
        out[col] = _fiberlyse_perf_get_bandpass_cached(self.res, low_hz, high_hz, fs, order=DEFAULT_BUTTER_ORDER)
    df = pd.DataFrame(out)
    meta = self._meta_df()
    extra = pd.DataFrame([('freq_acq_fps_hz', acq), ('freq_eff_fs_hz', fs), ('method', 'butter_sos_segmentwise_cached' if _HAVE_SCIPY_SIGNAL else 'fft_cached')], columns=['key', 'value'])
    meta2 = pd.concat([meta, extra], ignore_index=True)
    return {'freq_bands': df, 'meta': meta2}

ChannelTabsTk._draw_frequency = _fiberlyse_perf_draw_frequency
ChannelTabsTk._export_frequency = _fiberlyse_perf_export_frequency

# -- Batch average statistics cache -----------------------------------------
_FIBERLYSE_ORIGINAL_BATCH_GROUP_STATS_FOR_PERF = BatchAverageTk._group_stats
def _fiberlyse_perf_group_stats_cached(self, ids: List[str], mode: str, smooth_win: int, do_smooth: bool):
    signature = []
    for mid in list(ids):
        res = self.app.get_mouse_result(mid)
        if res is None:
            continue
        signature.append((
            str(mid),
            int(getattr(res, '_data_version', 0)),
            float(getattr(res, 'zf_interval_start_s', DEFAULT_ZF_INTERVAL_START_S)),
            float(getattr(res, 'zf_interval_end_s', DEFAULT_ZF_INTERVAL_END_S)),
        ))
    key = (tuple(signature), str(mode), int(smooth_win), bool(do_smooth), float(BATCH_HIDE_INITIAL_SECONDS))
    cache = getattr(self, '_stats_cache', None)
    if not isinstance(cache, dict):
        cache = {}
        self._stats_cache = cache
    if key not in cache:
        if len(cache) > 64:
            cache.clear()
        cache[key] = _FIBERLYSE_ORIGINAL_BATCH_GROUP_STATS_FOR_PERF(self, ids, mode, smooth_win, do_smooth)
    return cache[key]
BatchAverageTk._group_stats = _fiberlyse_perf_group_stats_cached

# -- Incremental updates for normalization tabs -----------------------------
def _fiberlyse_perf_draw_norm(self):
    ax = self.tab_norm.ax
    y, ylabel, label, title = self._get_norm_series()
    line = getattr(self, '_perf_norm_line', None)
    if line is not None and getattr(line, 'axes', None) is ax:
        try:
            line.set_data(self.res.t_exc, y)
            line.set_label(label)
            ax.set_title(title);ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel)
            ax.relim();ax.autoscale_view()
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            ax.legend(loc='best')
            self.tab_norm.redraw()
            return
        except Exception:
            pass
    ax.clear()
    self._perf_norm_line, = ax.plot(self.res.t_exc, y, label=label)
    ax.set_title(title);ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel);ax.legend(loc='best')
    self.tab_norm.redraw()

def _fiberlyse_perf_draw_norm_smooth(self):
    ax = self.tab_norm_smooth.ax
    y, ylabel, label, title = self._get_norm_series()
    w = int(getattr(self.res, 'smooth_window', DEFAULT_SMOOTH_WINDOW))
    y_s = get_smoothed_norm_array(self.res, self.norm_mode, window_size=w)
    raw_line = getattr(self, '_perf_norm_smooth_raw_line', None)
    smooth_line = getattr(self, '_perf_norm_smooth_line', None)
    if raw_line is not None and smooth_line is not None and getattr(raw_line, 'axes', None) is ax and getattr(smooth_line, 'axes', None) is ax:
        try:
            raw_line.set_data(self.res.t_exc, y)
            raw_line.set_label(f'{label} (raw)')
            smooth_line.set_data(self.res.t_exc, y_s)
            smooth_line.set_label(f'{ylabel} smoothed (win={w})')
            ax.set_title(f'{title} - smoothed');ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel)
            ax.relim();ax.autoscale_view()
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            ax.legend(loc='best')
        except Exception:
            ax.clear();raw_line = None;smooth_line = None
    if raw_line is None or smooth_line is None:
        ax.clear()
        self._perf_norm_smooth_raw_line, = ax.plot(self.res.t_exc, y, alpha=0.25, linewidth=1.0, zorder=1, label=f'{label} (raw)')
        self._perf_norm_smooth_line, = ax.plot(self.res.t_exc, y_s, linewidth=1.4, zorder=3, color=self._perf_norm_smooth_raw_line.get_color(), label=f'{ylabel} smoothed (win={w})')
        ax.set_title(f'{title} - smoothed');ax.set_xlabel('Time (s)');ax.set_ylabel(ylabel);ax.legend(loc='best')
    if self.parent_app is not None and hasattr(self.parent_app, '_results'):
        try:
            all_results = self.parent_app._results
            global_y_min = np.inf;global_y_max = -np.inf
            for _mid, res_other in all_results.items():
                y_s_other = get_smoothed_norm_array(res_other, self.norm_mode, window_size=w)
                finite_vals = y_s_other[np.isfinite(y_s_other)]
                if len(finite_vals) > 0:
                    global_y_min = min(global_y_min, np.nanmin(finite_vals));global_y_max = max(global_y_max, np.nanmax(finite_vals))
            if np.isfinite(global_y_min) and np.isfinite(global_y_max):
                padding = (global_y_max - global_y_min) * 0.05
                ax.set_ylim(global_y_min - padding, global_y_max + padding)
        except Exception:
            pass
    self.tab_norm_smooth.redraw()

ChannelTabsTk._draw_norm = _fiberlyse_perf_draw_norm
ChannelTabsTk._draw_norm_smooth = _fiberlyse_perf_draw_norm_smooth

# -- Asynchronous Excel writing ---------------------------------------------
_FIBERLYSE_ORIGINAL_EXPORT_EXCEL_FOR_PERF = PlotTabTk.export_excel
def _fiberlyse_perf_export_excel(self):
    suggested = f'{self.default_filename_prefix}_{self.tab_name}'.strip('_')
    path = filedialog.asksaveasfilename(title='Export data (Excel)', defaultextension='.xlsx', initialfile=f'{suggested}.xlsx', filetypes=[('Excel workbook', '*.xlsx'), ('All files', '*.*')])
    if not path:
        return
    try:
        if callable(self.export_provider):
            payload = self.export_provider()
        else:
            payload = self._payload_from_artists_fallback()
        if not isinstance(payload, dict) or not payload:
            raise ValueError('Export provider returned no data.')
    except Exception as e:
        messagebox.showerror('Export failed', f'Could not prepare Excel data:\n\n{e}')
        return
    app = None
    try:
        app = getattr(self.winfo_toplevel(), '_fiberlyse_app', None)
    except Exception:
        app = None
    def _work():
        self._write_xlsx_minimal(path, payload)
        return path
    def _apply(done_path):
        try:messagebox.showinfo('Export complete', f'Saved Excel file:\n{done_path}')
        except Exception:pass
        try:
            if app is not None and getattr(app, 'status', None) is not None:
                app.status.set(f'Export complete: {os.path.basename(done_path)}')
        except Exception:
            pass
    if app is not None and hasattr(app, '_fiberlyse_perf_submit_job'):
        app._fiberlyse_perf_submit_job('Exporting Excel workbook in background...', _work, _apply, kind=f'export:{id(self)}')
    else:
        def _thread_target():
            try:
                done = _work()
                self.after(0, lambda: _apply(done))
            except Exception as exc:
                self.after(0, lambda exc=exc: messagebox.showerror('Export failed', f'Could not export Excel data:\n\n{exc}'))
        threading.Thread(target=_thread_target, daemon=True).start()
PlotTabTk.export_excel = _fiberlyse_perf_export_excel

# -- Debounced queue filter --------------------------------------------------
_FIBERLYSE_QUEUE_TREE_REFRESH_IMMEDIATE_FOR_PERF = _fiberlyse_refresh_queue_tree
def _fiberlyse_perf_debounced_queue_filter(self):
    try:
        old = getattr(self, '_fiberlyse_queue_filter_after_id', None)
        if old is not None:
            self.root.after_cancel(old)
    except Exception:
        pass
    def _run():
        try:self._fiberlyse_queue_filter_after_id = None
        except Exception:pass
        _FIBERLYSE_QUEUE_TREE_REFRESH_IMMEDIATE_FOR_PERF(self)
    try:
        self._fiberlyse_queue_filter_after_id = self.root.after(150, _run)
    except Exception:
        _run()

# -- MainApp initialization and generation bookkeeping -----------------------
_FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_PERF_PATCH = MainAppTk.__init__
def _fiberlyse_perf_mainapp_init(self, *args, **kwargs):
    _FIBERLYSE_ORIGINAL_MAINAPPTK_INIT_FOR_PERF_PATCH(self, *args, **kwargs)
    try:
        self._fiberlyse_result_source_epoch = int(getattr(self, '_fiberlyse_result_source_epoch', 0))
        _fiberlyse_perf_init_background_jobs(self)
        self._fiberlyse_perf_submit_job = lambda status_text, work_fn, apply_fn, kind='pipeline': _fiberlyse_perf_submit_job(self, status_text, work_fn, apply_fn, kind=kind)
        self._fiberlyse_perf_poll_job_queue = lambda: _fiberlyse_perf_poll_job_queue(self)
        self._fiberlyse_perf_apply_pipeline_result = lambda payload: _fiberlyse_perf_apply_pipeline_result(self, payload)
        # Replace the immediate queue-search trace installed by the modern shell
        # with a debounced refresh.  Other direct queue refreshes remain immediate.
        var = getattr(self, 'var_queue_filter', None)
        if var is not None:
            try:
                for info in var.trace_info():
                    if info and str(info[0]).find('write') >= 0:
                        try:var.trace_remove('write', info[1])
                        except Exception:pass
            except Exception:
                pass
            try:var.trace_add('write', lambda *_args, app=self: _fiberlyse_perf_debounced_queue_filter(app))
            except Exception:pass
    except Exception as e:
        try:print(f'FiberLyse performance coordinator failed to initialize: {e}', file=sys.stderr)
        except Exception:pass

MainAppTk.__init__ = _fiberlyse_perf_mainapp_init
MainAppTk._fiberlyse_perf_init_background_jobs = _fiberlyse_perf_init_background_jobs
MainAppTk._fiberlyse_perf_submit_job = _fiberlyse_perf_submit_job
MainAppTk._fiberlyse_perf_poll_job_queue = _fiberlyse_perf_poll_job_queue
MainAppTk._fiberlyse_perf_apply_pipeline_result = _fiberlyse_perf_apply_pipeline_result
MainAppTk._fiberlyse_perf_schedule_pipeline_update = _fiberlyse_perf_schedule_pipeline_update
MainAppTk.on_artifact_enabled_toggled = _fiberlyse_perf_on_artifact_enabled_toggled
MainAppTk.on_interp_toggled = _fiberlyse_perf_on_interp_toggled
MainAppTk.apply_normalization = _fiberlyse_perf_apply_normalization

_FIBERLYSE_ORIGINAL_ON_ANALYSIS_FINISHED_FOR_PERF = MainAppTk.on_analysis_finished
def _fiberlyse_perf_on_analysis_finished(self, *args, **kwargs):
    try:self._fiberlyse_result_source_epoch = int(getattr(self, '_fiberlyse_result_source_epoch', 0)) + 1
    except Exception:pass
    return _FIBERLYSE_ORIGINAL_ON_ANALYSIS_FINISHED_FOR_PERF(self, *args, **kwargs)
MainAppTk.on_analysis_finished = _fiberlyse_perf_on_analysis_finished

_FIBERLYSE_ORIGINAL_CLEAR_CSVS_FOR_PERF = MainAppTk.clear_csvs
def _fiberlyse_perf_clear_csvs(self, *args, **kwargs):
    try:self._fiberlyse_result_source_epoch = int(getattr(self, '_fiberlyse_result_source_epoch', 0)) + 1
    except Exception:pass
    return _FIBERLYSE_ORIGINAL_CLEAR_CSVS_FOR_PERF(self, *args, **kwargs)
MainAppTk.clear_csvs = _fiberlyse_perf_clear_csvs

_FIBERLYSE_ORIGINAL_REMOVE_SELECTED_FILES_FOR_PERF = MainAppTk.remove_selected_files
def _fiberlyse_perf_remove_selected_files(self, *args, **kwargs):
    try:self._fiberlyse_result_source_epoch = int(getattr(self, '_fiberlyse_result_source_epoch', 0)) + 1
    except Exception:pass
    return _FIBERLYSE_ORIGINAL_REMOVE_SELECTED_FILES_FOR_PERF(self, *args, **kwargs)
MainAppTk.remove_selected_files = _fiberlyse_perf_remove_selected_files
# ---- End FiberLyse deep-research performance implementation ----

# ---- FiberLyse V22 stabilization bindings ---------------------------------
_FIBERLYSE_V22_PRE_BIND_INIT = MainAppTk.__init__
def _fiberlyse_v22_mainapp_init(self, *args, **kwargs):
    _FIBERLYSE_V22_PRE_BIND_INIT(self, *args, **kwargs)
    def schedule_artifact(_event=None):
        if getattr(self, '_results', None):
            try:self._fiberlyse_perf_schedule_pipeline_update(reason='artifact settings')
            except Exception:pass
    def schedule_smoothing(_event=None):
        if getattr(self, '_results', None):
            try:self._fiberlyse_perf_schedule_pipeline_update(reason='smoothing setting')
            except Exception:pass
    def schedule_acq(_event=None):
        if getattr(self, '_results', None):
            try:self._fiberlyse_perf_schedule_pipeline_update(reason='acquisition metadata')
            except Exception:pass
    for widget in [getattr(self, 'spin_factor', None), getattr(self, 'spin_pad', None)]:
        if widget is not None:
            try:widget.configure(command=schedule_artifact)
            except Exception:pass
            try:widget.bind('<Return>', schedule_artifact, add='+');widget.bind('<FocusOut>', schedule_artifact, add='+')
            except Exception:pass
    try:self.chk_shared.configure(command=schedule_artifact)
    except Exception:pass
    widget = getattr(self, 'spin_smooth_win', None)
    if widget is not None:
        try:widget.configure(command=schedule_smoothing)
        except Exception:pass
        try:widget.bind('<Return>', schedule_smoothing, add='+');widget.bind('<FocusOut>', schedule_smoothing, add='+')
        except Exception:pass
    widget = getattr(self, 'spin_acq_fps', None)
    if widget is not None:
        try:widget.configure(command=schedule_acq)
        except Exception:pass
        try:widget.bind('<Return>', schedule_acq, add='+');widget.bind('<FocusOut>', schedule_acq, add='+')
        except Exception:pass
MainAppTk.__init__ = _fiberlyse_v22_mainapp_init

# V22 invariant: frequency filtering uses measured per-channel sample timing.
# The user-entered acquisition FPS is retained/exported as metadata and fallback
# only if timestamps are insufficient.


# ---- FiberLyse V23 feature extension --------------------------------------
# V23 deliberately leaves the V22 scientific pipeline unchanged.  This layer
# adds GUI night mode, a centralized multi-file/multi-graph export workflow,
# and configurable Ctrl+I line/interval event annotations.
FIBERLYSE_VERSION = 'V23'

_FIBERLYSE_V23_THEME_LIGHT = {
    'bg': '#f6f7f9',
    'panel': '#ffffff',
    'fg': '#20252b',
    'muted': '#5f6b7a',
    'entry': '#ffffff',
    'border': '#c7ccd3',
    'select': '#dbeafe',
    'select_fg': '#111827',
    'tooltip_bg': '#222222',
    'tooltip_fg': '#ffffff',
    'plot_bg': '#ffffff',
    'plot_fg': '#20252b',
    'grid': '#d8dde4',
}
_FIBERLYSE_V23_THEME_DARK = {
    'bg': '#171a1f',
    'panel': '#22262d',
    'fg': '#e6edf3',
    'muted': '#a7b0bb',
    'entry': '#2b3038',
    'border': '#454c57',
    'select': '#344a67',
    'select_fg': '#f3f6f9',
    'tooltip_bg': '#0f1115',
    'tooltip_fg': '#f3f6f9',
    'plot_bg': '#20242b',
    'plot_fg': '#e6edf3',
    'grid': '#454c57',
}

def _fiberlyse_v23_theme_tokens(app=None, widget=None) -> Dict[str, str]:
    if app is None and widget is not None:
        try:app = getattr(widget.winfo_toplevel(), '_fiberlyse_app', None)
        except Exception:app = None
    dark = False
    try:dark = bool(app.var_night_mode.get())
    except Exception:
        try:dark = bool(getattr(app, '_fiberlyse_night_mode', False))
        except Exception:dark = False
    return dict(_FIBERLYSE_V23_THEME_DARK if dark else _FIBERLYSE_V23_THEME_LIGHT)

def _fiberlyse_v23_apply_ttk_theme(app) -> None:
    tok = _fiberlyse_v23_theme_tokens(app=app)
    try:
        style = ttk.Style(app.root)
        try:
            if 'clam' in style.theme_names():style.theme_use('clam')
        except Exception:pass
        app.root.configure(background=tok['bg'])
        style.configure('TFrame', background=tok['bg'])
        style.configure('Fiber.TFrame', background=tok['bg'])
        style.configure('Panel.TFrame', background=tok['panel'], relief='flat')
        style.configure('TLabel', background=tok['bg'], foreground=tok['fg'])
        style.configure('Panel.TLabel', background=tok['panel'], foreground=tok['fg'])
        style.configure('Title.TLabel', background=tok['panel'], foreground=tok['fg'], font=('TkDefaultFont', 11, 'bold'))
        style.configure('Hint.TLabel', background=tok['panel'], foreground=tok['muted'])
        style.configure('Tooltip.TLabel', background=tok['tooltip_bg'], foreground=tok['tooltip_fg'], relief='solid', borderwidth=1)
        style.configure('TButton', background=tok['panel'], foreground=tok['fg'], padding=(7, 3))
        style.map('TButton', background=[('active', tok['select']), ('pressed', tok['select'])], foreground=[('disabled', tok['muted'])])
        style.configure('Accent.TButton', background=tok['select'], foreground=tok['fg'], padding=(10, 4))
        style.map('Accent.TButton', background=[('active', tok['select']), ('pressed', tok['select'])])
        style.configure('TCheckbutton', background=tok['bg'], foreground=tok['fg'])
        style.map('TCheckbutton', background=[('active', tok['bg'])], foreground=[('disabled', tok['muted'])])
        style.configure('TRadiobutton', background=tok['bg'], foreground=tok['fg'])
        style.map('TRadiobutton', background=[('active', tok['bg'])], foreground=[('disabled', tok['muted'])])
        for sty in ['TEntry', 'TSpinbox', 'TCombobox']:
            style.configure(sty, fieldbackground=tok['entry'], background=tok['entry'], foreground=tok['fg'], arrowcolor=tok['fg'])
            style.map(sty, fieldbackground=[('readonly', tok['entry']), ('disabled', tok['panel'])], foreground=[('readonly', tok['fg']), ('disabled', tok['muted'])])
        style.configure('TNotebook', background=tok['bg'], borderwidth=0)
        style.configure('TNotebook.Tab', background=tok['panel'], foreground=tok['fg'], padding=(12, 5))
        style.map('TNotebook.Tab', background=[('selected', tok['select']), ('active', tok['entry'])], foreground=[('selected', tok['select_fg'])])
        style.configure('TLabelframe', background=tok['bg'], foreground=tok['fg'])
        style.configure('TLabelframe.Label', background=tok['bg'], foreground=tok['fg'], font=('TkDefaultFont', 10, 'bold'))
        style.configure('Queue.Treeview', rowheight=24, fieldbackground=tok['panel'], background=tok['panel'], foreground=tok['fg'])
        style.map('Queue.Treeview', background=[('selected', tok['select'])], foreground=[('selected', tok['select_fg'])])
        style.configure('Queue.Treeview.Heading', background=tok['entry'], foreground=tok['fg'], font=('TkDefaultFont', 9, 'bold'))
        style.map('Queue.Treeview.Heading', background=[('active', tok['select'])])
        style.configure('Treeview', fieldbackground=tok['panel'], background=tok['panel'], foreground=tok['fg'])
        style.map('Treeview', background=[('selected', tok['select'])], foreground=[('selected', tok['select_fg'])])
        style.configure('Horizontal.TProgressbar', background=tok['select'], troughcolor=tok['entry'])
    except Exception as e:
        try:print(f'V23 ttk theme apply failed: {e}', file=sys.stderr)
        except Exception:pass

def _fiberlyse_v23_apply_classic_widget_theme(widget, tok: Dict[str, str]) -> None:
    """Theme Tk widgets that do not obey ttk.Style (Text/Listbox/Canvas/toolbars)."""
    try:
        cls = str(widget.winfo_class())
    except Exception:
        cls = ''
    try:
        if cls in ('Listbox',):
            widget.configure(background=tok['panel'], foreground=tok['fg'], selectbackground=tok['select'], selectforeground=tok['select_fg'], highlightbackground=tok['border'], highlightcolor=tok['select'])
        elif cls in ('Text',):
            widget.configure(background=tok['panel'], foreground=tok['fg'], insertbackground=tok['fg'], selectbackground=tok['select'], selectforeground=tok['select_fg'], highlightbackground=tok['border'])
        elif cls in ('Canvas',):
            widget.configure(background=tok['panel'], highlightbackground=tok['border'])
        elif cls in ('Frame', 'Labelframe'):
            widget.configure(background=tok['panel'])
        elif cls in ('Label',):
            widget.configure(background=tok['panel'], foreground=tok['fg'])
        elif cls in ('Button', 'Checkbutton', 'Radiobutton'):
            widget.configure(background=tok['panel'], foreground=tok['fg'], activebackground=tok['select'], activeforeground=tok['select_fg'], highlightbackground=tok['border'])
        elif cls in ('Entry', 'Spinbox'):
            widget.configure(background=tok['entry'], foreground=tok['fg'], insertbackground=tok['fg'], highlightbackground=tok['border'])
    except Exception:
        pass
    try:
        for child in widget.winfo_children():
            _fiberlyse_v23_apply_classic_widget_theme(child, tok)
    except Exception:
        pass

def _fiberlyse_v23_apply_plot_theme(tab) -> None:
    app = None
    try:app = getattr(tab.winfo_toplevel(), '_fiberlyse_app', None)
    except Exception:app = None
    tok = _fiberlyse_v23_theme_tokens(app=app, widget=tab)
    try:tab.fig.patch.set_facecolor(tok['plot_bg'])
    except Exception:pass
    for ax in list(getattr(getattr(tab, 'fig', None), 'axes', []) or []):
        try:ax.set_facecolor(tok['plot_bg'])
        except Exception:pass
        try:
            for spine in ax.spines.values():spine.set_color(tok['plot_fg'])
        except Exception:pass
        try:ax.tick_params(axis='both', colors=tok['plot_fg'])
        except Exception:pass
        try:ax.xaxis.label.set_color(tok['plot_fg']);ax.yaxis.label.set_color(tok['plot_fg']);ax.title.set_color(tok['plot_fg'])
        except Exception:pass
        try:
            for txt in ax.texts:
                # Avoid forcing explicit user color overrides on data annotations.
                if mcolors.to_hex(txt.get_color(), keep_alpha=False) in ('#000000', '#20252b', '#e6edf3', '#ffffff'):
                    txt.set_color(tok['plot_fg'])
                try:
                    box = txt.get_bbox_patch()
                    if box is not None:
                        box.set_facecolor(tok['panel']);box.set_edgecolor(tok['border'])
                except Exception:pass
        except Exception:pass
        try:
            leg = ax.get_legend()
            if leg is not None:
                leg.get_frame().set_facecolor(tok['panel']);leg.get_frame().set_edgecolor(tok['border'])
                leg.get_frame().set_alpha(0.92)
                for txt in leg.get_texts():txt.set_color(tok['plot_fg'])
                try:leg.get_title().set_color(tok['plot_fg'])
                except Exception:pass
        except Exception:pass
    try:
        st = getattr(tab.fig, '_suptitle', None)
        if st is not None:st.set_color(tok['plot_fg'])
    except Exception:pass
    try:_fiberlyse_v23_apply_classic_widget_theme(tab.toolbar, tok)
    except Exception:pass
    try:tab.canvas_widget.configure(background=tok['plot_bg'], highlightbackground=tok['border'])
    except Exception:pass

def _fiberlyse_v23_iter_plot_tabs(app):
    seen = set()
    for cw in list(getattr(app, '_channel_widgets', {}).values()):
        for name in ['tab_raw', 'tab_art', 'tab_fit', 'tab_norm', 'tab_norm_smooth', 'tab_freq']:
            tab = getattr(cw, name, None)
            if tab is not None and id(tab) not in seen:
                seen.add(id(tab));yield tab
    for obj in [getattr(app, 'compare_widget', None), getattr(app, 'average_widget', None)]:
        tab = getattr(obj, 'plot', None) if obj is not None else None
        if tab is not None and id(tab) not in seen:
            seen.add(id(tab));yield tab

def _fiberlyse_v23_refresh_theme(app) -> None:
    try:app._fiberlyse_night_mode = bool(app.var_night_mode.get())
    except Exception:pass
    _fiberlyse_v23_apply_ttk_theme(app)
    tok = _fiberlyse_v23_theme_tokens(app=app)
    try:_fiberlyse_v23_apply_classic_widget_theme(app.root, tok)
    except Exception:pass
    try:
        if getattr(app, 'queue_tree', None) is not None:
            # Keep state colors legible in both modes.
            colors = [('queued', '#a7b0bb' if bool(app.var_night_mode.get()) else '#5f6b7a'), ('analyzing', '#f0b35a' if bool(app.var_night_mode.get()) else '#8a5a00'), ('analyzed', '#65d69a' if bool(app.var_night_mode.get()) else '#0f6b3f'), ('failed', '#ff7b7b' if bool(app.var_night_mode.get()) else '#a12020'), ('canceled', '#e5a85a' if bool(app.var_night_mode.get()) else '#7a4b00')]
            for tag, fg in colors:app.queue_tree.tag_configure(tag, foreground=fg)
    except Exception:pass
    for tab in _fiberlyse_v23_iter_plot_tabs(app):
        try:_fiberlyse_v23_apply_plot_theme(tab);tab.canvas.draw_idle()
        except Exception:pass
    try:app.status.set('Night mode enabled.' if app.var_night_mode.get() else 'Light mode enabled.')
    except Exception:pass

def _fiberlyse_v23_toggle_night_mode(app) -> None:
    _fiberlyse_v23_refresh_theme(app)

# Apply plot theme every time a plot schedules a redraw; plot-building functions
# frequently call ax.clear(), which resets face/text colors.
_FIBERLYSE_V23_REDRAW_BEFORE_THEME = PlotTabTk.redraw
def _fiberlyse_v23_redraw(self):
    try:_fiberlyse_v23_apply_plot_theme(self)
    except Exception:pass
    return _FIBERLYSE_V23_REDRAW_BEFORE_THEME(self)
PlotTabTk.redraw = _fiberlyse_v23_redraw

# ---- Configurable event annotations ---------------------------------------
_FIBERLYSE_V23_EVENT_LINESTYLES = {
    'Solid': '-',
    'Dashed': '--',
    'Dotted': ':',
    'Dash-dot': '-.',
}

def _fiberlyse_v23_event_color(value: str, fallback: str='#444444') -> str:
    try:return PlotTabTk._normalize_hex_color(value)
    except Exception:return fallback

def _fiberlyse_v23_apply_event_annotations(self) -> None:
    self._clear_time_marker_artists()
    marks = list(getattr(self, '_time_markers', []) or [])
    # Backwards compatibility with V22 marker dictionaries is intentionally kept.
    if not marks:
        for ax in list(getattr(self.fig, 'axes', [])):
            try:
                leg = ax.get_legend()
                if leg is not None:
                    leg.remove();ax.legend(loc='best') if ax.get_legend_handles_labels()[0] else None
            except Exception:pass
        return
    for ax in list(getattr(self.fig, 'axes', [])):
        if not self._axis_looks_like_time(ax):continue
        added_any = False
        for idx, m in enumerate(marks):
            kind = str(m.get('kind', 'line')).lower().strip()
            label = str(m.get('label') or '').strip()
            color = _fiberlyse_v23_event_color(str(m.get('color') or '#444444'))
            try:lw = max(0.25, float(m.get('linewidth', 1.4)))
            except Exception:lw = 1.4
            ls_name = str(m.get('linestyle', 'Dotted'))
            ls = _FIBERLYSE_V23_EVENT_LINESTYLES.get(ls_name, ':')
            if kind == 'interval':
                try:a = float(m.get('start'));b = float(m.get('end'))
                except Exception:continue
                if not np.isfinite(a) or not np.isfinite(b) or a == b:continue
                lo, hi = sorted([a, b])
                if not label:label = f'{lo:g}–{hi:g}s'
                try:alpha = min(0.95, max(0.02, float(m.get('alpha', 0.20))))
                except Exception:alpha = 0.20
                edge = bool(m.get('show_edges', True))
                try:
                    art = ax.axvspan(lo, hi, facecolor=color, edgecolor=(color if edge else 'none'), linewidth=(lw if edge else 0.0), linestyle=(ls if edge else '-'), alpha=alpha, label=label, zorder=0.8)
                    self._time_marker_artists.append(art);added_any = True
                except Exception:pass
            else:
                try:x = float(m.get('x'))
                except Exception:continue
                if not np.isfinite(x):continue
                if not label:label = f't={x:g}s'
                try:
                    art = ax.axvline(x, label=label, color=color, linestyle=ls, linewidth=lw, alpha=0.90, zorder=9)
                    self._time_marker_artists.append(art);added_any = True
                except Exception:pass
        if added_any:
            try:
                leg = ax.get_legend()
                if leg is not None:leg.remove()
                ax.legend(loc='best')
            except Exception:pass

PlotTabTk._apply_time_markers = _fiberlyse_v23_apply_event_annotations

def _fiberlyse_v23_add_event_annotation(self, event: Dict[str, Any]) -> None:
    cur = list(getattr(self, '_time_markers', []) or [])
    if len(cur) >= 4:
        try:messagebox.showwarning('Event annotations', 'You can add up to 4 event annotations per plot. Remove the last annotation with Ctrl+I then Ctrl+Backspace.')
        except Exception:pass
        return
    clean = dict(event or {})
    kind = str(clean.get('kind', 'line')).lower().strip()
    if kind not in ('line', 'interval'):kind = 'line'
    clean['kind'] = kind
    clean['color'] = _fiberlyse_v23_event_color(str(clean.get('color') or '#444444'))
    clean['linestyle'] = str(clean.get('linestyle') or 'Dotted')
    try:clean['linewidth'] = max(0.25, float(clean.get('linewidth', 1.4)))
    except Exception:clean['linewidth'] = 1.4
    if kind == 'interval':
        clean['start'] = float(clean['start']);clean['end'] = float(clean['end'])
        clean['alpha'] = min(0.95, max(0.02, float(clean.get('alpha', 0.20))))
        clean['show_edges'] = bool(clean.get('show_edges', True))
    else:
        clean['x'] = float(clean['x'])
    cur.append(clean);self._time_markers = cur;self.redraw()
PlotTabTk.add_event_annotation = _fiberlyse_v23_add_event_annotation

# Keep the legacy API for code paths that still call add_time_marker().
def _fiberlyse_v23_add_time_marker(self, x_s: float, label: Optional[str]=None) -> None:
    x = float(x_s)
    self.add_event_annotation({'kind':'line', 'x':x, 'label':(str(label).strip() if label else f't={x:g}s'), 'color':'#444444', 'linestyle':'Dotted', 'linewidth':1.4})
PlotTabTk.add_time_marker = _fiberlyse_v23_add_time_marker

def _fiberlyse_v23_show_event_dialog(app, tab: 'PlotTabTk') -> None:
    top = tk.Toplevel(app.root);top.title('Event annotation (Ctrl+I)');top.transient(app.root);top.grab_set()
    try:_fiberlyse_set_toplevel_pixel_geometry(top, app.root, 570, 560)
    except Exception:pass
    frm = ttk.Frame(top, padding=12);frm.pack(fill=tk.BOTH, expand=True);frm.columnconfigure(1, weight=1)
    ttk.Label(frm, text='Add an event to the currently active plot', font=('TkDefaultFont', 11, 'bold')).grid(row=0, column=0, columnspan=3, sticky='w')
    ttk.Label(frm, text='Choose a single vertical line or a shaded time interval. Appearance is stored with this plot.', justify='left').grid(row=1, column=0, columnspan=3, sticky='w', pady=(2, 10))
    kind_var = tk.StringVar(value=getattr(app, '_last_event_kind', 'Line'))
    ttk.Label(frm, text='Event type:').grid(row=2, column=0, sticky='w')
    kind_box = ttk.Combobox(frm, textvariable=kind_var, values=('Line', 'Interval'), state='readonly', width=15);kind_box.grid(row=2, column=1, sticky='w', pady=3)
    label_var = tk.StringVar(value='')
    ttk.Label(frm, text='Legend label:').grid(row=3, column=0, sticky='w');ttk.Entry(frm, textvariable=label_var, width=32).grid(row=3, column=1, columnspan=2, sticky='ew', pady=3)
    line_frame = ttk.LabelFrame(frm, text='Timing');line_frame.grid(row=4, column=0, columnspan=3, sticky='ew', pady=(8, 6));line_frame.columnconfigure(1, weight=1);line_frame.columnconfigure(3, weight=1)
    time_var = tk.StringVar(value=getattr(app, '_last_marker_time_str', ''))
    start_var = tk.StringVar(value=getattr(app, '_last_interval_start_str', ''))
    end_var = tk.StringVar(value=getattr(app, '_last_interval_end_str', ''))
    lbl_time = ttk.Label(line_frame, text='Time (s):');ent_time = ttk.Entry(line_frame, textvariable=time_var, width=16)
    lbl_start = ttk.Label(line_frame, text='Start (s):');ent_start = ttk.Entry(line_frame, textvariable=start_var, width=16)
    lbl_end = ttk.Label(line_frame, text='End (s):');ent_end = ttk.Entry(line_frame, textvariable=end_var, width=16)
    style_frame = ttk.LabelFrame(frm, text='Appearance');style_frame.grid(row=5, column=0, columnspan=3, sticky='ew', pady=(6, 6));style_frame.columnconfigure(1, weight=1)
    color_var = tk.StringVar(value=getattr(app, '_last_event_color', '#6b7280'))
    linestyle_var = tk.StringVar(value=getattr(app, '_last_event_linestyle', 'Dotted'))
    width_var = tk.StringVar(value=getattr(app, '_last_event_linewidth', '1.5'))
    alpha_var = tk.StringVar(value=getattr(app, '_last_event_alpha', '0.20'))
    edges_var = tk.BooleanVar(value=bool(getattr(app, '_last_event_edges', True)))
    ttk.Label(style_frame, text='Color (hex):').grid(row=0, column=0, sticky='w', padx=6, pady=4);ttk.Entry(style_frame, textvariable=color_var, width=14).grid(row=0, column=1, sticky='w', pady=4)
    ttk.Label(style_frame, text='Line style:').grid(row=1, column=0, sticky='w', padx=6, pady=4);ttk.Combobox(style_frame, textvariable=linestyle_var, values=tuple(_FIBERLYSE_V23_EVENT_LINESTYLES.keys()), state='readonly', width=14).grid(row=1, column=1, sticky='w', pady=4)
    ttk.Label(style_frame, text='Line width:').grid(row=2, column=0, sticky='w', padx=6, pady=4);ttk.Entry(style_frame, textvariable=width_var, width=14).grid(row=2, column=1, sticky='w', pady=4)
    lbl_alpha = ttk.Label(style_frame, text='Interval opacity:');ent_alpha = ttk.Entry(style_frame, textvariable=alpha_var, width=14);chk_edges = ttk.Checkbutton(style_frame, text='Show interval boundary lines', variable=edges_var)
    help_lbl = ttk.Label(frm, text='Ctrl+Backspace removes the most recently added annotation. You can still rename legend text by double-clicking it.', justify='left');help_lbl.grid(row=6, column=0, columnspan=3, sticky='w', pady=(6, 10))
    btn_row = ttk.Frame(frm);btn_row.grid(row=7, column=0, columnspan=3, sticky='ew');btn_row.columnconfigure(0, weight=1)
    def update_kind(*_):
        for w in [lbl_time, ent_time, lbl_start, ent_start, lbl_end, ent_end, lbl_alpha, ent_alpha, chk_edges]:
            try:w.grid_forget()
            except Exception:pass
        if kind_var.get() == 'Interval':
            lbl_start.grid(row=0, column=0, sticky='w', padx=6, pady=5);ent_start.grid(row=0, column=1, sticky='ew', padx=(0, 8), pady=5)
            lbl_end.grid(row=0, column=2, sticky='w', padx=6, pady=5);ent_end.grid(row=0, column=3, sticky='ew', padx=(0, 6), pady=5)
            lbl_alpha.grid(row=3, column=0, sticky='w', padx=6, pady=4);ent_alpha.grid(row=3, column=1, sticky='w', pady=4);chk_edges.grid(row=4, column=0, columnspan=2, sticky='w', padx=6, pady=4)
        else:
            lbl_time.grid(row=0, column=0, sticky='w', padx=6, pady=5);ent_time.grid(row=0, column=1, sticky='ew', padx=(0, 8), pady=5)
    kind_var.trace_add('write', update_kind);update_kind()
    def close():
        try:top.grab_release()
        except Exception:pass
        top.destroy()
    def remove_last(_e=None):
        tab.remove_last_time_marker();close();return 'break'
    def add_event(_e=None):
        try:
            color = PlotTabTk._normalize_hex_color(color_var.get())
            lw = float(width_var.get())
            if not np.isfinite(lw) or lw <= 0:raise ValueError('Line width must be a positive number.')
            event = {'kind':kind_var.get().lower(), 'label':label_var.get().strip(), 'color':color, 'linestyle':linestyle_var.get(), 'linewidth':lw}
            if kind_var.get() == 'Interval':
                a = float(start_var.get());b = float(end_var.get())
                if not np.isfinite(a) or not np.isfinite(b) or a == b:raise ValueError('Interval start and end must be different finite times.')
                alpha = float(alpha_var.get())
                if not np.isfinite(alpha) or not (0.02 <= alpha <= 0.95):raise ValueError('Interval opacity must be between 0.02 and 0.95.')
                event.update({'start':a, 'end':b, 'alpha':alpha, 'show_edges':bool(edges_var.get())})
                app._last_interval_start_str = start_var.get();app._last_interval_end_str = end_var.get();app._last_event_alpha = alpha_var.get();app._last_event_edges = bool(edges_var.get())
            else:
                x = float(time_var.get())
                if not np.isfinite(x):raise ValueError('Event time must be finite.')
                event['x'] = x;app._last_marker_time_str = time_var.get()
            app._last_event_kind = kind_var.get();app._last_event_color = color;app._last_event_linestyle = linestyle_var.get();app._last_event_linewidth = width_var.get()
            tab.add_event_annotation(event);close();return 'break'
        except Exception as e:
            messagebox.showerror('Invalid event annotation', str(e), parent=top);return 'break'
    ttk.Button(btn_row, text='Remove last', command=remove_last).grid(row=0, column=0, sticky='w')
    ttk.Button(btn_row, text='Cancel', command=close).grid(row=0, column=1, padx=(8, 8))
    ttk.Button(btn_row, text='Add event', command=add_event, style='Accent.TButton').grid(row=0, column=2)
    top.bind('<Return>', add_event);top.bind('<Escape>', lambda _e:(close(), 'break')[1]);top.bind('<Control-BackSpace>', remove_last);top.bind('<Control-Backspace>', remove_last)
    _fiberlyse_v23_apply_classic_widget_theme(top, _fiberlyse_v23_theme_tokens(app=app))
    try:(ent_start if kind_var.get() == 'Interval' else ent_time).focus_set()
    except Exception:pass
    app.root.wait_window(top)

def _fiberlyse_v23_show_time_marker_dialog(self, tab: 'PlotTabTk') -> None:
    _fiberlyse_v23_show_event_dialog(self, tab)
MainAppTk._show_time_marker_dialog = _fiberlyse_v23_show_time_marker_dialog

# ---- Centralized Export workflow ------------------------------------------
_FIBERLYSE_V23_EXPORT_GRAPHS = [
    ('raw', 'Raw'),
    ('artifact', 'Artifact remover'),
    ('fit', 'Control fit'),
    ('normalization', 'Normalization'),
    ('smoothed', 'Normalization - smoothed'),
    ('frequency', 'Frequency analysis'),
]

# Batch graphs combine recordings/channels across source files, so they are
# selected separately from the file x graph matrix below. Their export reflects
# the current selections/groups and options in the Batch compare/average tabs.
_FIBERLYSE_V23_BATCH_EXPORT_GRAPHS = [
    ('batch_compare', 'Batch compare'),
    ('batch_average', 'Batch average'),
]

def _fiberlyse_v23_safe_filename(text: str) -> str:
    s = re.sub(r'[^A-Za-z0-9._-]+', '_', str(text or '').strip()).strip('._')
    return s or 'FiberLyse'

def _fiberlyse_v23_graph_tab_and_draw(cw, graph_key: str):
    if graph_key == 'raw':cw._draw_raw();return cw.tab_raw
    if graph_key == 'artifact':cw._draw_artifact();return cw.tab_art
    if graph_key == 'fit':
        cw._draw_fit_and_attach_selector();cw._fit_initialized = True;return cw.tab_fit
    if graph_key == 'normalization':cw._draw_norm();return cw.tab_norm
    if graph_key == 'smoothed':cw._draw_norm_smooth();return cw.tab_norm_smooth
    if graph_key == 'frequency':cw._draw_frequency();cw._freq_drawn_version = int(getattr(cw.res, '_data_version', 0));return cw.tab_freq
    raise ValueError(f'Unknown graph type: {graph_key}')

def _fiberlyse_v23_save_payload_csvs(payload: Dict[str, pd.DataFrame], output_dir: str, base: str) -> List[str]:
    written = []
    for sheet, df in (payload or {}).items():
        safe_sheet = _fiberlyse_v23_safe_filename(sheet)
        path = os.path.join(output_dir, f'{base}__{safe_sheet}.csv')
        if not isinstance(df, pd.DataFrame):df = pd.DataFrame(df)
        df.to_csv(path, index=False);written.append(path)
    return written

def _fiberlyse_v23_export_batch_selected(app, batch_graph_keys: List[str], include_csv: bool, image_format: str, output_dir: str) -> Tuple[int, int]:
    """Export the current batch-analysis figures and their associated data.

    Batch compare/average are application-level graphs rather than belonging to
    one source file. The exported figure therefore represents exactly the
    current selection/group configuration and display options in those tabs.
    """
    image_format = str(image_format or 'PNG').lower()
    if image_format not in ('png', 'svg', 'pdf'):image_format = 'png'
    n_graph = 0;n_csv = 0
    for graph_key in list(batch_graph_keys or []):
        if graph_key == 'batch_compare':
            widget = getattr(app, 'compare_widget', None)
            if widget is None:
                try:app._refresh_visible_tabs();widget = getattr(app, 'compare_widget', None)
                except Exception:widget = None
            if widget is None:continue
            try:widget.refresh_available();widget.refresh_plot()
            except Exception:pass
            tab = widget.plot;base = 'batch_compare'
        elif graph_key == 'batch_average':
            widget = getattr(app, 'average_widget', None)
            if widget is None:
                try:app._refresh_visible_tabs();widget = getattr(app, 'average_widget', None)
                except Exception:widget = None
            if widget is None:continue
            try:widget.refresh_available();widget.refresh_plot()
            except Exception:pass
            tab = widget.plot;base = 'batch_average'
        else:
            continue
        _fiberlyse_v23_apply_plot_theme(tab)
        base = _fiberlyse_v23_safe_filename(base)
        graph_path = os.path.join(output_dir, f'{base}.{image_format}')
        tab.fig.savefig(graph_path, bbox_inches='tight', facecolor=tab.fig.get_facecolor())
        n_graph += 1
        if include_csv:
            payload = tab.export_provider() if callable(getattr(tab, 'export_provider', None)) else tab._payload_from_artists_fallback()
            n_csv += len(_fiberlyse_v23_save_payload_csvs(payload, output_dir, base))
        try:app.status.set(f'Exporting: {graph_key}...');app.root.update_idletasks()
        except Exception:pass
    return n_graph, n_csv

def _fiberlyse_v23_export_selected(app, file_keys: List[str], graph_keys: List[str], include_csv: bool, image_format: str, output_dir: str) -> Tuple[int, int]:
    image_format = str(image_format or 'PNG').lower()
    if image_format not in ('png', 'svg', 'pdf'):image_format = 'png'
    n_graph = 0;n_csv = 0
    original_active = getattr(app, '_active_file_key', None)
    original_sel = None
    try:original_sel = app.outer_tabs.select()
    except Exception:pass
    try:
        for file_key in file_keys:
            alias = getattr(app, '_file_alias_by_key', {}).get(file_key, file_key)
            mids = [mid for mid in list(getattr(app, '_mouse_order', [])) if str(mid).startswith(f'{file_key}:') and mid in (getattr(app, '_results', {}) or {})]
            if not mids:continue
            # Ensure hidden frames can host lazy-created channel widgets.
            for mid in mids:
                app._ensure_mouse_widget(mid)
                cw = getattr(app, '_channel_widgets', {}).get(mid)
                if cw is None:continue
                channel = app._mouse_channel_label(mid)
                for graph_key in graph_keys:
                    tab = _fiberlyse_v23_graph_tab_and_draw(cw, graph_key)
                    _fiberlyse_v23_apply_plot_theme(tab)
                    base = _fiberlyse_v23_safe_filename(f'{alias}_{channel}_{graph_key}')
                    graph_path = os.path.join(output_dir, f'{base}.{image_format}')
                    tab.fig.savefig(graph_path, bbox_inches='tight', facecolor=tab.fig.get_facecolor())
                    n_graph += 1
                    if include_csv:
                        payload = tab.export_provider() if callable(getattr(tab, 'export_provider', None)) else tab._payload_from_artists_fallback()
                        n_csv += len(_fiberlyse_v23_save_payload_csvs(payload, output_dir, base))
                    try:app.status.set(f'Exporting: {alias} / {channel} / {graph_key}...');app.root.update_idletasks()
                    except Exception:pass
    finally:
        try:app._active_file_key = original_active
        except Exception:pass
        try:
            if original_active is not None:app._refresh_visible_tabs()
            if original_sel:app.outer_tabs.select(original_sel)
        except Exception:pass
    return n_graph, n_csv

def _fiberlyse_v23_export_plan(app, plan: Dict[str, List[str]], include_csv: bool, image_format: str, output_dir: str) -> Tuple[int, int]:
    image_format = str(image_format or 'PNG').lower()
    if image_format not in ('png', 'svg', 'pdf'):image_format = 'png'
    n_graph = 0;n_csv = 0
    original_active = getattr(app, '_active_file_key', None)
    original_sel = None
    try:original_sel = app.outer_tabs.select()
    except Exception:pass
    try:
        for file_key, graph_keys in plan.items():
            if not graph_keys:continue
            alias = getattr(app, '_file_alias_by_key', {}).get(file_key, file_key)
            mids = [mid for mid in list(getattr(app, '_mouse_order', [])) if str(mid).startswith(f'{file_key}:') and mid in (getattr(app, '_results', {}) or {})]
            if not mids:continue
            for mid in mids:
                app._ensure_mouse_widget(mid)
                cw = getattr(app, '_channel_widgets', {}).get(mid)
                if cw is None:continue
                channel = app._mouse_channel_label(mid)
                for graph_key in graph_keys:
                    tab = _fiberlyse_v23_graph_tab_and_draw(cw, graph_key)
                    _fiberlyse_v23_apply_plot_theme(tab)
                    base = _fiberlyse_v23_safe_filename(f'{alias}_{channel}_{graph_key}')
                    graph_path = os.path.join(output_dir, f'{base}.{image_format}')
                    tab.fig.savefig(graph_path, bbox_inches='tight', facecolor=tab.fig.get_facecolor())
                    n_graph += 1
                    if include_csv:
                        payload = tab.export_provider() if callable(getattr(tab, 'export_provider', None)) else tab._payload_from_artists_fallback()
                        n_csv += len(_fiberlyse_v23_save_payload_csvs(payload, output_dir, base))
                    try:app.status.set(f'Exporting: {alias} / {channel} / {graph_key}...');app.root.update_idletasks()
                    except Exception:pass
    finally:
        try:app._active_file_key = original_active
        except Exception:pass
        try:
            if original_active is not None:app._refresh_visible_tabs()
            if original_sel:app.outer_tabs.select(original_sel)
        except Exception:pass
    return n_graph, n_csv

def _fiberlyse_v23_show_export_dialog(app, _source_tab=None) -> None:
    results = getattr(app, '_results', None) or {}
    if not results:
        messagebox.showinfo('Export', 'Run an analysis first. There are no analyzed graphs to export.', parent=app.root);return
    file_meta = []
    try:file_meta = _fiberlyse_file_meta_all(app)
    except Exception:
        for key in getattr(app, '_result_file_order', []):file_meta.append((key, app._file_alias_by_key.get(key, key), app._file_path_by_key.get(key, '')))
    file_meta = [(k,a,p) for (k,a,p) in file_meta if any(str(mid).startswith(f'{k}:') for mid in results)]
    if not file_meta:
        messagebox.showinfo('Export', 'No analyzed source files are available.', parent=app.root);return
    top = tk.Toplevel(app.root);top.title('Export');top.transient(app.root);top.grab_set()
    try:_fiberlyse_set_toplevel_pixel_geometry(top, app.root, 980, 640)
    except Exception:pass
    main = ttk.Frame(top, padding=12);main.pack(fill=tk.BOTH, expand=True);main.rowconfigure(2, weight=1);main.columnconfigure(0, weight=1)
    ttk.Label(main, text='Export graphs and optional CSV data', font=('TkDefaultFont', 12, 'bold')).grid(row=0, column=0, sticky='w')
    ttk.Label(main, text='Tick the exact file/graph combinations you want. Each checked cell exports that graph for every G channel found in that source file. Batch compare/average are selected separately below because they combine data across files.', justify='left', wraplength=900).grid(row=1, column=0, sticky='ew', pady=(2, 10))

    matrix_box = ttk.LabelFrame(main, text='Files × graphs');matrix_box.grid(row=2, column=0, sticky='nsew');matrix_box.rowconfigure(0, weight=1);matrix_box.columnconfigure(0, weight=1)
    canvas = tk.Canvas(matrix_box, borderwidth=0, highlightthickness=0)
    vsb = ttk.Scrollbar(matrix_box, orient='vertical', command=canvas.yview);hsb = ttk.Scrollbar(matrix_box, orient='horizontal', command=canvas.xview)
    canvas.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set);canvas.grid(row=0,column=0,sticky='nsew');vsb.grid(row=0,column=1,sticky='ns');hsb.grid(row=1,column=0,sticky='ew')
    inner = ttk.Frame(canvas);win = canvas.create_window((0,0),window=inner,anchor='nw')
    inner.bind('<Configure>', lambda _e: canvas.configure(scrollregion=canvas.bbox('all')), add='+')
    canvas.bind('<Configure>', lambda e: canvas.itemconfigure(win, height=max(1, e.height)) if inner.winfo_reqheight() < e.height else None, add='+')
    ttk.Label(inner,text='Source file',font=('TkDefaultFont',9,'bold')).grid(row=0,column=0,sticky='w',padx=(8,14),pady=6)
    for j, (_gkey, glabel) in enumerate(_FIBERLYSE_V23_EXPORT_GRAPHS, start=1):
        ttk.Label(inner,text=glabel,font=('TkDefaultFont',9,'bold'),justify='center',wraplength=110).grid(row=0,column=j,sticky='nsew',padx=6,pady=6)
    cell_vars: Dict[Tuple[str,str], Any] = {}
    for i, (fkey, alias, path) in enumerate(file_meta, start=1):
        label = f'{alias}\n{os.path.basename(path)}'
        ttk.Label(inner,text=label,justify='left',wraplength=220).grid(row=i,column=0,sticky='w',padx=(8,14),pady=5)
        for j, (gkey, _glabel) in enumerate(_FIBERLYSE_V23_EXPORT_GRAPHS, start=1):
            v=tk.BooleanVar(value=True);cell_vars[(fkey,gkey)]=v
            ttk.Checkbutton(inner,variable=v).grid(row=i,column=j,padx=12,pady=5)
    for c in range(1,1+len(_FIBERLYSE_V23_EXPORT_GRAPHS)):inner.columnconfigure(c,minsize=105)

    # Batch analysis graphs are global/current-state outputs, not tied to one
    # file row, so expose them in their own clearly labeled section.
    batch_vars: Dict[str, Any] = {}
    batch_box = ttk.LabelFrame(main, text='Batch analysis graphs (current selections/groups)');batch_box.grid(row=4,column=0,sticky='ew',pady=(10,0))
    ttk.Label(batch_box,text='These exports use the selections, groups, smoothing/SEM options, normalization mode, and other settings currently shown in the Batch tabs.',justify='left',wraplength=850).grid(row=0,column=0,columnspan=4,sticky='w',padx=8,pady=(7,4))
    for j,(gkey,glabel) in enumerate(_FIBERLYSE_V23_BATCH_EXPORT_GRAPHS):
        v=tk.BooleanVar(value=True);batch_vars[gkey]=v
        ttk.Checkbutton(batch_box,text=glabel,variable=v).grid(row=1,column=j,sticky='w',padx=(8,18),pady=(2,8))

    quick = ttk.Frame(main);quick.grid(row=3,column=0,sticky='ew',pady=(7,0))
    def set_all(value: bool):
        for v in cell_vars.values():v.set(bool(value))
        for v in batch_vars.values():v.set(bool(value))
    ttk.Button(quick,text='Select all',command=lambda:set_all(True)).pack(side=tk.LEFT)
    ttk.Button(quick,text='Clear all',command=lambda:set_all(False)).pack(side=tk.LEFT,padx=(6,0))
    # Convenient per-graph selection lets a user quickly select, e.g., only Fit
    # for every file and then add a few individual cells.
    ttk.Label(quick,text='Select graph for all files:').pack(side=tk.LEFT,padx=(16,6))
    graph_quick_var=tk.StringVar(value=_FIBERLYSE_V23_EXPORT_GRAPHS[0][1]);graph_quick=ttk.Combobox(quick,textvariable=graph_quick_var,values=tuple(x[1] for x in _FIBERLYSE_V23_EXPORT_GRAPHS),state='readonly',width=24);graph_quick.pack(side=tk.LEFT)
    def select_graph_column():
        label=graph_quick_var.get();key=next((k for k,l in _FIBERLYSE_V23_EXPORT_GRAPHS if l==label),None)
        if key:
            for fkey, _a, _p in file_meta:cell_vars[(fkey,key)].set(True)
    ttk.Button(quick,text='Select',command=select_graph_column).pack(side=tk.LEFT,padx=(6,0))

    opts = ttk.LabelFrame(main, text='Export options');opts.grid(row=5,column=0,sticky='ew',pady=(10,0));opts.columnconfigure(3, weight=1)
    include_csv_var = tk.BooleanVar(value=True);fmt_var = tk.StringVar(value='PNG')
    ttk.Checkbutton(opts,text='Also export CSV data related to each selected graph',variable=include_csv_var).grid(row=0,column=0,columnspan=2,sticky='w',padx=8,pady=7)
    ttk.Label(opts,text='Graph format:').grid(row=1,column=0,sticky='w',padx=8,pady=(0,8));ttk.Combobox(opts,textvariable=fmt_var,values=('PNG','SVG','PDF'),state='readonly',width=10).grid(row=1,column=1,sticky='w',pady=(0,8))
    ttk.Label(opts,text='CSV export writes one CSV for each data table associated with a graph (for example raw excitatory, raw isosbestic, and metadata).',justify='left',wraplength=850).grid(row=2,column=0,columnspan=4,sticky='ew',padx=8,pady=(0,8))
    btn = ttk.Frame(main);btn.grid(row=6,column=0,sticky='ew',pady=(12,0));btn.columnconfigure(0,weight=1)
    def close():
        try:top.grab_release()
        except Exception:pass
        top.destroy()
    def do_export():
        plan: Dict[str,List[str]] = {}
        for fkey, _a, _p in file_meta:
            chosen=[gkey for gkey,_label in _FIBERLYSE_V23_EXPORT_GRAPHS if bool(cell_vars[(fkey,gkey)].get())]
            if chosen:plan[fkey]=chosen
        batch_chosen=[gkey for gkey,_label in _FIBERLYSE_V23_BATCH_EXPORT_GRAPHS if bool(batch_vars[gkey].get())]
        if not plan and not batch_chosen:
            messagebox.showerror('Export', 'Select at least one file/graph combination or batch-analysis graph.', parent=top);return
        out = filedialog.askdirectory(title='Choose export folder', parent=top)
        if not out:return
        try:
            top.configure(cursor='watch');top.update_idletasks()
            n_graph, n_csv = _fiberlyse_v23_export_plan(app, plan, bool(include_csv_var.get()), fmt_var.get(), out) if plan else (0,0)
            b_graph, b_csv = _fiberlyse_v23_export_batch_selected(app, batch_chosen, bool(include_csv_var.get()), fmt_var.get(), out)
            n_graph += b_graph;n_csv += b_csv
            top.configure(cursor='')
            messagebox.showinfo('Export complete', f'Exported {n_graph} graph file(s)' + (f' and {n_csv} CSV file(s)' if include_csv_var.get() else '') + f'.\n\nFolder:\n{out}', parent=top)
            try:app.status.set(f'Export complete: {n_graph} graphs' + (f', {n_csv} CSVs' if include_csv_var.get() else ''))
            except Exception:pass
            close()
        except Exception as e:
            try:top.configure(cursor='')
            except Exception:pass
            messagebox.showerror('Export failed', f'Could not complete export:\n\n{e}', parent=top)
    ttk.Button(btn,text='Cancel',command=close).grid(row=0,column=1,padx=(0,8));ttk.Button(btn,text='Export selected',command=do_export,style='Accent.TButton').grid(row=0,column=2)
    top.bind('<Escape>', lambda _e:(close(),'break')[1])
    _fiberlyse_v23_apply_classic_widget_theme(top, _fiberlyse_v23_theme_tokens(app=app))
    app.root.wait_window(top)

MainAppTk.show_export_dialog = _fiberlyse_v23_show_export_dialog

def _fiberlyse_v23_export_excel_redirect(self):
    try:app = getattr(self.winfo_toplevel(), '_fiberlyse_app', None)
    except Exception:app = None
    if app is not None:return _fiberlyse_v23_show_export_dialog(app, self)
    return None
PlotTabTk.export_excel = _fiberlyse_v23_export_excel_redirect

# Reconfigure every PlotTab export button after the full V22 initializer chain has run.
_FIBERLYSE_V23_PLOTTAB_INIT_BEFORE_EXPORT = PlotTabTk.__init__
def _fiberlyse_v23_plottab_init(self, *args, **kwargs):
    _FIBERLYSE_V23_PLOTTAB_INIT_BEFORE_EXPORT(self, *args, **kwargs)
    try:
        self.export_btn.configure(text='Export…', command=lambda tab=self: _fiberlyse_v23_show_export_dialog(getattr(tab.winfo_toplevel(), '_fiberlyse_app', None), tab) if getattr(tab.winfo_toplevel(), '_fiberlyse_app', None) is not None else messagebox.showinfo('Export', 'The central export dialog is available from a FiberLyse analysis window.'))
    except Exception:pass
    try:_fiberlyse_v23_apply_plot_theme(self)
    except Exception:pass
PlotTabTk.__init__ = _fiberlyse_v23_plottab_init

# ---- Main V23 initialization ----------------------------------------------
_FIBERLYSE_V23_MAIN_INIT_BEFORE_FEATURES = MainAppTk.__init__
def _fiberlyse_v23_mainapp_init(self, *args, **kwargs):
    _FIBERLYSE_V23_MAIN_INIT_BEFORE_FEATURES(self, *args, **kwargs)
    try:self.root._fiberlyse_app = self
    except Exception:pass
    # Night mode is an explicit GUI preference. Start in light mode for backwards
    # familiarity and let the user switch instantly without rerunning analysis.
    self.var_night_mode = tk.BooleanVar(value=False)
    self._fiberlyse_night_mode = False
    try:
        top = self.btn_run.master
        self.chk_night_mode = ttk.Checkbutton(top, text='Night mode', variable=self.var_night_mode, command=lambda app=self:_fiberlyse_v23_toggle_night_mode(app))
        # Place after existing progress/cancel controls; grid adapts even if the
        # modern shell changes neighboring columns.
        self.chk_night_mode.grid(row=0, column=6, sticky='e', padx=(10, 0))
        _fiberlyse_attach_tooltip(self.chk_night_mode, 'Switch the FiberLyse GUI and Matplotlib plots between light and dark appearance. This does not change data or analysis.')
    except Exception:pass
    _fiberlyse_v23_apply_ttk_theme(self)
    try:_fiberlyse_v23_apply_classic_widget_theme(self.root, _fiberlyse_v23_theme_tokens(app=self))
    except Exception:pass
    # Avoid the legacy duplicate Ctrl+I bindings and register exactly one event
    # annotation handler.
    try:self.root.unbind_all('<Control-i>');self.root.unbind_all('<Control-I>')
    except Exception:pass
    try:self.root.bind_all('<Control-i>', self._on_ctrl_i, add='+');self.root.bind_all('<Control-I>', self._on_ctrl_i, add='+')
    except Exception:pass
    try:
        reg = getattr(self, 'shortcut_registry', None)
        if reg is not None and 'time_marker' in reg.actions:
            reg.actions['time_marker']['label'] = 'Add line / interval event annotation on active plot'
            reg.actions['time_marker']['callback'] = self._on_ctrl_i
    except Exception:pass
    try:self.status.set('Ready. V23: Night mode, centralized Export, and line/interval events available.')
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_v23_mainapp_init

# Re-theme plots created lazily after a night-mode switch.
_FIBERLYSE_V23_ENSURE_WIDGET_BEFORE_THEME = MainAppTk._ensure_mouse_widget
def _fiberlyse_v23_ensure_mouse_widget(self, mid: str) -> None:
    _FIBERLYSE_V23_ENSURE_WIDGET_BEFORE_THEME(self, mid)
    try:
        cw = getattr(self, '_channel_widgets', {}).get(mid)
        if cw is not None:
            for name in ['tab_raw','tab_art','tab_fit','tab_norm','tab_norm_smooth','tab_freq']:
                tab = getattr(cw, name, None)
                if tab is not None:_fiberlyse_v23_apply_plot_theme(tab)
    except Exception:pass
MainAppTk._ensure_mouse_widget = _fiberlyse_v23_ensure_mouse_widget

# Ensure batch tabs created by later UI actions also inherit the active theme.
_FIBERLYSE_V23_BATCH_COMPARE_INIT = BatchCompareTk.__init__
def _fiberlyse_v23_batch_compare_init(self, *args, **kwargs):
    _FIBERLYSE_V23_BATCH_COMPARE_INIT(self, *args, **kwargs)
    try:_fiberlyse_v23_apply_plot_theme(self.plot)
    except Exception:pass
BatchCompareTk.__init__ = _fiberlyse_v23_batch_compare_init
_FIBERLYSE_V23_BATCH_AVERAGE_INIT = BatchAverageTk.__init__
def _fiberlyse_v23_batch_average_init(self, *args, **kwargs):
    _FIBERLYSE_V23_BATCH_AVERAGE_INIT(self, *args, **kwargs)
    try:_fiberlyse_v23_apply_plot_theme(self.plot)
    except Exception:pass
BatchAverageTk.__init__ = _fiberlyse_v23_batch_average_init
# ---- End FiberLyse V23 feature extension ----------------------------------

# ---- FiberLyse V23.1 universal CSV setup extension ------------------------
# V23.1 keeps the V23 scientific pipeline intact and adds a beginner-facing
# import layer.  Hardware-specific CSVs are translated into the same internal
# ChannelResult representation before artifact removal/fitting/normalization.
import json as _fiberlyse_v231_json

FIBERLYSE_VERSION = 'V23.1'
_FIBERLYSE_V231_PROFILE_DIR = os.path.join(os.path.expanduser('~'), '.fiberlyse')
_FIBERLYSE_V231_PROFILE_PATH = os.path.join(_FIBERLYSE_V231_PROFILE_DIR, 'import_setup.json')
_FIBERLYSE_V231_LAYOUT_BUILTIN = 'neurophotometrics_bonsai'
_FIBERLYSE_V231_LAYOUT_SEPARATE = 'separate_columns'
_FIBERLYSE_V231_LAYOUT_INTERLEAVED = 'alternating_rows'
_FIBERLYSE_V231_TIME_UNITS = ('Seconds', 'Milliseconds', 'Microseconds')


def _fiberlyse_v231_builtin_profile() -> Dict[str, Any]:
    return {
        'profile_name': 'Neurophotometrics / Bonsai',
        'layout': _FIBERLYSE_V231_LAYOUT_BUILTIN,
        'exclude_initial_seconds': float(ANALYSIS_EXCLUDE_INITIAL_SECONDS),
        'remember': False,
    }


def _fiberlyse_v231_profile_label(profile: Optional[Dict[str, Any]]) -> str:
    p = profile or _fiberlyse_v231_builtin_profile()
    return str(p.get('profile_name') or 'CSV setup')


def _fiberlyse_v231_save_profile(profile: Dict[str, Any]) -> None:
    try:
        os.makedirs(_FIBERLYSE_V231_PROFILE_DIR, exist_ok=True)
        clean = dict(profile)
        clean['remember'] = True
        with open(_FIBERLYSE_V231_PROFILE_PATH, 'w', encoding='utf-8') as fh:
            _fiberlyse_v231_json.dump(clean, fh, indent=2, ensure_ascii=False)
    except Exception:
        pass


def _fiberlyse_v231_forget_saved_profile() -> None:
    try:
        if os.path.exists(_FIBERLYSE_V231_PROFILE_PATH):
            os.remove(_FIBERLYSE_V231_PROFILE_PATH)
    except Exception:
        pass


def _fiberlyse_v231_load_saved_profile() -> Optional[Dict[str, Any]]:
    try:
        if not os.path.isfile(_FIBERLYSE_V231_PROFILE_PATH):
            return None
        with open(_FIBERLYSE_V231_PROFILE_PATH, 'r', encoding='utf-8') as fh:
            p = _fiberlyse_v231_json.load(fh)
        if isinstance(p, dict) and p.get('layout') in {
            _FIBERLYSE_V231_LAYOUT_BUILTIN,
            _FIBERLYSE_V231_LAYOUT_SEPARATE,
            _FIBERLYSE_V231_LAYOUT_INTERLEAVED,
        }:
            return p
    except Exception:
        pass
    return None


def _fiberlyse_v231_time_scale(units: str) -> float:
    u = str(units or 'Seconds').strip().lower()
    if u.startswith('milli'):return 1.0 / 1000.0
    if u.startswith('micro'):return 1.0 / 1000000.0
    return 1.0


def _fiberlyse_v231_match_value(values, target: Any) -> np.ndarray:
    s = pd.Series(values)
    target_txt = str(target).strip()
    # Text comparison covers values such as "signal" / "reference".
    text_mask = s.astype(str).str.strip().to_numpy() == target_txt
    try:
        target_num = float(target_txt)
        nums = pd.to_numeric(s, errors='coerce').to_numpy(dtype=float)
        num_mask = np.isfinite(nums) & np.isclose(nums, target_num, rtol=0.0, atol=1e-12)
        return np.asarray(text_mask | num_mask, dtype=bool)
    except Exception:
        return np.asarray(text_mask, dtype=bool)


def _fiberlyse_v231_safe_channel_name(text: str, index: int) -> str:
    # Friendly internal names are intentionally hardware-neutral.
    t = re.sub(r'[^A-Za-z0-9_-]+', '_', str(text or '')).strip('_')
    return t[:40] if t else f'Channel{index}'


def _fiberlyse_v231_first_start_index(df: pd.DataFrame, profile: Dict[str, Any]) -> int:
    if str(profile.get('start_rule', 'first_row')) != 'column_value':
        return 0
    col = str(profile.get('start_column') or '')
    val = profile.get('start_value', '')
    if not col or col not in df.columns:
        raise ValueError("The CSV setup says recording start is marked by a column, but that column is missing.")
    idx = np.flatnonzero(_fiberlyse_v231_match_value(df[col], val))
    if idx.size == 0:
        raise ValueError(f"Could not find the recording-start value '{val}' in column '{col}'.")
    return int(idx[0])


def _fiberlyse_v231_prepare_time(df: pd.DataFrame, profile: Dict[str, Any]) -> Tuple[pd.DataFrame, np.ndarray, float]:
    time_col = str(profile.get('time_column') or '')
    if not time_col or time_col not in df.columns:
        raise ValueError("The column selected for time is missing from this CSV. Open 'CSV setup…' and check the mapping.")
    start_idx = _fiberlyse_v231_first_start_index(df, profile)
    df = df.iloc[start_idx:].reset_index(drop=True)
    try:
        raw_t = pd.to_numeric(df[time_col], errors='raise').to_numpy(dtype=float)
    except Exception as e:
        raise ValueError(f"The selected time column '{time_col}' must contain numbers: {e}")
    if raw_t.size < 2 or not np.all(np.isfinite(raw_t)):
        raise ValueError('The selected time column needs at least two finite numeric values.')
    raw_t = raw_t * _fiberlyse_v231_time_scale(str(profile.get('time_units', 'Seconds')))
    dt = np.diff(raw_t)
    if np.any(~np.isfinite(dt)) or np.any(dt <= 0):
        bad = int(np.sum(~np.isfinite(dt) | (dt <= 0)))
        raise ValueError(f'Time must increase from one row to the next; found {bad} duplicate/decreasing interval(s).')
    elapsed = raw_t - float(raw_t[0])
    try:exclude_s = max(0.0, float(profile.get('exclude_initial_seconds', 0.0)))
    except Exception:exclude_s = 0.0
    keep = np.isfinite(elapsed) & (elapsed >= exclude_s)
    if int(np.sum(keep)) < 4:
        raise ValueError(f'No usable recording remains after excluding the first {exclude_s:g} seconds.')
    return df.loc[keep].reset_index(drop=True), elapsed[keep], exclude_s


def _fiberlyse_v231_make_result(
    csv_path: str,
    channel_name: str,
    t_iso: np.ndarray,
    t_exc: np.ndarray,
    iso_raw: np.ndarray,
    exc_raw: np.ndarray,
    settings: Dict[str, Any],
    profile: Dict[str, Any],
    import_meta: Dict[str, Any],
) -> ChannelResult:
    t_iso = np.asarray(t_iso, dtype=float);t_exc = np.asarray(t_exc, dtype=float)
    iso_raw = np.asarray(iso_raw, dtype=float);exc_raw = np.asarray(exc_raw, dtype=float)
    if t_iso.size < 2 or t_exc.size < 2:
        raise ValueError(f'{channel_name}: not enough Signal/Reference samples after import.')
    measured_eff_fs = estimate_fs_from_t(t_exc)
    # V24.1: use the effective sampling rate measured from the study-signal
    # timestamps only. Original acquisition rate/downsampling history belongs in
    # the import/frequency metadata and cannot substitute for missing samples.
    supplied_acq = np.nan
    eff_fs = measured_eff_fs if np.isfinite(measured_eff_fs) and measured_eff_fs > 0 else np.nan
    acq_fps_val = np.nan
    artifact_enabled = bool(settings['artifact_enabled']);artifact_factor = float(settings['artifact_factor']);artifact_pad = int(settings['artifact_pad']);require_shared = bool(settings['require_shared']);align_mode = str(settings['align_mode']);use_linear_interp = bool(settings['use_interp'])
    if artifact_enabled:
        base_iso = detect_artifacts_by_derivative(t_iso, iso_raw, factor=artifact_factor, method='mad', pad=0)
        base_exc = detect_artifacts_by_derivative(t_exc, exc_raw, factor=artifact_factor, method='mad', pad=0)
        if require_shared:base_iso, base_exc = shared_artifacts_by_time(t_iso, base_iso, t_exc, base_exc)
        art_iso = expand_artifact_mask(base_iso, artifact_pad);art_exc = expand_artifact_mask(base_exc, artifact_pad)
    else:
        art_iso = np.zeros_like(iso_raw, dtype=bool);art_exc = np.zeros_like(exc_raw, dtype=bool)
    iso_clean_holes = remove_with_holes(iso_raw, art_iso);exc_clean_holes = remove_with_holes(exc_raw, art_exc)
    iso_clean_interp = linear_interpolate_by_time(t_iso, iso_clean_holes);exc_clean_interp = linear_interpolate_by_time(t_exc, exc_clean_holes)
    iso_on_exc_holes = align_iso_to_exc_no_interp(t_iso, iso_clean_holes, t_exc, mode=align_mode);iso_on_exc_interp = align_iso_to_exc_no_interp(t_iso, iso_clean_interp, t_exc, mode=align_mode)
    iso_on_exc_holes = np.asarray(iso_on_exc_holes, dtype=float);iso_on_exc_interp = np.asarray(iso_on_exc_interp, dtype=float);iso_on_exc_holes[np.asarray(art_exc, dtype=bool)] = np.nan
    if use_linear_interp:iso_clean = iso_clean_interp;exc_clean = exc_clean_interp;iso_on_exc = iso_on_exc_interp
    else:iso_clean = iso_clean_holes;exc_clean = exc_clean_holes;iso_on_exc = iso_on_exc_holes
    fit_windows = list(settings.get('fit_windows') or [(float(profile.get('exclude_initial_seconds', 0.0)), 6500.0)])
    res = ChannelResult(
        gcol=channel_name, source_path=os.path.abspath(str(csv_path)), analysis_exclude_initial_seconds=float(profile.get('exclude_initial_seconds', 0.0)),
        artifact_enabled=artifact_enabled, artifact_factor=artifact_factor, artifact_pad=artifact_pad, require_shared=require_shared, align_mode=align_mode,
        t_iso=t_iso, t_exc=t_exc, iso_raw=iso_raw, exc_raw=exc_raw, art_iso=np.asarray(art_iso, dtype=bool), art_exc=np.asarray(art_exc, dtype=bool),
        iso_clean_holes=iso_clean_holes, exc_clean_holes=exc_clean_holes, iso_clean_interp=iso_clean_interp, exc_clean_interp=exc_clean_interp,
        iso_on_exc_holes=iso_on_exc_holes, iso_on_exc_interp=iso_on_exc_interp, use_interpolation=use_linear_interp,
        iso_clean=iso_clean, exc_clean=exc_clean, iso_on_exc=iso_on_exc, windows=fit_windows,
        acq_fps_hz=float(acq_fps_val), eff_fs_hz=float(eff_fs), slope=np.nan, intercept=np.nan, r2=np.nan,
        fitted_iso_on_exc=np.full_like(t_exc, np.nan, dtype=float), residual=np.full_like(t_exc, np.nan, dtype=float), dF=np.full_like(t_exc, np.nan, dtype=float), dFF=np.full_like(t_exc, np.nan, dtype=float),
        slope_nointerp=np.nan, intercept_nointerp=np.nan, r2_nointerp=np.nan, fitted_iso_on_exc_nointerp=np.full_like(t_exc, np.nan, dtype=float), residual_nointerp=np.full_like(t_exc, np.nan, dtype=float), dF_nointerp=np.full_like(t_exc, np.nan, dtype=float), dFF_nointerp=np.full_like(t_exc, np.nan, dtype=float),
        smooth_window=max(1, int(settings['smooth_win'])), zf_interval_start_s=float(settings['interval_start']), zf_interval_end_s=float(settings['interval_end']),
        zF_global=np.full_like(t_exc, np.nan, dtype=float), zF_interval=np.full_like(t_exc, np.nan, dtype=float),
    )
    res.import_profile_name = _fiberlyse_v231_profile_label(profile)
    res.import_layout = str(profile.get('layout', ''))
    res.import_time_column = str(profile.get('time_column', ''))
    res.import_time_units = str(profile.get('time_units', 'Seconds'))
    for k, v in dict(import_meta or {}).items():
        try:setattr(res, str(k), v)
        except Exception:pass
    recompute_fit_and_downstream(res)
    return res


def _fiberlyse_v231_analyze_separate_columns(csv_path: str, profile: Dict[str, Any], settings: Dict[str, Any]) -> Dict[str, ChannelResult]:
    header = pd.read_csv(csv_path, nrows=0)
    pairs = list(profile.get('channel_pairs') or [])
    if not pairs:raise ValueError("No Signal/Reference column pair is defined. Open 'CSV setup…'.")
    required = {str(profile.get('time_column') or '')}
    if str(profile.get('start_rule', 'first_row')) == 'column_value':required.add(str(profile.get('start_column') or ''))
    for pair in pairs:
        required.add(str(pair.get('signal_column') or ''));required.add(str(pair.get('reference_column') or ''))
    missing = sorted(c for c in required if c and c not in header.columns)
    if missing:raise ValueError('This CSV does not match the active setup. Missing column(s): ' + ', '.join(missing))
    df = pd.read_csv(csv_path, usecols=[c for c in header.columns if c in required])
    df, elapsed, _exclude = _fiberlyse_v231_prepare_time(df, profile)
    results: Dict[str, ChannelResult] = {}
    for i, pair in enumerate(pairs, start=1):
        sig_col = str(pair.get('signal_column') or '');ref_col = str(pair.get('reference_column') or '')
        label = str(pair.get('label') or f'Channel{i}')
        sig = pd.to_numeric(df[sig_col], errors='coerce').to_numpy(dtype=float);ref = pd.to_numeric(df[ref_col], errors='coerce').to_numpy(dtype=float)
        key = _fiberlyse_v231_safe_channel_name(label, i)
        results[key] = _fiberlyse_v231_make_result(csv_path, key, elapsed, elapsed, ref, sig, settings, profile, {
            'import_signal_column': sig_col, 'import_reference_column': ref_col, 'import_channel_label': label,
        })
    return results


def _fiberlyse_v231_analyze_alternating_rows(csv_path: str, profile: Dict[str, Any], settings: Dict[str, Any]) -> Dict[str, ChannelResult]:
    header = pd.read_csv(csv_path, nrows=0)
    value_cols = [str(c) for c in list(profile.get('measurement_columns') or []) if str(c)]
    if not value_cols:raise ValueError("No fluorescence-value column is defined. Open 'CSV setup…'.")
    row_type_col = str(profile.get('row_type_column') or '')
    required = {str(profile.get('time_column') or ''), row_type_col, *value_cols}
    if str(profile.get('start_rule', 'first_row')) == 'column_value':required.add(str(profile.get('start_column') or ''))
    missing = sorted(c for c in required if c and c not in header.columns)
    if missing:raise ValueError('This CSV does not match the active setup. Missing column(s): ' + ', '.join(missing))
    df = pd.read_csv(csv_path, usecols=[c for c in header.columns if c in required])
    df, elapsed, _exclude = _fiberlyse_v231_prepare_time(df, profile)
    sig_mask = _fiberlyse_v231_match_value(df[row_type_col], profile.get('signal_row_value', ''))
    ref_mask = _fiberlyse_v231_match_value(df[row_type_col], profile.get('reference_row_value', ''))
    if int(np.sum(sig_mask)) < 2 or int(np.sum(ref_mask)) < 2:
        raise ValueError('The selected Signal/Reference row values did not find enough rows. Check the row-type column and values in CSV setup.')
    results: Dict[str, ChannelResult] = {}
    for i, col in enumerate(value_cols, start=1):
        vals = pd.to_numeric(df[col], errors='coerce').to_numpy(dtype=float)
        label = _fiberlyse_v231_safe_channel_name(col, i)
        results[label] = _fiberlyse_v231_make_result(csv_path, label, elapsed[ref_mask], elapsed[sig_mask], vals[ref_mask], vals[sig_mask], settings, profile, {
            'import_signal_column': col, 'import_reference_column': col, 'import_row_type_column': row_type_col,
            'import_signal_row_value': str(profile.get('signal_row_value', '')), 'import_reference_row_value': str(profile.get('reference_row_value', '')),
        })
    return results


def _fiberlyse_v231_analyze_file(csv_path: str, profile: Dict[str, Any], settings: Dict[str, Any]) -> Dict[str, ChannelResult]:
    layout = str((profile or {}).get('layout') or _FIBERLYSE_V231_LAYOUT_BUILTIN)
    if layout == _FIBERLYSE_V231_LAYOUT_BUILTIN:
        return analyze_csv(csv_path, artifact_enabled=settings['artifact_enabled'], artifact_factor=settings['artifact_factor'], artifact_method='mad', artifact_pad=settings['artifact_pad'], require_shared=settings['require_shared'], align_mode=settings['align_mode'], fit_windows=settings['fit_windows'], acq_fps_hz=settings['acq_fps'] if np.isfinite(settings['acq_fps']) and settings['acq_fps'] > 0 else None, smooth_window=settings['smooth_win'], zf_interval_start_s=settings['interval_start'], zf_interval_end_s=settings['interval_end'], use_linear_interp=settings['use_interp'])
    if layout == _FIBERLYSE_V231_LAYOUT_SEPARATE:
        return _fiberlyse_v231_analyze_separate_columns(csv_path, profile, settings)
    if layout == _FIBERLYSE_V231_LAYOUT_INTERLEAVED:
        return _fiberlyse_v231_analyze_alternating_rows(csv_path, profile, settings)
    raise ValueError('Unknown CSV setup. Open CSV setup and choose a supported layout.')


def _fiberlyse_v231_profile_summary(profile: Dict[str, Any]) -> str:
    layout = str(profile.get('layout', ''))
    if layout == _FIBERLYSE_V231_LAYOUT_BUILTIN:
        return 'Built-in Neurophotometrics/Bonsai mapping: SystemTimestamp = time, LedState = row type, G0/G1/... = fluorescence.'
    if layout == _FIBERLYSE_V231_LAYOUT_SEPARATE:
        pairs = list(profile.get('channel_pairs') or [])
        return f"Separate columns • Time: {profile.get('time_column','?')} • {len(pairs)} Signal/Reference pair(s) • first {float(profile.get('exclude_initial_seconds',0)):g} s excluded"
    cols = list(profile.get('measurement_columns') or [])
    return f"Alternating rows • Time: {profile.get('time_column','?')} • Row type: {profile.get('row_type_column','?')} • {len(cols)} fluorescence column(s) • first {float(profile.get('exclude_initial_seconds',0)):g} s excluded"


def _fiberlyse_v231_choose_example_csv(parent) -> Optional[str]:
    p = filedialog.askopenfilename(parent=parent, title='Choose one example CSV so FiberLyse can show its column names', filetypes=[('CSV files','*.csv'),('All files','*.*')])
    return str(p) if p else None


def _fiberlyse_v231_read_example(path: str) -> Tuple[List[str], pd.DataFrame]:
    df = pd.read_csv(path, nrows=8)
    return [str(c) for c in df.columns], df


def _fiberlyse_v231_preview_text(df: pd.DataFrame) -> str:
    try:return df.head(6).to_string(index=False, max_cols=12, max_colwidth=18)
    except Exception:return ''


def _fiberlyse_v231_suggest_column(columns: List[str], keywords: List[str], fallback_index: int=0, exclude: Optional[set]=None) -> str:
    exclude = set(exclude or set())
    lowered=[(c,str(c).lower()) for c in columns if c not in exclude]
    for kw in keywords:
        for c,lc in lowered:
            if kw in lc:return c
    remaining=[c for c in columns if c not in exclude]
    if not remaining:return ''
    return remaining[min(max(0,int(fallback_index)),len(remaining)-1)]


def _fiberlyse_v231_configure_custom(app, layout: str, initial_example: Optional[str]=None) -> Optional[Dict[str, Any]]:
    example_path = initial_example or _fiberlyse_v231_choose_example_csv(app.root)
    if not example_path:return None
    try:columns, preview_df = _fiberlyse_v231_read_example(example_path)
    except Exception as e:
        messagebox.showerror('Could not read CSV', f'FiberLyse could not read that example file:\n\n{e}', parent=app.root);return None
    top = tk.Toplevel(app.root);top.title('Define your CSV');top.transient(app.root);top.grab_set();top.geometry('960x760');top.minsize(820, 620)
    outer = ttk.Frame(top, padding=12);outer.pack(fill=tk.BOTH, expand=True)
    title = 'Tell FiberLyse what the columns mean'
    ttk.Label(outer, text=title, font=('TkDefaultFont', 14, 'bold')).pack(anchor='w')
    ttk.Label(outer, text='You do not need to know FiberLyse terminology. Choose the column that best matches each plain-language question below. Nothing is changed in your CSV.', wraplength=900, justify='left').pack(anchor='w', pady=(4,8))
    ttk.Label(outer, text=f'Example file: {os.path.basename(example_path)}').pack(anchor='w')
    preview_box = ttk.LabelFrame(outer, text='Small preview of your CSV');preview_box.pack(fill=tk.X, pady=(8,10))
    txt = tk.Text(preview_box, height=7, wrap='none');txt.pack(fill=tk.X, expand=True, padx=6, pady=6);txt.insert('1.0', _fiberlyse_v231_preview_text(preview_df));txt.configure(state='disabled')
    form = ttk.Frame(outer);form.pack(fill=tk.BOTH, expand=True)
    form.columnconfigure(1, weight=1)
    time_var = tk.StringVar(value=_fiberlyse_v231_suggest_column(columns, ['timestamp','time','clock','seconds','sec'], 0))
    units_var = tk.StringVar(value='Seconds')
    start_rule_var = tk.StringVar(value='first_row')
    start_col_var = tk.StringVar(value=_fiberlyse_v231_suggest_column(columns, ['start','ledstate','state','event','trigger'], 0))
    start_val_var = tk.StringVar(value='')
    exclude_var = tk.StringVar(value='0')
    profile_name_var = tk.StringVar(value='My CSV setup')

    row = 0
    ttk.Label(form, text='1. Which column tells FiberLyse when each measurement was recorded?', wraplength=430, justify='left').grid(row=row,column=0,sticky='w',pady=5)
    ttk.Combobox(form,textvariable=time_var,values=columns,state='readonly',width=34).grid(row=row,column=1,sticky='ew',padx=(12,0),pady=5);row+=1
    ttk.Label(form, text='What units are those time values in?').grid(row=row,column=0,sticky='w',pady=5)
    ttk.Combobox(form,textvariable=units_var,values=_FIBERLYSE_V231_TIME_UNITS,state='readonly',width=20).grid(row=row,column=1,sticky='w',padx=(12,0),pady=5);row+=1

    dynamic = ttk.LabelFrame(form, text='2. Where are the fluorescence values?');dynamic.grid(row=row,column=0,columnspan=2,sticky='nsew',pady=(10,6));dynamic.columnconfigure(1,weight=1);row+=1
    pair_vars: List[Tuple[tk.StringVar, tk.StringVar, tk.StringVar]] = []
    measurement_list = None
    rowtype_var = tk.StringVar(value=_fiberlyse_v231_suggest_column(columns, ['ledstate','state','wavelength','type','kind','channel'], 0))
    signal_code_var = tk.StringVar(value='2')
    reference_code_var = tk.StringVar(value='1')
    count_var = tk.IntVar(value=1)
    pair_holder = ttk.Frame(dynamic)

    def rebuild_pairs(*_args):
        for c in pair_holder.winfo_children():c.destroy()
        old = [(a.get(),b.get(),c.get()) for a,b,c in pair_vars];pair_vars.clear()
        n=max(1,min(8,int(count_var.get() or 1)))
        for i in range(n):
            prev = old[i] if i < len(old) else ('', '', f'Channel {i+1}')
            sig_guess=_fiberlyse_v231_suggest_column(columns,['470','465','signal','sensor','green','activity'],1,exclude={time_var.get()});ref_guess=_fiberlyse_v231_suggest_column(columns,['405','415','reference','control','isosbestic','iso'],2,exclude={time_var.get(),sig_guess});sv=tk.StringVar(value=prev[0] or sig_guess);rv=tk.StringVar(value=prev[1] or ref_guess);lv=tk.StringVar(value=prev[2] or f'Channel {i+1}');pair_vars.append((sv,rv,lv))
            ttk.Label(pair_holder,text=f'Recording channel {i+1}').grid(row=i,column=0,sticky='w',pady=3)
            ttk.Label(pair_holder,text='Signal you want to study:').grid(row=i,column=1,sticky='e',padx=(8,4));ttk.Combobox(pair_holder,textvariable=sv,values=columns,state='readonly',width=24).grid(row=i,column=2,sticky='ew')
            ttk.Label(pair_holder,text='Comparison/reference used for correction:').grid(row=i,column=3,sticky='e',padx=(8,4));ttk.Combobox(pair_holder,textvariable=rv,values=columns,state='readonly',width=24).grid(row=i,column=4,sticky='ew')
            ttk.Entry(pair_holder,textvariable=lv,width=14).grid(row=i,column=5,sticky='ew',padx=(8,0))
        pair_holder.columnconfigure(2,weight=1);pair_holder.columnconfigure(4,weight=1)

    if layout == _FIBERLYSE_V231_LAYOUT_SEPARATE:
        ttk.Label(dynamic,text='Use this when the biological Signal and the Reference/control already have their own columns in every row.',wraplength=820,justify='left').grid(row=0,column=0,columnspan=2,sticky='w',padx=6,pady=(6,4))
        countrow=ttk.Frame(dynamic);countrow.grid(row=1,column=0,columnspan=2,sticky='w',padx=6,pady=4);ttk.Label(countrow,text='How many Signal + Reference pairs are in each file?').pack(side=tk.LEFT);sp=ttk.Spinbox(countrow,from_=1,to=8,textvariable=count_var,width=5,command=rebuild_pairs);sp.pack(side=tk.LEFT,padx=8)
        pair_holder.grid(row=2,column=0,columnspan=2,sticky='ew',padx=6,pady=4);rebuild_pairs()
        try:sp.bind('<Return>',rebuild_pairs);sp.bind('<FocusOut>',rebuild_pairs)
        except Exception:pass
    else:
        ttk.Label(dynamic,text='Use this when Signal and Reference measurements take turns in the rows. One column says which kind of measurement each row is.',wraplength=820,justify='left').grid(row=0,column=0,columnspan=2,sticky='w',padx=6,pady=(6,4))
        ttk.Label(dynamic,text='Which column tells FiberLyse whether a row is Signal or Reference?').grid(row=1,column=0,sticky='w',padx=6,pady=4);ttk.Combobox(dynamic,textvariable=rowtype_var,values=columns,state='readonly',width=30).grid(row=1,column=1,sticky='ew',padx=6,pady=4)
        codeframe=ttk.Frame(dynamic);codeframe.grid(row=2,column=0,columnspan=2,sticky='w',padx=6,pady=4)
        ttk.Label(codeframe,text='Value meaning Signal:').pack(side=tk.LEFT);ttk.Entry(codeframe,textvariable=signal_code_var,width=10).pack(side=tk.LEFT,padx=(4,16));ttk.Label(codeframe,text='Value meaning Reference/control:').pack(side=tk.LEFT);ttk.Entry(codeframe,textvariable=reference_code_var,width=10).pack(side=tk.LEFT,padx=4)
        ttk.Label(dynamic,text='Which column(s) contain the recorded brightness/fluorescence values? Select one or more:',wraplength=430,justify='left').grid(row=3,column=0,sticky='nw',padx=6,pady=4)
        lf=ttk.Frame(dynamic);lf.grid(row=3,column=1,sticky='nsew',padx=6,pady=4);measurement_list=tk.Listbox(lf,selectmode=tk.EXTENDED,height=6,exportselection=False);sb=ttk.Scrollbar(lf,orient='vertical',command=measurement_list.yview);measurement_list.configure(yscrollcommand=sb.set);measurement_list.pack(side=tk.LEFT,fill=tk.BOTH,expand=True);sb.pack(side=tk.RIGHT,fill=tk.Y)
        for c in columns:measurement_list.insert(tk.END,c)
        # Select likely fluorescence columns but never silently commit them.
        for i,c in enumerate(columns):
            lc=c.lower()
            if re.match(r'^g\d+$', c, flags=re.I) or any(k in lc for k in ['fluor','signal','470','465']):measurement_list.selection_set(i)

    startbox=ttk.LabelFrame(form,text='3. Where should the recording begin?');startbox.grid(row=row,column=0,columnspan=2,sticky='ew',pady=(10,6));startbox.columnconfigure(1,weight=1);row+=1
    ttk.Radiobutton(startbox,text='The recording begins on the first row of the CSV',variable=start_rule_var,value='first_row').grid(row=0,column=0,columnspan=2,sticky='w',padx=6,pady=3)
    ttk.Radiobutton(startbox,text='The recording begins when a column reaches a particular value',variable=start_rule_var,value='column_value').grid(row=1,column=0,columnspan=2,sticky='w',padx=6,pady=3)
    ttk.Label(startbox,text='Column:').grid(row=2,column=0,sticky='e',padx=(20,4),pady=3);ttk.Combobox(startbox,textvariable=start_col_var,values=columns,state='readonly',width=28).grid(row=2,column=1,sticky='w',pady=3)
    ttk.Label(startbox,text='Value that means “recording starts”:').grid(row=3,column=0,sticky='e',padx=(20,4),pady=3);ttk.Entry(startbox,textvariable=start_val_var,width=18).grid(row=3,column=1,sticky='w',pady=3)
    ttk.Label(form,text='4. How much of the beginning should be ignored before any FiberLyse analysis?',wraplength=430,justify='left').grid(row=row,column=0,sticky='w',pady=5);exfr=ttk.Frame(form);exfr.grid(row=row,column=1,sticky='w',padx=(12,0),pady=5);ttk.Entry(exfr,textvariable=exclude_var,width=10).pack(side=tk.LEFT);ttk.Label(exfr,text='seconds').pack(side=tk.LEFT,padx=4);row+=1
    ttk.Label(form,text='Name this CSV setup:').grid(row=row,column=0,sticky='w',pady=5);ttk.Entry(form,textvariable=profile_name_var,width=32).grid(row=row,column=1,sticky='w',padx=(12,0),pady=5);row+=1
    remember_var=tk.BooleanVar(value=True);ttk.Checkbutton(form,text='Remember this setup next time FiberLyse starts',variable=remember_var).grid(row=row,column=0,columnspan=2,sticky='w',pady=(8,4));row+=1
    result: Dict[str, Any] = {}

    def save_and_close():
        try:
            exclude=max(0.0,float(exclude_var.get()))
        except Exception:
            messagebox.showerror('Check ignored time','Please enter the ignored beginning in seconds, for example 0 or 600.',parent=top);return
        p={'profile_name':profile_name_var.get().strip() or 'My CSV setup','layout':layout,'time_column':time_var.get(),'time_units':units_var.get(),'start_rule':start_rule_var.get(),'start_column':start_col_var.get(),'start_value':start_val_var.get(),'exclude_initial_seconds':exclude,'remember':bool(remember_var.get()),'example_file':example_path}
        if layout == _FIBERLYSE_V231_LAYOUT_SEPARATE:
            rebuild_pairs();pairs=[]
            for i,(sv,rv,lv) in enumerate(pair_vars,start=1):
                if not sv.get() or not rv.get():messagebox.showerror('Missing column',f'Choose both Signal and Reference/control columns for recording channel {i}.',parent=top);return
                if sv.get()==rv.get():messagebox.showerror('Check columns',f'Recording channel {i} uses the same column for Signal and Reference/control. Choose two different columns.',parent=top);return
                pairs.append({'signal_column':sv.get(),'reference_column':rv.get(),'label':lv.get().strip() or f'Channel {i}'})
            p['channel_pairs']=pairs
        else:
            p['row_type_column']=rowtype_var.get();p['signal_row_value']=signal_code_var.get().strip();p['reference_row_value']=reference_code_var.get().strip()
            sel=list(measurement_list.curselection()) if measurement_list is not None else []
            cols=[columns[int(i)] for i in sel]
            if not cols:messagebox.showerror('Missing fluorescence column','Select at least one column containing fluorescence values.',parent=top);return
            if not p['row_type_column'] or not p['signal_row_value'] or not p['reference_row_value']:messagebox.showerror('Missing row information','Choose the row-type column and enter the values meaning Signal and Reference/control.',parent=top);return
            p['measurement_columns']=cols
        # Validate against the example before accepting the setup.
        try:
            fake_settings=_fiberlyse_parse_analysis_settings(app) or {'artifact_enabled':False,'artifact_factor':11.9,'artifact_pad':1,'require_shared':True,'align_mode':DEFAULT_ALIGN_MODE,'acq_fps':40.0,'smooth_win':20,'mode':NORM_DFF,'interval_start':exclude,'interval_end':6500.0,'use_interp':True,'fit_windows':[(exclude,6500.0)]}
            fake_settings=dict(fake_settings);fake_settings['fit_windows']=[(exclude,6500.0)];fake_settings['interval_start']=max(exclude,float(fake_settings.get('interval_start',exclude)))
            test_results=_fiberlyse_v231_analyze_file(example_path,p,fake_settings)
            if not test_results:raise ValueError('No usable channels were produced.')
        except Exception as e:
            messagebox.showerror('This setup does not match the example CSV',f'FiberLyse tried the mapping before saving it and found a problem:\n\n{e}',parent=top);return
        result.update(p)
        top.destroy()
    btns=ttk.Frame(outer);btns.pack(fill=tk.X,pady=(10,0));ttk.Button(btns,text='Cancel',command=top.destroy).pack(side=tk.RIGHT);ttk.Button(btns,text='Use this CSV setup',command=save_and_close).pack(side=tk.RIGHT,padx=(0,8))
    _fiberlyse_v23_apply_classic_widget_theme(top,_fiberlyse_v23_theme_tokens(app=app))
    app.root.wait_window(top)
    return dict(result) if result else None


def _fiberlyse_v231_show_setup_chooser(app, first_run: bool=False) -> None:
    top=tk.Toplevel(app.root);top.title('FiberLyse CSV setup');top.transient(app.root);top.grab_set();top.geometry('700x500');top.resizable(True,True)
    frm=ttk.Frame(top,padding=16);frm.pack(fill=tk.BOTH,expand=True)
    ttk.Label(frm,text='How is your CSV organized?',font=('TkDefaultFont',15,'bold')).pack(anchor='w')
    ttk.Label(frm,text='Pick the description that looks most like your files. The wording below describes where the numbers are stored — you do not need to know fiber-photometry terminology.',wraplength=650,justify='left').pack(anchor='w',pady=(4,14))
    choice=tk.StringVar(value=_FIBERLYSE_V231_LAYOUT_BUILTIN)
    options=[
        (_FIBERLYSE_V231_LAYOUT_BUILTIN,'Neurophotometrics / Bonsai','Choose this if your CSV contains columns such as SystemTimestamp, LedState and G0/G1. FiberLyse already knows what those columns mean.'),
        (_FIBERLYSE_V231_LAYOUT_SEPARATE,'Main signal and comparison/reference are in separate columns','Choose this if each row already contains a time, a column for the signal you want to study, and a separate reference/control column.'),
        (_FIBERLYSE_V231_LAYOUT_INTERLEAVED,'Main signal and comparison/reference take turns in the rows','Choose this if one row is a Signal measurement, the next row is a Reference/control measurement, and another column tells you which kind each row is.'),
    ]
    for value,title,desc in options:
        box=ttk.Frame(frm);box.pack(fill=tk.X,pady=5);ttk.Radiobutton(box,text=title,variable=choice,value=value).pack(anchor='w');ttk.Label(box,text=desc,wraplength=610,justify='left',foreground='gray40').pack(anchor='w',padx=(26,0),pady=(1,0))
    remember_builtin=tk.BooleanVar(value=True)
    ttk.Checkbutton(frm,text='Remember this choice next time FiberLyse starts',variable=remember_builtin).pack(anchor='w',pady=(14,0))
    def apply_choice():
        c=choice.get()
        if c==_FIBERLYSE_V231_LAYOUT_BUILTIN:
            p=_fiberlyse_v231_builtin_profile();p['remember']=bool(remember_builtin.get());app.import_profile=p
            if p['remember']:_fiberlyse_v231_save_profile(p)
            else:_fiberlyse_v231_forget_saved_profile()
            _fiberlyse_v231_after_profile_changed(app);top.destroy();return
        # close chooser before opening detailed mapping to avoid nested grabs
        top.grab_release();top.destroy()
        p=_fiberlyse_v231_configure_custom(app,c)
        if p:
            app.import_profile=p
            if p.get('remember'):_fiberlyse_v231_save_profile(p)
            else:_fiberlyse_v231_forget_saved_profile()
            _fiberlyse_v231_after_profile_changed(app)
    btn=ttk.Frame(frm);btn.pack(side=tk.BOTTOM,fill=tk.X,pady=(16,0));ttk.Button(btn,text='Cancel',command=top.destroy).pack(side=tk.RIGHT);ttk.Button(btn,text='Continue',command=apply_choice).pack(side=tk.RIGHT,padx=(0,8))
    _fiberlyse_v23_apply_classic_widget_theme(top,_fiberlyse_v23_theme_tokens(app=app));app.root.wait_window(top)


def _fiberlyse_v231_after_profile_changed(app) -> None:
    p=getattr(app,'import_profile',None) or _fiberlyse_v231_builtin_profile();app.import_profile=p
    try:
        app.lbl_import_profile.configure(text=f"CSV setup: {_fiberlyse_v231_profile_label(p)}")
        _fiberlyse_attach_tooltip(app.lbl_import_profile,_fiberlyse_v231_profile_summary(p))
    except Exception:pass
    # Make the default baseline interval follow the amount intentionally removed
    # by the selected importer instead of assuming every hardware setup needs 600 s.
    try:
        exc=float(p.get('exclude_initial_seconds',0.0));app.var_interval_start.set(f'{exc:g}')
    except Exception:pass
    try:app.status.set('CSV setup active: '+_fiberlyse_v231_profile_summary(p))
    except Exception:pass


def _fiberlyse_v231_show_current_setup(app) -> None:
    p=getattr(app,'import_profile',None) or _fiberlyse_v231_builtin_profile()
    msg=_fiberlyse_v231_profile_summary(p)
    if p.get('layout')==_FIBERLYSE_V231_LAYOUT_SEPARATE:
        lines=[]
        for pair in p.get('channel_pairs',[]):lines.append(f"• {pair.get('label','Channel')}: Signal = {pair.get('signal_column')} | Reference = {pair.get('reference_column')}")
        if lines:msg += '\n\n'+'\n'.join(lines)
    elif p.get('layout')==_FIBERLYSE_V231_LAYOUT_INTERLEAVED:
        msg += f"\n\nSignal row value: {p.get('signal_row_value')}\nReference row value: {p.get('reference_row_value')}\nFluorescence columns: {', '.join(p.get('measurement_columns',[]))}"
    messagebox.showinfo('Current CSV setup',msg,parent=app.root)


# Override only the file-analysis dispatcher. All V23 queue/progress behavior is
# retained, but each source file is now read through the selected import profile.
def _fiberlyse_v231_run_analysis_paths(self, selected_paths: Optional[List[str]]=None) -> None:
    _fiberlyse_sync_queue_from_paths(self)
    all_paths=list(getattr(self,'csv_paths',[]) or [])
    if not all_paths:
        try:messagebox.showwarning('No files','Please add one or more CSV files first.')
        except Exception:pass
        return
    paths_to_analyze=[_fiberlyse_norm_path(p) for p in (selected_paths if selected_paths is not None else all_paths)];paths_to_analyze=[p for p in paths_to_analyze if p in set(all_paths)]
    if not paths_to_analyze:
        try:messagebox.showinfo('Analyze selected','No selected files are available to analyze.')
        except Exception:pass
        return
    settings=_fiberlyse_parse_analysis_settings(self)
    if settings is None:return
    profile=dict(getattr(self,'import_profile',None) or _fiberlyse_v231_builtin_profile())
    try:exclude=float(profile.get('exclude_initial_seconds',ANALYSIS_EXCLUDE_INITIAL_SECONDS if profile.get('layout')==_FIBERLYSE_V231_LAYOUT_BUILTIN else 0.0))
    except Exception:exclude=0.0
    settings=dict(settings);settings['fit_windows']=[(exclude,6500.0)]
    if float(settings.get('interval_start',exclude)) < exclude:settings['interval_start']=exclude
    recs_snapshot={p:dict(r) for p,r in _fiberlyse_existing_records(self).items()};file_meta_all=_fiberlyse_file_meta_all(self);mode=settings['mode'];self._analysis_cancel_requested=False;_fiberlyse_set_buttons_for_analysis(self,True)
    try:self.progress_bar.config(maximum=max(1,len(paths_to_analyze)));self.progress_var.set(0.0)
    except Exception:pass
    try:self.status.set(f"Analyzing {len(paths_to_analyze)} file(s) using CSV setup: {_fiberlyse_v231_profile_label(profile)}…")
    except Exception:pass
    for p in paths_to_analyze:_fiberlyse_mark_queue_status(self,p,'queued',channels=recs_snapshot.get(p,{}).get('channels',''),note='waiting')
    def worker():
        local_results=dict(getattr(self,'_results',{}) or {});local_display=dict(getattr(self,'_mouse_display',{}) or {});successes=0;failures=0;canceled=False
        for idx,path in enumerate(paths_to_analyze,start=1):
            if bool(getattr(self,'_analysis_cancel_requested',False)):
                canceled=True;self.root.after(0,lambda p=path:_fiberlyse_mark_queue_status(self,p,'canceled',note='not run'));continue
            rec=recs_snapshot.get(path,{});file_key=str(rec.get('key') or f'F{idx}');alias=str(rec.get('alias') or _basename_no_ext(path));self.root.after(0,lambda idx=idx,alias=alias:self.status.set(f'Analyzing {idx}/{len(paths_to_analyze)}: {alias}…'));self.root.after(0,lambda p=path:_fiberlyse_mark_queue_status(self,p,'analyzing',note='running'))
            try:
                per_file=_fiberlyse_v231_analyze_file(path,profile,settings);_fiberlyse_remove_results_for_file(local_results,local_display,file_key)
                for gcol in sorted(per_file.keys()):
                    mid=f'{file_key}:{gcol}';local_results[mid]=per_file[gcol];local_display[mid]=f'{alias}:{gcol}'
                successes+=1;self.root.after(0,lambda p=path,n=len(per_file):_fiberlyse_mark_queue_status(self,p,'analyzed',channels=n,note=''))
            except Exception as exc:
                _fiberlyse_remove_results_for_file(local_results,local_display,file_key);failures+=1;msg=str(exc);self.root.after(0,lambda p=path,msg=msg:_fiberlyse_mark_queue_status(self,p,'failed',channels='',note=msg))
            finally:self.root.after(0,lambda v=idx:self.progress_var.set(float(v)) if hasattr(self,'progress_var') else None)
        order=[]
        for _key,_alias,_path in file_meta_all:order.extend(sorted([mid for mid in local_results if str(mid).startswith(str(_key)+':')],key=lambda s:str(s).split(':',1)[-1]))
        def finish():
            _fiberlyse_set_buttons_for_analysis(self,False)
            try:self.btn_cancel.config(state=tk.DISABLED)
            except Exception:pass
            if local_results:
                self.on_analysis_finished(local_results,local_display,order,mode,file_meta_all);suffix=' Canceled.' if canceled else ''
                try:self.status.set(f"Done. Analyzed {successes} file(s); failed {failures}. CSV setup: {_fiberlyse_v231_profile_label(profile)}.{suffix}")
                except Exception:pass
            else:
                self._results={};self._mouse_display={};self._mouse_order=[]
                try:self.build_tabs(norm_mode=mode)
                except Exception:pass
                msg='No files were analyzed successfully.'+(' Analysis was canceled.' if canceled else '')
                try:self.status.set(msg)
                except Exception:pass
                if failures:
                    try:messagebox.showwarning('Analysis completed with errors',msg+'\n\nCheck the queue Notes column. If the column names differ from your setup, use CSV setup… and remap an example file.')
                    except Exception:pass
            _fiberlyse_refresh_queue_tree(self)
        self.root.after(0,finish)
    self._analysis_thread=threading.Thread(target=worker,daemon=True);self._analysis_thread.start()

MainAppTk.run_analysis = lambda self: _fiberlyse_v231_run_analysis_paths(self, selected_paths=None)
MainAppTk.run_analysis_selected = lambda self: _fiberlyse_v231_run_analysis_paths(self, selected_paths=_fiberlyse_selected_queue_paths(self)) if _fiberlyse_selected_queue_paths(self) else messagebox.showinfo('Analyze selected','Select one or more files in the queue first.')
MainAppTk.retry_failed_files = lambda self: _fiberlyse_v231_run_analysis_paths(self, selected_paths=[p for p,r in _fiberlyse_existing_records(self).items() if r.get('status')=='failed'])

# Extend export metadata with import provenance without changing graph data.
_FIBERLYSE_V231_META_BEFORE = ChannelTabsTk._meta_df
def _fiberlyse_v231_meta_df(self) -> pd.DataFrame:
    base=_FIBERLYSE_V231_META_BEFORE(self)
    r=self.res
    rows=[
        ('import_profile_name',getattr(r,'import_profile_name','Neurophotometrics / Bonsai')),
        ('import_layout',getattr(r,'import_layout',_FIBERLYSE_V231_LAYOUT_BUILTIN)),
        ('import_time_column',getattr(r,'import_time_column','SystemTimestamp')),
        ('import_time_units',getattr(r,'import_time_units','Seconds')),
        ('import_signal_column',getattr(r,'import_signal_column',getattr(r,'gcol',''))),
        ('import_reference_column',getattr(r,'import_reference_column',getattr(r,'gcol',''))),
    ]
    return pd.concat([base,pd.DataFrame(rows,columns=['key','value'])],ignore_index=True)
ChannelTabsTk._meta_df = _fiberlyse_v231_meta_df

# Main GUI integration: a visible CSV setup control plus a first-run chooser.
_FIBERLYSE_V231_MAIN_INIT_BEFORE = MainAppTk.__init__
def _fiberlyse_v231_main_init(self,*args,**kwargs):
    _FIBERLYSE_V231_MAIN_INIT_BEFORE(self,*args,**kwargs)
    saved=_fiberlyse_v231_load_saved_profile();self.import_profile=dict(saved or _fiberlyse_v231_builtin_profile())
    try:
        top=self.btn_run.master
        self.btn_csv_setup=ttk.Button(top,text='CSV setup…',command=lambda app=self:_fiberlyse_v231_show_setup_chooser(app,first_run=False));self.btn_csv_setup.grid(row=0,column=7,sticky='e',padx=(8,0))
        self.lbl_import_profile=ttk.Label(top,text='');self.lbl_import_profile.grid(row=1,column=6,columnspan=2,sticky='e',padx=(8,0),pady=(2,0))
        _fiberlyse_attach_tooltip(self.btn_csv_setup,'Tell FiberLyse what the columns in your CSV mean. Use this when changing acquisition hardware or CSV layout.')
    except Exception:pass
    _fiberlyse_v231_after_profile_changed(self)
    if saved is None:
        # Opening after Tk has drawn avoids a blank root window behind the modal.
        try:self.root.after(180,lambda app=self:_fiberlyse_v231_show_setup_chooser(app,first_run=True))
        except Exception:pass
    try:self.status.set('Ready. V23.1: choose or define a CSV setup, then analyze as normal.')
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_v231_main_init

# Responsive top-bar patch: keep CSV setup discoverable even on smaller screens.
try:
    _FIBERLYSE_V231_RESPONSIVE_BEFORE = MainAppTk._fiberlyse_apply_responsive_layout
except Exception:
    _FIBERLYSE_V231_RESPONSIVE_BEFORE = None
if _FIBERLYSE_V231_RESPONSIVE_BEFORE is not None:
    def _fiberlyse_v231_responsive(self,*args,**kwargs):
        out=_FIBERLYSE_V231_RESPONSIVE_BEFORE(self,*args,**kwargs)
        try:
            w=int(self.root.winfo_width())
            if w < 1050:
                self.btn_csv_setup.grid(row=2,column=3,columnspan=2,sticky='ew',padx=(4,0),pady=(2,0));self.lbl_import_profile.grid(row=3,column=0,columnspan=5,sticky='w',pady=(2,0))
            else:
                self.btn_csv_setup.grid(row=0,column=7,sticky='e',padx=(8,0),pady=0);self.lbl_import_profile.grid(row=1,column=6,columnspan=2,sticky='e',padx=(8,0),pady=(2,0))
        except Exception:pass
        return out
    MainAppTk._fiberlyse_apply_responsive_layout = _fiberlyse_v231_responsive


# The performance layer's resize callbacks call the module-level responsive
# function directly, so wrap that function as well (rather than only a class
# method) to keep the new CSV setup controls visible on laptop-sized screens.
_FIBERLYSE_V231_GLOBAL_RESPONSIVE_BEFORE = _fiberlyse_apply_responsive_layout
def _fiberlyse_apply_responsive_layout(self, force: bool=False) -> None:
    _FIBERLYSE_V231_GLOBAL_RESPONSIVE_BEFORE(self, force=force)
    try:w=int(self.root.winfo_width())
    except Exception:w=1200
    try:
        if w < 760:
            self.chk_night_mode.grid(row=3,column=0,columnspan=2,sticky='w',padx=(0,4),pady=(2,0))
            self.btn_csv_setup.grid(row=3,column=2,columnspan=3,sticky='ew',padx=(4,0),pady=(2,0))
            self.lbl_import_profile.grid(row=4,column=0,columnspan=5,sticky='w',pady=(2,0))
        elif w < 980:
            self.chk_night_mode.grid(row=2,column=0,columnspan=2,sticky='w',padx=(0,4),pady=(2,0))
            self.btn_csv_setup.grid(row=2,column=2,columnspan=3,sticky='ew',padx=(4,0),pady=(2,0))
            self.lbl_import_profile.grid(row=3,column=0,columnspan=5,sticky='w',pady=(2,0))
        else:
            self.chk_night_mode.grid(row=0,column=6,sticky='e',padx=(10,0),pady=0)
            self.btn_csv_setup.grid(row=0,column=7,sticky='e',padx=(8,0),pady=0)
            self.lbl_import_profile.grid(row=1,column=6,columnspan=2,sticky='e',padx=(8,0),pady=(2,0))
    except Exception:pass


# ---- V23.1 CSV format/dialect compatibility patch --------------------------
# Custom import profiles now retain how a CSV is encoded: what separates
# columns and which decimal mark is used.  This supports, for example,
# semicolon-separated European CSV files with decimal commas.
import csv as _fiberlyse_v231_csv
import io as _fiberlyse_v231_io

_FIBERLYSE_V231_SEPARATOR_CHOICES = (
    'Comma (,)',
    'Semicolon (;)',
    'Tab',
    'Pipe (|)',
)
_FIBERLYSE_V231_DECIMAL_CHOICES = (
    'Period (12.34)',
    'Comma (12,34)',
)


def _fiberlyse_v231_separator_from_choice(choice: str) -> str:
    c=str(choice or '').strip().lower()
    if 'semicolon' in c or c == ';':return ';'
    if 'tab' in c or c == '\\t':return '\t'
    if 'pipe' in c or c == '|':return '|'
    return ','


def _fiberlyse_v231_separator_choice(separator: str) -> str:
    s=str(separator or ',')
    if s == ';':return 'Semicolon (;)'
    if s == '\t':return 'Tab'
    if s == '|':return 'Pipe (|)'
    return 'Comma (,)'


def _fiberlyse_v231_decimal_from_choice(choice: str) -> str:
    c=str(choice or '').strip().lower()
    return ',' if ('comma' in c or c == ',') else '.'


def _fiberlyse_v231_decimal_choice(decimal_mark: str) -> str:
    return 'Comma (12,34)' if str(decimal_mark or '.') == ',' else 'Period (12.34)'


def _fiberlyse_v231_read_text_sample(path: str, max_bytes: int=160000) -> Tuple[str, str]:
    with open(path,'rb') as fh:raw=fh.read(max_bytes)
    last=None
    for enc in ('utf-8-sig','utf-8','cp1252','latin-1'):
        try:return raw.decode(enc), enc
        except Exception as e:last=e
    raise ValueError(f'Could not decode the CSV text: {last}')


def _fiberlyse_v231_detect_csv_format(path: str) -> Dict[str, str]:
    text, encoding = _fiberlyse_v231_read_text_sample(path)
    lines=[ln for ln in text.splitlines() if ln.strip()][:40]
    if not lines:raise ValueError('The selected CSV appears to be empty.')
    sample='\n'.join(lines)
    candidates=[',',';','\t','|']
    best_sep=',';best_score=-1e18
    for sep in candidates:
        try:
            rows=list(_fiberlyse_v231_csv.reader(_fiberlyse_v231_io.StringIO(sample),delimiter=sep))
        except Exception:
            continue
        widths=[len(r) for r in rows if r]
        if not widths:continue
        # Reward multiple, consistent columns and strongly reward a header whose
        # width agrees with the data rows.  This avoids mistaking decimal commas
        # for delimiters in semicolon-separated European CSV files.
        vals={w:widths.count(w) for w in set(widths)}
        mode_w=max(vals,key=vals.get);consistency=vals[mode_w]/float(len(widths))
        header_w=widths[0]
        score=(max(0,mode_w-1)*12.0)+(consistency*10.0)
        if header_w == mode_w:score += 12.0
        else:score -= 8.0*abs(header_w-mode_w)
        if mode_w <= 1:score -= 25.0
        if score > best_score:best_sep,best_score=sep,score
    # Inspect already-separated cells to determine whether numeric values use
    # a comma or period as their decimal mark.
    comma_hits=0;dot_hits=0
    num_comma=re.compile(r'^[-+]?\d+,\d+(?:[eE][-+]?\d+)?$')
    num_dot=re.compile(r'^[-+]?\d+\.\d+(?:[eE][-+]?\d+)?$')
    try:rows=list(_fiberlyse_v231_csv.reader(_fiberlyse_v231_io.StringIO(sample),delimiter=best_sep))
    except Exception:rows=[]
    for row in rows[1:]:
        for cell in row:
            v=str(cell).strip().strip('"').strip("'")
            if num_comma.match(v):comma_hits += 1
            if num_dot.match(v):dot_hits += 1
    decimal_mark=',' if comma_hits > dot_hits else '.'
    return {'csv_separator':best_sep,'decimal_mark':decimal_mark,'csv_encoding':encoding}


def _fiberlyse_v231_read_csv_with_format(path: str, profile: Optional[Dict[str,Any]]=None, **kwargs) -> pd.DataFrame:
    p=dict(profile or {})
    sep=str(p.get('csv_separator', ',') or ',')
    dec=str(p.get('decimal_mark', '.') or '.')
    enc=str(p.get('csv_encoding', 'utf-8-sig') or 'utf-8-sig')
    try:
        return pd.read_csv(path, sep=sep, decimal=dec, encoding=enc, **kwargs)
    except UnicodeDecodeError:
        # A saved profile can be moved between systems. If encoding is the only
        # mismatch, redetect it while preserving the user's separator/decimal.
        detected=_fiberlyse_v231_detect_csv_format(path)
        return pd.read_csv(path, sep=sep, decimal=dec, encoding=detected.get('csv_encoding','utf-8-sig'), **kwargs)


def _fiberlyse_v231_read_example(path: str, profile: Optional[Dict[str,Any]]=None) -> Tuple[List[str], pd.DataFrame]:
    p=dict(profile or _fiberlyse_v231_detect_csv_format(path))
    df=_fiberlyse_v231_read_csv_with_format(path,p,nrows=8)
    return [str(c) for c in df.columns],df


def _fiberlyse_v231_unique_column_values(path: str, column: str, profile: Dict[str,Any], max_rows: int=2000) -> List[str]:
    if not column:return []
    try:
        d=_fiberlyse_v231_read_csv_with_format(path,profile,usecols=[column],nrows=max_rows)
        vals=[]
        for v in d[column].dropna().tolist():
            if isinstance(v,(float,np.floating)) and np.isfinite(v) and float(v).is_integer():txt=str(int(v))
            else:txt=str(v).strip()
            if txt and txt not in vals:vals.append(txt)
            if len(vals)>=80:break
        return vals
    except Exception:return []


def _fiberlyse_v231_confirm_csv_format(app, path: str, detected: Dict[str,str]) -> Optional[Dict[str,str]]:
    top=tk.Toplevel(app.root);top.title('How is this CSV written?');top.transient(app.root);top.grab_set();top.resizable(True,True)
    try:
        sw=max(700,int(top.winfo_screenwidth()));sh=max(520,int(top.winfo_screenheight()));top.geometry(f'{min(820,sw-80)}x{min(620,sh-100)}')
    except Exception:top.geometry('820x600')
    frm=ttk.Frame(top,padding=14);frm.pack(fill=tk.BOTH,expand=True)
    ttk.Label(frm,text='First, let FiberLyse read the CSV correctly',font=('TkDefaultFont',14,'bold')).pack(anchor='w')
    ttk.Label(frm,text='Different machines save CSV files differently. FiberLyse has guessed the two settings below. If the preview looks like a normal table, you can simply press Next.',wraplength=760,justify='left').pack(anchor='w',pady=(4,10))
    sep_var=tk.StringVar(value=_fiberlyse_v231_separator_choice(detected.get('csv_separator',',')))
    dec_var=tk.StringVar(value=_fiberlyse_v231_decimal_choice(detected.get('decimal_mark','.')))
    detected_lbl=ttk.Label(frm,text='',justify='left');detected_lbl.pack(anchor='w',pady=(0,8))
    settings=ttk.LabelFrame(frm,text='CSV number and column format');settings.pack(fill=tk.X,pady=(0,8));settings.columnconfigure(1,weight=1)
    ttk.Label(settings,text='What separates one column from the next?').grid(row=0,column=0,sticky='w',padx=8,pady=6)
    sep_combo=ttk.Combobox(settings,textvariable=sep_var,values=_FIBERLYSE_V231_SEPARATOR_CHOICES,state='readonly',width=24);sep_combo.grid(row=0,column=1,sticky='w',padx=8,pady=6)
    ttk.Label(settings,text='How are decimal numbers written?').grid(row=1,column=0,sticky='w',padx=8,pady=6)
    dec_combo=ttk.Combobox(settings,textvariable=dec_var,values=_FIBERLYSE_V231_DECIMAL_CHOICES,state='readonly',width=24);dec_combo.grid(row=1,column=1,sticky='w',padx=8,pady=6)
    preview_box=ttk.LabelFrame(frm,text='Preview');preview_box.pack(fill=tk.BOTH,expand=True,pady=(4,8))
    txt=tk.Text(preview_box,height=10,wrap='none');ys=ttk.Scrollbar(preview_box,orient='vertical',command=txt.yview);xs=ttk.Scrollbar(preview_box,orient='horizontal',command=txt.xview);txt.configure(yscrollcommand=ys.set,xscrollcommand=xs.set);txt.grid(row=0,column=0,sticky='nsew',padx=(6,0),pady=(6,0));ys.grid(row=0,column=1,sticky='ns',pady=(6,0));xs.grid(row=1,column=0,sticky='ew',padx=(6,0),pady=(0,6));preview_box.rowconfigure(0,weight=1);preview_box.columnconfigure(0,weight=1)
    status=tk.StringVar(value='');ttk.Label(frm,textvariable=status,wraplength=760,justify='left').pack(anchor='w')
    current: Dict[str,str]={}
    def refresh(*_):
        prof={'csv_separator':_fiberlyse_v231_separator_from_choice(sep_var.get()),'decimal_mark':_fiberlyse_v231_decimal_from_choice(dec_var.get()),'csv_encoding':detected.get('csv_encoding','utf-8-sig')}
        try:
            cols,df=_fiberlyse_v231_read_example(path,prof)
            if len(cols)<2:raise ValueError('Only one column was found. Try a different column-separator option.')
            txt.configure(state='normal');txt.delete('1.0',tk.END);txt.insert('1.0',_fiberlyse_v231_preview_text(df));txt.configure(state='disabled')
            current.clear();current.update(prof)
            sep_show={'\t':'Tab'}.get(prof['csv_separator'],repr(prof['csv_separator']))
            dec_show='comma (12,34)' if prof['decimal_mark']==',' else 'period (12.34)'
            detected_lbl.configure(text=f"FiberLyse initially detected: columns separated by {_fiberlyse_v231_separator_choice(detected.get('csv_separator',','))}; decimals use {_fiberlyse_v231_decimal_choice(detected.get('decimal_mark','.'))}.")
            status.set(f'✓ Preview read successfully: {len(cols)} columns found. Current separator: {sep_show}; decimal style: {dec_show}.')
        except Exception as e:
            current.clear();txt.configure(state='normal');txt.delete('1.0',tk.END);txt.insert('1.0','Could not create a table with these settings.');txt.configure(state='disabled');status.set(f'⚠ {e}')
    sep_combo.bind('<<ComboboxSelected>>',refresh);dec_combo.bind('<<ComboboxSelected>>',refresh);refresh()
    result: Dict[str,str]={}
    def accept():
        refresh()
        if not current:
            messagebox.showerror('CSV format is not readable','Adjust the two CSV-format choices until the preview shows separate columns.',parent=top);return
        result.update(current);top.destroy()
    actions=ttk.Frame(frm);actions.pack(fill=tk.X,pady=(8,0));ttk.Button(actions,text='Cancel',command=top.destroy).pack(side=tk.RIGHT);ttk.Button(actions,text='Next: define columns →',command=accept).pack(side=tk.RIGHT,padx=(0,8))
    top.bind('<Return>',lambda _e:accept());top.bind('<Escape>',lambda _e:top.destroy())
    try:_fiberlyse_v23_apply_classic_widget_theme(top,_fiberlyse_v23_theme_tokens(app=app))
    except Exception:pass
    app.root.wait_window(top)
    return dict(result) if result else None


# Replace the custom analyzers so the real analysis uses exactly the same CSV
# separator/decimal settings used by the wizard preview.
def _fiberlyse_v231_analyze_separate_columns(csv_path: str, profile: Dict[str, Any], settings: Dict[str, Any]) -> Dict[str, ChannelResult]:
    header=_fiberlyse_v231_read_csv_with_format(csv_path,profile,nrows=0)
    pairs=list(profile.get('channel_pairs') or [])
    if not pairs:raise ValueError("No Signal/Reference column pair is defined. Open 'CSV setup…'.")
    required={str(profile.get('time_column') or '')}
    if str(profile.get('start_rule','first_row'))=='column_value':required.add(str(profile.get('start_column') or ''))
    for pair in pairs:required.add(str(pair.get('signal_column') or ''));required.add(str(pair.get('reference_column') or ''))
    missing=sorted(c for c in required if c and c not in header.columns)
    if missing:raise ValueError('This CSV does not match the active setup. Missing column(s): '+', '.join(missing))
    df=_fiberlyse_v231_read_csv_with_format(csv_path,profile,usecols=[c for c in header.columns if c in required])
    df,elapsed,_exclude=_fiberlyse_v231_prepare_time(df,profile)
    results: Dict[str,ChannelResult]={}
    for i,pair in enumerate(pairs,start=1):
        sig_col=str(pair.get('signal_column') or '');ref_col=str(pair.get('reference_column') or '');label=str(pair.get('label') or f'Channel{i}')
        sig=pd.to_numeric(df[sig_col],errors='coerce').to_numpy(dtype=float);ref=pd.to_numeric(df[ref_col],errors='coerce').to_numpy(dtype=float);key=_fiberlyse_v231_safe_channel_name(label,i)
        results[key]=_fiberlyse_v231_make_result(csv_path,key,elapsed,elapsed,ref,sig,settings,profile,{'import_signal_column':sig_col,'import_reference_column':ref_col,'import_channel_label':label,'import_csv_separator':str(profile.get('csv_separator',',')),'import_decimal_mark':str(profile.get('decimal_mark','.'))})
    return results


def _fiberlyse_v231_analyze_alternating_rows(csv_path: str, profile: Dict[str, Any], settings: Dict[str, Any]) -> Dict[str, ChannelResult]:
    header=_fiberlyse_v231_read_csv_with_format(csv_path,profile,nrows=0)
    value_cols=[str(c) for c in list(profile.get('measurement_columns') or []) if str(c)]
    if not value_cols:raise ValueError("No fluorescence-value column is defined. Open 'CSV setup…'.")
    row_type_col=str(profile.get('row_type_column') or '');required={str(profile.get('time_column') or ''),row_type_col,*value_cols}
    if str(profile.get('start_rule','first_row'))=='column_value':required.add(str(profile.get('start_column') or ''))
    missing=sorted(c for c in required if c and c not in header.columns)
    if missing:raise ValueError('This CSV does not match the active setup. Missing column(s): '+', '.join(missing))
    df=_fiberlyse_v231_read_csv_with_format(csv_path,profile,usecols=[c for c in header.columns if c in required]);df,elapsed,_exclude=_fiberlyse_v231_prepare_time(df,profile)
    sig_mask=_fiberlyse_v231_match_value(df[row_type_col],profile.get('signal_row_value',''));ref_mask=_fiberlyse_v231_match_value(df[row_type_col],profile.get('reference_row_value',''))
    if int(np.sum(sig_mask))<2 or int(np.sum(ref_mask))<2:raise ValueError('The selected Study/Reference row values did not find enough rows. Check the measurement-type column and values in CSV setup.')
    results: Dict[str,ChannelResult]={}
    for i,col in enumerate(value_cols,start=1):
        vals=pd.to_numeric(df[col],errors='coerce').to_numpy(dtype=float);label=_fiberlyse_v231_safe_channel_name(col,i)
        results[label]=_fiberlyse_v231_make_result(csv_path,label,elapsed[ref_mask],elapsed[sig_mask],vals[ref_mask],vals[sig_mask],settings,profile,{'import_signal_column':col,'import_reference_column':col,'import_row_type_column':row_type_col,'import_signal_row_value':str(profile.get('signal_row_value','')),'import_reference_row_value':str(profile.get('reference_row_value','')),'import_csv_separator':str(profile.get('csv_separator',',')),'import_decimal_mark':str(profile.get('decimal_mark','.'))})
    return results


def _fiberlyse_v231_profile_summary(profile: Dict[str, Any]) -> str:
    layout=str(profile.get('layout',''));sep=str(profile.get('csv_separator',','));dec=str(profile.get('decimal_mark','.'))
    fmt='' if layout==_FIBERLYSE_V231_LAYOUT_BUILTIN else f" • CSV: {_fiberlyse_v231_separator_choice(sep)}, {_fiberlyse_v231_decimal_choice(dec)}"
    if layout==_FIBERLYSE_V231_LAYOUT_BUILTIN:return 'Built-in Neurophotometrics/Bonsai mapping: SystemTimestamp = time, LedState = row type, G0/G1/... = fluorescence.'
    if layout==_FIBERLYSE_V231_LAYOUT_SEPARATE:
        pairs=list(profile.get('channel_pairs') or []);return f"Separate columns • Time: {profile.get('time_column','?')} • {len(pairs)} Signal/Reference pair(s) • first {float(profile.get('exclude_initial_seconds',0)):g} s excluded{fmt}"
    cols=list(profile.get('measurement_columns') or []);return f"Alternating rows • Time: {profile.get('time_column','?')} • Row type: {profile.get('row_type_column','?')} • {len(cols)} fluorescence column(s) • first {float(profile.get('exclude_initial_seconds',0)):g} s excluded{fmt}"


def _fiberlyse_v231_configure_custom(app, layout: str, initial_example: Optional[str]=None) -> Optional[Dict[str, Any]]:
    example_path=initial_example or _fiberlyse_v231_choose_example_csv(app.root)
    if not example_path:return None
    try:detected=_fiberlyse_v231_detect_csv_format(example_path)
    except Exception as e:messagebox.showerror('Could not inspect CSV',f'FiberLyse could not inspect that example file:\n\n{e}',parent=app.root);return None
    csvfmt=_fiberlyse_v231_confirm_csv_format(app,example_path,detected)
    if not csvfmt:return None
    try:columns,preview_df=_fiberlyse_v231_read_example(example_path,csvfmt)
    except Exception as e:messagebox.showerror('Could not read CSV',f'FiberLyse could not read that example file with the selected CSV-format settings:\n\n{e}',parent=app.root);return None
    top=tk.Toplevel(app.root);top.title('Define your CSV columns');top.transient(app.root);top.grab_set();top.resizable(True,True)
    try:
        sw=max(840,int(top.winfo_screenwidth()));sh=max(650,int(top.winfo_screenheight()));top.geometry(f'{min(1040,sw-80)}x{min(780,sh-100)}')
    except Exception:top.geometry('980x720')
    shell=ttk.Frame(top);shell.pack(fill=tk.BOTH,expand=True)
    # Fixed top action bar: always visible regardless of form height.
    topbar=ttk.Frame(shell,padding=(10,8));topbar.pack(side=tk.TOP,fill=tk.X)
    ttk.Label(topbar,text='Define what each column means',font=('TkDefaultFont',14,'bold')).pack(side=tk.LEFT)
    content_holder=ttk.Frame(shell);content_holder.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
    canvas=tk.Canvas(content_holder,borderwidth=0,highlightthickness=0);vs=ttk.Scrollbar(content_holder,orient='vertical',command=canvas.yview);canvas.configure(yscrollcommand=vs.set);canvas.pack(side=tk.LEFT,fill=tk.BOTH,expand=True);vs.pack(side=tk.RIGHT,fill=tk.Y)
    outer=ttk.Frame(canvas,padding=12);win=canvas.create_window((0,0),window=outer,anchor='nw')
    outer.bind('<Configure>',lambda _e:canvas.configure(scrollregion=canvas.bbox('all')));canvas.bind('<Configure>',lambda e:canvas.itemconfigure(win,width=max(1,e.width)))
    ttk.Label(outer,text='Choose the column that best matches each plain-language question below. Nothing is changed in your original CSV.',wraplength=900,justify='left').pack(anchor='w',pady=(0,6))
    ttk.Label(outer,text=f"Example file: {os.path.basename(example_path)} • {_fiberlyse_v231_separator_choice(csvfmt['csv_separator'])} • {_fiberlyse_v231_decimal_choice(csvfmt['decimal_mark'])}").pack(anchor='w')
    preview_box=ttk.LabelFrame(outer,text='Small preview of your CSV');preview_box.pack(fill=tk.X,pady=(8,10));txt=tk.Text(preview_box,height=7,wrap='none');txt.pack(fill=tk.X,expand=True,padx=6,pady=6);txt.insert('1.0',_fiberlyse_v231_preview_text(preview_df));txt.configure(state='disabled')
    form=ttk.Frame(outer);form.pack(fill=tk.BOTH,expand=True);form.columnconfigure(1,weight=1)
    time_var=tk.StringVar(value=_fiberlyse_v231_suggest_column(columns,['timestamp','time','clock','seconds','sec'],0));units_var=tk.StringVar(value='Seconds');start_rule_var=tk.StringVar(value='first_row');start_col_var=tk.StringVar(value=_fiberlyse_v231_suggest_column(columns,['start','ledstate','state','event','trigger'],0));start_val_var=tk.StringVar(value='');exclude_var=tk.StringVar(value='0');profile_name_var=tk.StringVar(value='My CSV setup')
    row=0
    ttk.Label(form,text='1. Which column tells FiberLyse when each measurement was recorded?',wraplength=430,justify='left').grid(row=row,column=0,sticky='w',pady=5);ttk.Combobox(form,textvariable=time_var,values=columns,state='readonly',width=34).grid(row=row,column=1,sticky='ew',padx=(12,0),pady=5);row+=1
    ttk.Label(form,text='What units are those time values in?').grid(row=row,column=0,sticky='w',pady=5);ttk.Combobox(form,textvariable=units_var,values=_FIBERLYSE_V231_TIME_UNITS,state='readonly',width=20).grid(row=row,column=1,sticky='w',padx=(12,0),pady=5);row+=1
    dynamic=ttk.LabelFrame(form,text='2. Where are the fluorescence values?');dynamic.grid(row=row,column=0,columnspan=2,sticky='nsew',pady=(10,6));dynamic.columnconfigure(1,weight=1);row+=1
    pair_vars: List[Tuple[tk.StringVar,tk.StringVar,tk.StringVar]]=[];measurement_list=None;rowtype_var=tk.StringVar(value=_fiberlyse_v231_suggest_column(columns,['ledstate','state','wavelength','type','kind','channel'],0));signal_code_var=tk.StringVar(value='2');reference_code_var=tk.StringVar(value='1');count_var=tk.IntVar(value=1);pair_holder=ttk.Frame(dynamic)
    def rebuild_pairs(*_):
        for c in pair_holder.winfo_children():c.destroy()
        old=[(a.get(),b.get(),c.get()) for a,b,c in pair_vars];pair_vars.clear();n=max(1,min(8,int(count_var.get() or 1)))
        for i in range(n):
            prev=old[i] if i<len(old) else ('','',f'Channel {i+1}');sig_guess=_fiberlyse_v231_suggest_column(columns,['470','465','signal','sensor','green','activity'],1,exclude={time_var.get()});ref_guess=_fiberlyse_v231_suggest_column(columns,['405','415','reference','control','isosbestic','iso'],2,exclude={time_var.get(),sig_guess});sv=tk.StringVar(value=prev[0] or sig_guess);rv=tk.StringVar(value=prev[1] or ref_guess);lv=tk.StringVar(value=prev[2] or f'Channel {i+1}');pair_vars.append((sv,rv,lv))
            ttk.Label(pair_holder,text=f'Recording channel {i+1}').grid(row=i,column=0,sticky='w',pady=3);ttk.Label(pair_holder,text='Signal you want to study:').grid(row=i,column=1,sticky='e',padx=(8,4));ttk.Combobox(pair_holder,textvariable=sv,values=columns,state='readonly',width=24).grid(row=i,column=2,sticky='ew');ttk.Label(pair_holder,text='Comparison/reference used for correction:').grid(row=i,column=3,sticky='e',padx=(8,4));ttk.Combobox(pair_holder,textvariable=rv,values=columns,state='readonly',width=24).grid(row=i,column=4,sticky='ew');ttk.Entry(pair_holder,textvariable=lv,width=14).grid(row=i,column=5,sticky='ew',padx=(8,0))
        pair_holder.columnconfigure(2,weight=1);pair_holder.columnconfigure(4,weight=1)
    rowtype_combo=None;signal_combo=None;reference_combo=None
    if layout==_FIBERLYSE_V231_LAYOUT_SEPARATE:
        ttk.Label(dynamic,text='Use this when the signal you want to study and its Reference/control already have their own columns in every row.',wraplength=820,justify='left').grid(row=0,column=0,columnspan=2,sticky='w',padx=6,pady=(6,4));countrow=ttk.Frame(dynamic);countrow.grid(row=1,column=0,columnspan=2,sticky='w',padx=6,pady=4);ttk.Label(countrow,text='How many Signal + Reference pairs are in each file?').pack(side=tk.LEFT);sp=ttk.Spinbox(countrow,from_=1,to=8,textvariable=count_var,width=5,command=rebuild_pairs);sp.pack(side=tk.LEFT,padx=8);pair_holder.grid(row=2,column=0,columnspan=2,sticky='ew',padx=6,pady=4);rebuild_pairs();sp.bind('<Return>',rebuild_pairs);sp.bind('<FocusOut>',rebuild_pairs)
    else:
        ttk.Label(dynamic,text='Use this when Study signal and Reference use the same fluorescence column(s), but different rows. A column such as LedState, LED, Wavelength or Type tells FiberLyse which kind of measurement each row contains.',wraplength=820,justify='left').grid(row=0,column=0,columnspan=2,sticky='w',padx=6,pady=(6,4))
        ttk.Label(dynamic,text='Which column says what kind of measurement each row contains?').grid(row=1,column=0,sticky='w',padx=6,pady=4);rowtype_combo=ttk.Combobox(dynamic,textvariable=rowtype_var,values=columns,state='readonly',width=30);rowtype_combo.grid(row=1,column=1,sticky='ew',padx=6,pady=4)
        codeframe=ttk.Frame(dynamic);codeframe.grid(row=2,column=0,columnspan=2,sticky='ew',padx=6,pady=4);ttk.Label(codeframe,text='Which value means Study signal?').grid(row=0,column=0,sticky='w');signal_combo=ttk.Combobox(codeframe,textvariable=signal_code_var,width=14);signal_combo.grid(row=0,column=1,sticky='w',padx=(5,18));ttk.Label(codeframe,text='Which value means Reference/control?').grid(row=0,column=2,sticky='w');reference_combo=ttk.Combobox(codeframe,textvariable=reference_code_var,width=14);reference_combo.grid(row=0,column=3,sticky='w',padx=5)
        found_values_var=tk.StringVar(value='');ttk.Label(dynamic,textvariable=found_values_var,wraplength=780,justify='left').grid(row=3,column=0,columnspan=2,sticky='w',padx=6,pady=(0,4))
        ttk.Label(dynamic,text='Which column(s) contain the recorded brightness/fluorescence values? Select one or more. For Neurophotometrics these are normally G0, G1, etc.:',wraplength=430,justify='left').grid(row=4,column=0,sticky='nw',padx=6,pady=4);lf=ttk.Frame(dynamic);lf.grid(row=4,column=1,sticky='nsew',padx=6,pady=4);measurement_list=tk.Listbox(lf,selectmode=tk.EXTENDED,height=6,exportselection=False);sb=ttk.Scrollbar(lf,orient='vertical',command=measurement_list.yview);measurement_list.configure(yscrollcommand=sb.set);measurement_list.pack(side=tk.LEFT,fill=tk.BOTH,expand=True);sb.pack(side=tk.RIGHT,fill=tk.Y)
        for c in columns:measurement_list.insert(tk.END,c)
        for i,c in enumerate(columns):
            lc=c.lower()
            if re.match(r'^g\d+$',c,flags=re.I) or any(k in lc for k in ['fluor','signal','470','465']):measurement_list.selection_set(i)
        def refresh_codes(*_):
            vals=_fiberlyse_v231_unique_column_values(example_path,rowtype_var.get(),csvfmt);signal_combo.configure(values=vals);reference_combo.configure(values=vals);found_values_var.set('Values FiberLyse found in this column: '+(', '.join(vals[:20]) if vals else 'none found'))
            lc=rowtype_var.get().lower()
            if 'ledstate' in lc:
                if '2' in vals:signal_code_var.set('2')
                if '1' in vals:reference_code_var.set('1')
        rowtype_combo.bind('<<ComboboxSelected>>',refresh_codes);refresh_codes()
    startbox=ttk.LabelFrame(form,text='3. Where should the recording begin?');startbox.grid(row=row,column=0,columnspan=2,sticky='ew',pady=(10,6));startbox.columnconfigure(1,weight=1);row+=1
    ttk.Radiobutton(startbox,text='The recording begins on the first row of the CSV',variable=start_rule_var,value='first_row').grid(row=0,column=0,columnspan=2,sticky='w',padx=6,pady=3);ttk.Radiobutton(startbox,text='The recording begins when a column reaches a particular value',variable=start_rule_var,value='column_value').grid(row=1,column=0,columnspan=2,sticky='w',padx=6,pady=3);ttk.Label(startbox,text='Column:').grid(row=2,column=0,sticky='e',padx=(20,4),pady=3);start_col_combo=ttk.Combobox(startbox,textvariable=start_col_var,values=columns,state='readonly',width=28);start_col_combo.grid(row=2,column=1,sticky='w',pady=3);ttk.Label(startbox,text='Value that means “recording starts”:').grid(row=3,column=0,sticky='e',padx=(20,4),pady=3);start_val_combo=ttk.Combobox(startbox,textvariable=start_val_var,width=18);start_val_combo.grid(row=3,column=1,sticky='w',pady=3)
    def refresh_start_values(*_):
        vals=_fiberlyse_v231_unique_column_values(example_path,start_col_var.get(),csvfmt);start_val_combo.configure(values=vals)
        if 'ledstate' in start_col_var.get().lower() and '7' in vals:start_val_var.set('7')
    start_col_combo.bind('<<ComboboxSelected>>',refresh_start_values);refresh_start_values()
    ttk.Label(form,text='4. How much of the beginning should be ignored before any FiberLyse analysis?',wraplength=430,justify='left').grid(row=row,column=0,sticky='w',pady=5);exfr=ttk.Frame(form);exfr.grid(row=row,column=1,sticky='w',padx=(12,0),pady=5);ttk.Entry(exfr,textvariable=exclude_var,width=10).pack(side=tk.LEFT);ttk.Label(exfr,text='seconds').pack(side=tk.LEFT,padx=4);row+=1
    ttk.Label(form,text='Name this CSV setup:').grid(row=row,column=0,sticky='w',pady=5);ttk.Entry(form,textvariable=profile_name_var,width=32).grid(row=row,column=1,sticky='w',padx=(12,0),pady=5);row+=1
    remember_var=tk.BooleanVar(value=True);ttk.Checkbutton(form,text='Remember this setup next time FiberLyse starts',variable=remember_var).grid(row=row,column=0,columnspan=2,sticky='w',pady=(8,4));row+=1
    result: Dict[str,Any]={}
    def save_and_close():
        try:exclude=max(0.0,float(exclude_var.get()))
        except Exception:messagebox.showerror('Check ignored time','Please enter the ignored beginning in seconds, for example 0 or 600.',parent=top);return
        p={'profile_name':profile_name_var.get().strip() or 'My CSV setup','layout':layout,'time_column':time_var.get(),'time_units':units_var.get(),'start_rule':start_rule_var.get(),'start_column':start_col_var.get(),'start_value':start_val_var.get(),'exclude_initial_seconds':exclude,'remember':bool(remember_var.get()),'example_file':example_path,**csvfmt}
        if layout==_FIBERLYSE_V231_LAYOUT_SEPARATE:
            rebuild_pairs();pairs=[]
            for i,(sv,rv,lv) in enumerate(pair_vars,start=1):
                if not sv.get() or not rv.get():messagebox.showerror('Missing column',f'Choose both Signal and Reference/control columns for recording channel {i}.',parent=top);return
                if sv.get()==rv.get():messagebox.showerror('Check columns',f'Recording channel {i} uses the same column for Signal and Reference/control.',parent=top);return
                pairs.append({'signal_column':sv.get(),'reference_column':rv.get(),'label':lv.get().strip() or f'Channel {i}'})
            p['channel_pairs']=pairs
        else:
            p['row_type_column']=rowtype_var.get();p['signal_row_value']=signal_code_var.get().strip();p['reference_row_value']=reference_code_var.get().strip();sel=list(measurement_list.curselection()) if measurement_list is not None else [];cols=[columns[int(i)] for i in sel]
            if not cols:messagebox.showerror('Missing fluorescence column','Select at least one column containing fluorescence values.',parent=top);return
            if not p['row_type_column'] or not p['signal_row_value'] or not p['reference_row_value']:messagebox.showerror('Missing row information','Choose the measurement-type column and the values meaning Study signal and Reference/control.',parent=top);return
            p['measurement_columns']=cols
        try:
            fake=_fiberlyse_parse_analysis_settings(app) or {'artifact_enabled':False,'artifact_factor':11.9,'artifact_pad':1,'require_shared':True,'align_mode':DEFAULT_ALIGN_MODE,'acq_fps':40.0,'smooth_win':20,'mode':NORM_DFF,'interval_start':exclude,'interval_end':6500.0,'use_interp':True,'fit_windows':[(exclude,6500.0)]};fake=dict(fake);fake['fit_windows']=[(exclude,6500.0)];fake['interval_start']=max(exclude,float(fake.get('interval_start',exclude)));test=_fiberlyse_v231_analyze_file(example_path,p,fake)
            if not test:raise ValueError('No usable channels were produced.')
        except Exception as e:messagebox.showerror('This setup does not match the example CSV',f'FiberLyse tried the mapping before saving it and found a problem:\n\n{e}',parent=top);return
        result.update(p);top.destroy()
    ttk.Button(topbar,text='Cancel',command=top.destroy).pack(side=tk.RIGHT);ttk.Button(topbar,text='Save & use this setup',command=save_and_close).pack(side=tk.RIGHT,padx=(0,8))
    bottom=ttk.Frame(outer);bottom.pack(fill=tk.X,pady=(12,4));ttk.Button(bottom,text='Cancel',command=top.destroy).pack(side=tk.RIGHT);ttk.Button(bottom,text='Save & use this setup',command=save_and_close).pack(side=tk.RIGHT,padx=(0,8))
    top.bind('<Control-Return>',lambda _e:save_and_close());top.bind('<Escape>',lambda _e:top.destroy())
    try:_fiberlyse_v23_apply_classic_widget_theme(top,_fiberlyse_v23_theme_tokens(app=app))
    except Exception:pass
    app.root.wait_window(top);return dict(result) if result else None


# Add the CSV-format choices to graph/data export metadata.
_FIBERLYSE_V231_DIALECT_META_BEFORE=ChannelTabsTk._meta_df
def _fiberlyse_v231_dialect_meta_df(self) -> pd.DataFrame:
    base=_FIBERLYSE_V231_DIALECT_META_BEFORE(self);r=self.res
    rows=[('import_csv_separator',getattr(r,'import_csv_separator',',')),('import_decimal_mark',getattr(r,'import_decimal_mark','.'))]
    return pd.concat([base,pd.DataFrame(rows,columns=['key','value'])],ignore_index=True)
ChannelTabsTk._meta_df=_fiberlyse_v231_dialect_meta_df

# ---- End V23.1 CSV format/dialect compatibility patch ----------------------

# ---- End FiberLyse V23.1 universal CSV setup extension --------------------



# ---- FiberLyse V24 universal frequency-analysis extension -----------------
# V24 keeps the V23.7 photometry preprocessing/fitting pipeline intact and
# replaces the final frequency-analysis layer with a timestamp-aware system.
# All filtering uses the CURRENT data contained in the loaded recording.
# Optional original-acquisition/downsampling fields are metadata only and are
# never used to pretend that discarded high-frequency samples still exist.
FIBERLYSE_VERSION = 'V24.1'

_FIBERLYSE_V24_FREQ_MODE_STANDARD = 'FiberLyse standard bands'
_FIBERLYSE_V24_FREQ_MODE_AUTO = 'Automatic bands for this recording'
_FIBERLYSE_V24_FREQ_MODE_CUSTOM = 'Custom bands'
_FIBERLYSE_V24_FREQ_MODES = [
    _FIBERLYSE_V24_FREQ_MODE_STANDARD,
    _FIBERLYSE_V24_FREQ_MODE_AUTO,
    _FIBERLYSE_V24_FREQ_MODE_CUSTOM,
]
_FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S = 60.0
_FIBERLYSE_V24_DEFAULT_PSD_OVERLAP_PCT = 50.0
_FIBERLYSE_V24_MAX_AUTO_BANDS = 8
_FIBERLYSE_V24_GAP_MULTIPLE = 1.75


def _fiberlyse_v24_float_or_nan(v: Any) -> float:
    try:
        f = float(v)
        return f if np.isfinite(f) else np.nan
    except Exception:
        return np.nan


def _fiberlyse_v24_default_custom_bands_text() -> str:
    bands = sorted([(float(a), float(b)) for a, b in FREQ_BANDS], key=lambda z: z[0])
    return '\n'.join(f'{a:g} - {b:g}' for a, b in bands)


def _fiberlyse_v24_init_frequency_settings(res: ChannelResult) -> None:
    if not hasattr(res, 'freq_band_mode'):
        res.freq_band_mode = _FIBERLYSE_V24_FREQ_MODE_STANDARD
    if not hasattr(res, 'freq_custom_bands'):
        res.freq_custom_bands = sorted([(float(a), float(b)) for a, b in FREQ_BANDS], key=lambda z: z[0])
    if not hasattr(res, 'freq_psd_window_s'):
        res.freq_psd_window_s = float(_FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S)
    if not hasattr(res, 'freq_psd_overlap_pct'):
        res.freq_psd_overlap_pct = float(_FIBERLYSE_V24_DEFAULT_PSD_OVERLAP_PCT)
    if not hasattr(res, 'freq_filter_order'):
        res.freq_filter_order = int(DEFAULT_BUTTER_ORDER)
    if not hasattr(res, 'freq_original_fs_hz'):
        res.freq_original_fs_hz = np.nan
    if not hasattr(res, 'freq_downsample_factor'):
        res.freq_downsample_factor = np.nan
    if not hasattr(res, 'freq_antialias_status'):
        res.freq_antialias_status = 'Unknown'


def _fiberlyse_v24_parse_custom_bands(text: str) -> List[Tuple[float, float]]:
    out: List[Tuple[float, float]] = []
    for raw in str(text or '').splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        # Accept: 0.1-0.3, 0.1 to 0.3, 0,1 - 0,3, or 0.1;0.3.
        nums = re.findall(r'[-+]?\d+(?:[\.,]\d+)?(?:[eE][-+]?\d+)?', line)
        if len(nums) < 2:
            raise ValueError(f"Could not read a low and high frequency from: '{line}'")
        try:
            low = float(nums[0].replace(',', '.'))
            high = float(nums[1].replace(',', '.'))
        except Exception:
            raise ValueError(f"Could not read numeric frequencies from: '{line}'")
        if not np.isfinite(low) or not np.isfinite(high) or low <= 0 or high <= low:
            raise ValueError(f"Frequency range must be positive and increasing: '{line}'")
        out.append((float(low), float(high)))
    if not out:
        raise ValueError('Enter at least one frequency range, for example: 0.15 - 0.30')
    # De-duplicate while retaining the user's order.
    dedup: List[Tuple[float, float]] = []
    seen = set()
    for b in out:
        key = (round(b[0], 12), round(b[1], 12))
        if key not in seen:
            seen.add(key); dedup.append(b)
    return dedup


def _fiberlyse_v24_time_segments(t: np.ndarray, y: np.ndarray, gap_multiple: float=_FIBERLYSE_V24_GAP_MULTIPLE) -> List[Tuple[int, int]]:
    """Finite contiguous segments, additionally split at large timestamp gaps."""
    t = np.asarray(t, dtype=float).reshape(-1); y = np.asarray(y, dtype=float).reshape(-1)
    n = min(t.size, y.size)
    if n < 2:
        return []
    t = t[:n]; y = y[:n]
    finite = np.isfinite(t) & np.isfinite(y)
    dt_all = np.diff(t[np.isfinite(t)])
    dt_all = dt_all[np.isfinite(dt_all) & (dt_all > 0)]
    med_dt = float(np.median(dt_all)) if dt_all.size else np.nan
    runs: List[Tuple[int, int]] = []
    start = None; prev = None
    for i in range(n):
        if not finite[i]:
            if start is not None and prev is not None and prev > start:
                runs.append((start, prev))
            start = prev = None
            continue
        if start is None:
            start = prev = i
            continue
        split = False
        try:
            d = float(t[i] - t[prev])
            if not np.isfinite(d) or d <= 0:
                split = True
            elif np.isfinite(med_dt) and med_dt > 0 and d > float(gap_multiple) * med_dt:
                split = True
        except Exception:
            split = True
        if split:
            if prev is not None and prev > start:
                runs.append((start, prev))
            start = i
        prev = i
    if start is not None and prev is not None and prev > start:
        runs.append((start, prev))
    return runs


def _fiberlyse_v24_frequency_diagnostics(t: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    t = np.asarray(t, dtype=float).reshape(-1); y = np.asarray(y, dtype=float).reshape(-1)
    n = min(t.size, y.size); t = t[:n]; y = y[:n]
    tf = t[np.isfinite(t)]
    dt = np.diff(tf) if tf.size >= 2 else np.array([], dtype=float)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    med_dt = float(np.median(dt)) if dt.size else np.nan
    fs = float(1.0 / med_dt) if np.isfinite(med_dt) and med_dt > 0 else np.nan
    mad_dt = float(np.median(np.abs(dt - med_dt))) if dt.size and np.isfinite(med_dt) else np.nan
    robust_jitter_pct = float(100.0 * 1.4826 * mad_dt / med_dt) if np.isfinite(mad_dt) and np.isfinite(med_dt) and med_dt > 0 else np.nan
    max_gap_multiple = float(np.max(dt) / med_dt) if dt.size and np.isfinite(med_dt) and med_dt > 0 else np.nan
    duration = float(np.nanmax(tf) - np.nanmin(tf)) if tf.size >= 2 else np.nan
    segs = _fiberlyse_v24_time_segments(t, y)
    seg_durations = []
    for i0, i1 in segs:
        try:
            d = float(t[i1] - t[i0])
            if np.isfinite(d) and d > 0:
                seg_durations.append(d)
        except Exception:
            pass
    longest = max(seg_durations) if seg_durations else np.nan
    return {
        'fs_hz': fs,
        'nyquist_hz': fs / 2.0 if np.isfinite(fs) else np.nan,
        'median_dt_s': med_dt,
        'robust_jitter_pct': robust_jitter_pct,
        'max_gap_multiple': max_gap_multiple,
        'duration_s': duration,
        'longest_continuous_s': longest,
        'n_segments': float(len(segs)),
    }


def _fiberlyse_v24_automatic_bands(t: np.ndarray, y: np.ndarray, fs: float, psd_window_s: float) -> List[Tuple[float, float]]:
    diag = _fiberlyse_v24_frequency_diagnostics(t, y)
    nyq = float(fs) / 2.0 if np.isfinite(fs) and fs > 0 else np.nan
    if not np.isfinite(nyq) or nyq <= 0:
        return []
    upper = 0.95 * nyq
    longest = float(diag.get('longest_continuous_s', np.nan))
    low_candidates = [0.001]
    if np.isfinite(psd_window_s) and psd_window_s > 0:
        low_candidates.append(1.0 / float(psd_window_s))
    if np.isfinite(longest) and longest > 0:
        low_candidates.append(3.0 / longest)
    low = max(low_candidates)
    if not np.isfinite(low) or low <= 0 or low >= upper:
        return []
    ratio = upper / low
    try:
        n_bands = int(np.ceil(np.log2(ratio)))
    except Exception:
        n_bands = 4
    n_bands = max(1, min(_FIBERLYSE_V24_MAX_AUTO_BANDS, n_bands))
    edges = np.geomspace(low, upper, n_bands + 1)
    bands: List[Tuple[float, float]] = []
    for a, b in zip(edges[:-1], edges[1:]):
        if np.isfinite(a) and np.isfinite(b) and b > a > 0:
            bands.append((float(a), float(b)))
    return bands


def _fiberlyse_v24_get_frequency_bands(res: ChannelResult, t: np.ndarray, y: np.ndarray, fs: float) -> List[Tuple[float, float]]:
    _fiberlyse_v24_init_frequency_settings(res)
    mode = str(getattr(res, 'freq_band_mode', _FIBERLYSE_V24_FREQ_MODE_STANDARD))
    if mode == _FIBERLYSE_V24_FREQ_MODE_AUTO:
        return _fiberlyse_v24_automatic_bands(t, y, fs, float(getattr(res, 'freq_psd_window_s', _FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S)))
    if mode == _FIBERLYSE_V24_FREQ_MODE_CUSTOM:
        return list(getattr(res, 'freq_custom_bands', []) or [])
    return sorted([(float(a), float(b)) for a, b in FREQ_BANDS], key=lambda z: z[0])


def _fiberlyse_v24_band_status(t: np.ndarray, y: np.ndarray, low_hz: float, high_hz: float, fs: float) -> Tuple[bool, str]:
    if not np.isfinite(fs) or fs <= 0:
        return False, 'Sampling rate could not be determined from timestamps.'
    nyq = fs / 2.0
    if not np.isfinite(low_hz) or not np.isfinite(high_hz) or low_hz <= 0 or high_hz <= low_hz:
        return False, 'Invalid frequency range.'
    # Strict Nyquist rule: do not silently clip an impossible requested band.
    if high_hz >= nyq:
        return False, f'Unavailable: upper edge {high_hz:g} Hz reaches/exceeds Nyquist ({nyq:.4g} Hz).'
    required_s = 3.0 / low_hz
    segs = _fiberlyse_v24_time_segments(t, y)
    usable = False
    best = 0.0
    for i0, i1 in segs:
        try:
            dur = float(t[i1] - t[i0]); best = max(best, dur)
            if dur >= required_s:
                usable = True; break
        except Exception:
            pass
    if not usable:
        return False, f'Unavailable: needs ≥3 cycles at {low_hz:g} Hz (~{required_s:.3g} s continuous); longest usable segment is {best:.3g} s.'
    return True, 'Available'


def _fiberlyse_v24_filter_band(res: ChannelResult, t: np.ndarray, y: np.ndarray, low_hz: float, high_hz: float, fs: float, order: int) -> np.ndarray:
    y = np.asarray(y, dtype=float); t = np.asarray(t, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)
    ok, _reason = _fiberlyse_v24_band_status(t, y, low_hz, high_hz, fs)
    if not ok:
        return out
    nyq = fs / 2.0
    if _HAVE_SCIPY_SIGNAL and butter is not None and sosfiltfilt is not None:
        try:
            sos = butter(int(order), [float(low_hz) / nyq, float(high_hz) / nyq], btype='band', output='sos')
        except Exception:
            return out
        for i0, i1 in _fiberlyse_v24_time_segments(t, y):
            seg = np.asarray(y[i0:i1 + 1], dtype=float)
            min_samples = max(8, int(np.ceil((3.0 * float(fs)) / float(low_hz))) + 1)
            if seg.size < min_samples:
                continue
            try:
                out[i0:i1 + 1] = sosfiltfilt(sos, seg)
            except Exception:
                try:out[i0:i1 + 1] = sosfiltfilt(sos, seg, padlen=0)
                except Exception:pass
    else:
        # Segment-wise FFT fallback; unlike the old fallback this can preserve
        # NaN/time-gap holes by processing each finite continuous segment alone.
        for i0, i1 in _fiberlyse_v24_time_segments(t, y):
            seg = np.asarray(y[i0:i1 + 1], dtype=float)
            min_samples = max(8, int(np.ceil((3.0 * float(fs)) / float(low_hz))) + 1)
            if seg.size < min_samples:
                continue
            try:out[i0:i1 + 1] = bandpass_fft_no_interp(seg, low_hz, high_hz, fs)
            except Exception:pass
    return out


def _fiberlyse_v24_compute_psd(t: np.ndarray, y: np.ndarray, fs: float, window_s: float, overlap_pct: float) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    t = np.asarray(t, dtype=float); y = np.asarray(y, dtype=float)
    info: Dict[str, Any] = {'method': 'unavailable', 'segment_duration_s': np.nan, 'window_s_requested': window_s, 'window_s_actual': np.nan, 'overlap_pct': overlap_pct}
    if not np.isfinite(fs) or fs <= 0:
        return np.array([], dtype=float), np.array([], dtype=float), info
    segs = _fiberlyse_v24_time_segments(t, y)
    if not segs:
        return np.array([], dtype=float), np.array([], dtype=float), info
    i0, i1 = max(segs, key=lambda ab: ab[1] - ab[0])
    x = np.asarray(y[i0:i1 + 1], dtype=float)
    if x.size < 8:
        return np.array([], dtype=float), np.array([], dtype=float), info
    try:info['segment_duration_s'] = float(t[i1] - t[i0])
    except Exception:pass
    requested_n = max(8, int(round(max(float(window_s), 1.0 / fs) * fs)))
    nperseg = min(int(x.size), requested_n)
    if nperseg < 8:
        return np.array([], dtype=float), np.array([], dtype=float), info
    overlap_pct = min(max(float(overlap_pct), 0.0), 95.0)
    noverlap = min(nperseg - 1, int(round(nperseg * overlap_pct / 100.0)))
    x = x - float(np.mean(x))
    welch_fn = getattr(_scipy_signal, 'welch', None) if _HAVE_SCIPY_SIGNAL else None
    if callable(welch_fn):
        try:
            f, pxx = welch_fn(x, fs=float(fs), window='hann', nperseg=nperseg, noverlap=noverlap, detrend='constant', scaling='density')
            info['method'] = 'Welch PSD (longest continuous segment)'
            info['window_s_actual'] = float(nperseg / fs)
            return np.asarray(f, dtype=float), np.asarray(pxx, dtype=float), info
        except Exception:
            pass
    # Portable fallback: one Hann-windowed periodogram from the first nperseg samples.
    xx = np.asarray(x[:nperseg], dtype=float)
    win = np.hanning(xx.size)
    xxw = xx * win
    X = np.fft.rfft(xxw)
    f = np.fft.rfftfreq(xxw.size, d=1.0 / float(fs))
    denom = float(fs) * float(np.sum(win ** 2))
    pxx = (np.abs(X) ** 2) / denom if denom > 0 else np.full_like(f, np.nan, dtype=float)
    if pxx.size > 2:
        pxx[1:-1] *= 2.0
    info['method'] = 'Hann periodogram fallback (longest continuous segment)'
    info['window_s_actual'] = float(nperseg / fs)
    return np.asarray(f, dtype=float), np.asarray(pxx, dtype=float), info


def _fiberlyse_v24_sampling_quality_text(diag: Dict[str, float]) -> str:
    jitter = float(diag.get('robust_jitter_pct', np.nan)); gap = float(diag.get('max_gap_multiple', np.nan))
    if not np.isfinite(jitter):
        return 'Sampling regularity unavailable'
    if jitter <= 1.0 and (not np.isfinite(gap) or gap < 1.5):
        quality = 'very regular'
    elif jitter <= 5.0 and (not np.isfinite(gap) or gap < 2.0):
        quality = 'mostly regular'
    else:
        quality = 'irregular — inspect timestamps/gaps'
    return f'{quality} (robust jitter {jitter:.3g}%; largest interval {gap:.3g}× median)' if np.isfinite(gap) else f'{quality} (robust jitter {jitter:.3g}%)'



def _fiberlyse_v24_get_band_cached(res: ChannelResult, t: np.ndarray, y: np.ndarray, low_hz: float, high_hz: float, fs: float, order: int) -> np.ndarray:
    try:version=int(getattr(res,'_data_version',0))
    except Exception:version=0
    key=(version,round(float(low_hz),12),round(float(high_hz),12),round(float(fs),12),int(order))
    cache=getattr(res,'_freq_v24_cache',None)
    if not isinstance(cache,dict):cache={};setattr(res,'_freq_v24_cache',cache)
    if key not in cache:
        if len(cache)>24:cache.clear()
        cache[key]=_fiberlyse_v24_filter_band(res,t,y,low_hz,high_hz,fs,order)
    return cache[key]


def _fiberlyse_v24_plot_decimate(t: np.ndarray, y: np.ndarray, max_points: int=40000) -> Tuple[np.ndarray,np.ndarray]:
    """Display-only decimation; full-resolution arrays remain used/exported.

    NaN boundaries and large timestamp gaps are explicitly retained so visual
    decimation does not silently draw through known missing regions.
    """
    t=np.asarray(t,dtype=float).reshape(-1);y=np.asarray(y,dtype=float).reshape(-1);n=min(t.size,y.size);t=t[:n];y=y[:n]
    if n<=max_points:return t,y
    stride=max(1,int(np.ceil(n/float(max_points))))
    keep=set(range(0,n,stride));keep.add(n-1)
    bad=np.flatnonzero(~np.isfinite(y)|~np.isfinite(t))
    for i in bad:
        for j in (i-1,i,i+1):
            if 0<=j<n:keep.add(int(j))
    tf=t[np.isfinite(t)];dt=np.diff(tf);dt=dt[np.isfinite(dt)&(dt>0)];med=float(np.median(dt)) if dt.size else np.nan
    if np.isfinite(med) and med>0:
        d=np.diff(t)
        gaps=np.flatnonzero(np.isfinite(d)&(d>_FIBERLYSE_V24_GAP_MULTIPLE*med))
        for i in gaps:
            keep.add(int(i));keep.add(int(i+1))
    idx=np.asarray(sorted(keep),dtype=int)
    return t[idx],y[idx]


def _fiberlyse_v24_draw_frequency(self):
    res = self.res; _fiberlyse_v24_init_frequency_settings(res)
    fig = self.tab_freq.fig; fig.clear()
    t = np.asarray(res.t_exc, dtype=float); dff = np.asarray(res.dFF_nointerp, dtype=float)
    diag = _fiberlyse_v24_frequency_diagnostics(t, dff)
    fs = float(diag.get('fs_hz', np.nan))
    if not np.isfinite(fs) or fs <= 0:
        fs = float(getattr(res, 'eff_fs_hz', np.nan))
    nyq = fs / 2.0 if np.isfinite(fs) and fs > 0 else np.nan
    bands = _fiberlyse_v24_get_frequency_bands(res, t, dff, fs)
    order = max(1, int(getattr(res, 'freq_filter_order', DEFAULT_BUTTER_ORDER)))
    psd_window_s = max(0.001, float(getattr(res, 'freq_psd_window_s', _FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S)))
    overlap = min(max(float(getattr(res, 'freq_psd_overlap_pct', _FIBERLYSE_V24_DEFAULT_PSD_OVERLAP_PCT)), 0.0), 95.0)
    f_psd, p_psd, psd_info = _fiberlyse_v24_compute_psd(t, dff, fs, psd_window_s, overlap)

    n_bands = max(1, len(bands)); n_rows_bands = int(np.ceil(n_bands / 2.0))
    gs = fig.add_gridspec(1 + n_rows_bands, 2, height_ratios=[1.35] + [1.0] * n_rows_bands)
    ax_psd = fig.add_subplot(gs[0, :]); ax_psd._fiberlyse_not_time_axis = True
    band_axes = []
    for i in range(n_rows_bands * 2):
        r = 1 + i // 2; c = i % 2
        ax = fig.add_subplot(gs[r, c]); band_axes.append(ax)
    for ax in band_axes[len(bands):]: ax.axis('off')

    if f_psd.size and p_psd.size:
        pos = (f_psd > 0) & np.isfinite(f_psd) & np.isfinite(p_psd) & (p_psd > 0)
        if np.any(pos):
            ax_psd.plot(f_psd[pos], p_psd[pos], linewidth=1.0, label='PSD')
            ax_psd.set_xscale('log'); ax_psd.set_yscale('log')
            ax_psd.legend(loc='best')
        else:
            ax_psd.text(0.5, 0.5, 'PSD could not be displayed.', ha='center', va='center', transform=ax_psd.transAxes)
    else:
        ax_psd.text(0.5, 0.5, 'PSD unavailable (not enough continuous data).', ha='center', va='center', transform=ax_psd.transAxes)
    for low_hz, high_hz in bands:
        ok, _ = _fiberlyse_v24_band_status(t, dff, low_hz, high_hz, fs)
        if ok:
            try:ax_psd.axvspan(low_hz, high_hz, alpha=0.08)
            except Exception:pass
    ax_psd.set_title(f'Frequency overview — {psd_info.get("method", "PSD")}')
    ax_psd.set_xlabel('Frequency (Hz)'); ax_psd.set_ylabel('Power spectral density')
    if np.isfinite(nyq) and nyq > 0:
        try:
            positive_f = f_psd[f_psd > 0] if f_psd.size else np.array([], dtype=float)
            xmin = float(np.min(positive_f)) if positive_f.size else max(nyq / 1000.0, 1e-4)
            ax_psd.set_xlim(max(xmin, 1e-6), nyq)
        except Exception:pass

    for ax, (low_hz, high_hz) in zip(band_axes, bands):
        label = f'{low_hz:.4g}–{high_hz:.4g} Hz'
        ok, reason = _fiberlyse_v24_band_status(t, dff, low_hz, high_hz, fs)
        if ok:
            y_band = _fiberlyse_v24_get_band_cached(res, t, dff, low_hz, high_hz, fs, order)
            if np.any(np.isfinite(y_band)):
                tp, yp = _fiberlyse_v24_plot_decimate(t, y_band, max_points=40000)
                ax.plot(tp, yp, linewidth=1.0, label=label); ax.legend(loc='best')
            else:
                ax.text(0.5, 0.5, 'Filter produced no finite samples.', ha='center', va='center', transform=ax.transAxes)
        else:
            ax.text(0.5, 0.5, reason, ha='center', va='center', wrap=True, transform=ax.transAxes)
        ax.set_title(label + ('' if ok else ' — unavailable'))
        ax.set_xlabel('Time (s)'); ax.set_ylabel('ΔF/F')

    mode = str(getattr(res, 'freq_band_mode', _FIBERLYSE_V24_FREQ_MODE_STANDARD))
    regularity = _fiberlyse_v24_sampling_quality_text(diag)
    orig_fs = _fiberlyse_v24_float_or_nan(getattr(res, 'freq_original_fs_hz', np.nan))
    ds = _fiberlyse_v24_float_or_nan(getattr(res, 'freq_downsample_factor', np.nan))
    current_line = f'Current file Fs={fs:.5g} Hz | Nyquist={nyq:.5g} Hz | {regularity}' if np.isfinite(fs) else 'Current file sampling rate unavailable'
    history = ''
    if np.isfinite(orig_fs) or np.isfinite(ds):
        parts=[]
        if np.isfinite(orig_fs): parts.append(f'original Fs={orig_fs:.5g} Hz')
        if np.isfinite(ds): parts.append(f'downsample ×{ds:.5g}')
        history = '\nRecorded history (metadata only): ' + ', '.join(parts) + '; current file limits are used for analysis.'
    fig.suptitle(f'{res.gcol} — Universal frequency analysis\n{current_line}\nBand mode: {mode}{history}', y=0.997, fontsize=10)
    fig.tight_layout(rect=(0, 0.0, 1, 0.91))
    self.tab_freq.redraw()


def _fiberlyse_v24_export_frequency(self) -> Dict[str, pd.DataFrame]:
    res = self.res; _fiberlyse_v24_init_frequency_settings(res)
    t = np.asarray(res.t_exc, dtype=float); dff = np.asarray(res.dFF_nointerp, dtype=float)
    diag = _fiberlyse_v24_frequency_diagnostics(t, dff)
    fs = float(diag.get('fs_hz', np.nan))
    if not np.isfinite(fs) or fs <= 0: fs = float(getattr(res, 'eff_fs_hz', np.nan))
    nyq = fs / 2.0 if np.isfinite(fs) and fs > 0 else np.nan
    bands = _fiberlyse_v24_get_frequency_bands(res, t, dff, fs)
    order = max(1, int(getattr(res, 'freq_filter_order', DEFAULT_BUTTER_ORDER)))
    out: Dict[str, Any] = {'t_s': t, 'dFF_nointerp': dff}
    status_rows = []
    for low_hz, high_hz in bands:
        ok, reason = _fiberlyse_v24_band_status(t, dff, low_hz, high_hz, fs)
        col = f'band_{low_hz:.6g}_{high_hz:.6g}_Hz'
        out[col] = _fiberlyse_v24_get_band_cached(res, t, dff, low_hz, high_hz, fs, order) if ok else np.full_like(dff, np.nan, dtype=float)
        status_rows.append({'low_hz': low_hz, 'high_hz': high_hz, 'available': bool(ok), 'status': reason})
    psd_window_s = max(0.001, float(getattr(res, 'freq_psd_window_s', _FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S)))
    overlap = min(max(float(getattr(res, 'freq_psd_overlap_pct', _FIBERLYSE_V24_DEFAULT_PSD_OVERLAP_PCT)), 0.0), 95.0)
    f_psd, p_psd, psd_info = _fiberlyse_v24_compute_psd(t, dff, fs, psd_window_s, overlap)
    meta = self._meta_df()
    extra_rows = [
        ('frequency_band_mode', str(getattr(res, 'freq_band_mode', _FIBERLYSE_V24_FREQ_MODE_STANDARD))),
        ('frequency_current_fs_hz', fs), ('frequency_nyquist_hz', nyq),
        ('frequency_median_dt_s', diag.get('median_dt_s', np.nan)),
        ('frequency_sampling_robust_jitter_pct', diag.get('robust_jitter_pct', np.nan)),
        ('frequency_largest_interval_multiple_of_median', diag.get('max_gap_multiple', np.nan)),
        ('frequency_recording_duration_s', diag.get('duration_s', np.nan)),
        ('frequency_longest_continuous_segment_s', diag.get('longest_continuous_s', np.nan)),
        ('frequency_filter_type', 'Butterworth band-pass (segment-wise)' if _HAVE_SCIPY_SIGNAL else 'FFT band-pass fallback (segment-wise)'),
        ('frequency_filter_order', order), ('frequency_psd_method', psd_info.get('method', 'unavailable')),
        ('frequency_psd_window_requested_s', psd_window_s), ('frequency_psd_window_actual_s', psd_info.get('window_s_actual', np.nan)),
        ('frequency_psd_overlap_pct', overlap),
        ('frequency_original_acquisition_fs_hz_metadata_only', getattr(res, 'freq_original_fs_hz', np.nan)),
        ('frequency_downsample_factor_metadata_only', getattr(res, 'freq_downsample_factor', np.nan)),
        ('frequency_antialias_before_downsampling', getattr(res, 'freq_antialias_status', 'Unknown')),
        ('frequency_note', 'All frequency analysis uses the current file timestamps/current Nyquist. Original acquisition metadata never restores frequencies discarded before export.'),
    ]
    meta2 = pd.concat([meta, pd.DataFrame(extra_rows, columns=['key', 'value'])], ignore_index=True)
    return {
        'freq_bands': pd.DataFrame(out),
        'freq_band_status': pd.DataFrame(status_rows),
        'frequency_psd': pd.DataFrame({'frequency_hz': f_psd, 'power_density': p_psd}),
        'meta': meta2,
    }


def _fiberlyse_v24_show_frequency_settings(self):
    res = self.res; _fiberlyse_v24_init_frequency_settings(res)
    parent = self.winfo_toplevel(); t = np.asarray(res.t_exc, dtype=float); y = np.asarray(res.dFF_nointerp, dtype=float)
    diag = _fiberlyse_v24_frequency_diagnostics(t, y)
    fs = float(diag.get('fs_hz', np.nan))
    if not np.isfinite(fs) or fs <= 0: fs = float(getattr(res, 'eff_fs_hz', np.nan))
    nyq = fs / 2.0 if np.isfinite(fs) else np.nan
    top = tk.Toplevel(parent); top.title('Frequency analysis settings'); top.transient(parent); top.grab_set(); top.resizable(True, True)
    try:
        sw=max(800,int(top.winfo_screenwidth())); sh=max(650,int(top.winfo_screenheight())); top.geometry(f'{min(900,sw-80)}x{min(780,sh-100)}')
    except Exception: top.geometry('900x760')
    outer = ttk.Frame(top, padding=12); outer.pack(fill=tk.BOTH, expand=True)
    ttk.Label(outer, text='Frequency analysis for this recording', font=('TkDefaultFont', 14, 'bold')).pack(anchor='w')
    ttk.Label(outer, text='FiberLyse uses the timestamps in the CURRENT file to determine sampling rate and Nyquist. Optional original-acquisition information below is stored only as metadata.', wraplength=830, justify='left').pack(anchor='w', pady=(4,8))
    info = ttk.LabelFrame(outer, text='What FiberLyse measured from this file'); info.pack(fill=tk.X, pady=(0,8))
    regtxt = _fiberlyse_v24_sampling_quality_text(diag)
    info_text = f'Current signal sampling rate: {fs:.6g} Hz\nNyquist frequency: {nyq:.6g} Hz\nMedian sample interval: {diag.get("median_dt_s", np.nan):.6g} s\nSampling regularity: {regtxt}\nRecording duration: {diag.get("duration_s", np.nan):.6g} s\nLongest continuous usable segment: {diag.get("longest_continuous_s", np.nan):.6g} s'
    ttk.Label(info, text=info_text, justify='left').pack(anchor='w', padx=8, pady=6)

    mode_var = tk.StringVar(value=str(getattr(res, 'freq_band_mode', _FIBERLYSE_V24_FREQ_MODE_STANDARD)))
    psd_window_var = tk.StringVar(value=f'{float(getattr(res,"freq_psd_window_s",_FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S)):g}')
    overlap_var = tk.StringVar(value=f'{float(getattr(res,"freq_psd_overlap_pct",_FIBERLYSE_V24_DEFAULT_PSD_OVERLAP_PCT)):g}')
    order_var = tk.StringVar(value=str(int(getattr(res, 'freq_filter_order', DEFAULT_BUTTER_ORDER))))
    orig_var = tk.StringVar(value='' if not np.isfinite(_fiberlyse_v24_float_or_nan(getattr(res,'freq_original_fs_hz',np.nan))) else f'{float(getattr(res,"freq_original_fs_hz")):g}')
    ds_var = tk.StringVar(value='' if not np.isfinite(_fiberlyse_v24_float_or_nan(getattr(res,'freq_downsample_factor',np.nan))) else f'{float(getattr(res,"freq_downsample_factor")):g}')
    aa_var = tk.StringVar(value=str(getattr(res, 'freq_antialias_status', 'Unknown')))
    apply_file_var = tk.BooleanVar(value=True)

    bands_box = ttk.LabelFrame(outer, text='1. Which frequency ranges should be shown?'); bands_box.pack(fill=tk.X, pady=(0,8)); bands_box.columnconfigure(1, weight=1)
    for r, mode in enumerate(_FIBERLYSE_V24_FREQ_MODES):
        ttk.Radiobutton(bands_box, text=mode, variable=mode_var, value=mode).grid(row=r, column=0, sticky='w', padx=8, pady=2)
    ttk.Label(bands_box, text='Standard bands stay identical between experiments; impossible bands are marked unavailable rather than clipped.\nAutomatic bands are log-spaced within this file’s measurable range.\nCustom bands let you define exact ranges.', justify='left', wraplength=500).grid(row=0, column=1, rowspan=3, sticky='w', padx=10, pady=3)
    custom_frame = ttk.LabelFrame(outer, text='Custom ranges (one per line, low - high in Hz)'); custom_frame.pack(fill=tk.BOTH, expand=True, pady=(0,8))
    custom_text = tk.Text(custom_frame, height=7, wrap='none'); custom_text.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
    cur_custom = list(getattr(res, 'freq_custom_bands', []) or [])
    custom_text.insert('1.0', '\n'.join(f'{a:g} - {b:g}' for a,b in cur_custom) if cur_custom else _fiberlyse_v24_default_custom_bands_text())

    psd_box = ttk.LabelFrame(outer, text='2. Spectrum and filter settings'); psd_box.pack(fill=tk.X, pady=(0,8)); psd_box.columnconfigure(1, weight=1)
    ttk.Label(psd_box, text='Spectrum window duration:').grid(row=0,column=0,sticky='w',padx=8,pady=4); ttk.Entry(psd_box,textvariable=psd_window_var,width=12).grid(row=0,column=1,sticky='w',pady=4); ttk.Label(psd_box,text='seconds').grid(row=0,column=2,sticky='w',padx=4)
    ttk.Label(psd_box, text='Spectrum overlap:').grid(row=1,column=0,sticky='w',padx=8,pady=4); ttk.Entry(psd_box,textvariable=overlap_var,width=12).grid(row=1,column=1,sticky='w',pady=4); ttk.Label(psd_box,text='%').grid(row=1,column=2,sticky='w',padx=4)
    ttk.Label(psd_box, text='Butterworth filter order:').grid(row=2,column=0,sticky='w',padx=8,pady=4); ttk.Spinbox(psd_box,from_=1,to=10,textvariable=order_var,width=10).grid(row=2,column=1,sticky='w',pady=4)
    ttk.Label(psd_box, text='Spectrum uses Welch PSD when SciPy is available. Window duration is expressed in seconds so it has the same meaning across hardware.', wraplength=520, justify='left').grid(row=0,column=3,rowspan=3,sticky='w',padx=12,pady=4)

    hist = ttk.LabelFrame(outer, text='3. Optional acquisition history (metadata only)'); hist.pack(fill=tk.X, pady=(0,8)); hist.columnconfigure(1,weight=1)
    ttk.Label(hist,text='Original acquisition sampling rate:').grid(row=0,column=0,sticky='w',padx=8,pady=4); ttk.Entry(hist,textvariable=orig_var,width=14).grid(row=0,column=1,sticky='w',pady=4); ttk.Label(hist,text='Hz').grid(row=0,column=2,sticky='w',padx=4)
    ttk.Label(hist,text='Downsampling factor:').grid(row=1,column=0,sticky='w',padx=8,pady=4); ttk.Entry(hist,textvariable=ds_var,width=14).grid(row=1,column=1,sticky='w',pady=4)
    ttk.Label(hist,text='Low-pass/anti-alias filter before downsampling?').grid(row=2,column=0,sticky='w',padx=8,pady=4); ttk.Combobox(hist,textvariable=aa_var,values=['Yes','No','Unknown'],state='readonly',width=12).grid(row=2,column=1,sticky='w',pady=4)
    ttk.Label(hist,text='These fields never change the current Nyquist limit. They help document how the file was produced and translate old sample-based settings.',wraplength=520,justify='left').grid(row=0,column=3,rowspan=3,sticky='w',padx=12,pady=4)
    ttk.Checkbutton(outer,text='Apply these frequency settings to all analyzed channels from the same source file',variable=apply_file_var).pack(anchor='w',pady=(0,8))

    def save_settings(_event=None):
        try:
            psd_s=float(psd_window_var.get()); ov=float(overlap_var.get()); order=int(float(order_var.get()))
            if not np.isfinite(psd_s) or psd_s <= 0: raise ValueError('Spectrum window duration must be a positive number of seconds.')
            if not np.isfinite(ov) or ov < 0 or ov > 95: raise ValueError('Spectrum overlap must be between 0 and 95%.')
            if order < 1 or order > 10: raise ValueError('Butterworth order must be between 1 and 10.')
            custom = _fiberlyse_v24_parse_custom_bands(custom_text.get('1.0', tk.END))
            orig = np.nan if not orig_var.get().strip() else float(orig_var.get())
            ds = np.nan if not ds_var.get().strip() else float(ds_var.get())
            if np.isfinite(orig) and orig <= 0: raise ValueError('Original acquisition sampling rate must be positive.')
            if np.isfinite(ds) and ds <= 0: raise ValueError('Downsampling factor must be positive.')
            # Help the user by deriving one history field from the other/current Fs.
            if np.isfinite(orig) and not np.isfinite(ds) and np.isfinite(fs) and fs > 0: ds = orig / fs
            elif np.isfinite(ds) and not np.isfinite(orig) and np.isfinite(fs) and fs > 0: orig = fs * ds
        except Exception as e:
            messagebox.showerror('Check frequency settings', str(e), parent=top); return 'break'
        targets=[res]
        if bool(apply_file_var.get()) and self.parent_app is not None and getattr(self.parent_app,'_results',None):
            src=os.path.abspath(str(getattr(res,'source_path','')))
            targets=[r for r in self.parent_app._results.values() if os.path.abspath(str(getattr(r,'source_path',''))) == src] or [res]
        for r in targets:
            _fiberlyse_v24_init_frequency_settings(r)
            r.freq_band_mode=mode_var.get(); r.freq_custom_bands=list(custom); r.freq_psd_window_s=float(psd_s); r.freq_psd_overlap_pct=float(ov); r.freq_filter_order=int(order); r.freq_original_fs_hz=float(orig) if np.isfinite(orig) else np.nan; r.freq_downsample_factor=float(ds) if np.isfinite(ds) else np.nan; r.freq_antialias_status=aa_var.get()
            try:r._freq_cache={}; r._freq_v24_cache={}
            except Exception:pass
        try:self._freq_drawn_version=None; self._draw_frequency(); self._freq_drawn_version=int(getattr(res,'_data_version',0))
        except Exception:pass
        try:top.grab_release()
        except Exception:pass
        top.destroy(); return 'break'
    def cancel(_event=None):
        try:top.grab_release()
        except Exception:pass
        top.destroy(); return 'break'
    actions=ttk.Frame(outer); actions.pack(fill=tk.X, pady=(4,0)); ttk.Button(actions,text='Cancel',command=cancel).pack(side=tk.RIGHT); ttk.Button(actions,text='Apply & redraw',command=save_settings).pack(side=tk.RIGHT,padx=(0,8))
    top.bind('<Control-Return>',save_settings); top.bind('<Escape>',cancel)
    try:_fiberlyse_v23_apply_classic_widget_theme(top,_fiberlyse_v23_theme_tokens(app=getattr(self,'parent_app',None)))
    except Exception:pass
    parent.wait_window(top)


# Do not let event-time annotations appear on the PSD frequency axis.
_FIBERLYSE_V24_AXIS_LOOKS_LIKE_TIME_BEFORE = PlotTabTk._axis_looks_like_time
def _fiberlyse_v24_axis_looks_like_time(ax) -> bool:
    if bool(getattr(ax, '_fiberlyse_not_time_axis', False)):
        return False
    return bool(_FIBERLYSE_V24_AXIS_LOOKS_LIKE_TIME_BEFORE(ax))
PlotTabTk._axis_looks_like_time = staticmethod(_fiberlyse_v24_axis_looks_like_time)

# Final authoritative V24 frequency methods (placed after all V23.7 overrides).
ChannelTabsTk._draw_frequency = _fiberlyse_v24_draw_frequency
ChannelTabsTk._export_frequency = _fiberlyse_v24_export_frequency
ChannelTabsTk.show_frequency_settings = _fiberlyse_v24_show_frequency_settings

# Add a visible settings button to every channel's frequency plot controls.
_FIBERLYSE_V24_CHANNELTABS_INIT_BEFORE = ChannelTabsTk.__init__
def _fiberlyse_v24_channeltabs_init(self, *args, **kwargs):
    _FIBERLYSE_V24_CHANNELTABS_INIT_BEFORE(self, *args, **kwargs)
    try:
        _fiberlyse_v24_init_frequency_settings(self.res)
        parent = self.tab_freq.save_btn.master
        self.freq_settings_btn = ttk.Button(parent, text='Frequency settings…', command=self.show_frequency_settings)
        self.freq_settings_btn.pack(side=tk.LEFT, padx=(0,8))
    except Exception:
        pass
ChannelTabsTk.__init__ = _fiberlyse_v24_channeltabs_init

# Update the application title/status without changing the V23.7 layout.
_FIBERLYSE_V24_MAINAPP_INIT_BEFORE = MainAppTk.__init__
def _fiberlyse_v24_mainapp_init(self, *args, **kwargs):
    _FIBERLYSE_V24_MAINAPP_INIT_BEFORE(self, *args, **kwargs)
    try:self.root.title('FiberLyse V24')
    except Exception:pass
    try:self.status.set('Ready. V24: universal CSV import + timestamp-aware adaptive frequency analysis.')
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_v24_mainapp_init

# Frequency metadata should also appear in the normal graph export metadata
# once settings have been configured, without altering the scientific arrays.
_FIBERLYSE_V24_META_BEFORE = ChannelTabsTk._meta_df
def _fiberlyse_v24_meta_df(self) -> pd.DataFrame:
    base = _FIBERLYSE_V24_META_BEFORE(self); r=self.res; _fiberlyse_v24_init_frequency_settings(r)
    diag = _fiberlyse_v24_frequency_diagnostics(np.asarray(r.t_exc,dtype=float), np.asarray(r.dFF_nointerp,dtype=float))
    rows=[
        ('frequency_current_fs_hz', diag.get('fs_hz', np.nan)),
        ('frequency_current_nyquist_hz', diag.get('nyquist_hz', np.nan)),
        ('frequency_band_mode', getattr(r,'freq_band_mode',_FIBERLYSE_V24_FREQ_MODE_STANDARD)),
        ('frequency_psd_window_s', getattr(r,'freq_psd_window_s',_FIBERLYSE_V24_DEFAULT_PSD_WINDOW_S)),
        ('frequency_psd_overlap_pct', getattr(r,'freq_psd_overlap_pct',_FIBERLYSE_V24_DEFAULT_PSD_OVERLAP_PCT)),
        ('frequency_filter_order', getattr(r,'freq_filter_order',DEFAULT_BUTTER_ORDER)),
        ('frequency_original_fs_hz_metadata_only', getattr(r,'freq_original_fs_hz',np.nan)),
        ('frequency_downsample_factor_metadata_only', getattr(r,'freq_downsample_factor',np.nan)),
        ('frequency_antialias_status', getattr(r,'freq_antialias_status','Unknown')),
    ]
    return pd.concat([base,pd.DataFrame(rows,columns=['key','value'])],ignore_index=True)
ChannelTabsTk._meta_df = _fiberlyse_v24_meta_df

# ---- End FiberLyse V24 universal frequency-analysis extension -------------

# ---- FiberLyse V24.1 timestamp-only sampling-rate UI extension ------------
# The editable Acq FPS control is removed from the visible GUI.  The legacy
# attribute remains internally so older extension code cannot fail, but V24.1
# settings readers always supply NaN and every importer derives eff_fs_hz from
# the actual study-signal timestamps.

_FIBERLYSE_V241_MAINAPP_INIT_BEFORE = MainAppTk.__init__
def _fiberlyse_v241_mainapp_init(self, *args, **kwargs):
    _FIBERLYSE_V241_MAINAPP_INIT_BEFORE(self, *args, **kwargs)
    try:self.root.title('FiberLyse V24.1')
    except Exception:pass
    try:
        # Keep the compatibility variable, but make it impossible for legacy
        # paths to interpret the old default 40 Hz as a user instruction.
        self.var_acq_fps.set('nan')
    except Exception:pass
    try:
        row2 = self.spin_smooth_win.master
        # Permanently remove the obsolete acquisition-FPS widgets.  V23.7's
        # responsive-layout callbacks can re-grid widgets after initialization,
        # so grid_remove() alone is not sufficient.  Destroying them means no
        # later resize/layout callback can make the control visible again.
        try:
            old_spin = getattr(self, 'spin_acq_fps', None)
            if old_spin is not None:
                old_spin.destroy()
        except Exception:pass
        try:
            old_lbl = _fiberlyse_find_label(row2, 'Acq FPS (Hz):')
            if old_lbl is not None:old_lbl.destroy()
        except Exception:pass
        try:
            old_note = _fiberlyse_find_label(row2, '(filter fs measured from timestamps)')
            if old_note is not None:old_note.destroy()
        except Exception:pass
        # Reclaim the space for smoothing and a read-only timestamp-derived rate.
        try:
            smooth_lbl = _fiberlyse_find_label(row2, 'Smooth win (samples):')
            if smooth_lbl is not None:smooth_lbl.grid(row=0, column=7, sticky='w')
        except Exception:pass
        try:self.spin_smooth_win.grid(row=0, column=8, padx=(6,12), sticky='w')
        except Exception:pass
        self.var_detected_signal_rate = tk.StringVar(value='Detected signal rate: after analysis')
        self.lbl_detected_signal_rate = ttk.Label(row2, textvariable=self.var_detected_signal_rate)
        self.lbl_detected_signal_rate.grid(row=0, column=9, columnspan=3, sticky='w', padx=(4,0))
    except Exception:pass
    try:self.status.set('Ready. V24.1: sampling rate is calculated automatically from timestamps.')
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_v241_mainapp_init


def _fiberlyse_v241_update_detected_rate_label(self) -> None:
    var = getattr(self, 'var_detected_signal_rate', None)
    if var is None:return
    results = dict(getattr(self, '_results', {}) or {})
    if not results:
        try:var.set('Detected signal rate: after analysis')
        except Exception:pass
        return
    # Prefer traces belonging to the currently selected source file.
    selected_key = None
    try:selected_key = self._get_selected_view_file_key()
    except Exception:selected_key = None
    chosen=[]
    for mid,res in results.items():
        if selected_key is not None:
            try:
                if self._mouse_file_key(mid) != selected_key:continue
            except Exception:pass
        try:
            fs = estimate_fs_from_t(np.asarray(res.t_exc, dtype=float))
        except Exception:
            fs = np.nan
        if np.isfinite(fs) and fs > 0:chosen.append(float(fs))
    if not chosen and selected_key is not None:
        for res in results.values():
            try:fs=estimate_fs_from_t(np.asarray(res.t_exc,dtype=float))
            except Exception:fs=np.nan
            if np.isfinite(fs) and fs>0:chosen.append(float(fs))
    if not chosen:
        text='Detected signal rate: unavailable (check timestamps)'
    else:
        lo=float(min(chosen));hi=float(max(chosen))
        if abs(hi-lo) <= max(1e-9, 0.001*max(abs(lo),abs(hi),1.0)):
            text=f'Detected signal rate: {float(np.median(chosen)):.5g} Hz (timestamps)'
        else:
            text=f'Detected signal rates: {lo:.5g}–{hi:.5g} Hz (timestamps)'
    try:var.set(text)
    except Exception:pass
MainAppTk._update_detected_signal_rate_label = _fiberlyse_v241_update_detected_rate_label

_FIBERLYSE_V241_ANALYSIS_FINISHED_BEFORE = MainAppTk.on_analysis_finished
def _fiberlyse_v241_on_analysis_finished(self, *args, **kwargs):
    out = _FIBERLYSE_V241_ANALYSIS_FINISHED_BEFORE(self, *args, **kwargs)
    try:self._update_detected_signal_rate_label()
    except Exception:pass
    return out
MainAppTk.on_analysis_finished = _fiberlyse_v241_on_analysis_finished

_FIBERLYSE_V241_VIEW_FILE_CHANGED_BEFORE = MainAppTk.on_view_file_changed
def _fiberlyse_v241_on_view_file_changed(self, *args, **kwargs):
    out = _FIBERLYSE_V241_VIEW_FILE_CHANGED_BEFORE(self, *args, **kwargs)
    try:self._update_detected_signal_rate_label()
    except Exception:pass
    return out
MainAppTk.on_view_file_changed = _fiberlyse_v241_on_view_file_changed

_FIBERLYSE_V241_CLEAR_CSVS_BEFORE = MainAppTk.clear_csvs
def _fiberlyse_v241_clear_csvs(self, *args, **kwargs):
    out = _FIBERLYSE_V241_CLEAR_CSVS_BEFORE(self, *args, **kwargs)
    try:
        if not getattr(self, '_results', None):self._update_detected_signal_rate_label()
    except Exception:pass
    return out
MainAppTk.clear_csvs = _fiberlyse_v241_clear_csvs

# Remove obsolete acquisition-FPS keys from standard metadata.  Original
# acquisition history remains available through V24's explicit
# frequency_original_fs_hz_metadata_only/downsample metadata fields.
_FIBERLYSE_V241_META_BEFORE = ChannelTabsTk._meta_df
def _fiberlyse_v241_meta_df(self) -> pd.DataFrame:
    df = _FIBERLYSE_V241_META_BEFORE(self)
    try:
        if isinstance(df, pd.DataFrame) and 'key' in df.columns:
            obsolete = {'acq_fps_user_hz', 'freq_acq_fps_hz'}
            df = df.loc[~df['key'].astype(str).isin(obsolete)].reset_index(drop=True)
    except Exception:pass
    return df
ChannelTabsTk._meta_df = _fiberlyse_v241_meta_df

# ---- End FiberLyse V24.1 timestamp-only sampling-rate UI extension --------


# ---- FiberLyse V25 visual design refresh ---------------------------------
# V25 is intentionally cosmetic.  It does not change import parsing,
# preprocessing, fitting, normalization, frequency calculations, AUC, batch
# calculations, or export payloads.  This final layer only standardizes visual
# hierarchy, spacing, control emphasis, and plot presentation across the many
# legacy GUI-extension layers above.

FIBERLYSE_VERSION = 'V25'

_FIBERLYSE_V25_THEME_LIGHT = {
    'bg': '#F3F5F7',
    'panel': '#FFFFFF',
    'fg': '#18212B',
    'muted': '#687789',
    'entry': '#F8FAFC',
    'border': '#D9E0E7',
    'select': '#DCE8FF',
    'select_fg': '#153B73',
    'tooltip_bg': '#18212B',
    'tooltip_fg': '#FFFFFF',
    'plot_bg': '#FFFFFF',
    'plot_fg': '#243140',
    'grid': '#C8D1DC',
    'primary': '#316FEA',
    'primary_hover': '#285FCC',
    'primary_fg': '#FFFFFF',
    'danger': '#C74752',
    'danger_hover': '#AA3741',
    'danger_fg': '#FFFFFF',
    'success': '#247A59',
    'warning': '#9A6419',
    'surface_alt': '#F7F9FB',
}

_FIBERLYSE_V25_THEME_DARK = {
    'bg': '#0F151C',
    'panel': '#18212B',
    'fg': '#E7EDF4',
    'muted': '#98A6B6',
    'entry': '#202B36',
    'border': '#2D3A47',
    'select': '#294B79',
    'select_fg': '#F5F8FC',
    'tooltip_bg': '#080C11',
    'tooltip_fg': '#F5F8FC',
    'plot_bg': '#151D26',
    'plot_fg': '#E2E9F0',
    'grid': '#44515F',
    'primary': '#5C8FF4',
    'primary_hover': '#6D9CFA',
    'primary_fg': '#08111E',
    'danger': '#E36B73',
    'danger_hover': '#F07B83',
    'danger_fg': '#160709',
    'success': '#5BC18D',
    'warning': '#E2A955',
    'surface_alt': '#1C2732',
}


def _fiberlyse_v25_is_dark(app=None) -> bool:
    try:
        return bool(app.var_night_mode.get())
    except Exception:
        try:return bool(getattr(app, '_fiberlyse_night_mode', False))
        except Exception:return False


def _fiberlyse_v25_tokens(app=None, widget=None) -> Dict[str, str]:
    if app is None and widget is not None:
        try:app = getattr(widget.winfo_toplevel(), '_fiberlyse_app', None)
        except Exception:app = None
    return dict(_FIBERLYSE_V25_THEME_DARK if _fiberlyse_v25_is_dark(app) else _FIBERLYSE_V25_THEME_LIGHT)


# Rebind the shared theme-token provider.  Existing V23+ dialog and Matplotlib
# theme code calls this function dynamically, so the modern V25 palette also
# propagates to old dialogs without duplicating them.
_fiberlyse_v23_theme_tokens = _fiberlyse_v25_tokens


def _fiberlyse_v25_pick_ui_font(root) -> str:
    """Prefer native modern UI fonts while remaining portable."""
    try:
        import tkinter.font as _tkfont_v25
        families = {str(x).lower(): str(x) for x in _tkfont_v25.families(root)}
        for candidate in ('Segoe UI', 'SF Pro Text', 'Inter', 'Helvetica Neue', 'Arial'):
            hit = families.get(candidate.lower())
            if hit:return hit
    except Exception:
        pass
    return 'TkDefaultFont'


def _fiberlyse_v25_apply_ttk_theme(app) -> None:
    tok = _fiberlyse_v25_tokens(app=app)
    try:
        style = ttk.Style(app.root)
        try:
            if 'clam' in style.theme_names():style.theme_use('clam')
        except Exception:pass
        font_family = _fiberlyse_v25_pick_ui_font(app.root)
        try:
            import tkinter.font as _tkfont_v25
            for named, size, weight in [
                ('TkDefaultFont', 10, 'normal'),
                ('TkTextFont', 10, 'normal'),
                ('TkMenuFont', 10, 'normal'),
                ('TkHeadingFont', 10, 'bold'),
                ('TkCaptionFont', 10, 'bold'),
            ]:
                try:_tkfont_v25.nametofont(named).configure(family=font_family, size=size, weight=weight)
                except Exception:pass
        except Exception:pass

        app.root.configure(background=tok['bg'])
        # Surfaces / cards
        style.configure('TFrame', background=tok['bg'])
        style.configure('Fiber.TFrame', background=tok['bg'])
        style.configure('V25.Content.TFrame', background=tok['bg'])
        style.configure('Panel.TFrame', background=tok['panel'], borderwidth=0, relief='flat')
        style.configure('V25.Card.TFrame', background=tok['panel'], bordercolor=tok['border'], borderwidth=1, relief='solid')
        style.configure('V25.Toolbar.TFrame', background=tok['panel'], bordercolor=tok['border'], borderwidth=1, relief='solid')
        style.configure('V25.PlotTools.TFrame', background=tok['panel'])
        # Labels
        style.configure('TLabel', background=tok['bg'], foreground=tok['fg'])
        style.configure('Panel.TLabel', background=tok['panel'], foreground=tok['fg'])
        style.configure('V25.Card.TLabel', background=tok['panel'], foreground=tok['fg'])
        style.configure('V25.Toolbar.TLabel', background=tok['panel'], foreground=tok['fg'])
        style.configure('Title.TLabel', background=tok['panel'], foreground=tok['fg'], font=(font_family, 12, 'bold'))
        style.configure('Hint.TLabel', background=tok['panel'], foreground=tok['muted'], font=(font_family, 9))
        style.configure('V25.Badge.TLabel', background=tok['surface_alt'], foreground=tok['muted'], padding=(8, 4), font=(font_family, 9))
        style.configure('Tooltip.TLabel', background=tok['tooltip_bg'], foreground=tok['tooltip_fg'], relief='solid', borderwidth=1, padding=(7,4))
        # Buttons: clear hierarchy with one strong primary action.
        style.configure('TButton', background=tok['surface_alt'], foreground=tok['fg'], bordercolor=tok['border'], borderwidth=1, relief='flat', padding=(11, 6), font=(font_family, 10))
        style.map('TButton',
                  background=[('active', tok['select']), ('pressed', tok['select']), ('disabled', tok['panel'])],
                  foreground=[('disabled', tok['muted'])],
                  bordercolor=[('focus', tok['primary']), ('active', tok['primary'])])
        style.configure('Primary.TButton', background=tok['primary'], foreground=tok['primary_fg'], bordercolor=tok['primary'], borderwidth=1, relief='flat', padding=(13, 7), font=(font_family, 10, 'bold'))
        style.map('Primary.TButton', background=[('active', tok['primary_hover']), ('pressed', tok['primary_hover']), ('disabled', tok['border'])], foreground=[('disabled', tok['muted'])], bordercolor=[('active', tok['primary_hover']), ('focus', tok['primary_hover'])])
        # Accent is kept as an alias because older FiberLyse code already uses it.
        style.configure('Accent.TButton', background=tok['primary'], foreground=tok['primary_fg'], bordercolor=tok['primary'], borderwidth=1, relief='flat', padding=(13, 7), font=(font_family, 10, 'bold'))
        style.map('Accent.TButton', background=[('active', tok['primary_hover']), ('pressed', tok['primary_hover']), ('disabled', tok['border'])], foreground=[('disabled', tok['muted'])])
        style.configure('Secondary.TButton', background=tok['panel'], foreground=tok['fg'], bordercolor=tok['border'], borderwidth=1, relief='flat', padding=(11, 6))
        style.map('Secondary.TButton', background=[('active', tok['surface_alt']), ('pressed', tok['select'])], bordercolor=[('active', tok['primary']), ('focus', tok['primary'])])
        style.configure('Quiet.TButton', background=tok['panel'], foreground=tok['muted'], bordercolor=tok['panel'], borderwidth=0, relief='flat', padding=(9, 6))
        style.map('Quiet.TButton', background=[('active', tok['surface_alt']), ('pressed', tok['select'])], foreground=[('active', tok['fg']), ('disabled', tok['muted'])])
        style.configure('Danger.TButton', background=tok['panel'], foreground=tok['danger'], bordercolor=tok['border'], borderwidth=1, relief='flat', padding=(11, 6))
        style.map('Danger.TButton', background=[('active', tok['danger']), ('pressed', tok['danger_hover'])], foreground=[('active', tok['danger_fg']), ('pressed', tok['danger_fg']), ('disabled', tok['muted'])], bordercolor=[('active', tok['danger']), ('focus', tok['danger'])])
        # Form controls
        style.configure('TCheckbutton', background=tok['bg'], foreground=tok['fg'], padding=(2,3))
        style.map('TCheckbutton', background=[('active', tok['bg'])], foreground=[('disabled', tok['muted'])])
        style.configure('V25.Card.TCheckbutton', background=tok['panel'], foreground=tok['fg'], padding=(2,3))
        style.map('V25.Card.TCheckbutton', background=[('active', tok['panel'])], foreground=[('disabled', tok['muted'])])
        style.configure('TRadiobutton', background=tok['bg'], foreground=tok['fg'], padding=(2,3))
        style.map('TRadiobutton', background=[('active', tok['bg'])], foreground=[('disabled', tok['muted'])])
        for sty in ('TEntry','TSpinbox','TCombobox'):
            style.configure(sty, fieldbackground=tok['entry'], background=tok['entry'], foreground=tok['fg'], insertcolor=tok['fg'], arrowcolor=tok['fg'], bordercolor=tok['border'], lightcolor=tok['border'], darkcolor=tok['border'], padding=5)
            style.map(sty, fieldbackground=[('readonly', tok['entry']), ('disabled', tok['panel'])], foreground=[('readonly', tok['fg']), ('disabled', tok['muted'])], bordercolor=[('focus', tok['primary'])])
        # Tabs / data tables
        style.configure('TNotebook', background=tok['bg'], borderwidth=0, tabmargins=(0, 0, 0, 0))
        style.configure('TNotebook.Tab', background=tok['bg'], foreground=tok['muted'], borderwidth=0, padding=(14, 8), font=(font_family, 10))
        style.map('TNotebook.Tab', background=[('selected', tok['panel']), ('active', tok['surface_alt'])], foreground=[('selected', tok['fg']), ('active', tok['fg'])], expand=[('selected', (0,0,0,1))])
        style.configure('Queue.Treeview', rowheight=30, fieldbackground=tok['panel'], background=tok['panel'], foreground=tok['fg'], bordercolor=tok['border'], relief='flat', font=(font_family, 10))
        style.map('Queue.Treeview', background=[('selected', tok['select'])], foreground=[('selected', tok['select_fg'])])
        style.configure('Queue.Treeview.Heading', background=tok['surface_alt'], foreground=tok['muted'], font=(font_family, 9, 'bold'), relief='flat', padding=(6,6), bordercolor=tok['border'])
        style.map('Queue.Treeview.Heading', background=[('active', tok['select'])], foreground=[('active', tok['fg'])])
        style.configure('Treeview', rowheight=28, fieldbackground=tok['panel'], background=tok['panel'], foreground=tok['fg'])
        style.map('Treeview', background=[('selected', tok['select'])], foreground=[('selected', tok['select_fg'])])
        # Frames / progress / separators
        style.configure('TLabelframe', background=tok['bg'], foreground=tok['fg'], bordercolor=tok['border'], borderwidth=1, relief='solid')
        style.configure('TLabelframe.Label', background=tok['bg'], foreground=tok['fg'], font=(font_family, 10, 'bold'))
        style.configure('Horizontal.TProgressbar', background=tok['primary'], troughcolor=tok['surface_alt'], bordercolor=tok['border'], lightcolor=tok['primary'], darkcolor=tok['primary'], thickness=8)
        style.configure('TPanedwindow', background=tok['bg'], sashwidth=6)
        style.configure('TSeparator', background=tok['border'])
    except Exception as e:
        try:print(f'V25 visual theme failed: {e}', file=sys.stderr)
        except Exception:pass


# Make the legacy refresh path apply V25 styles last, after all older ttk rules.
_FIBERLYSE_V25_TTK_THEME_BEFORE = _fiberlyse_v23_apply_ttk_theme
def _fiberlyse_v23_apply_ttk_theme(app) -> None:
    try:_FIBERLYSE_V25_TTK_THEME_BEFORE(app)
    except Exception:pass
    _fiberlyse_v25_apply_ttk_theme(app)


def _fiberlyse_v25_style_frame_tree(widget, frame_style: str, label_style: str, check_style: str) -> None:
    """Apply matching surface styles to descendants without reparenting them."""
    try:children = list(widget.winfo_children())
    except Exception:children = []
    for child in children:
        try:cls = str(child.winfo_class())
        except Exception:cls = ''
        try:
            if cls == 'TFrame':child.configure(style=frame_style)
            elif cls == 'TLabel':child.configure(style=label_style)
            elif cls == 'TCheckbutton':child.configure(style=check_style)
        except Exception:pass
        if cls in ('TFrame','Frame','TLabelframe','Labelframe'):
            _fiberlyse_v25_style_frame_tree(child, frame_style, label_style, check_style)


def _fiberlyse_v25_theme_classic_widgets(widget, tok: Dict[str,str]) -> None:
    """Modern flat treatment for classic Tk widgets and Matplotlib toolbar buttons."""
    try:cls = str(widget.winfo_class())
    except Exception:cls = ''
    try:
        if cls == 'Listbox':
            widget.configure(background=tok['panel'], foreground=tok['fg'], selectbackground=tok['select'], selectforeground=tok['select_fg'], relief='flat', borderwidth=1, highlightthickness=1, highlightbackground=tok['border'], highlightcolor=tok['primary'])
        elif cls == 'Text':
            widget.configure(background=tok['panel'], foreground=tok['fg'], insertbackground=tok['fg'], selectbackground=tok['select'], selectforeground=tok['select_fg'], relief='flat', borderwidth=1, highlightthickness=1, highlightbackground=tok['border'], highlightcolor=tok['primary'], padx=7, pady=6)
        elif cls == 'Canvas':
            widget.configure(background=tok['panel'], highlightbackground=tok['border'], highlightthickness=0)
        elif cls == 'Button':
            widget.configure(background=tok['panel'], foreground=tok['fg'], activebackground=tok['select'], activeforeground=tok['select_fg'], relief='flat', borderwidth=0, highlightthickness=0, padx=6, pady=4)
        elif cls in ('Checkbutton','Radiobutton'):
            widget.configure(background=tok['panel'], foreground=tok['fg'], activebackground=tok['panel'], activeforeground=tok['fg'], highlightthickness=0)
        elif cls in ('Entry','Spinbox'):
            widget.configure(background=tok['entry'], foreground=tok['fg'], insertbackground=tok['fg'], relief='flat', borderwidth=1, highlightthickness=1, highlightbackground=tok['border'], highlightcolor=tok['primary'])
    except Exception:pass
    try:
        for child in widget.winfo_children():_fiberlyse_v25_theme_classic_widgets(child, tok)
    except Exception:pass


# Extend the older recursive classic-widget theming so night-mode toggles also
# retain the modern border/relief treatment.
_FIBERLYSE_V25_CLASSIC_THEME_BEFORE = _fiberlyse_v23_apply_classic_widget_theme
def _fiberlyse_v23_apply_classic_widget_theme(widget, tok: Dict[str,str]) -> None:
    try:_FIBERLYSE_V25_CLASSIC_THEME_BEFORE(widget, tok)
    except Exception:pass
    _fiberlyse_v25_theme_classic_widgets(widget, tok)


def _fiberlyse_v25_apply_control_hierarchy(app) -> None:
    # Main toolbar
    for name, style_name in [
        ('btn_run','Primary.TButton'),
        ('btn_add','Secondary.TButton'),
        ('btn_clear','Quiet.TButton'),
        ('btn_cancel','Danger.TButton'),
        ('btn_csv_setup','Secondary.TButton'),
        ('btn_add_folder','Secondary.TButton'),
        ('btn_analyze_selected','Secondary.TButton'),
        ('btn_retry_failed','Secondary.TButton'),
        ('btn_remove_selected','Danger.TButton'),
        ('btn_open_folder','Quiet.TButton'),
        ('btn_shortcuts','Quiet.TButton'),
    ]:
        try:getattr(app, name).configure(style=style_name)
        except Exception:pass
    try:app.chk_night_mode.configure(style='V25.Card.TCheckbutton')
    except Exception:pass
    try:app.lbl_detected_signal_rate.configure(style='V25.Badge.TLabel')
    except Exception:pass
    try:app.lbl_import_profile.configure(style='V25.Toolbar.TLabel')
    except Exception:pass


def _fiberlyse_v25_polish_main_layout(app) -> None:
    """Use existing containers as modern cards; no controls are replaced."""
    try:
        top = app.btn_run.master
        top.configure(style='V25.Toolbar.TFrame', padding=(12,9))
        top.pack_configure(padx=14, pady=(12,6), fill=tk.X)
        _fiberlyse_v25_style_frame_tree(top, 'V25.Toolbar.TFrame', 'V25.Toolbar.TLabel', 'V25.Card.TCheckbutton')
    except Exception:pass
    try:
        analysis_row = app.spin_smooth_win.master
        analysis_row.configure(style='V25.Card.TFrame', padding=(12,8))
        analysis_row.pack_configure(padx=14, pady=(0,6), fill=tk.X)
        _fiberlyse_v25_style_frame_tree(analysis_row, 'V25.Card.TFrame', 'V25.Card.TLabel', 'V25.Card.TCheckbutton')
    except Exception:pass
    try:
        norm_row = app.cmb_norm.master
        norm_row.configure(style='V25.Card.TFrame', padding=(12,8))
        norm_row.pack_configure(padx=14, pady=(0,8), fill=tk.X)
        _fiberlyse_v25_style_frame_tree(norm_row, 'V25.Card.TFrame', 'V25.Card.TLabel', 'V25.Card.TCheckbutton')
    except Exception:pass
    try:
        app.workspace_pane.pack_configure(padx=14, pady=(0,10), fill=tk.BOTH, expand=True)
        app.queue_panel.configure(style='V25.Card.TFrame', padding=(12,12))
        app.workspace_holder.configure(style='V25.Content.TFrame')
        _fiberlyse_v25_style_frame_tree(app.queue_panel, 'Panel.TFrame', 'Panel.TLabel', 'V25.Card.TCheckbutton')
        app.root.after_idle(lambda: app.workspace_pane.sashpos(0, 315) if app.workspace_pane.winfo_width() > 850 else None)
    except Exception:pass
    try:
        app.queue_tree.column('#0', width=170, minwidth=110, stretch=True)
        app.queue_tree.column('status', width=82, minwidth=76, stretch=False)
        app.queue_tree.column('channels', width=46, minwidth=40, stretch=False, anchor='center')
        app.queue_tree.column('note', width=120, minwidth=70, stretch=True)
    except Exception:pass
    try:
        # Status bar remains visually separate but understated.
        for child in app.root.winfo_children():
            try:
                if str(child.winfo_class()) == 'TLabel' and str(child.cget('textvariable')) == str(app.status):
                    child.configure(style='Hint.TLabel', padding=(10,5), relief='flat')
            except Exception:pass
    except Exception:pass
    _fiberlyse_v25_apply_control_hierarchy(app)


def _fiberlyse_v25_polish_plot(tab) -> None:
    tok = _fiberlyse_v25_tokens(widget=tab)
    try:tab.fig.patch.set_facecolor(tok['bg'])
    except Exception:pass
    for ax in list(getattr(tab.fig, 'axes', []) or []):
        try:ax.set_facecolor(tok['plot_bg'])
        except Exception:pass
        try:
            ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
            for side in ('left','bottom'):
                ax.spines[side].set_color(tok['border']);ax.spines[side].set_linewidth(0.8)
        except Exception:pass
        try:ax.tick_params(axis='both', colors=tok['muted'], labelcolor=tok['muted'], length=3.5, width=0.8)
        except Exception:pass
        try:
            ax.xaxis.label.set_color(tok['plot_fg']);ax.yaxis.label.set_color(tok['plot_fg']);ax.title.set_color(tok['plot_fg'])
        except Exception:pass
        try:
            ax.set_axisbelow(True);ax.grid(True, color=tok['grid'], alpha=0.20, linewidth=0.65)
        except Exception:pass
        try:
            leg=ax.get_legend()
            if leg is not None:
                frame=leg.get_frame();frame.set_facecolor(tok['plot_bg']);frame.set_edgecolor(tok['border']);frame.set_alpha(0.92)
                for txt in leg.get_texts():txt.set_color(tok['plot_fg'])
                try:leg.get_title().set_color(tok['plot_fg'])
                except Exception:pass
        except Exception:pass
    try:
        tools = tab.save_btn.master
        tools.configure(style='V25.PlotTools.TFrame', padding=(8,6))
        tab.export_btn.configure(style='Primary.TButton')
        tab.save_btn.configure(style='Secondary.TButton')
        for name in ('auc_btn','axis_interval_btn'):
            btn=getattr(tab,name,None)
            if btn is not None:btn.configure(style='Quiet.TButton')
    except Exception:pass
    try:_fiberlyse_v25_theme_classic_widgets(tab.toolbar, tok)
    except Exception:pass


# Plot widgets are created lazily, so polish every new PlotTab after the entire
# inherited initializer chain has completed.
_FIBERLYSE_V25_PLOTTAB_INIT_BEFORE = PlotTabTk.__init__
def _fiberlyse_v25_plottab_init(self, *args, **kwargs):
    _FIBERLYSE_V25_PLOTTAB_INIT_BEFORE(self, *args, **kwargs)
    try:_fiberlyse_v25_polish_plot(self)
    except Exception:pass
PlotTabTk.__init__ = _fiberlyse_v25_plottab_init

_FIBERLYSE_V25_REDRAW_BEFORE = PlotTabTk.redraw
def _fiberlyse_v25_redraw(self):
    out = _FIBERLYSE_V25_REDRAW_BEFORE(self)
    try:_fiberlyse_v25_polish_plot(self);self.canvas.draw_idle()
    except Exception:pass
    return out
PlotTabTk.redraw = _fiberlyse_v25_redraw


# Frequency settings is added after PlotTab initialization, so style that
# channel-level button at the end of ChannelTabs construction.
_FIBERLYSE_V25_CHANNELTABS_INIT_BEFORE = ChannelTabsTk.__init__
def _fiberlyse_v25_channeltabs_init(self, *args, **kwargs):
    _FIBERLYSE_V25_CHANNELTABS_INIT_BEFORE(self, *args, **kwargs)
    try:self.freq_settings_btn.configure(style='Quiet.TButton')
    except Exception:pass
    try:
        for tab in (self.tab_raw,self.tab_art,self.tab_fit,self.tab_norm,self.tab_norm_smooth,self.tab_freq):
            _fiberlyse_v25_polish_plot(tab)
    except Exception:pass
ChannelTabsTk.__init__ = _fiberlyse_v25_channeltabs_init


def _fiberlyse_v25_theme_toplevel(event, app) -> None:
    try:w = event.widget
    except Exception:return
    try:
        if not isinstance(w, tk.Toplevel):return
    except Exception:return
    tok = _fiberlyse_v25_tokens(app=app)
    try:w.configure(background=tok['bg'])
    except Exception:pass
    try:w.after_idle(lambda win=w,t=tok:_fiberlyse_v25_theme_classic_widgets(win,t))
    except Exception:pass


_FIBERLYSE_V25_MAIN_INIT_BEFORE = MainAppTk.__init__
def _fiberlyse_v25_mainapp_init(self, *args, **kwargs):
    _FIBERLYSE_V25_MAIN_INIT_BEFORE(self, *args, **kwargs)
    try:self.root.title('FiberLyse V25')
    except Exception:pass
    try:_fiberlyse_v25_apply_ttk_theme(self)
    except Exception:pass
    try:_fiberlyse_v25_polish_main_layout(self)
    except Exception:pass
    try:_fiberlyse_v25_theme_classic_widgets(self.root, _fiberlyse_v25_tokens(app=self))
    except Exception:pass
    try:self.root.bind_all('<Map>', lambda e,app=self:_fiberlyse_v25_theme_toplevel(e,app), add='+')
    except Exception:pass
    # Re-polish after responsive callbacks have made their first layout pass.
    try:self.root.after(350, lambda app=self:(_fiberlyse_v25_apply_ttk_theme(app), _fiberlyse_v25_polish_main_layout(app)))
    except Exception:pass
MainAppTk.__init__ = _fiberlyse_v25_mainapp_init


# Keep visual hierarchy after theme toggles.  Scientific state is untouched.
_FIBERLYSE_V25_REFRESH_THEME_BEFORE = _fiberlyse_v23_refresh_theme
def _fiberlyse_v23_refresh_theme(app) -> None:
    _FIBERLYSE_V25_REFRESH_THEME_BEFORE(app)
    try:_fiberlyse_v25_apply_ttk_theme(app);_fiberlyse_v25_polish_main_layout(app)
    except Exception:pass
    try:_fiberlyse_v25_theme_classic_widgets(app.root, _fiberlyse_v25_tokens(app=app))
    except Exception:pass
    for tab in _fiberlyse_v23_iter_plot_tabs(app):
        try:_fiberlyse_v25_polish_plot(tab);tab.canvas.draw_idle()
        except Exception:pass

# ---- End FiberLyse V25 visual design refresh -----------------------------

if __name__ == '__main__':main()
