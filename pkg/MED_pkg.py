"""
MED_pkg — Movement Element Decomposition Package
=============================================

Overview
--------
This module implements the Movement Element Decomposition (MED) algorithm, a method
for decomposing continuous kinematic time-series data into discrete movement
elements (MEs). Each element is a bounded interval of movement (motor primitive),
identified by its velocity profile and validated against configurable thresholds
for displacement, duration, and velocity.

The algorithm is designed for positional data in one, two, or three spatial
dimensions, and produces per-dimension and aggregate features for each detected
ME.

Pipeline Summary
----------------
Raw position data undergoes the following sequential processing steps:

    1. Unit conversion  — input coordinates are normalised to meters.
    2. Low-pass filtering — a zero-phase Butterworth filter (default: 4th order,
       10 Hz cutoff) removes high-frequency noise while preserving movement
       kinematics. Border samples affected by filter edge effects are discarded
       unless trim_border is set to False.
    3. Velocity computation — instantaneous velocity is derived via first-order
       finite differences of the filtered position signal.
    4. Segmentation (segment_MED) — velocity peaks and zero-crossings are located
       per dimension to define candidate movement elements.
    5. Element validation — candidates failing the minimum displacement (minD),
       duration (minT), or velocity (minV) criteria are rejected.
    6. Feature extraction (analyze_elements_MED) — for each valid element,
       scalar features are computed: displacement D [m], mean velocity V [m/s],
       duration T [s], similarity index W from the Hoff bell-shaped profile,
       coefficient of determination R², and peak count P.
    7. Scaling analysis (scaling_MED) — a power-law relationship V = K · D^alpha
       is estimated from the log-log regression of V on D per dimension.

Public API
----------
MED(*args, trim_border=True)
    Main entry point. Accepts position data and processing parameters, and returns
    a structured DataFrame of per-dimension movement statistics.

segment_MED(t, r, v, min_D, min_T, min_V)
    Identifies the frame indices delimiting each movement element in a 1-D
    velocity signal.

analyze_elements_MED(t, r, v, ME)
    Computes kinematic descriptors (D, V, T, W, R², P) for an array of
    pre-segmented movement elements.

scaling_MED(D_dict, V_dict, plot_name=None)
    Estimates power-law scaling exponents (alpha), coefficients (K), and
    goodness-of-fit (R²) from displacement–velocity pairs per dimension.
    Optionally saves a log-log scatter plot.

Dependencies
------------
numpy, scipy (signal), pandas, matplotlib

Notes
-----
- The Hoff velocity profile used for the similarity index (W) is the
  minimum-jerk bell-shaped curve: H(t) = V_mean · 30 · (t⁴ − 2t³ + t²),
  where t ∈ [0, 1] is normalised time.
- The border trimming applied after filtering removes `borderEffect` frames from
  each end, where borderEffect = round(2 · FPS / lp · order / 4). Disable
  with trim_border=False.
- This module is a Python port of an equivalent MATLAB implementation; numerical
  equivalence with the MATLAB version has been verified for the filtering step
  (padlen matched to MATLAB's default).

How to cite
-----------
DOI: https://doi.org/10.5281/zenodo.20510859

"""

import numpy as np
from scipy.signal import find_peaks
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as s


def MED(*args, trim_border=True):
    """
    Main entry point of the MED algorithm.

    Positional arguments (all after FPS are optional; pass '' or empty array to
    use the default):
      args[0]  movementData  — (N × dim) position array [m, cm, or mm]
      args[1]  FPS           — capture frequency in frames per second
      args[2]  unit          — 'm' (default) | 'cm' | 'mm'
      args[3]  limits        — np.array([minD, minT, minV]); default [0.003, 0.1, 0.01]
      args[4]  filter        — np.array([lp, order]); default [10, 4]
      args[5]  outputNames   — list of variable names to include in the output
      args[6]  plotScalingName — file path for the scaling plot, or '' to skip

    Keyword argument:
      trim_border  bool, default True; set False to retain the filter
                   edge-effect samples at both ends of the signal.
    """
    numberOfArgs = len(args)

    movementData = args[0]
    FPS = args[1]

    unit = "m"
    minD, minT, minV = 0.003, 0.1, 0.01
    lp, order = 10, 4
    outputNames = ["scaling", "D", "V", "T", "N", "Nt", "W", "R2", "P",
                   "D_all", "V_all", "T_all", "W_all", "R2_all", "P_all", "timeSeries", "ME"]
    plotScalingName = np.array([])

    # Conditional for if args is >= 3, setting units
    if numberOfArgs >= 3:
        if args[2] != "":
            unit = args[2]

    # Conditional for if args is >= 4, setting limits for identifying valid movement elements
    # minD is the minimum displacement [m], minT is the minimum duration [s], minV is the minimum velocity [m/s]
    if numberOfArgs >= 4:
        if args[3].size != 0:
            limits = args[3]
            minD, minT, minV = limits[0], limits[1], limits[2]

    # Conditional for if args is >= 5; lp is low-pass cutoff [Hz], order is filter order
    if numberOfArgs >= 5:
        if args[4].size != 0:
            filter = args[4]
            lp, order = filter[0], filter[1]

    if numberOfArgs >= 6:
        outputNames = args[5]

    # Conditional for if args is >= 7, checks if user specified plotScalingName
    if numberOfArgs >= 7:
        if args[6] != "":
            plotScalingName = args[6]

    # Convert units to meters if necessary
    if unit == "mm":
        movementData = movementData / 1000
    elif unit == "cm":
        movementData = movementData / 100

    # Apply Butterworth low-pass filter to position data
    b, a = s.butter(order, (2 * lp) / FPS, btype='low', analog=False, output='ba')
    movementData = s.filtfilt(b, a, movementData, axis=0, padlen=3 * (max(len(b), len(a)) - 1))

    borderEffect = int(round(2 * FPS / lp * order / 4))
    if trim_border:
        movementData = movementData[borderEffect:-borderEffect, :]

    # Calculate velocity and time
    v = np.diff(movementData, axis=0) * FPS
    movementData = movementData[:-1, :]
    t = np.arange(1, v.shape[0] + 1) / FPS

    # Determine dimensionality and prepare labels
    nDim = movementData.shape[1]

    if nDim == 1:
        dim = ["all"]
    elif nDim == 2:
        dim = ["x", "y", "all"]
    elif nDim == 3:
        dim = ["x", "y", "z", "all"]
    else:
        raise ValueError("nDim must be 1, 2, or 3.")

    outputVar = pd.DataFrame(
        index=dim,
        columns=["scaling", "D", "V", "T", "N", "Nt", "W", "R2", "P",
                 "D_all", "V_all", "T_all", "W_all", "R2_all", "P_all", "timeSeries", "ME"])

    # Find frames delimiting valid movement elements per dimension
    for i in range(nDim):
        outputVar.at[dim[i], "timeSeries"] = {"r": movementData[:, i], "v": v[:, i], "t": t}
        values = segment_MED(outputVar.loc[dim[i], "timeSeries"]["t"],
                             outputVar.loc[dim[i], "timeSeries"]["r"],
                             outputVar.loc[dim[i], "timeSeries"]["v"],
                             minD, minT, minV)
        if values.size != 0:
            outputVar.at[dim[i], "ME"] = values
        else:
            outputVar.at[dim[i], "ME"] = None

    cols_array = ["D_all", "V_all", "T_all", "W_all", "R2_all", "P_all"]
    for c in cols_array:
        if c not in outputVar.columns:
            outputVar[c] = None
        outputVar[c] = outputVar[c].astype(object)

    # Analyze movement elements per dimension (iterate backwards to allow safe deletion)
    for i in range(len(dim) - 2, -1, -1):
        it = dim[i]
        me_value = outputVar.loc[it, "ME"]

        # Drop dimension if no valid segments were found
        if me_value is None:
            outputVar.drop(index=it, inplace=True)
            del dim[i]
            nDim -= 1
            continue

        outputVar.at[it, "N"] = 0
        outputVar.at[it, "Nt"] = 0

        # Compute kinematic descriptors for each movement element
        (valueNIt, valueNtIt, valueDallIt, valueVallIt, valueTallIt,
         valueWallIt, valueR2allIt, valuePallIt) = analyze_elements_MED(
            outputVar.loc[it, "timeSeries"]["t"],
            outputVar.loc[it, "timeSeries"]["r"],
            outputVar.loc[it, "timeSeries"]["v"],
            outputVar.loc[it, "ME"])

        # Store per-element descriptor arrays
        outputVar.at[it, "N"] = valueNIt
        outputVar.at[it, "Nt"] = valueNtIt
        outputVar.at[it, "D_all"] = np.atleast_1d(valueDallIt).tolist()
        outputVar.at[it, "V_all"] = np.atleast_1d(valueVallIt).tolist()
        outputVar.at[it, "T_all"] = np.atleast_1d(valueTallIt).tolist()
        outputVar.at[it, "W_all"] = np.atleast_1d(valueWallIt).tolist()
        outputVar.at[it, "R2_all"] = np.atleast_1d(valueR2allIt).tolist()
        outputVar.at[it, "P_all"] = np.atleast_1d(valuePallIt).tolist()

    if nDim == 0:
        return pd.DataFrame(np.nan, index=[0], columns=outputNames)

    # Initialize combined results across all dimensions
    outputVar.at["all", "D_all"] = []
    outputVar.at["all", "V_all"] = []
    outputVar.at["all", "T_all"] = []
    outputVar.at["all", "W_all"] = []
    outputVar.at["all", "R2_all"] = []
    outputVar.at["all", "P_all"] = []
    outputVar.at["all", "N"] = 0
    outputVar.at["all", "Nt"] = 0

    for i in range(nDim):
        it = dim[i]
        outputVar.at["all", "D_all"].append(outputVar.at[it, "D_all"])
        outputVar.at["all", "V_all"].append(outputVar.at[it, "V_all"])
        outputVar.at["all", "N"] = outputVar["N"]["all"] + outputVar["N"][it]
        outputVar.at["all", "T_all"].append(outputVar.at[it, "T_all"])
        outputVar.at["all", "W_all"].append(outputVar.at[it, "W_all"])
        outputVar.at["all", "R2_all"].append(outputVar.at[it, "R2_all"])
        outputVar.at["all", "P_all"].append(outputVar.at[it, "P_all"])

    # Weighted average movement element rate
    try:
        outputVar.at["all", "Nt"] = (outputVar["N"]["all"] /
                                     (outputVar.loc[dim[i], "timeSeries"]["t"][-1] -
                                      outputVar.loc[dim[i], "timeSeries"]["t"][0]))
    except Exception:
        outputVar.at["all", "Nt"] = 0

    def safe_concat(data_list):
        valid_items = [x for x in data_list if x is not None]
        if valid_items:
            return np.concatenate(valid_items).tolist()
        else:
            return []

    for col in cols_array:
        outputVar.at["all", col] = safe_concat(outputVar.at["all", col])

    # Compute scaling law from displacement and velocity
    if "scaling" in outputVar.columns:
        alpha, K, R2 = scaling_MED(outputVar.loc[:, "D_all"], outputVar.loc[:, "V_all"], plotScalingName)
        for i in range(nDim + 1):
            outputVar.at[dim[i], "scaling"] = {"alpha": alpha[dim[i]], "K": K[dim[i]], "R2_alpha": R2[dim[i]]}

    # Compute average values per dimension
    for i in range(nDim + 1):
        it = dim[i]
        outputVar.loc[it, "D"] = np.nanmean(outputVar.loc[it, "D_all"])
        outputVar.loc[it, "V"] = np.nanmean(outputVar.loc[it, "V_all"])
        outputVar.loc[it, "T"] = np.nanmean(outputVar.loc[it, "T_all"])
        outputVar.loc[it, "W"] = np.nanmean(outputVar.loc[it, "W_all"])
        outputVar.loc[it, "R2"] = np.nanmean(outputVar.loc[it, "R2_all"])
        outputVar.loc[it, "P"] = np.nanmean(outputVar.loc[it, "P_all"])

    # Populate output with selected fields
    outputFinal = pd.DataFrame(index=dim, columns=outputNames)
    for i in range(len(outputFinal.columns)):
        varName = outputFinal.columns[i]
        outputFinal.loc[:, varName] = outputVar.loc[:, varName]

    return outputFinal


def segment_MED(t, r, v, min_D, min_T, min_V):
    """
    Identifies frame indices delimiting each movement element in a 1-D signal.

    Parameters
    ----------
    t     : array — time vector
    r     : array — 1-D position
    v     : array — 1-D velocity
    min_D : float — minimum displacement threshold [m]
    min_T : float — minimum duration threshold [s]
    min_V : float — velocity threshold for noise classification [m/s]

    Returns
    -------
    ME : (N × 2) array of start and end frame indices, or empty array.
    """
    v = v.ravel()

    # Find peaks
    pk_pos_index, _ = find_peaks(v)
    pk_neg_index, _ = find_peaks(-v)

    # Combine and sort the indices of all critical points
    pk_index = np.sort(np.concatenate([pk_pos_index, pk_neg_index]))

    if pk_index.size == 0:
        return np.array([])

    # Classify peaks based on their velocity relative to the min_V threshold
    pk_velocities = v[pk_index]
    pk_class = np.zeros_like(pk_index, dtype=int)
    pk_class[pk_velocities > min_V] = 1
    pk_class[pk_velocities < -min_V] = -1

    # Early exit if there are too few critical points or no significant peaks
    if len(pk_class) < 3 or np.sum(np.abs(pk_class)) == 0:
        return np.array([])

    # Find points that cross zero velocity between critical points
    buffer_index = []
    buffer_class = []

    # HEAD — before the first peak: find zero-crossing
    if pk_class[0] != 0 and pk_index[0] > 0:
        segment = np.abs(v[0 : pk_index[0]])
        if segment.size > 0:
            loc = np.argmin(segment)
            min_val = segment[loc]
            if min_val <= min_V:
                buffer_index.append(loc)
                buffer_class.append(0)

    # BODY — loop through all peaks
    for i in range(len(pk_class)):
        # Add the current peak to the buffer
        buffer_index.append(pk_index[i])
        buffer_class.append(pk_class[i])

        # If sign changes between adjacent peaks, insert a zero-crossing between them
        if i < len(pk_class) - 1:
            if abs(pk_class[i] - pk_class[i + 1]) == 2:
                start_slice = pk_index[i]
                end_slice = pk_index[i + 1] + 1
                segment = np.abs(v[start_slice : end_slice])
                if segment.size > 0:
                    relative_loc = np.argmin(segment)
                    buffer_index.append(start_slice + relative_loc)
                    buffer_class.append(0)

    # TAIL — after the last peak: find zero-crossing
    if pk_class[-1] != 0 and pk_index[-1] < len(v) - 1:
        start_slice = pk_index[-1] + 1
        segment = np.abs(v[start_slice :])
        if segment.size > 0:
            relative_loc = np.argmin(segment)
            min_val = segment[relative_loc]
            if min_val <= min_V:
                buffer_index.append(start_slice + relative_loc)
                buffer_class.append(0)

    # Finalize variables
    pk_index = np.array(buffer_index)
    pk_class = np.abs(np.array(buffer_class))

    # Creates a vector indicating transitions from zero to non-zero peak (element boundaries)
    seg_class = np.diff(pk_class)
    # seg_i: beginning of each element
    seg_i = pk_index[:-1][seg_class == 1]
    # seg_f: end of each element
    seg_f = pk_index[1:][seg_class == -1]

    # Defensive checks for segments (empty or scalar)
    if seg_i.size == 0 or seg_f.size == 0:
        return np.array([])

    if seg_i.size == 1 and seg_f.size == 1:
        if seg_i[0] >= seg_f[0]:
            return np.array([])

    # If the movement ends with peak velocity with no 0 after the last peak
    if pk_class[-1] == 1:
        seg_i = seg_i[:-1]
    # If the movement starts with peak velocity with no 0 before the first peak
    if pk_class[0] == 1:
        seg_f = seg_f[1:]

    # Handles the case where acceleration causes velocity to jump from <-minV to >+minV in one frame
    if seg_f.size > seg_i.size:
        seg_f = seg_f[: seg_i.size]
    elif seg_f.size < seg_i.size:
        seg_i = seg_i[: seg_f.size]

    # Exclude segments where start equals end
    exclude = (seg_f == seg_i)
    seg_i = seg_i[~exclude]
    seg_f = seg_f[~exclude]

    ME = np.column_stack((seg_i, seg_f))

    if ME.shape[0] == 0:
        return np.array([])

    # Select elements that pass the displacement and duration filters
    valid_mask = np.ones(ME.shape[0], dtype=bool)
    for i in range(ME.shape[0]):
        start_idx, end_idx = ME[i, 0], ME[i, 1]
        # Guard against out-of-bounds indices
        if end_idx >= len(r) or end_idx >= len(t):
            valid_mask[i] = False
            continue
        displacement = np.abs(r[end_idx] - r[start_idx])
        duration = t[end_idx] - t[start_idx]
        if displacement < min_D or duration < min_T:
            valid_mask[i] = False

    return ME[valid_mask]


def analyze_elements_MED(t, r, v, ME):
    """
    Computes kinematic descriptors for each pre-segmented movement element.

    Parameters
    ----------
    t  : array — time vector
    r  : array — 1-D position
    v  : array — 1-D velocity
    ME : (N × 2) array — start and end frame indices of each element

    Returns
    -------
    N, Nt, D, V, T, W, R2, P
    """
    N = ME.shape[0]
    if t[-1] != t[0]:
        Nt = N / (t[-1] - t[0])
    else:
        Nt = 0

    D = np.zeros(N)
    V = np.zeros(N)
    T = np.zeros(N)
    W = np.zeros(N)
    R2 = np.zeros(N)
    P = np.zeros(N)

    # Loop through every movement element
    for i in range(N):
        idx_range = range(ME[i, 0], ME[i, 1] + 1)
        r_i = r[ME[i, 0] : ME[i, 1] + 1]
        v_i = v[idx_range]
        t_i = t[idx_range]

        D[i] = np.abs(r_i[-1] - r_i[0])
        T[i] = t_i[-1] - t_i[0]
        V[i] = np.abs(np.mean(v_i))

        t_interp = np.linspace(0, 1, len(v_i))
        Hoff = V[i] * 30 * ((t_interp**4) - 2*(t_interp**3) + (t_interp**2))

        if V[i] != 0:
            W[i] = np.std(Hoff - v_i, ddof=1) / V[i]
        else:
            W[i] = 0
        if len(v_i) > 1:
            corr = np.corrcoef(Hoff, v_i)[0, 1]
        else:
            corr = 0
        if not np.isnan(corr):
            R2[i] = corr**2
        else:
            R2[i] = 0
        P[i] = len(find_peaks(np.abs(v_i))[0])

    return N, Nt, D, V, T, W, R2, P


def scaling_MED(D_dict, V_dict, plot_name=None):
    """
    Estimates power-law scaling parameters V = K · D^alpha per dimension.

    Parameters
    ----------
    D_dict    : dict — displacement arrays keyed by dimension label
    V_dict    : dict — mean velocity arrays keyed by dimension label
    plot_name : str or None — file path to save log-log scatter plot

    Returns
    -------
    alpha, K, R2_vals : dicts keyed by dimension label
    """
    # Initialize dictionaries for scaling exponent (alpha), coefficient (K), and R^2
    alpha = {}
    K = {}
    R2_vals = {}

    n = len(D_dict)  # Number of dimensions

    # Determine subplot grid size based on number of dimensions
    if n == 1:
        nrows, ncols = 1, 1
    elif n == 3:
        nrows, ncols = 1, 3
    else:
        nrows, ncols = 2, 2

    plt.figure(figsize=(10, 5))

    # Iterate over each dimension
    for i, dim in enumerate(D_dict.keys()):
        # Take logarithm of displacement and velocity arrays
        logD = np.log(D_dict[dim])
        logV = np.log(V_dict[dim])
        # Fit a line to log-log data (power-law: V = K * D^alpha)
        p = np.polyfit(logD, logV, 1)

        alpha[dim] = p[0]                        # Scaling exponent
        K[dim] = np.exp(p[1])                    # Scaling coefficient
        corr = np.corrcoef(logD, logV)[0, 1]     # Correlation coefficient
        R2_vals[dim] = corr**2                   # R^2 value

        # If plot_name is provided, generate and plot the fit
        if plot_name:
            x_fit = np.linspace(np.min(logD), np.max(logD), 100)
            y_fit = np.polyval(p, x_fit)
            plt.subplot(nrows, ncols, i + 1)
            plt.loglog(D_dict[dim], V_dict[dim], 'o')                    # Data points
            plt.loglog(np.exp(x_fit), np.exp(y_fit), '-', linewidth=1.5) # Fitted line
            plt.title(dim)
            plt.xlabel("Displacement (m)")
            plt.ylabel("Mean Velocity (m/s)")

    # Save the plot if a filename is given
    if plot_name:
        plt.suptitle("Scaling of Mean Velocity (m/s) vs Displacement (m)")
        plt.savefig(plot_name)
        plt.close()

    return alpha, K, R2_vals
