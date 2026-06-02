# MED_pkg — Movement Element Decomposition Package
# =============================================
#
# Overview
# --------
# This file implements the Movement Element Decomposition (MED) algorithm, a method
# for decomposing continuous kinematic time-series data into discrete movement
# elements (MEs). Each element is a bounded interval of movement (motor primitive),
# identified by its velocity profile and validated against configurable thresholds
# for displacement, duration, and velocity.
#
# The algorithm supports positional data in one, two, or three spatial dimensions
# and produces per-dimension and aggregate features for each detected ME.
#
# Pipeline Summary
# ----------------
# Raw position data undergoes the following sequential processing steps:
#
#   1. Unit conversion    — input coordinates are normalised to meters.
#   2. Low-pass filtering — a zero-phase Butterworth filter (default: 4th order,
#      10 Hz cutoff) removes high-frequency noise while preserving movement
#      kinematics. Border samples affected by filter edge effects are discarded
#      unless trim_border is set to FALSE.
#   3. Velocity computation — instantaneous velocity is derived via first-order
#      finite differences of the filtered position signal.
#   4. Segmentation (segment_MED) — velocity peaks and zero-crossings are located
#      per dimension to define candidate movement elements.
#   5. Element validation — candidates failing minimum displacement (min_D),
#      duration (min_T), or velocity (min_V) criteria are rejected.
#   6. Feature extraction (analyze_elements_MED) — for each valid element, scalar
#      features are computed: displacement D [m], mean velocity V [m/s], duration
#      T [s], similarity index W relative to the Hoff bell-shaped profile,
#      coefficient of determination R², and peak count P.
#   7. Scaling analysis (scaling_MED) — a power-law relationship V = K · D^alpha
#      is estimated via log-log regression.
#
# Public API
# ----------
# MED(movementData, FPS, ...)
#   Main entry point. Returns a named list of per-dimension movement statistics.
#
# segment_MED(t, r, v, min_D, min_T, min_V)
#   Identifies the frame indices delimiting each movement element in a 1-D
#   velocity signal.
#
# analyze_elements_MED(t, r, v, ME)
#   Computes kinematic features (D, V, T, W, R², P) for an array of
#   pre-segmented movement elements.
#
# scaling_MED(D, V, name=NULL)
#   Estimates power-law scaling exponents (alpha), coefficients (K), and
#   goodness-of-fit (R²) from displacement–velocity pairs per dimension.
#   Optionally saves log-log scatter plots via ggplot2.
#
# Dependencies
# ------------
# gsignal, signal, pracma, ggplot2
#
# Notes
# -----
# - The Hoff velocity profile used for similarity index (W) is the
#   minimum-jerk bell-shaped curve: H(t) = V_mean · 30 · (t⁴ − 2t³ + t²),
#   where t ∈ [0, 1] is normalised time.
# - borderEffect = round(2 · FPS / lp · order / 4) frames are removed from
#   each end of the signal after filtering to suppress edge artefacts. This
#   behaviour can be disabled by passing trim_border = FALSE.
# - Note: the variables N and T shadow base R objects (NULL and TRUE alias).
#   They are used here to maintain naming consistency with the MATLAB and
#   Python implementations.
#
# How to cite
# -----------
# DOI: https://doi.org/10.5281/zenodo.20510859
#
# Usage example
# -------------
# source("MED_pkg.R")
# movementData <- as.matrix(read.csv("your_data.csv"))
# FPS    <- 100
# result <- MED(movementData, FPS)

library(gsignal)
library(signal)
library(pracma)
library(ggplot2)


segment_MED <- function(t, r, v, min_D, min_T, min_V) {
    # Identifies frame indices delimiting each movement element in a 1-D signal.
    #
    # Parameters:
    #   t     — time vector
    #   r     — 1-D position vector
    #   v     — 1-D velocity vector
    #   min_D — minimum displacement threshold [m]
    #   min_T — minimum duration threshold [s]
    #   min_V — velocity threshold for noise classification [m/s]
    #
    # Returns:
    #   Matrix with two columns (start, end frame indices), or NULL.

    # Find peaks
    pos_peaks <- pracma::findpeaks(v)
    pk_pos_index <- pos_peaks[, 2]

    neg_peaks <- pracma::findpeaks(-v)
    pk_neg_index <- neg_peaks[, 2]

    # Combine and sort the indices of all critical points
    pk_index <- sort(c(pk_pos_index, pk_neg_index))

    # Classify peaks based on their velocity relative to the min_V threshold
    pk_class <- rep(NA, length(pk_index))
    pk_class[v[pk_index] >  min_V] <-  1
    pk_class[v[pk_index] < -min_V] <- -1
    pk_class[abs(v[pk_index]) <= min_V] <- 0

    # Early exit if there are too few critical points or no significant peaks
    if (length(pk_class) < 3 || sum(abs(pk_class), na.rm = TRUE) == 0) {
        return(NULL)
    }

    # Find points that cross zero velocity between critical points
    buffer_index <- rep(0, length(pk_index) * 2)
    buffer_class <- rep(0, length(pk_class) * 2)
    j <- 1

    # HEAD — before the first peak: find zero-crossing
    if (pk_class[1] != 0 && pk_index[1] > 1) {
        range_vals <- abs(v[1 : (pk_index[1] - 1)])
        loc <- which.min(range_vals)
        min_val <- range_vals[loc]
        if (min_val <= min_V) {
            buffer_index[j] <- loc
            buffer_class[j] <- 0
            j <- j + 1
        }
    }

    # BODY — loop through all peaks
    for (i in 1 : length(pk_class)) {
        # Add the current peak to the buffer
        buffer_index[j] <- pk_index[i]
        buffer_class[j] <- pk_class[i]
        j <- j + 1

        # If sign changes between adjacent peaks, insert a zero-crossing between them
        if (i < length(pk_class)) {
            if (abs(pk_class[i] - pk_class[i + 1]) == 2 && (pk_index[i + 1] - pk_index[i]) > 1) {
                range_idxs <- pk_index[i] : pk_index[i + 1]
                current_segment <- abs(v[range_idxs])
                loc <- which.min(current_segment)
                buffer_index[j] <- range_idxs[loc]
                buffer_class[j] <- 0
                j <- j + 1
            }
        }
    }

    # TAIL — after the last peak: find zero-crossing
    last_idx <- length(pk_index)
    if (pk_class[last_idx] != 0 && pk_index[last_idx] < length(v)) {
        range_idxs <- (pk_index[last_idx] + 1) : length(v)
        current_segment <- abs(v[range_idxs])
        loc <- which.min(current_segment)
        min_val <- current_segment[loc]
        if (min_val <= min_V) {
            buffer_index[j] <- range_idxs[loc]
            buffer_class[j] <- 0
            j <- j + 1
        }
    }

    # Finalize variables
    pk_class <- abs(buffer_class[1 : (j - 1)])
    pk_index <- buffer_index[1 : (j - 1)]

    # Creates a vector indicating transitions from zero to non-zero peak (element boundaries)
    seg_class <- diff(pk_class)
    # seg_i: beginning of each element
    seg_i <- pk_index[which(seg_class == 1)]
    # seg_f: end of each element
    seg_f <- pk_index[which(seg_class == -1) + 1]

    # Defensive checks for segments (empty or scalar)
    if (length(seg_i) == 0 || length(seg_f) == 0) {
        return(NULL)
    }

    if (length(seg_i) == 1 && length(seg_f) == 1) {
        if (seg_i[1] >= seg_f[1]) {
            return(NULL)
        }
    }

    # If the movement ends with peak velocity with no 0 after the last peak
    if (pk_class[length(pk_class)] == 1) {
        seg_i <- seg_i[1 : (length(seg_i) - 1)]
    }
    # If the movement starts with peak velocity with no 0 before the first peak
    if (pk_class[1] == 1) {
        seg_f <- seg_f[2 : length(seg_f)]
    }

    # Handles the case where acceleration causes velocity to jump from <-min_V to >+min_V in one frame
    if (length(seg_f) > length(seg_i)) {
        seg_f <- seg_f[1 : length(seg_i)]
    } else if (length(seg_f) < length(seg_i)) {
        seg_i <- seg_i[1 : length(seg_f)]
    }

    # Exclude segments where start equals end
    exclude <- (seg_f == seg_i)
    seg_i <- seg_i[!exclude]
    seg_f <- seg_f[!exclude]

    ME <- cbind(seg_i, seg_f)

    # Select elements that pass the displacement and duration filters
    keep <- rep(TRUE, nrow(ME))
    for (i in 1 : nrow(ME)) {
        displacement <- abs(r[ME[i, 2]] - r[ME[i, 1]])
        duration <- t[ME[i, 2]] - t[ME[i, 1]]
        if (displacement < min_D || duration < min_T) {
            keep[i] <- FALSE
        }
    }

    ME <- ME[keep, , drop = FALSE]
    return(ME)
}


analyze_elements_MED <- function(t, r, v, ME) {
    # Computes kinematic features for each pre-segmented movement element.
    #
    # Parameters:
    #   t  — time vector
    #   r  — 1-D position vector
    #   v  — 1-D velocity vector
    #   ME — matrix with two columns (start, end frame indices)
    #
    # Returns a named list: N, Nt, D, V, T, W, R2, peak.

    N <- nrow(ME)
    Nt <- N / (t[length(t)] - t[1])
    D <- numeric(N)
    V <- numeric(N)
    T <- numeric(N)
    W <- numeric(N)
    R2 <- numeric(N)
    peak <- numeric(N)

    # Loop through every movement element
    for (i in 1 : N) {
        v_i <- v[ME[i, 1] : ME[i, 2]]
        t_i <- t[ME[i, 1] : ME[i, 2]]
        D[i] <- abs(r[ME[i, 2]] - r[ME[i, 1]])
        T[i] <- t_i[length(t_i)] - t_i[1]
        V[i] <- abs(mean(v_i))
        t_Hoff <- seq(0, 1, length.out = length(t_i))
        Hoff <- V[i] * 30 * ((t_Hoff^4) - 2*(t_Hoff^3) + (t_Hoff^2))
        dvHoff <- Hoff - v_i
        W[i] <- sd(dvHoff) / abs(V[i])
        R2[i] <- cor(Hoff, v_i)^2
        peak[i] <- nrow(pracma::findpeaks(abs(v_i)))
    }

    return(list(N = N, Nt = Nt, D = D, V = V, T = T, W = W, R2 = R2, peak = peak))
}


scaling_MED <- function(D, V, name = NULL) {
    # Estimates power-law scaling parameters V = K · D^alpha per dimension.
    #
    # Parameters:
    #   D    — named list of displacement vectors
    #   V    — named list of mean velocity vectors
    #   name — file path prefix for ggplot2 PNG outputs, or NULL to skip plots
    #
    # Returns: named list(alpha, K, R2), or invisibly returns plot list if name
    #          is provided.

    # Number of dimensions
    dim <- names(D)

    # Initialize lists for scaling exponent (alpha), coefficient (K), and R^2
    alpha <- list()
    K <- list()
    R2 <- list()

    if (!is.null(name)) {
        plot_list <- list()
    } else {
        plot_list <- NULL
    }

    # Iterate over each dimension
    for (i in seq_along(dim)) {
        it <- dim[i]
        # Take logarithm of displacement and velocity arrays
        logD <- log(D[[it]])
        logV <- log(V[[it]])
        # Fit a line to log-log data (power-law: V = K * D^alpha)
        fit <- lm(logV ~ logD)
        p <- coef(fit)

        alpha[[it]] <- unname(p[2])         # Scaling exponent
        K[[it]] <- exp(unname(p[1]))        # Scaling coefficient
        R2[[it]] <- summary(fit)$r.squared  # R^2 value

        # If name is provided, generate and plot the fit
        if (!is.null(name)) {
            df <- data.frame(D = D[[it]], V = V[[it]])
            x_fit <- seq(min(logD), max(logD), length.out = 100)
            y_fit <- p[1] + p[2] * x_fit
            fit_df <- data.frame(D = exp(x_fit), V = exp(y_fit))

            plt <- ggplot(df, aes(x = D, y = V)) +
                geom_point(color = "#0072BD", size = 2) +
                geom_line(data = fit_df, aes(x = D, y = V), color = "firebrick", linewidth = 1.2) +
                scale_x_log10() +
                scale_y_log10() +
                labs(title = paste("Scaling:", it),
                     x = "Displacement (m)",
                     y = "Mean Velocity (m/s)") +
                theme_minimal()

            plot_list[[it]] <- plt
            # Save the plot if a filename is given
            ggsave(filename = paste0(name, "_", it, ".png"), plot = plt, width = 6, height = 4)
        }
    }

    if (!is.null(name)) {
        return(invisible(plot_list))
    }
    return(list(alpha = alpha, K = K, R2 = R2))
}


MED <- function(movementData, FPS, unit = 'm', limits = c(0.003, 0.1, 0.01),
                filter = c(10, 4), outputVar = NULL, scaling_plot_name = NULL,
                t = NULL, trim_border = TRUE) {
    # Main entry point of the MED algorithm.
    #
    # Parameters:
    #   movementData      — (N × dim) position matrix [m, cm, or mm]
    #   FPS               — capture frequency in frames per second
    #   unit              — 'm' (default) | 'cm' | 'mm'
    #   limits            — c(min_D, min_T, min_V); defaults c(0.003, 0.1, 0.01)
    #   filter            — c(lp, order); defaults c(10, 4)
    #   outputVar         — character vector of requested output fields
    #   scaling_plot_name — file path prefix for scaling PNG plots, or NULL
    #   t                 — custom time vector; if NULL, derived from FPS
    #   trim_border       — logical, default TRUE; set FALSE to retain the
    #                       filter edge-effect samples at both ends of the signal

    if (is.null(outputVar)) {
        outputVar <- c("scaling", "D", "V", "T", "N", "Nt", "W", "R2", "P",
                       "D_all", "V_all", "T_all", "W_all", "R2_all", "P_all", "timeSeries", "ME")
    }

    # Conditional for if args is >= 3, setting units
    if (is.null(unit) || unit == "") {
        unit <- "m"
    }

    # Conditional for if args is >= 4, setting limits for identifying valid movement elements
    # min_D is the minimum displacement [m], min_T is the minimum duration [s], min_V is the minimum velocity [m/s]
    if (is.null(limits) || length(limits) != 3) {
        limits <- c(0.003, 0.1, 0.01)
    }

    # Conditional for if args is >= 5; lp is low-pass cutoff [Hz], order is filter order
    if (is.null(filter) || length(filter) != 2) {
        filter <- c(10, 4)
    }

    min_D <- limits[1]; min_T <- limits[2]; min_V <- limits[3]
    lp <- filter[1]; order <- filter[2]

    # Convert units to meters if necessary
    if (unit == "mm") {
        movementData <- movementData / 1000
    }
    if (unit == "cm") {
        movementData <- movementData / 100
    }

    # Apply Butterworth low-pass filter to position data
    bf <- gsignal::butter(order, (2 * lp) / FPS)
    movementData <- apply(movementData, 2, function(col) gsignal::filtfilt(bf$b, bf$a, col))

    borderEffect <- round(2 * FPS / lp * order / 4)
    if (trim_border) {
        movementData <- movementData[(borderEffect + 1) : (nrow(movementData) - borderEffect), , drop = FALSE]
    }

    # Calculate velocity and time
    v <- apply(movementData, 2, function(col) diff(col) * FPS)
    movementData <- movementData[1 : (nrow(movementData) - 1), , drop = FALSE]

    if (is.null(t)) {
        t <- seq(1, nrow(v)) / FPS
    }

    timeSeries <- list(r = movementData, v = v, t = t)

    # Determine dimensionality and prepare labels
    nDim <- ncol(movementData)

    if (nDim == 1) {
        dimLabels <- "all"
    } else {
        dimLabels <- c("x", "y", "z")[1 : nDim]
    }

    keepMask <- rep(FALSE, nDim)
    ME <- list()

    # Find frames delimiting valid movement elements per dimension
    for (i in 1 : nDim) {
        lbl <- dimLabels[i]
        currentME <- segment_MED(timeSeries$t, timeSeries$r[, i], timeSeries$v[, i], min_D, min_T, min_V)
        ME[[lbl]] <- currentME
        if (!is.null(currentME) && nrow(currentME) > 0) {
            keepMask[i] <- TRUE
        }
    }

    timeSeries$r <- timeSeries$r[, keepMask, drop = FALSE]
    timeSeries$v <- timeSeries$v[, keepMask, drop = FALSE]

    validLabels <- dimLabels[keepMask]
    nValid <- length(validLabels)

    if (nValid == 0) {
        return(NULL)
    }

    N <- list(); Nt <- list()
    D_all <- list(); V_all <- list(); T_all <- list()
    W_all <- list(); R2_all <- list(); P_all <- list()

    # Analyze movement elements per dimension
    for (i in 1 : nValid) {
        lbl <- validLabels[i]
        res <- analyze_elements_MED(timeSeries$t, timeSeries$r[, i], timeSeries$v[, i], ME[[lbl]])
        N[[lbl]] <- res$N
        Nt[[lbl]] <- res$Nt
        D_all[[lbl]] <- res$D
        V_all[[lbl]] <- res$V
        T_all[[lbl]] <- res$T
        W_all[[lbl]] <- res$W
        R2_all[[lbl]] <- res$R2
        P_all[[lbl]] <- res$peak
    }

    # Initialize combined results across all dimensions
    if (nDim > 1 && nValid > 0) {
        combine_dims <- function(dataList) unlist(dataList[validLabels])
        D_all$all <- combine_dims(D_all)
        V_all$all <- combine_dims(V_all)
        T_all$all <- combine_dims(T_all)
        W_all$all <- combine_dims(W_all)
        R2_all$all <- combine_dims(R2_all)
        P_all$all <- combine_dims(P_all)
        N$all <- sum(unlist(N[validLabels]))
        # Weighted average movement element rate
        Nt$all <- N$all / (t[length(t)] - t[1])
    }

    if (!"all" %in% dimLabels) {
        dimLabels <- c(dimLabels, "all")
    }

    # Compute scaling law from displacement and velocity
    if ("scaling" %in% outputVar) {
        scaling <- scaling_MED(D_all, V_all, scaling_plot_name)
    }

    D <- list(); V <- list(); T <- list(); W <- list(); R2 <- list(); P <- list()

    # Compute average values per dimension
    for (it in dimLabels) {
        D[[it]] <- mean(D_all[[it]])
        V[[it]] <- mean(V_all[[it]])
        T[[it]] <- mean(T_all[[it]])
        W[[it]] <- mean(W_all[[it]])
        R2[[it]] <- mean(R2_all[[it]])
        P[[it]] <- mean(P_all[[it]])
    }

    # Populate output with selected fields
    output <- list()
    for (varName in outputVar) {
        if (exists(varName, inherits = FALSE)) {
            output[[varName]] <- get(varName, inherits = FALSE)
        } else {
            output[[varName]] <- NULL
        }
    }
    return(output)
}
