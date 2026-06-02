# MED — Movement Element Decomposition Package

Studying human motor control requires metrics that are both sensitive to underlying physiological phenomena and interpretable in clinical and research contexts. The Movement Element Decomposition (MED) method was developed to address this need: grounded in motor-primitive and movement-optimization theories, it segments a continuous kinematic recording into discrete movement elements (MEs) and extracts features that are robust and directly linked to known mechanisms of motor control.

This repository provides a unified, open-source implementation of MED across the three principal scientific computing platforms (MATLAB, Python, and R) in a single package file per language. The toolbox automates low-pass filtering, velocity-based segmentation, and the extraction of three categories of features per movement element: global features (displacement, velocity, duration, movement element count, and rate), shape features (similarity to the minimum-jerk bell profile and peak count), and scaling features (the power-law exponent and coefficient of the velocity–displacement relationship).

---

## Repository Structure

```
MED_pkg/
├── pkg/         ← package source (one file per platform)
│   ├── MED_pkg.m
│   ├── MED_pkg.py
│   └── MED_pkg.R
│
├── main/        ← dataset-specific processing scripts
│   ├── main_adult.m / .py / .R
│   ├── main_children.m / .py / .R
│   └── main_upper.m / .py / .R
│
├── analysis/    ← cross-platform validation scripts
│   ├── comparison_table.m
│   ├── corr_variables.m
│   ├── figure_features.py
│   └── figure_cross_plataform.py
│
├── data/        ← place datasets here (not tracked by git)
│   ├── dataset_Adult_Gait/
│   ├── CP child gait data/td/
│   └── complex-upper-limb-movements-1.0.0/
│
├── output/      ← generated automatically (not tracked by git)
│   ├── adult/
│   ├── children/
│   ├── upper/
│   └── analysis/
│       └── figures/
│
└── README.md
```

---

## Validated Datasets

The toolbox was tested on three publicly available datasets covering different movement modalities and recording technologies.

| Dataset | Modality | Sensor | Folder | Citation |
|---|---|---|---|---|
| Adult Gait | 3-D walking gait | IMU | data/dataset_Adult_Gait/ | Voisard et al., 2025 |
| Child Gait (TD) | 3-D walking gait | Optical tracker | data/CP child gait data/td/ | Meyns, 2022 |
| Upper Limb | 3-D reaching movements | Optical tracker | data/complex-upper-limb-movements-1.0.0/ | Miranda et al., 2018 |

Download each dataset and place it in the corresponding subfolder before running the main scripts.

---

## Dependencies

### MATLAB
- Signal Processing Toolbox (`butter`, `filtfilt`)
- Statistics and Machine Learning Toolbox (`corrcoef`, `findpeaks`)

### Python
```bash
pip install numpy scipy pandas matplotlib seaborn
```

### R
```r
install.packages(c("gsignal", "signal", "pracma", "ggplot2"))
```

---

## Usage

### Step 1 — Place datasets

Download the datasets listed above and place them under `data/` as shown in the repository structure.

### Step 2 — Run a main script

Each main script processes one dataset and writes a CSV to `output/<dataset>/`. All nine scripts (three datasets × three languages) must be run to fully populate `output/` before the analysis scripts can be executed.

**MATLAB** — open MATLAB with the repository root as the working directory:
```matlab
cd path/to/MED_pkg
run('main/main_adult.m')
```

**Python** — run from the repository root:
```bash
cd path/to/MED_pkg
python main/main_adult.py
```

**R** — set the working directory to the repository root before sourcing:
```r
setwd("path/to/MED_pkg")
source("main/main_adult.R")
```

### Step 3 — Run the analysis scripts

With all nine output CSVs in place, generate the cross-platform comparison table first (its output is the input to the bar plot script):

```matlab
run('analysis/comparison_table.m')
```

Then generate all figures:

```matlab
run('analysis/corr_variables.m')
```

```bash
python analysis/figure_features.py
python analysis/figure_cross_plataform.py
```

All figures are written to `output/analysis/figures/`.

---

## API

The main function has the same interface across all three languages. In MATLAB it is named `MED_pkg` (required to match the filename); in Python and R it is `MED`.

### Signature

**MATLAB**
```matlab
output = MED_pkg(movementData, FPS)
output = MED_pkg(movementData, FPS, unit, limits, filter, outputVar, ...
                 scaling_plot_name, t, trim_border)
```

**Python**
```python
output = MED(movementData, FPS)
output = MED(movementData, FPS, unit, limits, filter, outputNames,
             plotScalingName, trim_border=True)
```

**R**
```r
output <- MED(movementData, FPS)
output <- MED(movementData, FPS, unit, limits, filter, outputVar,
              scaling_plot_name, t, trim_border)
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `movementData` | — | N × dim position matrix |
| `FPS` | — | Capture frequency [Hz] |
| `unit` | `'m'` | Length unit of the input data: `'m'` \| `'cm'` \| `'mm'` |
| `limits` | `[0.003, 0.1, 0.01]` | Segmentation thresholds `[min_D (m), min_T (s), min_V (m/s)]` |
| `filter` | `[10, 4]` | Low-pass filter settings `[cutoff Hz, Butterworth order]` |
| `trim_border` | `true` | Discard filter edge-effect samples after filtering |

### Output Fields

Each field is a struct / dict / named list indexed by spatial dimension (`x`, `y`, `z`, `all`):

| Field | Category | Description |
|---|---|---|
| `D` | Global | Mean displacement per ME [m] |
| `V` | Global | Mean velocity per ME [m/s] |
| `T` | Global | Mean duration per ME [s] |
| `N` | Global | Total count of MEs |
| `Nt` | Global | Rate of MEs [ME/s] |
| `W` | Shape | Mean Hoff waveform deviation index |
| `R2` | Shape | Mean R² of ME profile vs. Hoff curve |
| `P` | Shape | Mean velocity peak count per ME |
| `scaling` | Scaling | Power-law fit: `alpha` (exponent), `K` (coefficient), `R2_alpha` (fit R²) |
| `timeSeries` | Internal | Filtered position `r`, velocity `v`, and time vector `t` |
| `ME` | Internal | Start and end frame indices of each detected ME |

---

## Output File Format

Each main script produces a single CSV file with one row per trial. Column names follow the pattern `<variable>_<dimension>`:

```
file, ind, task, trial,
D_all, D_x, D_y, D_z,
V_all, V_x, V_y, V_z,
T_all, N_all, Nt_all, W_all, R2_all, P_all,
alpha_all, K_all, R2_alpha_all,
...  (repeated for _x, _y, _z where applicable)
```

The dimensions present depend on the dataset: the upper-limb dataset produces `_x` and `_y`; the gait datasets produce `_x`, `_y`, and `_z`.

---

## Processing Pipeline

Each call to the main function executes the following steps in sequence:

1. **Unit conversion** — input coordinates are normalized to meters.
2. **Low-pass filtering** — a zero-phase Butterworth filter suppresses high-frequency noise. Edge-effect samples are discarded unless `trim_border` is set to false.
3. **Velocity computation** — instantaneous velocity is obtained by first-order finite differencing of the filtered position.
4. **Segmentation** — velocity peaks and zero-crossings define candidate movement element boundaries per spatial dimension.
5. **Validation** — candidates failing the minimum displacement, duration, or velocity thresholds are discarded.
6. **Feature extraction** — global, shape, and scaling features are computed for each retained element.
7. **Scaling fit** — a power-law V = K · D^α is estimated by log-log linear regression across all elements.

---

## Notes on the Hoff Profile

The shape features W and R² quantify how closely each ME's velocity profile resembles the minimum-jerk bell-shaped curve, defined as:

```
H(t) = V_mean · 30 · (t⁴ − 2t³ + t²),   t ∈ [0, 1]
```

where t is normalized time, and V_mean is the mean velocity of the element.

---

## How to Cite

If you use this package in your research, please cite it as follows:

**Plain Text Citation:**
Silva, M. S., Qadri, L., Vijay, N., Miranda, J. G. V., & Daneault, J.-F. (2026). MED — Movement Element Decomposition Package. Zenodo. https://doi.org/10.5281/zenodo.20510859

**BibTeX:**
```bash
@software{med_package,
author       = {Silva, Mateus Souza and Qadri, Leila and Vijay, Neil and Miranda, José Garcia Vivas and Daneault, Jean-Francois},
title        = {MED — Movement Element Decomposition Package},
year         = {2026},
publisher    = {Zenodo},
doi          = {10.5281/zenodo.20510859},
url          = {https://doi.org/10.5281/zenodo.20510859}
}
```
