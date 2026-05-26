import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
from sklearn.metrics import pairwise_distances
from mpl_toolkits.mplot3d import Axes3D
from skfda.representation.grid import FDataGrid
from skfda.preprocessing.dim_reduction.feature_extraction import FPCA
from matplotlib.colors import to_rgba
from skfda.preprocessing.smoothing import BasisSmoother
from skfda.representation.basis import BSplineBasis
from scipy.stats import norm
from scipy.stats import f_oneway
from scipy.stats import ttest_ind
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from itertools import combinations
from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm
from patsy import dmatrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_curve,
    auc,
    confusion_matrix,
    ConfusionMatrixDisplay
)
from sklearn.discriminant_analysis import (
    LinearDiscriminantAnalysis,
    QuadraticDiscriminantAnalysis
)
from sklearn.tree import DecisionTreeClassifier

# ------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------

n_per_group = 30

diagnostics = pd.read_excel("Diagnostics.xlsx")

selected_ids = (
    diagnostics[
        diagnostics["Rhythm"].isin(["AFIB", "SB", "SR"])
    ]
    .sort_values(["Rhythm", "FileName"])
    .groupby("Rhythm")
    .head(n_per_group)["FileName"]
    .tolist()
)

len(selected_ids)
selected_ids[:5]

files = glob.glob("ECGData/*.csv")

files_sample = [
    f for f in files
    if os.path.splitext(os.path.basename(f))[0] in selected_ids
]


def read_file(file):
    df = pd.read_csv(file)

    if "II" in df.columns:
        column_to_leave = df["II"].values
    else:
        column_to_leave = df.iloc[:, 1].values

    file_name = os.path.basename(file)
    file_id = os.path.splitext(file_name)[0]

    return np.concatenate(([file_id], column_to_leave))


result_list = [read_file(f) for f in files_sample]

df = pd.DataFrame(result_list)
df = df.rename(columns={0: "ID"})

kovariantes = pd.read_excel("Diagnostics.xlsx")

kovariantes_sample = kovariantes[
    kovariantes["FileName"].isin(df["ID"])
]

df.to_csv("ECGData/ECG_combined_true.csv", index=False)
df = pd.read_csv("ECGData/ECG_combined_true.csv")

# ------------------------------------------------------------------------
# 2. FUNCTIONAL DATA OBJECT
# ------------------------------------------------------------------------

time = np.linspace(0, 10, 5000)

df_be_id = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
df_scaled = df_be_id / 1000

data_matrix = df_scaled.to_numpy(dtype=float)

fdata_obj_scaled = data_matrix

plt.figure()

for i in range(min(10, data_matrix.shape[0])):
    plt.plot(time, data_matrix[i])

plt.title("ECG Lead II — raw signals")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (mV)")
plt.show()

# ------------------------------------------------------------------------
# 3. SMOOTHING PARAMETER SELECTION — GCV curve
# ------------------------------------------------------------------------

Y = fdata_obj_scaled
t = time

n, m = Y.shape

# B-spline basis
nbasis = 500
degree = 3

knots = np.linspace(t.min(), t.max(), nbasis - 2)

knots_full = np.concatenate((
    np.repeat(t.min(), degree),
    knots,
    np.repeat(t.max(), degree)
))

def bspline_basis_matrix(x):
    B = np.zeros((len(x), nbasis))
    for i in range(nbasis):
        c = np.zeros(nbasis)
        c[i] = 1
        B[:, i] = BSpline(knots_full, c, degree)(x)
    return B

B = bspline_basis_matrix(t)

lambda_fine = 10 ** np.arange(-8, -1.5, 0.5)

results = []

I = np.eye(nbasis)

for lam in lambda_fine:

    BtB = B.T @ B
    BtY = B.T @ Y.T

    coef = np.linalg.solve(BtB + lam * I, BtY)

    Yhat = (B @ coef).T

    resid = Y - Yhat
    rss = np.sum(resid ** 2)

    # smoother matrix
    S = B @ np.linalg.inv(BtB + lam * I) @ B.T
    df = np.trace(S)

    gcv = rss / ((n * m - df) ** 2)

    results.append((np.log10(lam), gcv))

cv = pd.DataFrame(results, columns=["log_lambda", "GCV"])

plt.figure(figsize=(7,4))
plt.plot(cv["log_lambda"], cv["GCV"], linewidth=2)
plt.axvline(-4, color="red", linestyle="--")

plt.xlabel("log10(lambda)")
plt.ylabel("GCV")
plt.title("GCV vs smoothing parameter (R-like)")

plt.show()

# ------------------------------------------------------------------------
# 4. FINAL SMOOTHED OBJECT
#    nbasis=500, lambda=1e-04, penalty order=1 (penalise 1st derivative)
# ------------------------------------------------------------------------

# -----------------------------
# INPUT
# -----------------------------
Y = fdata_obj_scaled
t = time

n, m = Y.shape

# -----------------------------
# BASIS (create.bspline.basis)
# -----------------------------
nbasis = 500
degree = 3

knots = np.linspace(t.min(), t.max(), nbasis - degree - 1)

knots_full = np.concatenate((
    np.repeat(t.min(), degree),
    knots,
    np.repeat(t.max(), degree)
))

# -----------------------------
# B-spline design matrix
# -----------------------------
def bspline_design(x):
    B = np.zeros((len(x), nbasis))
    for i in range(nbasis):
        c = np.zeros(nbasis)
        c[i] = 1
        B[:, i] = BSpline(knots_full, c, degree)(x)
    return B

B = bspline_design(t)

# -----------------------------
# fdPar equivalent: lambda = 1e-4, Lfd = 1
# (first derivative penalty → roughness penalty approximation)
# -----------------------------
lam = 1e-4
I = np.eye(nbasis)

BtB = B.T @ B
BtY = B.T @ Y.T

coef = np.linalg.solve(BtB + lam * I, BtY)

# -----------------------------
# smoothed.fd
# -----------------------------
Y_smooth = (B @ coef).T

fd_final = Y_smooth

# ------------------------------------------------------------------------
# 5. VALIDATION PLOTS
# ------------------------------------------------------------------------

# --- Plot A: Raw vs smoothed for 3 patients ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, i in zip(axes, [0, 44, 79]):
    ax.plot(time, fdata_obj_scaled[i], color="black", linewidth=1)
    ax.plot(time, fd_final[i], color="red", linewidth=2)
    ax.set_title(f"Patient {i+1}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")
    ax.legend(["Raw", "Smoothed"])

plt.tight_layout()
plt.show()

# --- Single patient ---
plt.figure()
plt.plot(time, fdata_obj_scaled[0], color="black", linewidth=1)
plt.plot(time, fd_final[0], color="red", linewidth=2)
plt.title("Patient 1")
plt.xlabel("Time (s)")
plt.ylabel("mV")
plt.legend(["Raw", "Smoothed"])
plt.show()


# --- Plot C: Lambda comparison ---

def smooth_lambda(lam):
    BtB = B.T @ B
    BtY = B.T @ fdata_obj_scaled.T
    coef = np.linalg.solve(BtB + lam * np.eye(B.shape[1]), BtY)
    return (B @ coef).T


fig, axes = plt.subplots(1, 3, figsize=(15, 4))

lambdas = [1e-8, 1e-4, 1]
titles = ["lambda = 1e-08", "lambda = 1e-04", "lambda = 1"]
colors = ["red", "green", "blue"]

for ax, lam, title, col in zip(axes, lambdas, titles, colors):
    sm = smooth_lambda(lam)
    i = 0

    ax.plot(time, fdata_obj_scaled[i], color="grey")
    ax.plot(time, sm[i], color=col, linewidth=2)
    ax.set_title(title)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")

plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------
# OUTLIERS
# ------------------------------------------------------------------------

# -----------------------------
# 1. EVALUATE SMOOTHED CURVES
# -----------------------------
tmp = fd_final
tmp = np.asarray(tmp)

# ensure correct orientation
if tmp.shape[0] != data_matrix.shape[0]:
    tmp = tmp.T

# --------------------------------------------------------
# 2. BAND DEPTH (Rough Python equivalent via centrality)
# --------------------------------------------------------
# approximation: mean distance from all curves
dist = pairwise_distances(tmp)
bd = -np.mean(dist, axis=1)

plt.figure()
plt.plot(bd)
plt.title("Band Depth (approximated)")
plt.xlabel("Curve index")
plt.ylabel("Depth (higher = more central)")
plt.show()

# --------------------------------------------------------
# 3. MODIFIED BAND DEPTH (robust version)
# --------------------------------------------------------
# approximation: median-based depth
mbd = -np.median(dist, axis=1)

plt.figure()
plt.plot(mbd)
plt.title("Modified Band Depth (approximated)")
plt.xlabel("Curve index")
plt.ylabel("MBD Value")
plt.show()

# --------------------------------------------------------
# 4. FUNCTIONAL BOXPLOT (MBD STYLE)
# --------------------------------------------------------
center_idx = np.argmax(mbd)
center_curve = tmp[center_idx]

dist_to_center = np.linalg.norm(tmp - center_curve, axis=1)
threshold = np.percentile(dist_to_center, 95)

outliers = np.where(dist_to_center > threshold)[0]

plt.figure()
for i in range(min(20, len(tmp))):
    plt.plot(time, tmp[i], color="lightgrey")

for i in outliers:
    plt.plot(time, tmp[i], color="red")

plt.plot(time, center_curve, color="black", linewidth=2)

plt.title("Functional Boxplot (approx MBD)")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (mV)")
plt.show()

# --------------------------------------------------------
# 5. MUOD (approximation: amplitude + shape + magnitude)
# --------------------------------------------------------

mean_curve = np.mean(tmp, axis=0)

# amplitude deviation
amp_dev = np.abs(np.max(tmp, axis=1) - np.max(mean_curve))

# shape deviation
shape_dev = np.linalg.norm(tmp - mean_curve, axis=1)

# magnitude deviation
mag_dev = np.linalg.norm(tmp, axis=1)

out_amp = np.where(amp_dev > np.percentile(amp_dev, 95))[0]
out_shape = np.where(shape_dev > np.percentile(shape_dev, 95))[0]
out_mag = np.where(mag_dev > np.percentile(mag_dev, 95))[0]

print("MUOD OUTLIERS")
print("Amplitude:", out_amp)
print("Shape:", out_shape)
print("Magnitude:", out_mag)

#----------------------------------------
# EDA
#----------------------------------------

# SECTION 1 — mean + SD bands

mean_beat = np.mean(fd_final, axis=0)
stddev_beat = np.std(fd_final, axis=0)

plt.figure(figsize=(10, 5))

# all curves
for i in range(fd_final.shape[0]):
    plt.plot(time, fd_final[i], color="black", alpha=0.15, linewidth=0.8)

# mean and bands
plt.plot(time, mean_beat, color="red", linewidth=3, label="Mean")
plt.plot(time, mean_beat + 3 * stddev_beat, color="blue", linewidth=2, label="Mean ± 3xSD")
plt.plot(time, mean_beat - 3 * stddev_beat, color="blue", linewidth=2)

plt.axhline(0, color="grey", linestyle="--", linewidth=1)

plt.title("Smoothed ECG signals")
plt.xlabel("Time (s)")
plt.ylabel("mV")

plt.legend(frameon=True)
plt.show()

# SECTION 2 — Covariance function (Python equivalent)

# compute covariance over time points
ecgvar_mat = np.cov(fd_final, rowvar=False)

# subsample
t_sub = time[::5]
idx = np.arange(0, len(time), 5)
ecgvar_sub = ecgvar_mat[np.ix_(idx, idx)]

# -------------------------
# 3D surface plot
# -------------------------
X, Y = np.meshgrid(t_sub, t_sub)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, ecgvar_sub, cmap="viridis")

ax.set_title("Covariance surface — ECG beats")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Time (s)")
ax.set_zlabel("Covariance (mV²)")

plt.show()

# -------------------------
# contour plot
# -------------------------
plt.figure(figsize=(6, 5))
plt.contour(t_sub, t_sub, ecgvar_sub)
plt.title("Covariance contour — ECG beats")
plt.xlabel("Time (s)")
plt.ylabel("Time (s)")
plt.show()

# -------------------------
# heatmap
# -------------------------
plt.figure(figsize=(6, 5))
plt.imshow(
    ecgvar_sub,
    extent=[t_sub.min(), t_sub.max(), t_sub.min(), t_sub.max()],
    origin="lower",
    aspect="auto"
)

plt.title("Covariance heatmap — ECG beats")
plt.xlabel("Time (s)")
plt.ylabel("Time (s)")
plt.colorbar(label="Covariance")

plt.show()

import numpy as np
import matplotlib.pyplot as plt

# SECTION 3 — Derivatives (numerical approximation)

dt = time[1] - time[0]

# first and second derivatives
deriv1 = np.gradient(fd_final, dt, axis=1)
deriv2 = np.gradient(deriv1, dt, axis=1)

mean_d1 = np.mean(deriv1, axis=0)
mean_d2 = np.mean(deriv2, axis=0)

# -------------------------
# Plot derivatives
# -------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# First derivative
for i in range(fd_final.shape[0]):
    axes[0].plot(time, deriv1[i], color="black", alpha=0.15, linewidth=0.8)

axes[0].plot(time, mean_d1, color="red", linewidth=3, linestyle="--")
axes[0].axhline(0, color="grey", linestyle="--", linewidth=1)

axes[0].set_title("First derivative (slope)")
axes[0].set_xlabel("Time (s)")
axes[0].set_ylabel("mV/s")

axes[0].legend(["Mean slope"], frameon=True)

# Second derivative
for i in range(fd_final.shape[0]):
    axes[1].plot(time, deriv2[i], color="black", alpha=0.15, linewidth=0.8)

axes[1].plot(time, mean_d2, color="blue", linewidth=3, linestyle="--")
axes[1].axhline(0, color="grey", linestyle="--", linewidth=1)

axes[1].set_title("Second derivative (acceleration)")
axes[1].set_xlabel("Time (s)")
axes[1].set_ylabel("mV/s²")

axes[1].legend(["Mean acceleration"], frameon=True)

plt.tight_layout()
plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# SECTION 4 — Centrality measures by rhythm group

rhythms_v = kovariantes_sample["Rhythm"].values

# top groups
rhythm_tab = pd.Series(rhythms_v).value_counts()
top_groups = rhythm_tab.index[:min(4, len(rhythm_tab))]

print("Top rhythm groups:", list(top_groups))
print(rhythm_tab.loc[top_groups])

# ------------------------------------------------------------
# functional data container
# ------------------------------------------------------------
eval_grid_full = np.linspace(0, 10, fd_final.shape[1])

full_mat = fd_final.copy()

# group extractor
def get_group(r):
    return full_mat[rhythms_v == r]

def trimmed_mean(x, trim=0.15):
    low = int(trim * x.shape[0])
    high = x.shape[0] - low
    return np.mean(np.sort(x, axis=0)[low:high], axis=0)

# ------------------------------------------------------------
# plot centrality
# ------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes = axes.flatten()

for ax, r in zip(axes, top_groups):
    grp = get_group(r)
    n_r = grp.shape[0]

    mean_curve = np.mean(grp, axis=0)
    trim_curve = trimmed_mean(grp, 0.15)
    median_curve = np.median(grp, axis=0)

    ax.plot(eval_grid_full, mean_curve, label="mean")
    ax.plot(eval_grid_full, trim_curve, label="trimmed mean")
    ax.plot(eval_grid_full, median_curve, label="median")

    ax.set_title(f"Centrality — {r} (n={n_r})")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")
    ax.legend()

plt.tight_layout()
plt.show()

# SECTION 5 — Dispersion measures by rhythm group

rhythms_v = kovariantes_sample["Rhythm"].values
top_groups = pd.Series(rhythms_v).value_counts().index[:4]

def get_group(r):
    return fd_final[rhythms_v == r]

def trimmed_variance(x, trim=0.15):
    n = x.shape[0]
    k = int(trim * n)
    x_sorted = np.sort(x, axis=0)
    return np.var(x_sorted[k:n-k], axis=0)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for ax, r in zip(axes, top_groups):
    grp = get_group(r)
    n_r = grp.shape[0]

    var_curve = np.var(grp, axis=0)
    trim_var = trimmed_variance(grp, 0.15)

    ax.plot(time, var_curve, label="var")
    ax.plot(time, trim_var, label="trim var")

    ax.set_title(f"Dispersion — {r} (n={n_r})")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Variance (mV²)")
    ax.legend()

plt.tight_layout()
plt.show()

# SECTION 6 — Group mean curves overlaid

rhythms_v = kovariantes_sample["Rhythm"].values
unique_r = np.unique(rhythms_v)

# color map
cmap = plt.get_cmap("rainbow", len(unique_r))
cols_map = {r: cmap(i) for i, r in enumerate(unique_r)}

mean_beat = np.mean(fd_final, axis=0)

plt.figure(figsize=(10, 6))

# overall mean
plt.plot(time, mean_beat, color="black", linewidth=2, linestyle="--", label="Overall")

# group means
for r in unique_r:
    idx = np.where(rhythms_v == r)[0]
    if len(idx) >= 2:
        grp_mean = np.mean(fd_final[idx], axis=0)
        plt.plot(time, grp_mean, color=cols_map[r], linewidth=2, label=str(r))

plt.axhline(0, color="grey", linestyle="--", linewidth=1)

plt.title("Mean ECG beat by rhythm group")
plt.xlabel("Time (s)")
plt.ylabel("mV")
plt.legend()
plt.show()

# ------------------------------------------------------------------------
# 7. FUNCTIONAL PCA (CLEAN)
# ------------------------------------------------------------------------

fd = FDataGrid(
    data_matrix,
    grid_points=time
)

fpca = FPCA(n_components=5)
scores = fpca.fit_transform(fd)

harmonics = fpca.components_
mean_fn = fpca.mean_

var_explained = fpca.explained_variance_ratio_
print("Variance explained:", var_explained)

# ------------------------------------------------------------------------
# 1. Scree plot
# ------------------------------------------------------------------------

plt.figure()
plt.plot(range(1, len(var_explained) + 1), var_explained, marker="o")
plt.title("Scree Plot — FPCA")
plt.xlabel("PC")
plt.ylabel("Variance explained")
plt.show()

# ------------------------------------------------------------------------
# 2. Harmonics (FIXED SHAPE)
# ------------------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for i, ax in enumerate(axes):
    comp = harmonics[i].data_matrix.squeeze()
    ax.plot(time, comp)
    ax.set_title(f"PC{i+1}")
    ax.axhline(0, color="black", linewidth=0.5)

plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------
# 3. Mean ± variation (clean R-style)
# ------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

mean_vals = mean_fn.data_matrix.squeeze()

for i in range(2):
    comp = harmonics[i].data_matrix.squeeze()
    sd = np.sqrt(fpca.explained_variance_[i])

    ax = axes[i]
    ax.plot(time, mean_vals, color="black", linewidth=2)
    ax.plot(time, mean_vals + 2 * sd * comp, "b--")
    ax.plot(time, mean_vals - 2 * sd * comp, "r--")

    ax.set_title(f"PC{i+1} variation")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")

plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------
# 4. SCORE PLOT (FIXED COLORS)
# ------------------------------------------------------------------------

rhythms_v = np.array(rhythms_v)
unique = np.unique(rhythms_v)
color_map = {k: plt.cm.tab10(i % 10) for i, k in enumerate(unique)}

plt.figure()

for i in range(scores.shape[0]):
    r = rhythms_v[i]
    plt.scatter(
        scores[i, 0],
        scores[i, 1],
        color=color_map[r],
        s=20
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("FPCA Scores by Rhythm")
plt.show()

# ------------------------------------------------------------------------
# 5. VARIMAX
# ------------------------------------------------------------------------

def varimax(Phi, gamma=1.0, q=20, tol=1e-6):
    p, k = Phi.shape
    R = np.eye(k)
    d = 0

    for _ in range(q):
        d_old = d
        Lambda = Phi @ R
        u, s, vh = np.linalg.svd(
            Phi.T @ (Lambda ** 3 - (gamma / p) * Lambda @ np.diag(np.diag(Lambda.T @ Lambda)))
        )
        R = u @ vh
        d = np.sum(s)
        if d_old != 0 and d / d_old < 1 + tol:
            break

    return Phi @ R

# rotated scores
scores_rot = varimax(scores)

# use same harmonics
harm_rot = harmonics

# ------------------------------------------------------------------------
# 6. CLEAN 4-PANEL ROTATED PLOTS (READABLE)
# ------------------------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

for i, ax in enumerate(axes.flat):

    comp = harm_rot[i].data_matrix.squeeze()
    sd = np.sqrt(fpca.explained_variance_[i])

    ax.plot(time, mean_vals, color="black", linewidth=2)

    ax.plot(time, mean_vals + 2 * sd * comp, "r--")
    ax.plot(time, mean_vals - 2 * sd * comp, "r--")

    ax.set_title(f"Rotated PC{i+1}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")

plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------
# REGISTRATION
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Evaluate smoothed signals on grid
# ------------------------------------------------------------------------

eval_grid = time
smooth_mat = fd_final

n_patients = smooth_mat.shape[0]

# ------------------------------------------------------------------------
# QRS detection
# ------------------------------------------------------------------------

search_idx = np.where(eval_grid <= 1)[0]

qrs_times = np.zeros(n_patients)

for i in range(n_patients):
    seg = smooth_mat[i, search_idx]
    max_idx = np.argmax(seg)
    qrs_times[i] = eval_grid[search_idx[max_idx]]

# ------------------------------------------------------------------------
# Beat window definition
# ------------------------------------------------------------------------

window_before = 0.05
window_after = 0.40

beat_grid = np.linspace(0, window_before + window_after, 250)

# valid beats
valid = (qrs_times >= window_before) & (qrs_times <= (10 - window_after))

print("Valid patients:", np.sum(valid), "/", n_patients)

valid_idx = np.where(valid)[0]
n_valid = len(valid_idx)

# ------------------------------------------------------------------------
# Build aligned beat matrix
# ------------------------------------------------------------------------

beat_matrix = np.zeros((n_valid, len(beat_grid)))

for j, i in enumerate(valid_idx):

    t_window = np.linspace(
        qrs_times[i] - window_before,
        qrs_times[i] + window_after,
        len(beat_grid)
    )

    beat_matrix[j, :] = np.interp(t_window, eval_grid, smooth_mat[i, :])

# ------------------------------------------------------------------------
# Covariates alignment
# ------------------------------------------------------------------------

df_valid = kovariantes.iloc[valid_idx].copy()
rhythms_v = kovariantes.iloc[valid_idx]["Rhythm"].values

unique_rhythms = np.unique(rhythms_v)
cols_map = {r: plt.cm.tab10(k % 10) for k, r in enumerate(unique_rhythms)}

# ------------------------------------------------------------------------
# Plot: aligned beats
# ------------------------------------------------------------------------

plt.figure(figsize=(10, 4))

for j in range(n_valid):
    r = rhythms_v[j]
    plt.plot(
        beat_grid,
        beat_matrix[j],
        color=to_rgba(cols_map[r], 0.35),
        linewidth=1
    )

plt.title("Single beats aligned at QRS (by rhythm)")
plt.xlabel("Time relative to QRS (s)")
plt.ylabel("mV")
plt.ylim(-0.4, 1.2)
plt.show()

# ------------------------------------------------------------------------
# Functional representation
# ------------------------------------------------------------------------

basis = BSplineBasis(n_basis=100, domain_range=(0, 0.45))

fd_beats = FDataGrid(
    data_matrix=beat_matrix,
    grid_points=beat_grid
)

plt.figure()
plt.plot(beat_grid, beat_matrix.mean(axis=0), linewidth=2)
plt.title("Overall mean beat")
plt.xlabel("Time (s)")
plt.ylabel("mV")
plt.show()

# ------------------------------------------------------------------------
# Subsample plots (raw vs smoothed vs registered)
# ------------------------------------------------------------------------

plot_idx = np.arange(0, 5000, 5)
plot_time = time[plot_idx]

plot_raw = smooth_mat[:, plot_idx]

plt.figure(figsize=(15, 4))

plt.subplot(1, 3, 1)
for i in range(min(50, plot_raw.shape[0])):
    plt.plot(plot_time, plot_raw[i], alpha=0.3)
plt.title("Raw ECG")

plt.subplot(1, 3, 2)
for i in range(min(50, smooth_mat.shape[0])):
    plt.plot(plot_time, smooth_mat[i, plot_idx], alpha=0.3)
plt.title("Smoothed ECG")

plt.subplot(1, 3, 3)
for i in range(min(n_valid, 50)):
    plt.plot(beat_grid, beat_matrix[i], alpha=0.3)
plt.title("Registered beats")

plt.tight_layout()
plt.show()

###################
# EDA on registered beats
##################

# mean + sd (across curves)
mean_reg = fd_beats.data_matrix.mean(axis=0).squeeze()
std_reg  = fd_beats.data_matrix.std(axis=0).squeeze()

plt.figure(figsize=(10, 4))

# all curves
for i in range(fd_beats.data_matrix.shape[0]):
    plt.plot(beat_grid, fd_beats.data_matrix[i].squeeze(),
             color="black", alpha=0.15, linewidth=0.8)

# mean + envelopes
plt.plot(beat_grid, mean_reg, color="red", linewidth=3)
plt.plot(beat_grid, mean_reg + 3*std_reg, color="blue", linewidth=2)
plt.plot(beat_grid, mean_reg - 3*std_reg, color="blue", linewidth=2)

plt.axhline(0, color="grey", linestyle="--", linewidth=1)

plt.title("Registered ECG beats")
plt.xlabel("Time (s)")
plt.ylabel("mV")

plt.show()

X = fd_beats.data_matrix.squeeze()

ecgvar_mat = np.cov(X, rowvar=False)

t_sub = beat_grid[::2]
idx = np.arange(0, len(beat_grid), 2)
ecgvar_sub = ecgvar_mat[np.ix_(idx, idx)]

Xg, Yg = np.meshgrid(t_sub, t_sub)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(Xg, Yg, ecgvar_sub, cmap="viridis")

ax.set_title("Covariance surface — registered beats")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Time (s)")
ax.set_zlabel("Covariance (mV²)")
plt.show()

# contour
plt.figure()
plt.contour(t_sub, t_sub, ecgvar_sub)
plt.title("Covariance contour — registered beats")
plt.show()

# heatmap
plt.figure()
plt.imshow(ecgvar_sub, origin="lower", aspect="auto")
plt.title("Covariance heatmap — registered beats")
plt.colorbar()
plt.show()

dt = beat_grid[1] - beat_grid[0]

deriv1 = np.gradient(fd_beats.data_matrix, dt, axis=1)
deriv2 = np.gradient(deriv1, dt, axis=1)

mean_d1 = deriv1.mean(axis=0).squeeze()
mean_d2 = deriv2.mean(axis=0).squeeze()

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# FIRST DERIVATIVE
for i in range(fd_beats.data_matrix.shape[0]):
    axes[0].plot(beat_grid, deriv1[i], color="black", alpha=0.15)

axes[0].plot(beat_grid, mean_d1, color="red", linewidth=2)
axes[0].axhline(0, color="grey", linestyle="--")
axes[0].set_title("First derivative (slope)")
axes[0].set_ylabel("mV/s")

# SECOND DERIVATIVE
for i in range(fd_beats.data_matrix.shape[0]):
    axes[1].plot(beat_grid, deriv2[i], color="black", alpha=0.15)

axes[1].plot(beat_grid, mean_d2, color="blue", linewidth=2)
axes[1].axhline(0, color="grey", linestyle="--")
axes[1].set_title("Second derivative (acceleration)")
axes[1].set_ylabel("mV/s²")

plt.show()

df_valid = kovariantes.iloc[valid_idx].copy()
rhythms_v = kovariantes.iloc[valid_idx]["Rhythm"].values

# top groups
from collections import Counter
top_groups = [k for k, _ in Counter(rhythms_v).most_common(4)]

def get_group(r):
    return fd_beats.data_matrix[rhythms_v == r]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for ax, r in zip(axes, top_groups):
    grp = get_group(r)

    mean_curve = grp.mean(axis=0)
    median_curve = np.median(grp, axis=0)

    trim = int(0.15 * grp.shape[0])
    sorted_grp = np.sort(grp, axis=0)
    trimmed_mean = sorted_grp[trim:-trim].mean(axis=0)

    ax.plot(beat_grid, mean_curve, label="mean")
    ax.plot(beat_grid, trimmed_mean, label="trimmed mean")
    ax.plot(beat_grid, median_curve, label="median")

    ax.set_title(f"{r}")
    ax.legend()

plt.tight_layout()
plt.show()

def trimmed_var(x, trim=0.15):
    k = int(trim * x.shape[0])
    xs = np.sort(x, axis=0)[k:-k]
    return xs.var(axis=0)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for ax, r in zip(axes, top_groups):
    grp = fd_beats.data_matrix[rhythms_v == r]

    ax.plot(beat_grid, grp.var(axis=0), label="var")
    ax.plot(beat_grid, trimmed_var(grp), label="trim var")

    ax.set_title(f"{r}")
    ax.legend()

plt.tight_layout()
plt.show()

unique_r = np.unique(rhythms_v)
colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_r)))
color_map = dict(zip(unique_r, colors))

overall = fd_beats.data_matrix.mean(axis=0)

plt.figure(figsize=(10, 5))
plt.plot(beat_grid, overall, "k--", linewidth=2, label="Overall")

for r in unique_r:
    grp = fd_beats.data_matrix[rhythms_v == r]
    if grp.shape[0] >= 2:
        plt.plot(beat_grid, grp.mean(axis=0),
                 color=color_map[r], label=r)

plt.legend()
plt.title("Mean ECG beat by rhythm group")
plt.show()

fpca = FPCA(n_components=5)
scores = fpca.fit_transform(fd_beats)

harmonics = fpca.components_
mean_fn = fpca.mean_
var_exp = fpca.explained_variance_ratio_

print(var_exp)

# scree
plt.plot(var_exp, marker="o")
plt.title("Scree plot")
plt.show()

# harmonics
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for i in range(3):
    comp = harmonics[i].data_matrix.squeeze()
    axes[i].plot(beat_grid, comp)
    axes[i].set_title(f"PC{i+1}")

plt.show()

import numpy as np
import matplotlib.pyplot as plt

rhythms_v = np.array(rhythms_v)

unique = np.unique(rhythms_v)
color_map = {k: plt.cm.tab10(i % 10) for i, k in enumerate(unique)}

plt.figure()

for i in range(scores.shape[0]):
    plt.scatter(
        scores[i, 0],
        scores[i, 1],
        color=color_map[rhythms_v[i]],
        s=20
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("FPCA Scores by Rhythm (Registered Beats)")
plt.show()

mean_vals = fd_beats.data_matrix.mean(axis=0).squeeze()

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for i in range(2):
    comp = harmonics[i].data_matrix.squeeze()
    sd = np.sqrt(fpca.explained_variance_[i])

    ax = axes[i]
    ax.plot(beat_grid, mean_vals, color="black", linewidth=2)

    ax.plot(beat_grid, mean_vals + 2 * sd * comp, "b--")
    ax.plot(beat_grid, mean_vals - 2 * sd * comp, "r--")

    ax.set_title(f"PC{i+1} variation (registered)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")

plt.tight_layout()
plt.show()

def varimax(Phi, gamma=1.0, q=20, tol=1e-6):
    p, k = Phi.shape
    R = np.eye(k)
    d = 0

    for _ in range(q):
        d_old = d
        Lambda = Phi @ R

        u, s, vh = np.linalg.svd(
            Phi.T @ (
                Lambda**3 - (gamma / p) * Lambda @ np.diag(np.diag(Lambda.T @ Lambda))
            )
        )

        R = u @ vh
        d = np.sum(s)

        if d_old != 0 and d / d_old < 1 + tol:
            break

    return Phi @ R

scores_rot = varimax(scores)
harm_rot = harmonics

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

for i, ax in enumerate(axes.flat):

    comp = harm_rot[i].data_matrix.squeeze()
    sd = np.sqrt(fpca.explained_variance_[i])

    ax.plot(beat_grid, mean_vals, color="black", linewidth=2)

    ax.plot(beat_grid, mean_vals + 2 * sd * comp, "r--")
    ax.plot(beat_grid, mean_vals - 2 * sd * comp, "r--")

    ax.set_title(f"Rotated PC{i+1} (registered)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("mV")

plt.tight_layout()
plt.show()

# ------------------------------------------------------------------------
# HYPOTHESIS TESTING
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------

X = fd_beats.data_matrix.squeeze()

time_grid = beat_grid

# align covariates
kovariantes_valid = kovariantes.iloc[valid_idx].copy()

# ------------------------------------------------------------------------
# GENDER GROUPS
# ------------------------------------------------------------------------

male_idx = kovariantes_valid["Gender"] == "MALE"
female_idx = kovariantes_valid["Gender"] == "FEMALE"

fd_beats_m = X[male_idx]
fd_beats_f = X[female_idx]

# ------------------------------------------------------------------------
# PLOT — Mean ECG by gender
# ------------------------------------------------------------------------

overall_mean = X.mean(axis=0)
male_mean = fd_beats_m.mean(axis=0)
female_mean = fd_beats_f.mean(axis=0)

plt.figure(figsize=(10,5))

plt.plot(time_grid, overall_mean,
         color="black", linestyle="--", linewidth=2,
         label="Overall")

plt.plot(time_grid, female_mean,
         color="red", linewidth=2,
         label="Female")

plt.plot(time_grid, male_mean,
         color="blue", linewidth=2,
         label="Male")

plt.axhline(0, color="grey", linestyle="--")

plt.title("Mean ECG beat by gender")
plt.xlabel("Time (s)")
plt.ylabel("mV")

plt.legend()
plt.show()

# ------------------------------------------------------------------------
# POINTWISE TWO-SAMPLE Z TEST
# ------------------------------------------------------------------------

def Ztwosample(x, y, t_seq, alpha=0.05):

    n1 = x.shape[0]
    n2 = y.shape[0]

    mean1 = x.mean(axis=0)
    mean2 = y.mean(axis=0)

    var1 = x.var(axis=0, ddof=1)
    var2 = y.var(axis=0, ddof=1)

    zvals = (mean1 - mean2) / np.sqrt(var1/n1 + var2/n2)

    pvals = 2 * (1 - norm.cdf(np.abs(zvals)))

    significant = pvals < alpha

    plt.figure(figsize=(10,4))

    plt.plot(t_seq, pvals)

    plt.axhline(alpha,
                color="red",
                linestyle="--")

    plt.title("Pointwise Z-test p-values")
    plt.xlabel("Time")
    plt.ylabel("p-value")

    plt.show()

    return {
        "z": zvals,
        "pvalue": pvals,
        "significant": significant
    }


stat = Ztwosample(
    fd_beats_m,
    fd_beats_f,
    time_grid
)

# ------------------------------------------------------------------------
# L2 TWO-SAMPLE TEST
# ------------------------------------------------------------------------

def L2_stat_twosample(
        x,
        y,
        t_seq,
        method=1,
        replications=500
):

    mean1 = x.mean(axis=0)
    mean2 = y.mean(axis=0)

    dt = t_seq[1] - t_seq[0]

    stat = np.sum((mean1 - mean2)**2) * dt

    # asymptotic
    if method == 1:

        return {
            "statistic": stat,
            "pvalue": None
        }

    # permutation
    if method == 2:

        combined = np.vstack([x, y])

        n1 = x.shape[0]
        n2 = y.shape[0]

        perm_stats = []

        for _ in range(replications):

            idx = np.random.permutation(n1 + n2)

            g1 = combined[idx[:n1]]
            g2 = combined[idx[n1:]]

            m1 = g1.mean(axis=0)
            m2 = g2.mean(axis=0)

            s = np.sum((m1 - m2)**2) * dt

            perm_stats.append(s)

        perm_stats = np.array(perm_stats)

        pval = np.mean(perm_stats >= stat)

        plt.figure(figsize=(6,4))

        plt.hist(perm_stats, bins=30)

        plt.axvline(stat,
                    color="red",
                    linewidth=2)

        plt.title("Permutation distribution — L2 statistic")

        plt.show()

        return {
            "statistic": stat,
            "pvalue": pval,
            "perm_stats": perm_stats
        }


stat_l2_1 = L2_stat_twosample(
    fd_beats_m,
    fd_beats_f,
    time_grid,
    method=1
)

print(stat_l2_1)

stat_l2_2 = L2_stat_twosample(
    fd_beats_m,
    fd_beats_f,
    time_grid,
    method=2,
    replications=500
)

print(stat_l2_2["pvalue"])

# ------------------------------------------------------------------------
# FUNCTIONAL F-TYPE TEST
# ------------------------------------------------------------------------

def F_stat_twosample(
        x,
        y,
        t_seq,
        method=1,
        replications=500
):

    mean1 = x.mean(axis=0)
    mean2 = y.mean(axis=0)

    var1 = x.var(axis=0, ddof=1)
    var2 = y.var(axis=0, ddof=1)

    dt = t_seq[1] - t_seq[0]

    numerator = np.sum((mean1 - mean2)**2) * dt

    denominator = (
        np.sum(var1) + np.sum(var2)
    ) * dt

    stat = numerator / denominator

    if method == 1:

        return {
            "statistic": stat,
            "pvalue": None
        }

    if method == 2:

        combined = np.vstack([x, y])

        n1 = x.shape[0]
        n2 = y.shape[0]

        perm_stats = []

        for _ in range(replications):

            idx = np.random.permutation(n1 + n2)

            g1 = combined[idx[:n1]]
            g2 = combined[idx[n1:]]

            m1 = g1.mean(axis=0)
            m2 = g2.mean(axis=0)

            v1 = g1.var(axis=0, ddof=1)
            v2 = g2.var(axis=0, ddof=1)

            num = np.sum((m1 - m2)**2) * dt

            den = (np.sum(v1) + np.sum(v2)) * dt

            perm_stats.append(num / den)

        perm_stats = np.array(perm_stats)

        pval = np.mean(perm_stats >= stat)

        plt.figure(figsize=(6,4))

        plt.hist(perm_stats, bins=30)

        plt.axvline(stat,
                    color="red",
                    linewidth=2)

        plt.title("Permutation distribution — F statistic")

        plt.show()

        return {
            "statistic": stat,
            "pvalue": pval,
            "perm_stats": perm_stats
        }


stat_f_1 = F_stat_twosample(
    fd_beats_m,
    fd_beats_f,
    time_grid,
    method=1
)

print(stat_f_1)

stat_f_2 = F_stat_twosample(
    fd_beats_m,
    fd_beats_f,
    time_grid,
    method=2,
    replications=500
)

print(stat_f_2["pvalue"])

# ------------------------------------------------------------------------
# PERMUTATION T TEST
# ------------------------------------------------------------------------

def permutation_t_test(x, y, n_perm=1000):

    obs = np.mean(x) - np.mean(y)

    combined = np.concatenate([x.flatten(), y.flatten()])

    n1 = x.size

    perm_stats = []

    for _ in range(n_perm):

        perm = np.random.permutation(combined)

        g1 = perm[:n1]
        g2 = perm[n1:]

        perm_stats.append(np.mean(g1) - np.mean(g2))

    perm_stats = np.array(perm_stats)

    pval = np.mean(np.abs(perm_stats) >= np.abs(obs))

    return {
        "statistic": obs,
        "pvalue": pval
    }


stat_perm = permutation_t_test(
    fd_beats_m,
    fd_beats_f
)

print(stat_perm)

# ------------------------------------------------------------------------
# RESULT INTERPRETATION
# ------------------------------------------------------------------------

print("No significant gender differences if p-values > 0.05")

# ------------------------------------------------------------------------
# 3-GROUP ANALYSIS — RHYTHM GROUPS
# ------------------------------------------------------------------------

sb_idx = kovariantes_valid["Rhythm"] == "SB"
afib_idx = kovariantes_valid["Rhythm"] == "AFIB"
sr_idx = kovariantes_valid["Rhythm"] == "SR"

fd_beats_sb = X[sb_idx]
fd_beats_afib = X[afib_idx]
fd_beats_sr = X[sr_idx]

# ------------------------------------------------------------------------
# PLOT — Rhythm group means
# ------------------------------------------------------------------------

plt.figure(figsize=(10,5))

plt.plot(time_grid,
         overall_mean,
         "k--",
         linewidth=2,
         label="Overall")

plt.plot(time_grid,
         fd_beats_sb.mean(axis=0),
         color="green",
         linewidth=2,
         label="SB")

plt.plot(time_grid,
         fd_beats_afib.mean(axis=0),
         color="blue",
         linewidth=2,
         label="AFIB")

plt.plot(time_grid,
         fd_beats_sr.mean(axis=0),
         color="red",
         linewidth=2,
         label="SR")

plt.axhline(0,
            color="grey",
            linestyle="--")

plt.title("Mean ECG beat by rhythm groups")

plt.xlabel("Time (s)")
plt.ylabel("mV")

plt.legend()
plt.show()

# ------------------------------------------------------------------------
# POINTWISE FUNCTIONAL ANOVA
# ------------------------------------------------------------------------

def fANOVA_pointwise(data, groups, t_seq, alpha=0.05):

    n_time = data.shape[1]

    pvals = np.zeros(n_time)

    unique_groups = np.unique(groups)

    group_means = {}

    for g in unique_groups:
        group_means[g] = data[groups == g].mean(axis=0)

    overall_mean = data.mean(axis=0)

    # ------------------------------------
    # pointwise ANOVA
    # ------------------------------------

    for i in range(n_time):

        samples = [
            data[groups == g, i]
            for g in unique_groups
        ]

        _, p = f_oneway(*samples)

        pvals[i] = p

    # ------------------------------------
    # PLOT p-values
    # ------------------------------------

    plt.figure(figsize=(10,4))

    plt.plot(t_seq, pvals)

    plt.axhline(alpha,
                color="red",
                linestyle="--")

    plt.title("Pointwise ANOVA p-values")

    plt.xlabel("Time")
    plt.ylabel("p-value")

    plt.show()

    # ------------------------------------
    # PLOT group means
    # ------------------------------------

    plt.figure(figsize=(10,5))

    plt.plot(t_seq,
             overall_mean,
             color="black",
             linewidth=2,
             label="Overall")

    for g in unique_groups:

        plt.plot(
            t_seq,
            group_means[g],
            linewidth=2,
            label=g
        )

    plt.title("Group means")

    plt.xlabel("Time")
    plt.ylabel("Mean")

    plt.legend()

    plt.show()

    return {
        "pvalues": pvals,
        "group_means": group_means,
        "overall_mean": overall_mean
    }


fanova_results = fANOVA_pointwise(
    X,
    kovariantes_valid["Rhythm"].values,
    time_grid
)

# ------------------------------------------------------------------------
# GLOBAL FUNCTIONAL ANOVA
# ------------------------------------------------------------------------

group_labels = kovariantes_valid["Rhythm"].values

groups_unique = np.unique(group_labels)

samples = [
    X[group_labels == g].mean(axis=1)
    for g in groups_unique
]

F_stat, pval = f_oneway(*samples)

print("GLOBAL FUNCTIONAL ANOVA")
print("F statistic:", F_stat)
print("p-value:", pval)

def pointwise_fanova(data, groups, t_seq, alpha=0.05):
    HAS_TUKEY = True
    data = np.asarray(data)
    groups = np.asarray(groups)

    unique_groups = np.unique(groups)
    k = len(unique_groups)
    n_time = data.shape[1]

    pvals = np.zeros(n_time)
    group_means = np.zeros((n_time, k))

    pair_labels = [
        ("SR", "SB"),
        ("SR", "AFIB"),
        ("SB", "AFIB")
    ]

    tukey_mat = np.zeros((n_time, len(pair_labels)))

    # --------------------------
    # MAIN LOOP (POINTWISE TEST)
    # --------------------------
    for t in range(n_time):

        samples = [data[groups == g, t] for g in unique_groups]

        # ANOVA
        pvals[t] = f_oneway(*samples).pvalue

        # group means
        group_means[t, :] = [np.mean(s) for s in samples]

        # Tukey HSD (per time point)
        if HAS_TUKEY:
            tmp_df = pd.DataFrame({
                "value": data[:, t],
                "group": groups
            })

            tukey = pairwise_tukeyhsd(
                endog=tmp_df["value"],
                groups=tmp_df["group"],
                alpha=alpha
            )

            # clean extraction from Tukey summary table
            res = pd.DataFrame(
                data=tukey.summary().data[1:],  # skip header
                columns=tukey.summary().data[0]
            )

            for i, (g1, g2) in enumerate(pair_labels):

                match = res[
                    ((res["group1"] == g1) & (res["group2"] == g2)) |
                    ((res["group1"] == g2) & (res["group2"] == g1))
                    ]

                if len(match) > 0:
                    tukey_mat[t, i] = float(match["p-adj"].values[0])
                else:
                    tukey_mat[t, i] = np.nan

    # overall mean
    overall_mean = np.mean(data, axis=0)

    plt.figure(figsize=(10, 8))

    # ANOVA p-values
    plt.subplot(2, 1, 1)
    plt.plot(t_seq, pvals, color="black")
    plt.axhline(alpha, color="blue", linestyle="--")
    plt.title("Pointwise ANOVA p-values")
    plt.ylabel("p-value")

    # group means
    plt.subplot(2, 1, 2)
    plt.plot(t_seq, overall_mean, color="black", linewidth=2, label="Overall mean")

    for i, g in enumerate(unique_groups):
        plt.plot(t_seq, group_means[:, i], label=str(g))

    plt.title("Group means over time")
    plt.xlabel("Time")
    plt.ylabel("Mean")
    plt.legend()

    plt.tight_layout()
    plt.show()

    # --------------------------
    # TUKEY PLOTS
    # --------------------------
    if HAS_TUKEY:

        plt.figure(figsize=(10, 6))

        for i, (g1, g2) in enumerate(pair_labels):
            plt.plot(t_seq, tukey_mat[:, i], label=f"{g1} vs {g2}")

        plt.axhline(alpha, color="blue", linestyle="--")
        plt.title("Tukey HSD adjusted p-values (pairwise)")
        plt.xlabel("Time")
        plt.ylabel("p-value")
        plt.legend()
        plt.tight_layout()
        plt.show()

    else:
        print("Skipping plots")

    return {
        "p_values": pvals,
        "group_means": group_means,
        "overall_mean": overall_mean,
        "tukey": tukey_mat if HAS_TUKEY else None
    }

data = np.asarray(fd_beats.data_matrix)

groups = np.asarray(kovariantes_valid["Rhythm"])

if data.ndim == 3 and data.shape[2] == 1:
    data = data.squeeze(axis=2)

result = pointwise_fanova(
    data=data,
    groups=groups,
    t_seq=beat_grid
)

# ----------------------
# Functional regression
# -----------------------
from pygam import LinearGAM, s
beat_mat = np.asarray(fd_beats.data_matrix).squeeze()

dat = pd.DataFrame({
    "Rhythm": pd.Categorical(
        kovariantes_valid["Rhythm"],
        categories=["AFIB", "SR", "SB"]
    ),
    "Age": pd.to_numeric(kovariantes_valid["PatientAge"]),
    "HR": pd.to_numeric(kovariantes_valid["VentricularRate"]),
    "RAxis": pd.to_numeric(kovariantes_valid["RAxis"]),
    "TAxis": pd.to_numeric(kovariantes_valid["TAxis"])
})
#----------------
# Design matrices
#----------------
n_subjects = beat_mat.shape[0]
n_time = beat_mat.shape[1]

long_df = pd.DataFrame({
    "ECG": beat_mat.flatten(),

    "time": np.tile(
        beat_grid,
        n_subjects
    ),

    "Rhythm": np.repeat(
        dat["Rhythm"].values,
        n_time
    ),

    "Age": np.repeat(
        dat["Age"].values,
        n_time
    ),

    "HR": np.repeat(
        dat["HR"].values,
        n_time
    ),

    "RAxis": np.repeat(
        dat["RAxis"].values,
        n_time
    ),

    "TAxis": np.repeat(
        dat["TAxis"].values,
        n_time
    )
})

long_df = pd.get_dummies(
    long_df,
    columns=["Rhythm"],
    drop_first=True
)

print(long_df.head())

X = long_df[[
    "time",
    "Rhythm_SR",
    "Rhythm_SB",
    "Age",
    "HR",
    "RAxis",
    "TAxis"
]].values.astype(float)

y = long_df["ECG"].values.astype(float)


from pygam import LinearGAM, s, l

gam = LinearGAM(

    # baseline smooth over time
    s(0, n_splines=10) +

    # linear parametric effects
    l(1) +   # Rhythm_SR
    l(2) +   # Rhythm_SB
    l(3) +   # Age
    l(4) +   # HR
    l(5) +   # RAxis
    l(6) +   # TAxis

    # varying coefficient smooths
    s(0, by=1, n_splines=10) +
    s(0, by=2, n_splines=10) +

    s(0, by=3, n_splines=10) +
    s(0, by=4, n_splines=10) +
    s(0, by=5, n_splines=10) +
    s(0, by=6, n_splines=10)

).fit(X, y)

# --------------
# Model summary
# --------------

print(gam.summary())

# ------------------------
# Pseudo R2 table
# ------------------------
r2_table = pd.DataFrame({
    "Model": ["GAM full model"],
    "Pseudo_R2": [gam.statistics_["pseudo_r2"]["explained_deviance"]]
})

print("\nMODEL FIT")
print(r2_table.round(4))

titles = [
    "Intercept - baseline beat shape",
    "SR effect (vs AFIB)",
    "SB effect (vs AFIB)",
    "Age effect (per year)",
    "Heart rate effect (per bpm)",
    "R-axis effect (per degree)",
    "T-axis effect (per degree)"
]

fig, axes = plt.subplots(3, 3, figsize=(15, 12))

axes = axes.flatten()

for term in range(7):

    ax = axes[term]

    XX = gam.generate_X_grid(term=term)

    pdep = gam.partial_dependence(
        term=term,
        X=XX
    )

    ax.plot(
        XX[:, 0],
        pdep,
        linewidth=2
    )

    ax.axhline(
        0,
        linestyle="--",
        color="grey"
    )

    ax.set_title(titles[term])

    ax.set_xlabel("Time since QRS (s)")
    ax.set_ylabel(r"$\beta(t)$")

# hide unused panels
for j in range(7, 9):
    axes[j].axis("off")

plt.tight_layout()
plt.show()

# ---------------
# Observed vs fitted
# --------------
pred = gam.predict(X)

long_df["pred"] = pred

pred_matrix = pred.reshape(n_subjects, n_time)

observed_mean = beat_mat.mean(axis=0)
predicted_mean = pred_matrix.mean(axis=0)

plt.figure(figsize=(10, 5))

plt.plot(
    beat_grid,
    observed_mean,
    color="black",
    linewidth=2,
    label="Observed"
)

plt.plot(
    beat_grid,
    predicted_mean,
    color="red",
    linestyle="--",
    linewidth=2,
    label="Predicted"
)

plt.xlabel("Time since QRS (s)")
plt.ylabel("mV")

plt.title("Observed vs predicted mean ECG")

plt.legend()

plt.show()

# -------------------
# Sensitivity analysis
# --------------------
def fit_sensitivity_model(k):

    gam_k = LinearGAM(

        s(0, n_splines=k) +

        l(1) +
        l(2) +
        l(3) +
        l(4) +
        l(5) +
        l(6) +

        s(0, by=1, n_splines=k) +
        s(0, by=2, n_splines=k) +

        s(0, by=3, n_splines=k) +
        s(0, by=4, n_splines=k) +
        s(0, by=5, n_splines=k) +
        s(0, by=6, n_splines=k)

    ).fit(X, y)

    return gam_k

fit_k5 = fit_sensitivity_model(5)
fit_k10 = fit_sensitivity_model(10)
fit_k15 = fit_sensitivity_model(15)

sensitivity = pd.DataFrame({
    "k": [5, 10, 15],

    "Pseudo_R2": [
        fit_k5.statistics_["pseudo_r2"]["explained_deviance"],
        fit_k10.statistics_["pseudo_r2"]["explained_deviance"],
        fit_k15.statistics_["pseudo_r2"]["explained_deviance"]
    ],

    "EDF": [
        fit_k5.statistics_["edof"],
        fit_k10.statistics_["edof"],
        fit_k15.statistics_["edof"]
    ]
})

print("\nSENSITIVITY ANALYSIS")
print(sensitivity.round(4))
# ------------------------------------------------------------------------
# CLASSIFICATION
# ------------------------------------------------------------------------
# RESPONSE VARIABLE
# ------------------------------------------------------------------------

kovariantes_valid["bad_rhythm"] = np.where(
    kovariantes_valid["Rhythm"] == "AFIB",
    1,
    0
)

y = kovariantes_valid["bad_rhythm"].values

# ------------------------------------------------------------------------
# FUNCTIONAL DATA
# ------------------------------------------------------------------------

X_fd = np.asarray(
    fd_beats.data_matrix
).squeeze()

# ------------------------------------------------------------------------
# COVARIATES
# ------------------------------------------------------------------------

cov_df = pd.DataFrame({
    "Age": pd.to_numeric(
        kovariantes_valid["PatientAge"],
        errors="coerce"
    ),

    "Gender": kovariantes_valid["Gender"],

    "VenticularRate": pd.to_numeric(
        kovariantes_valid["VentricularRate"],
        errors="coerce"
    ),

    "Arti": pd.to_numeric(
        kovariantes_valid["AtrialRate"],
        errors="coerce"
    )
})

# ------------------------------------------------------------------------
# TRAIN / TEST SPLIT
# ------------------------------------------------------------------------

indices = np.arange(len(y))

train_idx, test_idx = train_test_split(
    indices,
    test_size=0.20,
    stratify=y,
    random_state=123,
    shuffle=True
)

X_fd_train = X_fd[train_idx]
X_fd_test  = X_fd[test_idx]

y_train = y[train_idx]
y_test  = y[test_idx]

cov_train = cov_df.iloc[train_idx].reset_index(drop=True)
cov_test  = cov_df.iloc[test_idx].reset_index(drop=True)

# ------------------------------------------------------------------------
# FPCA
# ------------------------------------------------------------------------

fd_train = FDataGrid(
    data_matrix=X_fd_train,
    grid_points=beat_grid
)

fd_test = FDataGrid(
    data_matrix=X_fd_test,
    grid_points=beat_grid
)

fpca = FPCA(n_components=3)

scores_train = fpca.fit_transform(fd_train)
scores_test = fpca.transform(fd_test)

# ------------------------------------------------------------------------
# BUILD FEATURE MATRICES
# ------------------------------------------------------------------------

scores_train_df = pd.DataFrame(
    scores_train,
    columns=["PC1", "PC2", "PC3"]
)

scores_test_df = pd.DataFrame(
    scores_test,
    columns=["PC1", "PC2", "PC3"]
)

X_train = pd.concat(
    [
        scores_train_df,
        cov_train.reset_index(drop=True)
    ],
    axis=1
)

X_test = pd.concat(
    [
        scores_test_df,
        cov_test.reset_index(drop=True)
    ],
    axis=1
)

# ------------------------------------------------------------------------
# ENCODE GENDER
# ------------------------------------------------------------------------

X_train = pd.get_dummies(
    X_train,
    columns=["Gender"],
    drop_first=True
)

X_test = pd.get_dummies(
    X_test,
    columns=["Gender"],
    drop_first=True
)

# ensure same columns
X_test = X_test.reindex(
    columns=X_train.columns,
    fill_value=0
)

# ------------------------------------------------------------------------
# STANDARDIZATION
# ------------------------------------------------------------------------

scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# ------------------------------------------------------------------------
# MODELS
# ------------------------------------------------------------------------

models = {

    "LDA": LinearDiscriminantAnalysis(),

    "QDA": QuadraticDiscriminantAnalysis(),

    "RPART": DecisionTreeClassifier(
        random_state=123
    )
}

# ------------------------------------------------------------------------
# ROC PLOT
# ------------------------------------------------------------------------

plt.figure(figsize=(8, 6))

colors = {
    "LDA": "red",
    "QDA": "blue",
    "RPART": "green"
}

results = {}
auc_values = {}

for model_name, model in models.items():

    # ------------------------------------------------------------
    # FIT
    # ------------------------------------------------------------

    model.fit(
        X_train_scaled,
        y_train
    )

    # ------------------------------------------------------------
    # TRAIN ROC
    # ------------------------------------------------------------

    train_prob = model.predict_proba(
        X_train_scaled
    )[:, 1]

    fpr_train, tpr_train, thresh_train = roc_curve(
        y_train,
        train_prob,
        pos_label=1
    )

    # Youden index (same logic as R pROC::coords(...,"best"))
    best_idx = np.argmax(
        tpr_train - fpr_train
    )

    best_threshold = thresh_train[best_idx]

    # ------------------------------------------------------------
    # TEST ROC
    # ------------------------------------------------------------

    test_prob = model.predict_proba(
        X_test_scaled
    )[:, 1]

    fpr_test, tpr_test, thresh_test = roc_curve(
        y_test,
        test_prob,
        pos_label=1
    )

    auc_val = auc(
        fpr_test,
        tpr_test
    )

    auc_values[model_name] = auc_val

    # ------------------------------------------------------------
    # ROC LINE
    # ------------------------------------------------------------

    plt.plot(
        fpr_test,
        tpr_test,
        color=colors[model_name],
        linewidth=2,
        label=f"{model_name} AUC={auc_val:.3f}"
    )

    # ------------------------------------------------------------
    # PREDICTED CLASSES
    # ------------------------------------------------------------

    pred_class = (
        test_prob >= best_threshold
    ).astype(int)

    # ------------------------------------------------------------
    # CONFUSION MATRIX
    # ------------------------------------------------------------

    cm = confusion_matrix(
        y_test,
        pred_class,
        labels=[0, 1]
    )

    results[model_name] = cm

    # ------------------------------------------------------------
    # PRINT RESULTS
    # ------------------------------------------------------------

    print("\n")
    print("MODEL:", model_name)
    print("AUC:", round(auc_val, 3))
    print(cm)

# diagonal
plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--",
    color="black"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")

plt.title("ROC curves")

plt.legend()

plt.show()

# ------------------------------------------------------------------------
# CONFUSION MATRIX — QDA
# ------------------------------------------------------------------------

cm_qda = results["QDA"]

fig, ax = plt.subplots(figsize=(6, 5))

disp = ConfusionMatrixDisplay(
    confusion_matrix=cm_qda,
    display_labels=["0", "1"]
)

disp.plot(
    cmap="Blues",
    ax=ax,
    colorbar=False
)

plt.title("Confusion Matrix — QDA model")

plt.show()

# ------------------------------------------------------------------------
# SUMMARY TABLE
# ------------------------------------------------------------------------

summary_df = pd.DataFrame({
    "Model": list(auc_values.keys()),
    "AUC": list(auc_values.values())
})

print("\nAUC SUMMARY")
print(summary_df.round(3))