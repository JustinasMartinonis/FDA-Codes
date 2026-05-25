library(dplyr)
library(fda.usc)
library(fda)

# ------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------
# selected_ids <- c(
#   "MUSE_20180113_171753_06000","MUSE_20180113_173729_65000",
#   "MUSE_20180712_154333_20000","MUSE_20180209_125550_60000",
#   "MUSE_20180114_070437_19000","MUSE_20180113_123828_45000",
#   "MUSE_20180210_120030_82000","MUSE_20180114_075642_73000",
#   "MUSE_20180115_124238_55000","MUSE_20180118_131051_42000",
#   "MUSE_20180210_131911_95000","MUSE_20180113_183530_81000",
#   "MUSE_20180115_133409_97000","MUSE_20180209_173515_96000",
#   "MUSE_20180209_174606_15000","MUSE_20180712_151719_95000",
#   "MUSE_20180114_065222_94000","MUSE_20180116_174129_66000",
#   "MUSE_20180113_175435_01000","MUSE_20180210_125415_39000",
#   "MUSE_20180712_151007_06000","MUSE_20180114_070350_28000",
#   "MUSE_20180113_115609_72000","MUSE_20180112_120916_00000",
#   "MUSE_20180114_172151_31000","MUSE_20180210_120343_67000",
#   "MUSE_20180119_175654_25000","MUSE_20180712_151031_25000",
#   "MUSE_20180209_172644_81000","MUSE_20180112_134110_84000",
#   "MUSE_20180118_170505_57000","MUSE_20180114_070704_83000",
#   "MUSE_20180114_073446_22000","MUSE_20180120_122519_07000",
#   "MUSE_20180114_130642_01000","MUSE_20180209_121003_67000",
#   "MUSE_20180112_072012_70000","MUSE_20180210_122813_91000",
#   "MUSE_20180209_114343_23000","MUSE_20180116_175157_61000",
#   "MUSE_20180209_130455_68000","MUSE_20180113_115416_77000",
#   "MUSE_20180114_134109_01000","MUSE_20180114_130026_36000",
#   "MUSE_20180114_131946_20000","MUSE_20180114_065438_26000",
#   "MUSE_20180112_073032_10000","MUSE_20180209_122213_15000",
#   "MUSE_20180111_160516_89000","MUSE_20180209_122642_22000",
#   "MUSE_20180118_180008_12000","MUSE_20180118_130122_04000",
#   "MUSE_20180118_174318_44000","MUSE_20180113_173559_20000",
#   "MUSE_20180118_131033_64000","MUSE_20180210_120643_69000",
#   "MUSE_20180209_132616_53000","MUSE_20180209_123437_34000",
#   "MUSE_20180712_151159_30000","MUSE_20180116_173530_20000",
#   "MUSE_20180712_152526_02000","MUSE_20180115_122241_82000",
#   "MUSE_20180113_131856_90000","MUSE_20180116_173550_05000",
#   "MUSE_20180112_124429_51000","MUSE_20180118_124235_14000",
#   "MUSE_20180116_122520_49000","MUSE_20180209_174351_56000",
#   "MUSE_20180118_171246_91000","MUSE_20180712_160516_95000",
#   "MUSE_20180210_124948_52000","MUSE_20180114_075504_92000",
#   "MUSE_20180114_130503_62000","MUSE_20180210_115910_62000",
#   "MUSE_20180114_073806_97000","MUSE_20180116_125007_72000",
#   "MUSE_20180712_151453_25000","MUSE_20180210_130027_42000",
#   "MUSE_20180115_121245_81000","MUSE_20180114_124111_74000",
#   "MUSE_20180119_180226_47000","MUSE_20180712_151423_33000",
#   "MUSE_20180116_134008_99000","MUSE_20180210_115410_26000",
#   "MUSE_20180119_180035_02000","MUSE_20180209_175337_71000",
#   "MUSE_20180113_073219_67000","MUSE_20180712_155838_47000",
#   "MUSE_20180115_123540_50000","MUSE_20180113_135514_17000",
#   "MUSE_20180114_131255_99000","MUSE_20180712_154636_61000",
#   "MUSE_20180112_123428_11000","MUSE_20180116_171345_23000",
#   "MUSE_20180115_132130_20000","MUSE_20180118_120822_25000",
#   "MUSE_20180209_171615_97000","MUSE_20180118_174246_35000",
#   "MUSE_20180209_173904_56000","MUSE_20180115_123740_94000"
# )

library(readxl)
library(dplyr)

n_per_group <- 30

diagnostics <- read_excel("Diagnostics.xlsx")

selected_ids <- diagnostics %>%
  filter(Rhythm %in% c("AFIB", "SB", "SR")) %>%
  arrange(Rhythm, FileName) %>%
  group_by(Rhythm) %>%
  slice_head(n = n_per_group) %>%
  ungroup() %>%
  pull(FileName)

length(selected_ids)
head(selected_ids)



files <- list.files("ECGData/", full.names = TRUE)

# Keep only files that match selected IDs
files_sample <- files[
  tools::file_path_sans_ext(basename(files)) %in% selected_ids
]

#set.seed(123)
#files_sample <- sample(files, 100)

read_file <- function(file) {
  df <- read.csv(file)
  if ("II" %in% colnames(df)) {
    column_to_leave <- df[["II"]]
  } else {
    column_to_leave <- df[[2]]
  }
  file_name <- tail(strsplit(file, "/|\\\\")[[1]], 1)
  file_id   <- strsplit(file_name, "\\.")[[1]][1]
  c(ID = file_id, t(column_to_leave))
}

result_list <- lapply(files_sample, read_file)
df          <- as.data.frame(do.call(rbind, result_list))

# Load covariates (rhythm labels, demographics)
kovariantes        <- readxl::read_excel("Diagnostics.xlsx")

kovariantes_sample <- kovariantes[kovariantes$FileName %in% df$ID, ]

write.csv(df, "ECGData/ECG_combined_true.csv", row.names = FALSE)
df <- read.csv("ECGData/ECG_combined_true.csv")

# ------------------------------------------------------------------------
# 2. FUNCTIONAL DATA OBJECT
# ------------------------------------------------------------------------

time      <- seq(0, 10, length.out = 500 * 10)   # 5000 points over 10s

df_be_id        <- df[, 2:ncol(df)]
df_be_id[]      <- lapply(df_be_id, as.numeric)

df_scaled       <- df_be_id / 1000                # microvolts -> millivolts
fdata_obj_scaled <- fdata(df_scaled, argvals = time)

plot(fdata_obj_scaled,
     main = "ECG Lead II — 90 patients (raw)",
     xlab = "Time (s)", ylab = "Amplitude (mV)")

# ------------------------------------------------------------------------
# 3. SMOOTHING PARAMETER SELECTION — GCV curve
# ------------------------------------------------------------------------

lambda_fine  <- 10^seq(-8, -2, by = 0.5)
nb_fixed     <- 500
cv_lambda    <- data.frame()

for (l in lambda_fine) {
  basis     <- create.bspline.basis(rangeval = c(0, 10), nbasis = nb_fixed)
  fdPar_obj <- fdPar(basis, Lfdobj = int2Lfd(1), lambda = l)
  sm        <- smooth.basis(fdata_obj_scaled$argvals,
                            t(fdata_obj_scaled$data), fdPar_obj)
  cv_lambda <- rbind(cv_lambda, data.frame(lambda = l, GCV = mean(sm$gcv)))
}

# Plot GCV curve — ieškome taško nuo kada pradeda kilti

plot(log10(cv_lambda$lambda), cv_lambda$GCV, type = "l", lwd = 2,
     xlab = "log10(lambda)", ylab = "Mean GCV",
     main = "GCV vs Smoothing Parameter (nbasis = 500)")
abline(v = log10(1e-04), col = "red", lty = 2, lwd = 1.5)
legend("topleft", legend = "Selected: lambda = 1e-04", col = "red", lty = 2)

# from log10(lambda) = -3 onward GCV climbs steeply. 
# This confirms that lambda above ~0.001 is oversmoothing.
# from -8 all the way to -4 the GCV barely changes. 
# This means lambda below 1e-04 makes essentially no difference 

# GCV difference nb 500 from 700 is negligible (0.00200 vs 0.00083 — both tiny)

# ------------------------------------------------------------------------
# 4. FINAL SMOOTHED OBJECT
#    nbasis=500, lambda=1e-04, penalty order=1 (penalise 1st derivative)
# ------------------------------------------------------------------------

basis     <- create.bspline.basis(rangeval = c(0, 10), nbasis = 500)
fdPar_obj <- fdPar(basis, Lfdobj = int2Lfd(1), lambda = 1e-04)
smoothed  <- smooth.basis(fdata_obj_scaled$argvals,
                          t(fdata_obj_scaled$data), fdPar_obj)
fd_final  <- smoothed$fd


# ------------------------------------------------------------------------
# 5. VALIDATION PLOTS
# ------------------------------------------------------------------------

# --- Plot A: Raw vs smoothed for 3 patients ---
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (i in c(1, 45,80)) {
  plot(fdata_obj_scaled[i],
       main = paste("Patient", i),
       xlab = "Time (s)", ylab = "mV", col = "black")
  lines(smoothed$fd[i], col = "red", lwd = 2)
  legend("topright", legend = c("Raw", "Smoothed"),
         col = c("black", "red"), lty = 1, cex = 0.7)
}
par(mfrow = c(1, 1))

# vienam pacientui raw vs smoothed
plot(fdata_obj_scaled[1],
     main = paste("Patient", 1),
     xlab = "Time (s)", ylab = "mV", col = "black")
lines(smoothed$fd[1], col = "red", lwd = 2)
legend("topright", legend = c("Raw", "Smoothed"),
       col = c("black", "red"), lty = 1, cex = 0.7)


# --- Plot C: Lambda comparison 

par(mfrow = c(1, 3))
i     <- 1
basis <- create.bspline.basis(c(0, 10), nbasis = 500)

# too wiggly
sm_under <- smooth.basis(fdata_obj_scaled$argvals, t(fdata_obj_scaled$data),
                         fdPar(basis, int2Lfd(1), lambda = 1e-08))
plot(fdata_obj_scaled[i], main = "lambda = 1e-08",
     xlab = "Time (s)", ylab = "mV", col = "grey70")
lines(sm_under$fd[i], col = "red", lwd = 2)

# optimal
sm_opt <- smooth.basis(fdata_obj_scaled$argvals, t(fdata_obj_scaled$data),
                       fdPar(basis, int2Lfd(1), lambda = 0.01))
plot(fdata_obj_scaled[i], main = "lambda = 1e-04",
     xlab = "Time (s)", ylab = "mV", col = "grey70")
lines(sm_opt$fd[i], col = "darkgreen", lwd = 2)

# oversmoothed
sm_over <- smooth.basis(fdata_obj_scaled$argvals, t(fdata_obj_scaled$data),
                        fdPar(basis, int2Lfd(1), lambda = 1))
plot(fdata_obj_scaled[i], main = "lambda = 1",
     xlab = "Time (s)", ylab = "mV", col = "grey70")
lines(sm_over$fd[i], col = "blue", lwd = 2)

par(mfrow = c(1, 1))

# ------------------------------------------------------------------------
# OUTLIERS
# ------------------------------------------------------------------------

library(fdaoutlier)
library(rainbow)

# Evaluate smoothed curves on time grid
eval_grid <- fdata_obj_scaled$argvals
tmp <- eval.fd(eval_grid, fd_final)

# --- Method 1: Band Depth ---
bd <- band_depth(dt = t(tmp))
plot(bd, type = "l", main = "Band Depth", 
     ylab = "BD Value", xlab = "Curve index")

mbd <- modified_band_depth(t(tmp))
plot(mbd, type = "l", main = "Modified Band Depth",
     ylab = "MBD Value", xlab = "Curve index")


# --- Method 3: Functional boxplot ---
fbplot(tmp, method = "MBD", main = "Functional Boxplot with MBD method", xlab = "Time (5000 = 10s)",
       ylab = "Amplitude (mV)")


# --- Method 4: MUOD (shape / amplitude / magnitude outliers) ---
m <- muod(t(tmp), cut_method = "boxplot")
m$outliers  # check all three types
# 72, 17, 28, 6 - kur kelis kartus pasikartoja

#suspicious <- c(23, 59,  3, 84,  9, 86, 92,  4, 10, 21) # rankiniu būdu atrinkti id
#patient_ids <- df$ID[suspicious]
#print(patient_ids)

#par(mfrow = c(1, 2))
#plot(fd_final[92], 
 #    main = paste("Curve", 92, ". Magnitude outlier"),
 #   xlab = "Time (s)", ylab = "mV", col = "red", lwd = 2)
#plot(fd_final[85], 
 #    main = paste("Curve", 85, ". Central curve"),
  #   xlab = "Time (s)", ylab = "mV", col = "red", lwd = 2)

# Plot them against the rest
#par(mfrow = c(2, 2))
#for(i in suspicious){
 # plot(fdata_obj_scaled[i], 
  #     main = paste("Curve", i, "— ID"),
   #    xlab = "Time (s)", ylab = "mV", col = "red", lwd = 2)
#}
#par(mfrow = c(1,1))


###################
# EDA
##################
library(fda.usc)
library(fields)

# SECTION 1 — skaidrese nenaudojom
# ------------------------------------------------------------------------
mean_beat   <- mean.fd(fd_final)
stddev_beat <- std.fd(fd_final)

plot(fd_final,
     main = "Smoothed ECG signals",
     xlab = "Time (s)", ylab = "mV",
     col  = rgb(0, 0, 0, 0.15), lwd = 0.8)
abline(h = 0, col = "grey80", lty = 3)
lines(mean_beat,                col = 2, lwd = 3, lty = 1)
lines(mean_beat + 3*stddev_beat,  col = 4, lwd = 2, lty = 1)
lines(mean_beat - 3*stddev_beat,  col = 4, lwd = 2, lty = 1)
legend("topright", lty = c(2, 2), col = c(2, 4), lwd = c(3, 2),
       legend = c("Mean", "Mean ± 3xSD"),
       box.col = "white", cex = 0.85)

# SECTION 2 — Bivariate covariance function
# ------------------------------------------------------------------------
library(fields)
ecgvar.bifd <- var.fd(fd_final)

# Subsample for computational efficiency
t_sub       <- time[seq(1, 500, by = 5)]   # 50 time points
ecgvar_mat  <- eval.bifd(t_sub, t_sub, ecgvar.bifd)
#ecgvar_mat  <- eval.bifd(time, time, ecgvar.bifd)

# Figure 1 — 3D perspective surface
persp(t_sub, t_sub, ecgvar_mat,
      theta = -45, phi = 25, r = 3, expand = 0.5,
      ticktype = "detailed", col = "lightblue",
      xlab = "Time (s)", ylab = "Time (s)",
      zlab = "Covariance (mV²)",
      main = "Covariance surface — ECG beats")

# Figure 2 — contour
contour(t_sub, t_sub, ecgvar_mat,
        xlab = "Time (s)", ylab = "Time (s)",
        main = "Covariance contour — ECG beats")

# Figure 3 — heatmap with contour overlay - SKAIDRESE
image.plot(t_sub, t_sub, ecgvar_mat,
           xlab = "Time (s)", ylab = "Time (s)",
           main = "Covariance heatmap — ECG beats")
contour(t_sub, t_sub, ecgvar_mat,
        col = "white", add = TRUE)

# SECTION 3 — Derivatives: slope and acceleration - SKAIDRESE
# ------------------------------------------------------------------------

deriv1 <- deriv.fd(fd_final, 1)   # 1st derivative = rate of change
deriv2 <- deriv.fd(fd_final, 2)   # 2nd derivative = acceleration

opar <- par(mfrow = c(1, 2))

plot(deriv1,
     main = "First derivative (slope)",
     xlab = "Time (s)", ylab = "mV/s",
     col  = rgb(0, 0, 0, 0.15), lwd = 0.8)
lines(mean.fd(deriv1), col = "red", lwd = 3, lty = 2)
abline(h = 0, col = "grey70", lty = 3)
legend("topright", legend = "Mean slope",
       col = "red", lty = 2, lwd = 2,
       box.col = "white", cex = 0.8)

plot(deriv2,
     main = "Second derivative (acceleration)",
     xlab = "Time (s)", ylab = "mV/s²",
     col  = rgb(0, 0, 0, 0.15), lwd = 0.8)
lines(mean.fd(deriv2), col = "blue", lwd = 3, lty = 2)
abline(h = 0, col = "grey70", lty = 3)
legend("topright", legend = "Mean acceleration",
       col = "blue", lty = 2, lwd = 2,
       box.col = "white", cex = 0.8)

par(opar)

# SECTION 4 — Centrality measures by rhythm group
# ------------------------------------------------------------------------

# Use largest rhythm groups
rhythms_v <- kovariantes_sample$Rhythm
rhythm_tab <- sort(table(rhythms_v), decreasing = TRUE)
top_groups <- names(rhythm_tab[1:min(4, length(rhythm_tab))])
cat("Top rhythm groups:", top_groups, "\n")
print(rhythm_tab[top_groups])

# Convert fd_final to fdata for fda.usc centrality/dispersion functions
eval_grid_full <- seq(0, 10, length.out = 500)
full_mat       <- t(eval.fd(eval_grid_full, fd_final))
fdata_full     <- fdata(full_mat, argvals = eval_grid_full)

get_fdata_group <- function(r) {
  fdata_full[which(rhythms_v == r), ]
}

opar <- par(mfrow = c(2, 2))

for(r in top_groups){
  grp  <- get_fdata_group(r)
  n_r  <- sum(rhythms_v == r)
  
  plot(func.mean(grp),
       main = paste0("Centrality — ", r, " (n=", n_r, ")"),
       xlab = "Time (s)", ylab = "mV",
       ylim = c(-0.2, 0.6))
  abline(h = 0, col = "grey80", lty = 3)
  legend("topright", cex = 0.65, box.col = "white",
         lty = 1:5, col = 1:5,
         legend = c("mean", "trim.mode", "trim.RP",
                    "median.mode", "median.RP"))
  lines(func.trim.mode(grp, trim = 0.15), col = 2, lty = 2)
  lines(func.trim.RP(grp,   trim = 0.15), col = 3, lty = 3)
  lines(func.med.mode(grp,  trim = 0.15), col = 4, lty = 4)
  lines(func.med.RP(grp,    trim = 0.15), col = 5, lty = 5)
}

par(opar)
# SECTION 5 — Dispersion measures by rhythm group
# ------------------------------------------------------------------------

opar <- par(mfrow = c(2, 2))

for(r in top_groups){
  grp <- get_fdata_group(r)
  n_r <- sum(rhythms_v == r)
  
  plot(func.var(grp),
       main = paste0("Dispersion — ", r, " (n=", n_r, ")"),
       xlab = "Time (s)", ylab = "Variance (mV²)")
  legend("topright", cex = 0.65, box.col = "white",
         lty = 1:3, col = 1:3,
         legend = c("var", "trimvar.mode", "trimvar.RP"))
  lines(func.trimvar.mode(grp, trim = 0.15), col = 2, lty = 2)
  lines(func.trimvar.RP(grp,   trim = 0.15), col = 3, lty = 3)
}

par(opar)
# SECTION 6 — Group mean curves overlaid
# ------------------------------------------------------------------------

unique_r <- unique(rhythms_v)
cols_map <- setNames(rainbow(length(unique_r)), unique_r)

plot(mean_beat,
     main = "Mean ECG beat by rhythm group",
     xlab = "Time (s)", ylab = "mV",
     col = "black", lwd = 2, lty = 2,
     ylim = c(-0.2, 0.6))
abline(h = 0, col = "grey80", lty = 3)

for(r in unique_r){
  idx <- which(rhythms_v == r)
  if(length(idx) >= 2)
    lines(mean.fd(fd_final[idx]), col = cols_map[r], lwd = 2)
}

legend("topright",
       legend = c("Overall", unique_r),
       col    = c("black", cols_map),
       lty    = c(2, rep(1, length(unique_r))),
       lwd = 2, cex = 0.65, box.col = "white")

# SECTION 7 - Functional Principal Component Analysis
# ------------------------------------------------------------------------

# Run FPCA from FD
pca_res <- pca.fd(fd_final, nharm = 5)

# Variance explained
var_explained <- pca_res$varprop
print(var_explained)

# Scree plot
plot(var_explained, type = "b", pch = 16,
     xlab = "Principal Component",
     ylab = "Variance Explained",
     main = "Scree Plot — FPCA")

# Plot first 3 harmonics (eigenfunctions)
par(mfrow = c(1, 3))
plot(pca_res$harmonics[1], main = "PC1")
plot(pca_res$harmonics[2], main = "PC2")
plot(pca_res$harmonics[3], main = "PC3")
par(mfrow = c(1,1))

# Mean ± variation for PC1 and PC2
mean_fd <- mean.fd(fd_final)

par(mfrow = c(1, 2))

# PC1 variation
plot(mean_fd, lwd = 2, main = "PC1 variation",
     xlab = "Time (s)", ylab = "mV")
lines(mean_fd + 2 * sqrt(pca_res$values[1]) * pca_res$harmonics[1],
      col = "blue", lty = 2)
lines(mean_fd - 2 * sqrt(pca_res$values[1]) * pca_res$harmonics[1],
      col = "red", lty = 2)
legend("topright", legend = c("Mean", "+2SD", "-2SD"),
       col = c("black", "blue", "red"), lty = c(1,2,2), cex = 0.8)

# PC2 variation
plot(mean_fd, lwd = 2, main = "PC2 variation",
     xlab = "Time (s)", ylab = "mV")
lines(mean_fd + 2 * sqrt(pca_res$values[2]) * pca_res$harmonics[2],
      col = "blue", lty = 2)
lines(mean_fd - 2 * sqrt(pca_res$values[2]) * pca_res$harmonics[2],
      col = "red", lty = 2)

par(mfrow = c(1,1))

# PCA scores
scores <- pca_res$scores

# Score plot
plot(scores[,1], scores[,2],
     xlab = "PC1", ylab = "PC2",
     main = "FPCA Scores",
     pch = 16, col = "blue")

# Score plot colored by rhythm (if available)
if (exists("rhythms_v")) {
  cols <- as.numeric(as.factor(rhythms_v))
  
  plot(scores[,1], scores[,2],
       col = cols, pch = 16,
       xlab = "PC1", ylab = "PC2",
       main = "FPCA Scores by Rhythm")
  
  legend("topright",
         legend = unique(rhythms_v),
         col = unique(cols),
         pch = 16, cex = 0.7)
}

# Rotated PCA plot - same as without rotation
#pdf("plots/PCA_rotated_plot.pdf", width = 8, height = 10)

pca_res <- pca.fd(fd_final, nharm = 5, centerfns = TRUE)
pca_rot <- varmx.pca.fd(pca_res)

mean_fd <- mean.fd(fd_final)
time_grid <- seq(0, 0.45, length.out = 250)  # same as beat_grid

par(mfrow = c(2,2), mar = c(4,4,3,1))

for(i in 1:4){
  harmon <- pca_rot$harmonics[i]
  sd_i   <- sqrt(pca_rot$values[i])
  var_pct <- round(100 * pca_rot$varprop[i], 1)
  
  # Evaluate numerically
  mean_vals <- eval.fd(time_grid, mean_fd)
  harm_vals <- eval.fd(time_grid, harmon)
  
  # Determine ylim
  ylim_vals <- range(c(mean_vals, 
                       mean_vals + 2*sd_i*harm_vals, 
                       mean_vals - 2*sd_i*harm_vals))
  
  # Plot
  plot(time_grid, mean_vals, type = "l", lwd = 2, ylim = ylim_vals,
       xlab = "Time (s)", ylab = "Amplitude (mV)",
       main = paste0("Rotated PCA function ", i,
                     " (Percentage of variability ", var_pct, ")"))
  
  # Add ±SD lines
  lines(time_grid, mean_vals + 1*sd_i*harm_vals, col = "blue", lty = 2)
  lines(time_grid, mean_vals - 1*sd_i*harm_vals, col = "blue", lty = 2)
  lines(time_grid, mean_vals + 2*sd_i*harm_vals, col = "red", lty = 3)
  lines(time_grid, mean_vals - 2*sd_i*harm_vals, col = "red", lty = 3)
}

par(mfrow = c(1,1))

# ------------------------------------------------------------------------
# REGISTRATION
# ------------------------------------------------------------------------
eval_grid_reg <- time
smooth_mat    <- eval.fd(eval_grid_reg, fd_final)
search_idx    <- which(eval_grid_reg <= 1)

qrs_times <- numeric(90)
for(i in 1:90){
  seg           <- smooth_mat[search_idx, i]
  qrs_times[i]  <- eval_grid_reg[search_idx[which.max(seg)]]
}
# Reduce pre-peak window so early-beat patients are included
window_before <- 0.05
window_after  <- 0.40   # slightly longer to keep full T wave
beat_grid     <- seq(0, window_before + window_after, length.out = 250)

# Now all patients with QRS > 0.05s are valid
valid    <- qrs_times >= window_before & qrs_times <= (10 - window_after)
cat("Valid patients:", sum(valid), "/ 99\n")

valid_idx         <- which(valid)
n_valid           <- sum(valid)
beat_matrix       <- matrix(NA, nrow = n_valid, ncol = 250)

for(j in seq_along(valid_idx)){
  i               <- valid_idx[j]
  t_window        <- seq(qrs_times[i] - window_before,
                         qrs_times[i] + window_after,
                         length.out = 250)
  beat_matrix[j,] <- eval.fd(t_window, fd_final[i])
}

df_valid          <- df[valid_idx, ]
kovariantes_valid <- kovariantes[kovariantes$FileName %in% df_valid$ID, ]
print(table(kovariantes_valid$Rhythm))

# Keep only the largest rhythm groups (≥10 patients)
n_valid           <- nrow(beat_matrix)
cat("Patients after group filter:", n_valid, "\n")

fdata_beats <- fdata(beat_matrix, argvals = beat_grid)

# Plot all beats coloured by rhythm group
rhythms_v <- kovariantes_valid$Rhythm
cols_map   <- setNames(rainbow(length(unique(rhythms_v))), unique(rhythms_v))

plot(fdata_beats[1,],
     main = "Single beat — aligned at QRS, coloured by rhythm",
     xlab = "Time relative to QRS peak (s)", ylab = "mV",
     col = "white", ylim = c(-0.4, 1.2))

for(j in 1:n_valid){
  r <- kovariantes_valid$Rhythm[j]
  lines(beat_grid, beat_matrix[j,],
        col = adjustcolor(cols_map[r], alpha.f = 0.4), lwd = 1)
}
legend("topright", legend = names(cols_map),
       col = cols_map, lty = 1, lwd = 2, cex = 0.7, box.col = "white")


# Smooth the beat matrix into a new fd object for FPCA
basis_beat <- create.bspline.basis(rangeval = c(0, 0.45), nbasis = 100)
fdPar_beat <- fdPar(basis_beat, Lfdobj = int2Lfd(2), lambda = 1e-04)
smoothed_beats <- smooth.basis(beat_grid, t(beat_matrix), fdPar_beat)
fd_beats <- smoothed_beats$fd

# Quick mean per rhythm group to confirm registration worked
plot(mean.fd(fd_beats), main = "Overall mean beat",
     xlab = "Time (s)", ylab = "mV", lwd = 2)


#################
# Subsample for cleaner plots — every 5th point is enough visually
plot_idx   <- seq(1, 5000, by = 5)  # 1000 points
plot_time  <- time[plot_idx]
plot_raw   <- t(fdata_obj_scaled$data)[plot_idx, ]
plot_smooth <- eval.fd(plot_time, fd_final)

# --- Panel 1: Raw ---
graphics::matplot(plot_time, plot_raw, type = "l", lty = 1, lwd = 0.5,
        col = rainbow(99, alpha = 0.5),
        main = "Raw ECG — 100 patients",
        xlab = "Time (s)", ylab = "mV",
        ylim = c(-0.5, 1.5))

# --- Panel 2: After B-spline smoothing ---
graphics::matplot(plot_time, plot_smooth, type = "l", lty = 1, lwd = 0.5,
        col = rainbow(99, alpha = 0.5),
        main = "After B-spline smoothing\n(nbasis=500, λ=1e-04)",
        xlab = "Time (s)", ylab = "mV",
        ylim = c(-0.5, 1.5))

# --- Panel 3: After registration ---
graphics::matplot(
  beat_grid, t(beat_matrix),
  type = "l", lty = 1, lwd = 0.6,
  col = rainbow(n_valid, alpha = 0.5),
  main = "After registration\n(single beat, QRS aligned)",
  xlab = "Time relative to QRS (s)", ylab = "mV",
  ylim = c(-0.4, 1.2)
)


###################
# EDA on registered beats
##################
library(fda.usc)
library(fields)

# SECTION 1 — Mean ± SD envelope
# ------------------------------------------------------------------------
mean_reg   <- mean.fd(fd_beats)
stddev_reg <- std.fd(fd_beats)

plot(fd_beats,
     main = "Registered ECG beats",
     xlab = "Time (s)", ylab = "mV",
     col  = rgb(0, 0, 0, 0.15), lwd = 0.8)
abline(h = 0, col = "grey80", lty = 3)
lines(mean_reg,                col = 2, lwd = 3, lty = 1)
lines(mean_reg + 3*stddev_reg, col = 4, lwd = 2, lty = 1)
lines(mean_reg - 3*stddev_reg, col = 4, lwd = 2, lty = 1)
legend("topright", lty = c(1, 1), col = c(2, 4), lwd = c(3, 2),
       legend = c("Mean", "Mean ± 3×SD"),
       box.col = "white", cex = 0.85)

# SECTION 2 — Bivariate covariance function
# ------------------------------------------------------------------------
ecgvar.bifd <- var.fd(fd_beats)

t_sub      <- seq(0, 0.45, length.out = 100)
ecgvar_mat <- eval.bifd(t_sub, t_sub, ecgvar.bifd)

# 3D perspective surface
persp(t_sub, t_sub, ecgvar_mat,
      theta = -45, phi = 25, r = 3, expand = 0.5,
      ticktype = "detailed", col = "lightblue",
      xlab = "Time (s)", ylab = "Time (s)",
      zlab = "Covariance (mV²)",
      main = "Covariance surface — registered beats")

# Contour
contour(t_sub, t_sub, ecgvar_mat,
        xlab = "Time (s)", ylab = "Time (s)",
        main = "Covariance contour — registered beats")

# Heatmap with contour overlay
image.plot(t_sub, t_sub, ecgvar_mat,
           xlab = "Time (s)", ylab = "Time (s)",
           main = "Covariance heatmap — registered beats")
contour(t_sub, t_sub, ecgvar_mat,
        col = "white", add = TRUE)

# SECTION 3 — Derivatives: slope and acceleration
# ------------------------------------------------------------------------
deriv1 <- deriv.fd(fd_beats, 1)
deriv2 <- deriv.fd(fd_beats, 2)

opar <- par(mfrow = c(1, 2))

plot(deriv1,
     main = "First derivative (slope)",
     xlab = "Time (s)", ylab = "mV/s",
     col  = rgb(0, 0, 0, 0.15), lwd = 0.8)
lines(mean.fd(deriv1), col = "red", lwd = 3, lty = 2)
abline(h = 0, col = "grey70", lty = 3)
legend("topright", legend = "Mean slope",
       col = "red", lty = 2, lwd = 2,
       box.col = "white", cex = 0.8)

plot(deriv2,
     main = "Second derivative (acceleration)",
     xlab = "Time (s)", ylab = "mV/s²",
     col  = rgb(0, 0, 0, 0.15), lwd = 0.8)
lines(mean.fd(deriv2), col = "blue", lwd = 3, lty = 2)
abline(h = 0, col = "grey70", lty = 3)
legend("topright", legend = "Mean acceleration",
       col = "blue", lty = 2, lwd = 2,
       box.col = "white", cex = 0.8)

par(opar)

# SECTION 4 — Centrality measures by rhythm group
# ------------------------------------------------------------------------

# Use largest rhythm groups
rhythms_v <- kovariantes_valid$Rhythm
rhythm_tab <- sort(table(rhythms_v), decreasing = TRUE)
top_groups <- names(rhythm_tab[1:min(4, length(rhythm_tab))])
cat("Top rhythm groups:", top_groups, "\n")
print(rhythm_tab[top_groups])

get_fdata_group <- function(r) {
  fdata_reg[which(rhythms_v == r), ]
}

# Convert registered fd to fdata for fda.usc centrality/dispersion functions
eval_grid_beats <- seq(0, 0.45, length.out = 250)
reg_mat         <- t(eval.fd(eval_grid_beats, fd_beats))
fdata_reg       <- fdata(reg_mat, argvals = eval_grid_beats)

opar <- par(mfrow = c(2, 2))

for(r in top_groups){
  grp  <- get_fdata_group(r)
  n_r  <- sum(rhythms_v == r)
  
  plot(func.mean(grp),
       main = paste0("Centrality — ", r, " (n=", n_r, ")"),
       xlab = "Time (s)", ylab = "mV",
       ylim = c(-0.2, 0.6))
  abline(h = 0, col = "grey80", lty = 3)
  legend("topright", cex = 0.65, box.col = "white",
         lty = 1:5, col = 1:5,
         legend = c("mean", "trim.mode", "trim.RP",
                    "median.mode", "median.RP"))
  lines(func.trim.mode(grp, trim = 0.15), col = 2, lty = 2)
  lines(func.trim.RP(grp,   trim = 0.15), col = 3, lty = 3)
  lines(func.med.mode(grp,  trim = 0.15), col = 4, lty = 4)
  lines(func.med.RP(grp,    trim = 0.15), col = 5, lty = 5)
}

par(opar)
# SECTION 5 — Dispersion measures by rhythm group
# ------------------------------------------------------------------------

opar <- par(mfrow = c(2, 2))

for(r in top_groups){
  grp <- get_fdata_group(r)
  n_r <- sum(rhythms_v == r)
  
  plot(func.var(grp),
       main = paste0("Dispersion — ", r, " (n=", n_r, ")"),
       xlab = "Time (s)", ylab = "Variance (mV²)")
  legend("topright", cex = 0.65, box.col = "white",
         lty = 1:3, col = 1:3,
         legend = c("var", "trimvar.mode", "trimvar.RP"))
  lines(func.trimvar.mode(grp, trim = 0.15), col = 2, lty = 2)
  lines(func.trimvar.RP(grp,   trim = 0.15), col = 3, lty = 3)
}

par(opar)
# SECTION 6 — Group mean curves overlaid
# ------------------------------------------------------------------------

unique_r <- unique(rhythms_v)
cols_map <- setNames(rainbow(length(unique_r)), unique_r)

plot(mean.fd(fd_beats),
     main = "Mean ECG beat by rhythm group",
     xlab = "Time (s)", ylab = "mV",
     col = "black", lwd = 2, lty = 2,
     ylim = c(-0.2, 0.6))
abline(h = 0, col = "grey80", lty = 3)

for(r in unique_r){
  idx <- which(rhythms_v == r)
  if(length(idx) >= 2)
    lines(mean.fd(fd_beats[idx]), col = cols_map[r], lwd = 2)
}

legend("topright",
       legend = c("Overall", unique_r),
       col    = c("black", cols_map),
       lty    = c(2, rep(1, length(unique_r))),
       lwd = 2, cex = 0.65, box.col = "white")

# SECTION 7 - Functional Principal Component Analysis
# ------------------------------------------------------------------------

# Run FPCA from FD
pca_res <- pca.fd(fd_beats, nharm = 5)

# Variance explained
var_explained <- pca_res$varprop
print(var_explained)

# Scree plot
plot(var_explained, type = "b", pch = 16,
     xlab = "Principal Component",
     ylab = "Variance Explained",
     main = "Scree Plot — FPCA")

# Plot first 3 harmonics (eigenfunctions)
par(mfrow = c(1, 3))
plot(pca_res$harmonics[1], main = "PC1")
plot(pca_res$harmonics[2], main = "PC2")
plot(pca_res$harmonics[3], main = "PC3")
par(mfrow = c(1,1))

# Mean ± variation for PC1 and PC2
mean_fd <- mean.fd(fd_beats)

par(mfrow = c(1, 2))

# PC1 variation
plot(mean_fd, lwd = 2, main = "PC1 variation",
     xlab = "Time (s)", ylab = "mV")
lines(mean_fd + 2 * sqrt(pca_res$values[1]) * pca_res$harmonics[1],
      col = "blue", lty = 2)
lines(mean_fd - 2 * sqrt(pca_res$values[1]) * pca_res$harmonics[1],
      col = "red", lty = 2)
legend("topright", legend = c("Mean", "+2SD", "-2SD"),
       col = c("black", "blue", "red"), lty = c(1,2,2), cex = 0.8)

# PC2 variation
plot(mean_fd, lwd = 2, main = "PC2 variation",
     xlab = "Time (s)", ylab = "mV")
lines(mean_fd + 2 * sqrt(pca_res$values[2]) * pca_res$harmonics[2],
      col = "blue", lty = 2)
lines(mean_fd - 2 * sqrt(pca_res$values[2]) * pca_res$harmonics[2],
      col = "red", lty = 2)

par(mfrow = c(1,1))

# PCA scores
scores <- pca_res$scores

# Score plot
plot(scores[,1], scores[,2],
     xlab = "PC1", ylab = "PC2",
     main = "FPCA Scores",
     pch = 16, col = "blue")

# Score plot colored by rhythm (if available)
if (exists("rhythms_v")) {
  cols <- as.numeric(as.factor(rhythms_v))
  
  plot(scores[,1], scores[,2],
       col = cols, pch = 16,
       xlab = "PC1", ylab = "PC2",
       main = "FPCA Scores by Rhythm")
  
  legend("topright",
         legend = unique(rhythms_v),
         col = unique(cols),
         pch = 16, cex = 0.7)
}

# Rotated PCA plot - same as without rotation
#pdf("plots/PCA_rotated_plot.pdf", width = 8, height = 10)

pca_res <- pca.fd(fd_beats, nharm = 5, centerfns = TRUE)
pca_rot <- varmx.pca.fd(pca_res)

mean_fd <- mean.fd(fd_beats)
time_grid <- seq(0, 0.45, length.out = 250)  # same as beat_grid

par(mfrow = c(2,2), mar = c(4,4,3,1))

for(i in 1:4){
  harmon <- pca_rot$harmonics[i]
  sd_i   <- sqrt(pca_rot$values[i])
  var_pct <- round(100 * pca_rot$varprop[i], 1)
  
  # Evaluate numerically
  mean_vals <- eval.fd(time_grid, mean_fd)
  harm_vals <- eval.fd(time_grid, harmon)
  
  # Determine ylim
  ylim_vals <- range(c(mean_vals, 
                       mean_vals + 2*sd_i*harm_vals, 
                       mean_vals - 2*sd_i*harm_vals))
  
  # Plot
  plot(time_grid, mean_vals, type = "l", lwd = 2, ylim = ylim_vals,
       xlab = "Time (s)", ylab = "Amplitude (mV)",
       main = paste0("Rotated PCA function ", i,
                     " (Percentage of variability ", var_pct, ")"))
  
  # Add ±SD lines
  lines(time_grid, mean_vals + 1*sd_i*harm_vals, col = "blue", lty = 2)
  lines(time_grid, mean_vals - 1*sd_i*harm_vals, col = "blue", lty = 2)
  lines(time_grid, mean_vals + 2*sd_i*harm_vals, col = "red", lty = 3)
  lines(time_grid, mean_vals - 2*sd_i*harm_vals, col = "red", lty = 3)
}

for(i in 1:4){
  harmon <- pca_res$harmonics[i]
  sd_i   <- sqrt(pca_res$values[i])
  var_pct <- round(100 * pca_res$varprop[i], 1)
  
  # Evaluate numerically
  mean_vals <- eval.fd(time_grid, mean_fd)
  harm_vals <- eval.fd(time_grid, harmon)
  
  # Determine ylim
  ylim_vals <- range(c(mean_vals, 
                       mean_vals + 2*sd_i*harm_vals, 
                       mean_vals - 2*sd_i*harm_vals))
  
  # Plot
  plot(time_grid, mean_vals, type = "l", lwd = 2, ylim = ylim_vals,
       xlab = "Time (s)", ylab = "Amplitude (mV)",
       main = paste0("PCA function ", i,
                     " (Percentage of variability ", var_pct, ")"))
  
  # Add ±SD lines
  lines(time_grid, mean_vals + 1*sd_i*harm_vals, col = "blue", lty = 2)
  lines(time_grid, mean_vals - 1*sd_i*harm_vals, col = "blue", lty = 2)
  lines(time_grid, mean_vals + 2*sd_i*harm_vals, col = "red", lty = 3)
  lines(time_grid, mean_vals - 2*sd_i*harm_vals, col = "red", lty = 3)
}
par(mfrow = c(1,1))

#-------------------
# Hypothesis testing
#-------------------

kovariantes_valid <- kovariantes[
  match(df_valid$ID, kovariantes$FileName),
]
all(kovariantes_valid$FileName == df_valid$ID)
male_idx   <- kovariantes_valid$Gender == "MALE"
female_idx <- kovariantes_valid$Gender == "FEMALE"
fd_beats_m <- fd_beats[male_idx]
fd_beats_f <- fd_beats[female_idx]
plot(mean.fd(fd_beats),
     main = "Mean ECG beat by gender",
     xlab = "Time (s)", ylab = "mV",
     col = "black", lwd = 2, lty = 2,
     ylim = c(-0.1, 0.3))
abline(h = 0, col = "grey80", lty = 3)
lines(mean.fd(fd_beats_f), col = "red", lwd = 2)
lines(mean.fd(fd_beats_m), col = "blue", lwd = 2)
legend("topright",
       legend = c("Overall", "Female", "Male"),
       col = c("black", "red", "blue"),
       lty = c(2, 1, 1),
       lwd = 2, cex = 0.8)
time_grid <-seq(0, 0.45, length.out = 250)  # same as beat_grid
# TWO-SAMPLE - POINTWISE
source("Ztwosample.R")
# Do ECG for the first half second is the same among men and women
# H0: mu(male) = mu(female)
# H1: mu(male) != mu(female)
stat <- Ztwosample(x=fd_beats_m, y=fd_beats_f, t.seq = time_grid)
stat
source("L2stattwosample.R")
source("trace.R")
stat_l2.1 <- L2.stat.twosample(x=fd_beats_m, y=fd_beats_f, t.seq = time_grid, method=1)
stat_l2.1
stat_l2.2 <- L2.stat.twosample(x=fd_beats_m, y=fd_beats_f, t.seq = time_grid, method=2, replications=500)
stat_l2.1$pvalue
source("Fstattwosample.R")
stat_f.1 <- F.stat.twosample(x=fd_beats_m, y=fd_beats_f, t.seq = time_grid, method=1)
stat_f.1
stat_f.2 <- F.stat.twosample(x=fd_beats_m, y=fd_beats_f, t.seq = time_grid, method=2, replications=500)
stat_f.2$pvalue
stat_perm <- tperm.fd(fd_beats_m,fd_beats_f)
stat_perm
#nebuvo skirtumo visur p-value > 0.05
# 3 SAMPLE - RITMAS
sb_idx   <- kovariantes_valid$Rhythm == "SB"
atir_idx <- kovariantes_valid$Rhythm == "AFIB"
sr_idx <- kovariantes_valid$Rhythm == "SR"
fd_beats_sb <- fd_beats[sb_idx]
fd_beats_atir <- fd_beats[atir_idx]
fd_beats_sr <- fd_beats[sr_idx]
plot(mean.fd(fd_beats),
     main = "Mean ECG beat by rhythm groups",
     xlab = "Time (s)", ylab = "mV",
     col = "black", lwd = 2, lty = 2, ylim=c(-0.1,0.3))
abline(h = 0, col = "grey80", lty = 3)
lines(mean.fd(fd_beats_sb), col = "green", lwd = 2)
lines(mean.fd(fd_beats_atir), col = "blue", lwd = 2)
lines(mean.fd(fd_beats_sr), col = "red", lwd = 2)
legend("topright",
       legend = c("Overall", "SB", "AFIB", "SR"),
       col = c("black","green", "blue", "red"),
       lty = c(2, 1, 1),
       lwd = 2, cex = 0.8)
library(fdANOVA)
###########################################################
#### Pointwise fANOVA function                         ####
###########################################################
fANOVA.pointwise <- function(data, groups, t.seq, alpha=0.05) {
  # data is matrix with time in rows and variables in columns
  # group is a list names separating columns into different groups, a factor
  # time scale for measures
  n <- nrow(data)
  pvals <- numeric(n)
  lv <- levels(groups)
  k <- length(lv)
  mean.p <- matrix(NA, ncol=k, nrow=n)
  perm <- factorial(k)/(factorial(2)*(factorial(k-2)))
  Tukey.posthoc <- matrix(NA, ncol=perm, nrow=n)
  for(i in 1:n) {
    dt <- data.frame((data[i,]), groups)
    names(dt) <- c("values", "groups")
    av <- aov(values~groups, data = dt)
    pvals[i] <- summary(av)[[1]]["Pr(>F)"][1,1]
    mean.p[i,]  <- as.matrix((dt %>% group_by(groups) %>% summarise(mean(values)))[,2])
    colnames(Tukey.posthoc) <- rownames(TukeyHSD(av)$groups)
    Tukey.posthoc[i,] <- TukeyHSD(av)$groups[,4]
  }
  
  overall_mean <- apply(data, 1, mean)
  
  opar1 <- par(mfrow=c(2,1))
  
  plot(t.seq, pvals, type="l", main = "Pointwise ANOVA p-values",
       xlab = "Time", ylab="p-value", ylim=c(0,1))
  lines(t.seq, rep(0.05, n), col="blue", lty=2)
  
  mn <- min(mean.p, overall_mean)
  mx <- max(mean.p, overall_mean)
  
  plot(t.seq, overall_mean, type = "l", main = "Group means",
       xlab = "Time", ylab = "Mean", ylim = c(mn-0.05, mx+0.05))
  for(i in 1:k) {
    lines(t.seq, mean.p[,i], col=i+1, lty=i+1)
  }
  
  legend("topright", legend=c("Overall", lv), lty=1:(k+1), col=1:(k+1), title="Group")
  
  par(opar1)
  
  
  opar2 <- par(mfrow=c(1,1), ask = TRUE)
  
  for(i in 1:perm) {
    plot(t.seq, Tukey.posthoc[,i], type="l", main = paste("Tukey HSD p-values", rownames(TukeyHSD(av)$groups)[i]),
         xlab = "Time", ylab = "p-value", ylim = c(0,1))
    lines(t.seq, rep(0.05, n), col="blue", lty=2)
  }
  
  par(opar2)
  
  return(list(p.values=pvals, TukeyHSD=Tukey.posthoc, gr.means = mean.p, overal.mean=overall_mean))
}
# ---------------------------------------
# Funtion on scalar regression using pffr
# ---------------------------------------

fd_beats.eval <- eval.fd(time_grid, fd_beats)
fANOVA.pointwise(data=fd_beats.eval, groups=as.factor(kovariantes_valid$Rhythm),
                 t.seq=time_grid, alpha=0.05)
(fanova.fd <- fanova.tests(x = fd_beats.eval, group.label = as.factor(kovariantes_valid$Rhythm),
                           parallel = TRUE))

#-----------------------------
# Funtion on scalar regression using pffr
#-----------------------------

library(refund)

beat_mat <- t(eval.fd(beat_grid, fd_beats))

dat <- data.frame(
  Rhythm = factor(kovariantes_valid$Rhythm, levels = c("AFIB", "SR", "SB")),
  Age    = as.numeric(kovariantes_valid$PatientAge),
  HR     = as.numeric(kovariantes_valid$VentricularRate),
  RAxis  = as.numeric(kovariantes_valid$RAxis),
  TAxis  = as.numeric(kovariantes_valid$TAxis)
)
dat$ECG <- beat_mat

m1 <- pffr(ECG ~ Rhythm, yind = beat_grid, data = dat)
m2 <- pffr(ECG ~ Rhythm + Age + HR, yind = beat_grid, data = dat)
m3 <- pffr(ECG ~ Rhythm + Age + HR + RAxis + TAxis, yind = beat_grid, data = dat)

r2_table <- data.frame(
  Model = c("M1: Rhythm",
            "M2: + Age + HR",
            "M3: + RAxis + TAxis"),
  R2_adj  = c(summary(m1)$r.sq, summary(m2)$r.sq, summary(m3)$r.sq),
  DevExpl = c(summary(m1)$dev.expl, summary(m2)$dev.expl, summary(m3)$dev.expl)
)
print(r2_table, digits = 4)
fosr_full <- m3
summary(fosr_full)

par(mfrow = c(3, 3), mar = c(4, 4.5, 2.5, 1))

plot(fosr_full, select = 1, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[0](t)),
     main = "Intercept - baseline beat shape")

plot(fosr_full, select = 2, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[AFIB](t)),
     main = "AFIB (reference, shrunk to 0)")
abline(h = 0, lty = 3, col = "grey60")

plot(fosr_full, select = 3, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[SR](t)),
     main = "SR effect (vs AFIB)")
abline(h = 0, lty = 3, col = "grey60")

plot(fosr_full, select = 4, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[SB](t)),
     main = "SB effect (vs AFIB)")
abline(h = 0, lty = 3, col = "grey60")

plot(fosr_full, select = 5, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[Age](t)),
     main = "Age effect (per year)")
abline(h = 0, lty = 3, col = "grey60")

plot(fosr_full, select = 6, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[HR](t)),
     main = "Heart rate effect (per bpm)")
abline(h = 0, lty = 3, col = "grey60")

plot(fosr_full, select = 7, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[RAxis](t)),
     main = "R-axis effect (per degree)")
abline(h = 0, lty = 3, col = "grey60")

plot(fosr_full, select = 8, scale = 0,
     xlab = "Time since QRS (s)", ylab = expression(beta[TAxis](t)),
     main = "T-axis effect (per degree)")
abline(h = 0, lty = 3, col = "grey60")

par(mfrow = c(1, 1))

# Sensitivity analysis
fit_k5  <- pffr(ECG ~ Rhythm + Age + HR + RAxis + TAxis,
                yind = beat_grid, data = dat,
                bs.yindex = list(bs = "ps", k = 5,  m = c(2, 1)))

fit_k10 <- pffr(ECG ~ Rhythm + Age + HR + RAxis + TAxis,
                yind = beat_grid, data = dat,
                bs.yindex = list(bs = "ps", k = 10, m = c(2, 1)))

fit_k15 <- pffr(ECG ~ Rhythm + Age + HR + RAxis + TAxis,
                yind = beat_grid, data = dat,
                bs.yindex = list(bs = "ps", k = 15, m = c(2, 1)))

sensitivity <- data.frame(k = c(5, 10, 15), 
                          R2_adj = c(summary(fit_k5)$r.sq, summary(fit_k10)$r.sq,
                                     summary(fit_k15)$r.sq),
                          EDF_TAxis = c(summary(fit_k5)$s.table["TAxis(yindex)","edf"],
                                        summary(fit_k10)$s.table["TAxis(yindex)","edf"],
                                        summary(fit_k15)$s.table["TAxis(yindex)","edf"]),
                          EDF_SR = c(summary(fit_k5)$s.table["RhythmSR(yindex)","edf"],
                                     summary(fit_k10)$s.table["RhythmSR(yindex)","edf"],
                                     summary(fit_k15)$s.table["RhythmSR(yindex)","edf"]))
print(sensitivity, digits = 4)

