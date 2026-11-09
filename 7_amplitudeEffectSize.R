library(dplyr)
library(here)

rm(list = ls())
setwd(here())
set.seed(42)


###############################################################################
# INDIVIDUAL 30+15 AMPLITUDE EFFECT SIZE vs A SMOOTHED-PERMUTATION NULL
#
# Quantifies the rhythm at the level of EFFECT SIZE, without AIC / model
# selection. For each individual we fit the combined 30+15 cosinor and read off
# the amplitude of each component (30-day lunar, 15-day semilunar). We ask whether
# those amplitudes are larger than what the SAME smoothing produces from a series
# with no real lunar relationship.
#
# Pipeline (identical for observed and null):
#   sort by DayinCycle -> SMA(SMA_N) smoothing -> fit
#     onset_smoothed ~ Weekend + cos/sin(30) + cos/sin(15)   (raw OLS)
#   amp30 = sqrt(c30^2 + s30^2); amp15 = sqrt(c15^2 + s15^2)
#
# NULL: permute the RAW onset (breaking any onset-lunar relationship), then apply
# the IDENTICAL smoothing + fit. Because the null goes through the same smoothing,
# it captures the SPURIOUS amplitude that a moving average produces from noise
# (smoothing reddens the spectrum and creates non-zero 15/30-day amplitude even in
# random data). Comparing observed vs this null therefore ISOLATES the real signal
# from the smoothing artifact -- which is why smoothing is legitimate here, unlike
# in AIC model selection (see 6_modelSelection_null.R).
#
# Two aggregate readouts per component:
#   (a) fraction of individuals whose observed amplitude exceeds their own null
#       95th percentile (presence, smoothing-calibrated);
#   (b) effect size: mean observed amplitude vs the permutation null distribution
#       of the mean amplitude (ratio + p).
#
# Output: results/amplitude_effectSize.csv
###############################################################################

SMA_N      <- 3      # smoothing window. 3 preserves the 15-day band much better than 7
                     # (a 7-day MA strongly attenuates the semilunar component); set to
                     # 7 to compare. Window 15 or 30 would NULL the 15/30-day signal.
n_perm     <- 500    # permutations
min_nights <- 14


# --- Trailing moving average; valid from index SMA_N onward ---
sma <- function(y, n) {
  if (n <= 1) return(y)
  out <- rep(NA_real_, length(y))
  for (i in n:length(y)) out[i] <- mean(y[(i - n + 1):i])
  out
}

# --- Design: intercept (+ optional Weekend) + 30 and 15 day harmonics ---
buildAmpDesign <- function(t, wk = NULL) {
  cols <- list(`(Intercept)` = rep(1, length(t)))
  if (!is.null(wk)) cols$wk <- wk
  cols$c30 <- cos(2*pi*t/30); cols$s30 <- sin(2*pi*t/30)
  cols$c15 <- cos(2*pi*t/15); cols$s15 <- sin(2*pi*t/15)
  do.call(cbind, cols)
}

# --- Component amplitudes from one (possibly permuted) raw, sorted onset vector ---
ampOf <- function(yraw, prep) {
  ys <- sma(yraw, SMA_N)[prep$keep]
  b  <- lm.fit(prep$X, ys)$coefficients
  c(amp30 = sqrt(b["c30"]^2 + b["s30"]^2),
    amp15 = sqrt(b["c15"]^2 + b["s15"]^2))
}

# --- Per-individual preparation (design fixed under permutation; only y permutes) ---
prepIndiv <- function(t, y, wk) {
  o    <- order(t)
  t <- t[o]; y <- y[o]; if (!is.null(wk)) wk <- wk[o]
  keep <- SMA_N:length(y)
  list(yraw = y, keep = keep,
       X = buildAmpDesign(t[keep], if (is.null(wk)) NULL else wk[keep]))
}


# --- Observed + smoothed-permutation null for one population ---
runPopulation <- function(df, id_col, t_col, y_col, label, use_weekday = TRUE) {
  df     <- df[!is.na(df[[y_col]]) & !is.na(df[[t_col]]), ]
  has_wk <- use_weekday && "typeOfDay" %in% names(df)

  cache <- list()
  for (id in unique(df[[id_col]])) {
    sub <- df[df[[id_col]] == id, ]
    if (nrow(sub) < min_nights) next
    t  <- as.numeric(sub[[t_col]]); y <- as.numeric(sub[[y_col]])
    wk <- if (has_wk) as.numeric(sub$typeOfDay == "Weekend") else NULL
    cache[[as.character(id)]] <- prepIndiv(t, y, wk)
  }
  n_ind <- length(cache)
  if (n_ind < 3) { cat(label, ": too few individuals\n"); return(NULL) }

  obs <- vapply(cache, function(e) ampOf(e$yraw, e), numeric(2))   # [2 x n_ind] (amp30, amp15)
  null30 <- matrix(NA_real_, n_perm, n_ind)
  null15 <- matrix(NA_real_, n_perm, n_ind)
  for (p in seq_len(n_perm)) {
    a <- vapply(cache, function(e) ampOf(sample(e$yraw), e), numeric(2))
    null30[p, ] <- a[1, ]; null15[p, ] <- a[2, ]
  }

  out <- list()
  for (comp in c("30", "15")) {
    obs_v  <- if (comp == "30") obs[1, ] else obs[2, ]
    null_m <- if (comp == "30") null30 else null15
    # (a) per-individual: observed > own null 95th percentile
    thr95   <- apply(null_m, 2, quantile, 0.95, na.rm = TRUE)
    frac_sig <- mean(obs_v > thr95)
    # (b) effect size on the cohort mean amplitude
    obs_mean   <- mean(obs_v)
    null_means <- rowMeans(null_m, na.rm = TRUE)
    p_agg      <- mean(null_means >= obs_mean)
    ratio      <- obs_mean / mean(null_means)
    cat(sprintf("  %-7s comp-%-2s | obs mean amp=%.3f  null mean=%.3f  ratio=%.2fx  p=%.3f | %2.0f%% individuals > null95\n",
                label, comp, obs_mean, mean(null_means), ratio, p_agg, 100 * frac_sig))
    out[[comp]] <- data.frame(population = label, component = paste0(comp, "-day"),
                              n = n_ind, obs_mean_amp = obs_mean,
                              null_mean_amp = mean(null_means), effect_ratio = ratio,
                              p_aggregate = p_agg, frac_indiv_above_null95 = frac_sig)
  }
  do.call(rbind, out)
}


# --- Data ---
seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", stringsAsFactors = FALSE)
seattle <- seattle[seattle$year %in% c("SS2020", "AW2020", "SS2024"), ]

tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", stringsAsFactors = FALSE)
tobas <- tobas[tobas$year != "Spring2017", ]

titi <- read.csv("data/Data_S3_all_onsetData_titi.csv", stringsAsFactors = FALSE)

cat("Amplitude effect size vs smoothed-permutation null (SMA=", SMA_N,
    ", ", n_perm, " permutations)\n\n", sep = "")

res <- rbind(
  runPopulation(seattle, "ID_year", "DayinCycle", "Start_Time", "Seattle", use_weekday = TRUE),
  runPopulation(tobas,   "ID_year", "DayinCycle", "Start_Time", "Tobas",   use_weekday = TRUE),
  runPopulation(titi,    "id2",     "DayinCycle", "time",       "Titi",    use_weekday = FALSE)
)

write.csv(res, "results/amplitude_effectSize.csv", row.names = FALSE)
cat("\nSaved: results/amplitude_effectSize.csv\n")
