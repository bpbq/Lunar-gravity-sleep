###############################################################################
# 5b_spectral_LombScargle_standalone.R   (STANDALONE - safe to delete)
#
# Independent, data-driven spectral check for the ~30 and ~15 day periodicities
# in sleep onset (Start_Time), at the INDIVIDUAL level, WITHOUT assuming the
# periods a priori (no 30-day folding).
#
# Method:
#   - Lomb-Scargle periodogram per individual (handles gaps / uneven sampling,
#     the correct tool for short irregular biological series).
#   - Linear detrend; optional 7-day moving average (SMA7) to suppress the strong
#     weekly/weekend social rhythm. IMPORTANT: when SMA7 is used, the AR(1)
#     red-noise surrogates pass through the SAME filter, so significance is not
#     an artifact of spectral reddening.
#   - Spectra interpolated onto a common period grid and averaged across
#     individuals; significance = exceedance of an AR(1) red-noise null.
#
# WHAT THIS SCRIPT SHOWS (see Claude verification, June 2026):
#   At the individual level these short series (25-80 nights) are UNDERPOWERED:
#   with a fair (identically-filtered) red-noise null, the 30 and 15 day peaks
#   sit at ~1.0x the 95% threshold, i.e. essentially AT the noise level. The
#   raw (unsmoothed) spectra are flat. This is expected: frequency resolution is
#   ~1/T, and a strong 7-day social rhythm competes with the lunar bands. The
#   appropriately-powered evidence for 30 & 15 day periodicity is the
#   cosinor + AIC model comparison and the cohort-level analyses, not this.
#
# Run independently of the rest of the pipeline. No external packages beyond
# base R + dplyr + lubridate (Lomb-Scargle implemented inline).
###############################################################################

library(dplyr)
library(lubridate)
set.seed(42)

## ---- Lomb-Scargle periodogram (normalized), base R --------------------------
lomb_scargle <- function(t, y, freqs) {
  y <- y - mean(y)
  vary <- sum(y^2) / length(y)
  if (vary < 1e-12) return(rep(NA_real_, length(freqs)))
  P <- numeric(length(freqs))
  for (k in seq_along(freqs)) {
    w <- 2 * pi * freqs[k]
    tau <- atan2(sum(sin(2 * w * t)), sum(cos(2 * w * t))) / (2 * w)
    wt <- w * (t - tau)
    c_ <- cos(wt); s_ <- sin(wt)
    P[k] <- 0.5 * ((sum(y * c_)^2) / sum(c_^2) + (sum(y * s_)^2) / sum(s_^2)) / vary
  }
  P
}

sma7 <- function(x, n = 7) as.numeric(stats::filter(x, rep(1/n, n), sides = 2))

## ---- pipeline: detrend (-> SMA7) -> Lomb-Scargle ----------------------------
pipeline <- function(tdays, y, freqs, smooth) {
  ok <- !is.na(y); tdays <- tdays[ok]; y <- y[ok]
  if (length(y) < 24) return(NULL)
  yd <- residuals(lm(y ~ tdays))
  td <- tdays
  if (smooth) {
    ys <- sma7(yd); keep <- !is.na(ys); td <- td[keep]; yd <- ys[keep]
    if (length(yd) < 14) return(NULL)
  }
  if (sd(yd) < 1e-9) return(NULL)
  lomb_scargle(td, yd, freqs)
}

ar1_surrogate <- function(yd) {
  n <- length(yd)
  phi <- acf(yd, lag.max = 1, plot = FALSE)$acf[2]
  if (is.na(phi)) phi <- 0
  phi <- max(min(phi, 0.99), -0.99)
  sig <- sqrt(max(var(yd) * (1 - phi^2), .Machine$double.eps))
  as.numeric(arima.sim(list(ar = phi), n = n, sd = sig))
}

period_grid <- seq(5, 65, by = 0.5)
freqs       <- 1 / period_grid
N_SURR      <- 200

analyze <- function(df, var = "Start_Time", smooth = FALSE) {
  specs <- list(); nulls <- list()
  for (idv in unique(df$ID_year)) {
    g <- df[df$ID_year == idv, ]
    g <- g[order(g$date), ]
    tdays <- as.numeric(g$date - min(g$date))
    y <- g[[var]]
    keep <- !duplicated(tdays); tdays <- tdays[keep]; y <- y[keep]
    obs <- pipeline(tdays, y, freqs, smooth)
    if (is.null(obs)) next
    specs[[idv]] <- obs
    # AR(1) fit on the RAW detrended series; surrogate passes through SAME filter
    ok <- !is.na(y); tr <- tdays[ok]; yr <- y[ok]
    yraw_dt <- residuals(lm(yr ~ tr))
    ns <- matrix(NA_real_, N_SURR, length(freqs))
    for (k in seq_len(N_SURR)) {
      sur <- ar1_surrogate(yraw_dt); td <- tr; s_ <- sur
      if (smooth) { s2 <- sma7(s_); kk <- !is.na(s2); td <- td[kk]; s_ <- s2[kk]
                    if (length(s_) < 14) next }
      ns[k, ] <- lomb_scargle(td, s_, freqs)
    }
    nulls[[idv]] <- ns
  }
  spec_mat   <- do.call(rbind, specs)
  mean_power <- colMeans(spec_mat, na.rm = TRUE)
  null_means <- matrix(NA_real_, N_SURR, length(freqs))
  ids <- names(nulls)
  for (k in seq_len(N_SURR))
    null_means[k, ] <- rowMeans(sapply(ids, function(i) nulls[[i]][k, ]), na.rm = TRUE)
  thr95 <- apply(null_means, 2, quantile, 0.95, na.rm = TRUE)
  thr99 <- apply(null_means, 2, quantile, 0.99, na.rm = TRUE)
  list(n = length(specs), period = period_grid,
       power = mean_power, thr95 = thr95, thr99 = thr99)
}

report <- function(tag, r) {
  cat(sprintf("\n=== %s | n=%d individuals ===\n", tag, r$n))
  for (P in c(30, 15)) {
    i <- which.min(abs(r$period - P))
    sig <- if (r$power[i] > r$thr99[i]) "SIG@99" else
           if (r$power[i] > r$thr95[i]) "SIG@95" else "n.s."
    cat(sprintf("  %dd: power=%.4f null95=%.4f ratio=%.2f %s\n",
                P, r$power[i], r$thr95[i], r$power[i]/r$thr95[i], sig))
  }
}

## ---- Load data (same cohorts as the main analysis) --------------------------
seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv")
seattle <- seattle[seattle$year %in% c("SS2020","AW2020","SS2024"), ]
seattle$date <- as.Date(seattle$night)                       # ISO yyyy-mm-dd

tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv")
tobas <- tobas[tobas$year != "Spring2017", ]
tobas$date <- as.Date(tobas$night, format = "%m/%d/%Y")

## ---- Run: raw vs 7-day-smoothed, both populations ---------------------------
res <- list()
for (nm in c("Seattle","Tobas")) {
  df <- if (nm == "Seattle") seattle else tobas
  for (sm in c(FALSE, TRUE)) {
    r <- analyze(df, "Start_Time", smooth = sm)
    res[[paste(nm, ifelse(sm,"SMA7","raw"))]] <- r
    report(paste(nm, ifelse(sm,"(7-day smoothed)","(raw)")), r)
  }
}

## ---- Plot -------------------------------------------------------------------
png("Figures/Individual_LombScargle_spectra.png", width = 1700, height = 1100, res = 160)
op <- par(mfrow = c(2,2), mar = c(4,4,3,1))
for (key in names(res)) {
  r <- res[[key]]; m <- r$period >= 5 & r$period <= 45
  plot(r$period[m], r$power[m], type = "l", lwd = 2,
       col = ifelse(grepl("Seattle", key), "#D55E00", "#0072B2"),
       xlab = "Period (days)", ylab = "mean LS power",
       main = sprintf("%s (n=%d)", key, r$n))
  lines(r$period[m], r$thr95[m], lty = 2, col = "gray40")
  abline(v = c(30,15), lty = c(2,3))
}
par(op); dev.off()
cat("\nSaved Figures/Individual_LombScargle_spectra.png\n")
cat("\nNOTE: individual short series are underpowered; see header comment.\n")
