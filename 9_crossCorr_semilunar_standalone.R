###############################################################################
# 9_crossCorr_semilunar_standalone.R   
#
# Tests for a ~15-day (semilunar) periodicity IN THE CROSS-CORRELOGRAM between
# sleep onset and lunisolar gravity, as requested by the Reviewer.
#
# Rationale: if both sleep onset and gravity carry a ~15-day component, the
# sleep-vs-gravity cross-correlation function (CCF) should itself oscillate with
# a ~15-day period (a secondary peak ~15 lags from the main one). We quantify
# this by fitting a 30 + 15-day harmonic (cosinor) model TO THE CROSS-
# CORRELOGRAM and comparing the 15-day amplitude (A15) with the 30-day
# amplitude (A30). Significance of A15 is assessed with a circular-shift null.
#
# Two complementary views:
#   (1) Cohort-folded CCF  : CCF of the by-night-in-cycle mean profiles
#                            (same object as in the main cross-correlation
#                            figures; high SNR), lags -15..14.
#   (2) Individual-level   : CCF of each participant's raw daily series over a
#                            wide lag window (+/- MAXLAG), averaged across
#                            individuals, with a circular-shift null.
#
# Expected result (see Python verification, June 2026):
#   Toba/Qom: clear ~15-day component in the CCF (A15 ~ A30).
#   Seattle : CCF dominated by the 30-day component (A15 << A30),
#             consistent with the weaker 15-day gravimetric forcing there.
#
# Dependencies: base R + dplyr + lubridate.  Run independently.
###############################################################################

library(dplyr)
library(lubridate)
set.seed(42)

MAXLAG <- 22                       # individual-level lag window (days)
SMOOTH <- TRUE                     # apply 7-day SMA to the sleep series at the
                                   # individual level (TRUE = consistent with Fig. 4;
                                   # set FALSE for the raw-series version)
GRAV   <- "data/Data_S4_meanDaily_Gravitational_variation.csv"

# 7-day centered simple moving average (same smoothing used for Fig. 4)
sma7 <- function(x, n = 7) as.numeric(stats::filter(x, rep(1/n, n), sides = 2))

## ---- load gravity (daily series per location) -------------------------------
g <- read.csv(GRAV)
g$date <- as.Date(sprintf("%04d-%02d-%02d", g$Year, g$Month, g$Day))
# base R (avoids dplyr::select being masked by other packages)
prep_grav <- function(g, cityname) {
  sub <- g[g$city == cityname, c("date", "mean_Total_C")]
  sub <- sub[!duplicated(sub$date), ]
  names(sub)[2] <- "grav"
  sub
}
gravS <- prep_grav(g, "Seattle")
gravT <- prep_grav(g, "Ingeniero Juarez")

## ---- load sleep -------------------------------------------------------------
seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv")
seattle <- seattle[seattle$year %in% c("SS2020","AW2020","SS2024"), ]
seattle$date <- as.Date(seattle$night)

tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv")
tobas <- tobas[tobas$year != "Spring2017", ]
tobas$date <- as.Date(tobas$night, format = "%m/%d/%Y")

## ---- helper: normalized CCF at integer lags ---------------------------------
ccf_at <- function(a, b, lags) {
  out <- numeric(length(lags))
  for (i in seq_along(lags)) {
    L <- lags[i]
    if (L > 0)      { x <- a[(L+1):length(a)];        y <- b[1:(length(b)-L)] }
    else if (L < 0) { x <- a[1:(length(a)+L)];        y <- b[(1-L):length(b)] }
    else            { x <- a;                          y <- b }
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 15) { out[i] <- NA; next }
    xx <- x[ok] - mean(x[ok]); yy <- y[ok] - mean(y[ok])
    d  <- sqrt(sum(xx^2) * sum(yy^2))
    out[i] <- if (d > 0) sum(xx * yy) / d else NA
  }
  out
}

## ---- helper: fit 30+15-day cosinor to a CCF curve ---------------------------
# Returns A30, A15, ratio = A15/A30, and the lm object.
harmonic_fit <- function(lag, ccf) {
  d <- data.frame(lag = lag, ccf = ccf); d <- d[is.finite(d$ccf), ]
  w <- 2*pi/30; w2 <- 2*pi/15
  m <- lm(ccf ~ cos(w*lag) + sin(w*lag) + cos(w2*lag) + sin(w2*lag), data = d)
  cf <- coef(m)
  A30 <- sqrt(cf[2]^2 + cf[3]^2); A15 <- sqrt(cf[4]^2 + cf[5]^2)
  list(A30 = unname(A30), A15 = unname(A15), ratio = unname(A15/A30), model = m)
}

###############################################################################
## (1) COHORT-FOLDED CCF + harmonic decomposition
###############################################################################
cat("\n================ COHORT-FOLDED CROSS-CORRELOGRAM ================\n")
lags_fold <- -15:14
summaryTab <- data.frame()         # collects A30/A15/ratio/p for the output table

N_FOLD_SHUF <- 1000                # permutations for the cohort-level A15 null

# circular CCF of two centered 30-pt profiles -> reordered to lags -15..14
circ_ccf <- function(sp, gp) {
  d <- sqrt(sum(sp^2) * sum(gp^2))
  c0 <- sapply(0:29, function(L) sum(sp * gp[((seq_len(30)-1-L) %% 30)+1]) / d)
  c(c0[16:30], c0[1:15])
}

# per-cohort centered by-night-in-cycle profiles (sleep + gravity)
folded_profiles <- function(df, grav, var = "Start_Time") {
  df <- merge(df, grav, by = "date", all.x = TRUE)
  prof <- list()
  for (c in unique(df$year)) {
    gc <- df[df$year == c, ]
    sp <- tapply(gc[[var]], gc$DayinCycle, mean, na.rm = TRUE)
    gp <- tapply(gc$grav,   gc$DayinCycle, mean, na.rm = TRUE)
    sp <- approx(as.numeric(names(sp)), sp, 0:29, rule = 2)$y
    gp <- approx(as.numeric(names(gp)), gp, 0:29, rule = 2)$y
    prof[[as.character(c)]] <- list(sp = sp - mean(sp), gp = gp - mean(gp))
  }
  prof
}

for (pop in c("Seattle","Tobas")) {
  df <- if (pop == "Seattle") seattle else tobas
  gr <- if (pop == "Seattle") gravS  else gravT
  prof <- folded_profiles(df, gr)
  ccfs <- t(sapply(prof, function(p) circ_ccf(p$sp, p$gp)))
  m  <- colMeans(ccfs, na.rm = TRUE)
  hf <- harmonic_fit(lags_fold, m)

  # NULL: circularly shift each cohort's sleep profile, recompute the mean CCF,
  # refit, and collect A30 and A15. p = P(null A >= observed A). Shifting the
  # profile destroys the sleep-gravity phase alignment but preserves the
  # sleep waveform's own amplitudes.
  a30_null <- numeric(N_FOLD_SHUF); a15_null <- numeric(N_FOLD_SHUF)
  for (k in seq_len(N_FOLD_SHUF)) {
    cc_k <- t(sapply(prof, function(p) {
      s <- sample(1:29, 1)
      sp_sh <- p$sp[((seq_len(30) - 1 + s) %% 30) + 1]
      circ_ccf(sp_sh, p$gp)
    }))
    h <- harmonic_fit(lags_fold, colMeans(cc_k, na.rm = TRUE))
    a30_null[k] <- h$A30; a15_null[k] <- h$A15
  }
  p30 <- (1 + sum(a30_null >= hf$A30)) / (1 + N_FOLD_SHUF)
  p15 <- (1 + sum(a15_null >= hf$A15)) / (1 + N_FOLD_SHUF)

  cat(sprintf("%-8s cohorts=%d | A30=%.3f p[A30]=%.3f | A15=%.3f p[A15]=%.3f | A15/A30=%.2f\n",
              pop, length(prof), hf$A30, p30, hf$A15, p15, hf$ratio))
  assign(paste0("folded_", pop), list(lag = lags_fold, mean = m, hf = hf, p30 = p30, p15 = p15))
  summaryTab <- rbind(summaryTab, data.frame(
    level = "cohort-folded", population = pop, n_units = length(prof),
    A30 = hf$A30, A15 = hf$A15, ratio = hf$ratio, p_A30 = p30, p_A15 = p15))
}

###############################################################################
## (2) INDIVIDUAL-LEVEL extended-lag CCF + phase-coherence null (A30 and A15)
###############################################################################
cat("\n================ INDIVIDUAL-LEVEL CROSS-CORRELOGRAM ================\n")
lags_ind <- -MAXLAG:MAXLAG
N_SHUF   <- 200

# Collect each participant's aligned daily series (raw + gravity).
indiv_series <- function(df, grav, var = "Start_Time") {
  df <- aggregate(df[[var]], by = list(ID_year = df$ID_year, date = df$date),
                  FUN = function(z) mean(z, na.rm = TRUE))
  names(df)[3] <- "v"
  df <- merge(df, grav, by = "date", all.x = TRUE)
  series <- list()
  for (idv in unique(df$ID_year)) {
    gc <- df[df$ID_year == idv, ]; gc <- gc[order(gc$date), ]
    full <- seq(min(gc$date), max(gc$date), by = "day")
    yraw <- approx(gc$date, gc$v,    full, rule = 1)$y
    gr   <- approx(gc$date, gc$grav, full, rule = 1)$y
    chk  <- if (SMOOTH) sma7(yraw) else yraw
    if (sum(is.finite(chk)) < 24) next
    series[[idv]] <- list(yraw = yraw, gr = gr)
  }
  series
}

# Mean cross-correlogram across individuals. If shift = TRUE, each individual's
# series is circularly shifted by an INDEPENDENT random amount before its CCF
# is computed: this is the phase-coherence null. It preserves every series' own
# autocorrelation/spectral content (so each individual's A30/A15 are unchanged)
# but destroys the across-individual phase alignment, so the harmonic amplitudes
# of the MEAN CCF collapse unless individuals share a consistent phase.
mean_ccf_indiv <- function(series, shift = FALSE) {
  mats <- sapply(series, function(p) {
    y <- p$yraw
    if (shift) { n <- length(y); y <- y[((seq_len(n) - 1 + sample(5:(n-5), 1)) %% n) + 1] }
    if (SMOOTH) y <- sma7(y)
    ccf_at(y, p$gr, lags_ind)
  })
  rowMeans(mats, na.rm = TRUE)
}

for (pop in c("Seattle","Tobas")) {
  df <- if (pop == "Seattle") seattle else tobas
  gr <- if (pop == "Seattle") gravS  else gravT
  series <- indiv_series(df, gr)
  obs <- mean_ccf_indiv(series, shift = FALSE)
  hf  <- harmonic_fit(lags_ind, obs)

  # phase-coherence null on the MEAN CCF, for BOTH the 30- and 15-day components
  a30_null <- numeric(N_SHUF); a15_null <- numeric(N_SHUF)
  for (k in seq_len(N_SHUF)) {
    h <- harmonic_fit(lags_ind, mean_ccf_indiv(series, shift = TRUE))
    a30_null[k] <- h$A30; a15_null[k] <- h$A15
  }
  p30 <- (1 + sum(a30_null >= hf$A30)) / (1 + N_SHUF)
  p15 <- (1 + sum(a15_null >= hf$A15)) / (1 + N_SHUF)

  cat(sprintf("%-8s n=%d | A30=%.3f p[A30]=%.3f | A15=%.3f p[A15]=%.3f | A15/A30=%.2f\n",
              pop, length(series), hf$A30, p30, hf$A15, p15, hf$ratio))
  assign(paste0("indiv_", pop), list(lag = lags_ind, mean = obs, hf = hf,
                                     p30 = p30, p15 = p15, n = length(series)))
  summaryTab <- rbind(summaryTab, data.frame(
    level = sprintf("individual (%s)", if (SMOOTH) "SMA7" else "raw"),
    population = pop, n_units = length(series),
    A30 = hf$A30, A15 = hf$A15, ratio = hf$ratio, p_A30 = p30, p_A15 = p15))
}

## ---- write the summary table (cohort + individual) --------------------------
dir.create("results", showWarnings = FALSE)
write.csv(summaryTab, "results/CrossCorr_semilunar_harmonics.csv", row.names = FALSE)
cat("\nSaved results/CrossCorr_semilunar_harmonics.csv\n")
print(summaryTab)

###############################################################################
## PLOTS
###############################################################################
dir.create("Figures", showWarnings = FALSE)

png("Figures/CrossCorr_semilunar_folded.png", width = 1700, height = 750, res = 160)
par(mfrow = c(1,2), mar = c(4,4,3,1))
for (pop in c("Seattle","Tobas")) {
  o <- get(paste0("folded_", pop)); col <- if (pop=="Seattle") "#D55E00" else "#0072B2"
  plot(o$lag, o$mean, type = "l", lwd = 2, col = col,
       xlab = "lag (days)", ylab = "cohort-folded CCF (onset vs gravity)",
       main = sprintf("%s  A15/A30=%.2f  p[A30]=%.3f  p[A15]=%.3f",
                      pop, o$hf$ratio, o$p30, o$p15))
  abline(h = 0, col = "gray50"); abline(v = c(-15,0,15), lty = c(2,3,2), col = c("green","black","green"))
  lines(o$lag, predict(o$hf$model), col = "black", lty = 2)   # 30+15 fit
}
dev.off(); cat("\nSaved Figures/CrossCorr_semilunar_folded.png\n")

png("Figures/CrossCorr_semilunar_individual.png", width = 1700, height = 750, res = 160)
par(mfrow = c(1,2), mar = c(4,4,3,1))
for (pop in c("Seattle","Tobas")) {
  o <- get(paste0("indiv_", pop)); col <- if (pop=="Seattle") "#D55E00" else "#0072B2"
  plot(o$lag, o$mean, type = "l", lwd = 2, col = col,
       xlab = "lag (days)", ylab = "mean individual CCF (onset vs gravity)",
       main = sprintf("%s  A15/A30=%.2f  p[A30]=%.3f  p[A15]=%.3f",
                      pop, o$hf$ratio, o$p30, o$p15))
  abline(h = 0, col = "gray50"); abline(v = c(-15,0,15), lty = c(2,3,2), col = c("green","black","green"))
}
dev.off(); cat("Saved Figures/CrossCorr_semilunar_individual.png\n")

cat("\n=== Done. See results/CrossCorr_semilunar_harmonics.csv for A30/A15, A15/A30 and p-values. ===\n")
