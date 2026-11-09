library(dplyr)
library(here)

rm(list = ls())
setwd(here())
set.seed(42)


###############################################################################
# INDIVIDUAL COSINOR MODEL SELECTION vs A NO-RHYTHM NULL
#
# Demonstrates that the periodic coupling with the lunar cycle is identifiable
# beyond chance.
#
# Per individual we fit FOUR models of onset as a function of DayinCycle (raw
# data; the weekly social cycle is controlled with a Weekend covariate, the SAME
# way the cohort models in 2_analysis.R do it -- NOT by smoothing, which biases
# AIC model selection, see note at the bottom):
#   null  : onset ~ Weekend                              (no lunar rhythm)
#   30    : onset ~ Weekend + cos/sin(30)                (lunar)
#   15    : onset ~ Weekend + cos/sin(15)                (semilunar)
#   30+15 : onset ~ Weekend + cos/sin(30) + cos/sin(15)  (combined)
# and select the model with the lowest AIC (and AICc).
#
# Two statistics are reported:
#   (a) the selection breakdown (% null / 30 / 15 / 30+15);
#   (b) the fraction of individuals selecting ANY rhythm model (i.e. not null) =
#       "has a lunar/semilunar rhythm".
#
# NULL: within each individual we randomly PERMUTE the onset values (breaking any
# onset-lunar-phase relationship), refit and re-select. This gives the chance /
# overfitting rate at which a rhythm model beats null. p = proportion of null
# "any-rhythm" fractions >= the observed one. (Circular shift would NOT work:
# it only rotates phase and leaves a sinusoid as the same sinusoid; circular
# shift is the correct null for the cross-correlation of two signals, see
# 3_crossCor, not for the periodicity of a single signal.)
#
# Output: results/modelSelection_null.csv
###############################################################################

n_perm     <- 500    # permutations
min_nights <- 14     # need enough nights for a 6-parameter (30+15 + Weekend) model


# --- Design matrices for the four models (intercept + optional Weekend) ---
buildModels <- function(t, wk = NULL) {
  base <- if (is.null(wk)) cbind(`(Intercept)` = rep(1, length(t)))
          else             cbind(`(Intercept)` = 1, wk = wk)
  list(
    null    = base,
    `30`    = cbind(base, c30 = cos(2*pi*t/30), s30 = sin(2*pi*t/30)),
    `15`    = cbind(base, c15 = cos(2*pi*t/15), s15 = sin(2*pi*t/15)),
    `30+15` = cbind(base, c30 = cos(2*pi*t/30), s30 = sin(2*pi*t/30),
                          c15 = cos(2*pi*t/15), s15 = sin(2*pi*t/15))
  )
}

# Relative (A)ICc from residual sum of squares (Gaussian OLS); constants drop out
ic_value <- function(rss, n, p, crit) {
  K <- p + 1                                        # +1 for the variance
  a <- n * log(max(rss, 1e-12) / n) + 2 * K
  if (crit == "aicc") a <- a + (2 * K * (K + 1)) / max(n - K - 1, 1)
  a
}

selectModel <- function(y, designs, crit) {
  n <- length(y)
  a <- vapply(designs, function(X) ic_value(sum(lm.fit(X, y)$residuals^2), n, ncol(X), crit),
              numeric(1))
  names(a)[which.min(a)]
}


# --- Observed selection + permutation null for one population ---
runPopulation <- function(df, id_col, t_col, y_col, label, use_weekday = TRUE) {
  df     <- df[!is.na(df[[y_col]]) & !is.na(df[[t_col]]), ]
  has_wk <- use_weekday && "typeOfDay" %in% names(df)

  cache <- list()                                   # per individual: y + 4 designs (fixed under permutation)
  for (id in unique(df[[id_col]])) {
    sub <- df[df[[id_col]] == id, ]
    if (nrow(sub) < min_nights) next
    t  <- as.numeric(sub[[t_col]]); y <- as.numeric(sub[[y_col]])
    wk <- if (has_wk) as.numeric(sub$typeOfDay == "Weekend") else NULL
    cache[[as.character(id)]] <- list(y = y, designs = buildModels(t, wk))
  }
  n_ind <- length(cache)
  cat(label, " (n=", n_ind, ", weekday covariate=", has_wk, "):\n", sep = "")

  out <- list()
  for (crit in c("aic", "aicc")) {
    obs     <- vapply(cache, function(e) selectModel(e$y, e$designs, crit), character(1))
    f       <- table(factor(obs, levels = c("null", "30", "15", "30+15"))) / n_ind
    any_obs <- mean(obs != "null")

    f_null <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
      sel <- vapply(cache, function(e) selectModel(sample(e$y), e$designs, crit), character(1))
      f_null[p] <- mean(sel != "null")
    }
    pval <- mean(f_null >= any_obs)

    # "selected:" (4 categories, sum to 100%) vs "ANY rhythm" (obs = 100% - null
    # category; chance = permutation expectation, a different quantity)
    cat(sprintf("  %-4s | selected: null=%2.0f%% 30=%2.0f%% 15=%2.0f%% 30+15=%2.0f%% | ANY rhythm: obs=%2.0f%% chance=%2.0f%% (95%%=%2.0f%%) p=%.3f\n",
                crit, 100*f["null"], 100*f["30"], 100*f["15"], 100*f["30+15"],
                100*any_obs, 100*mean(f_null), 100*quantile(f_null, 0.95), pval))

    out[[crit]] <- data.frame(
      population = label, criterion = crit, n = n_ind,
      frac_null = as.numeric(f["null"]), frac_30 = as.numeric(f["30"]),
      frac_15 = as.numeric(f["15"]), frac_3015 = as.numeric(f["30+15"]),
      any_rhythm_obs = any_obs, any_rhythm_null_mean = mean(f_null),
      any_rhythm_null_q95 = as.numeric(quantile(f_null, 0.95)), p_any_rhythm = pval)
  }
  do.call(rbind, out)
}


# --- Data ---
seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", stringsAsFactors = FALSE)
seattle <- seattle[seattle$year %in% c("SS2020", "AW2020", "SS2024"), ]

tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", stringsAsFactors = FALSE)
tobas <- tobas[tobas$year != "Spring2017", ]

titi <- read.csv("data/Data_S3_all_onsetData_titi.csv", stringsAsFactors = FALSE)

cat("Individual cosinor model selection vs no-rhythm null (", n_perm, " permutations)\n\n", sep = "")

res <- rbind(
  runPopulation(seattle, "ID_year", "DayinCycle", "Start_Time", "Seattle onset", use_weekday = TRUE),
  runPopulation(tobas,   "ID_year", "DayinCycle", "Start_Time", "Tobas onset",   use_weekday = TRUE),
  # Titi monkeys: no weekly social cycle, so no Weekend covariate.
  runPopulation(titi,    "id2",     "DayinCycle", "time",       "Titi onset",    use_weekday = FALSE)
)

write.csv(res, "results/modelSelection_null.csv", row.names = FALSE)
cat("\nSaved: results/modelSelection_null.csv\n")

# NOTE ON SMOOTHING: do NOT pre-smooth the series before this selection. A moving
# average induces autocorrelation and deflates the residual variance, miscalibrating
# the AIC and inflating the selection of the more complex models -- under a permuted
# (no-rhythm) null, a 7-day SMA makes 30+15 win ~77% of the time. The weekly cycle
# is therefore controlled with the Weekend covariate (as in 2_analysis.R), on raw data.
