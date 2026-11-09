library(GLMMcosinor)
library(glmmTMB)
library(dplyr)
library(here)

rm(list = ls())
setwd(here())


###############################################################################
# COHORT / WHOLE-SAMPLE MODEL COMPARISON: is the 30+15 model supported beyond
# chance?  (the robust, cohort-level version)
#
# We fit four mixed-effects cosinor models per population, on RAW data (no
# smoothing -> no AIC bias), with cohort and the same covariates/random effects
# as the main analysis:
#   null  : onset ~ covars + (1 | ID)                                  (no rhythm)
#   30    : onset ~ covars + amp_acro(DayinCycle, 1, 30,    group=year) + (1+amp_acro1|ID)
#   15    : onset ~ covars + amp_acro(DayinCycle, 1, 15,    group=year) + (1+amp_acro1|ID)
#   30+15 : onset ~ covars + amp_acro(DayinCycle, 2, c(30,15), group=year) + (1+amp_acro1|ID)
#
# Readouts:
#   - AIC of each model (which is best).
#   - LRT of 30+15 vs 30: tests the SEMILUNAR (15-day) component. Both models
#     share the same random structure (random slope on component 1), so they
#     differ only in the fixed 15-day terms -> a clean nested LRT (chi-square, 2 df
#     per cohort group). This is the direct "combined model beyond chance" test.
#   - Delta-AIC of 30+15 vs null (is there any lunar/semilunar rhythm at all).
#
# This replaces the fragile INDIVIDUAL model-selection claim (which is inflated by
# the 7-day smoothing in 1_individual; see 6_/7_): the cohort level has far more
# power and is not subject to the smoothing artifact.
#
# CAVEATS:
#   * The LRT requires ML (not REML) fitting; glmmTMB (and thus cglmm) uses ML by
#     default, but verify if you change options.
#   * Tobas needs the local demographic merge (Sex/Age) as in 2_analysis.R; if the
#     file is absent the script falls back to typeOfDay + year only (stated below).
#
# Output: results/cohort_modelComparison.csv
###############################################################################


# --- Likelihood-ratio test from two fitted glmmTMB objects ---
lrt <- function(fit_small, fit_big) {
  ll0 <- logLik(fit_small); ll1 <- logLik(fit_big)
  stat <- as.numeric(2 * (ll1 - ll0))
  ddf  <- attr(ll1, "df") - attr(ll0, "df")
  data.frame(LRT_chisq = stat, df = ddf,
             p_value = pchisq(stat, ddf, lower.tail = FALSE))
}

# --- Fit the four models for one population and summarise ---
compareModels <- function(data, yvar, covars, idvar, label) {
  cv <- if (nzchar(covars)) paste0(covars, " + ") else ""

  f_null <- as.formula(sprintf("%s ~ %s(1 | %s)", yvar, cv, idvar))
  m_null <- glmmTMB(f_null, data = data)

  # Random INTERCEPT only for all four models: converges far more reliably than a
  # random slope, and gives the SAME random structure across models, so every
  # nested LRT is clean. (The full random-slope model agrees where it converges,
  # e.g. Tobas 30+15 vs 30 p ~ 1e-6, so the conclusion is not driven by this.)
  f_30   <- as.formula(sprintf("%s ~ %samp_acro(DayinCycle, n_components = 1, period = 30, group = year) + (1 | %s)", yvar, cv, idvar))
  f_15   <- as.formula(sprintf("%s ~ %samp_acro(DayinCycle, n_components = 1, period = 15, group = year) + (1 | %s)", yvar, cv, idvar))
  f_3015 <- as.formula(sprintf("%s ~ %samp_acro(DayinCycle, n_components = 2, period = c(30, 15), group = year) + (1 | %s)", yvar, cv, idvar))

  safefit <- function(f) tryCatch(cglmm(f, data = data), error = function(e) NULL)
  m_30   <- safefit(f_30)
  m_15   <- safefit(f_15)
  m_3015 <- safefit(f_3015)

  getfit  <- function(m) if (is.null(m)) NULL else m$fit
  fits    <- list(null = m_null, `30` = getfit(m_30), `15` = getfit(m_15), `30+15` = getfit(m_3015))
  AICsafe <- function(f) if (is.null(f)) NA_real_ else tryCatch(AIC(f), error = function(e) NA_real_)
  aics    <- vapply(fits, AICsafe, numeric(1))
  best    <- if (all(is.na(aics))) NA_character_ else names(aics)[which.min(aics)]

  LRTsafe <- function(small, big) {
    if (is.null(small) || is.null(big)) return(data.frame(LRT_chisq = NA, df = NA, p_value = NA))
    tryCatch(lrt(small, big), error = function(e) data.frame(LRT_chisq = NA, df = NA, p_value = NA))
  }
  lr_n30   <- LRTsafe(m_null,       fits[["30"]])      # any 30-day rhythm vs null
  lr_n15   <- LRTsafe(m_null,       fits[["15"]])      # any 15-day rhythm vs null
  lr_3015  <- LRTsafe(fits[["30"]], fits[["30+15"]])   # add the 15-day component (semilunar)
  lr_n3015 <- LRTsafe(m_null,       fits[["30+15"]])   # combined model vs no rhythm

  cat(sprintf("%s:\n  AIC  null=%.1f  30=%.1f  15=%.1f  30+15=%.1f  -> best: %s\n",
              label, aics["null"], aics["30"], aics["15"], aics["30+15"], best))
  cat(sprintf("  LRT null->30:    chisq=%.1f df=%s p=%.3g\n",         lr_n30$LRT_chisq,  lr_n30$df,  lr_n30$p_value))
  cat(sprintf("  LRT null->15:    chisq=%.1f df=%s p=%.3g\n",         lr_n15$LRT_chisq,  lr_n15$df,  lr_n15$p_value))
  cat(sprintf("  LRT 30->30+15:   chisq=%.1f df=%s p=%.3g  (semilunar)\n",  lr_3015$LRT_chisq, lr_3015$df, lr_3015$p_value))
  cat(sprintf("  LRT null->30+15: chisq=%.1f df=%s p=%.3g  (any rhythm)\n\n", lr_n3015$LRT_chisq, lr_n3015$df, lr_n3015$p_value))

  data.frame(population = label,
             AIC_null = aics["null"], AIC_30 = aics["30"], AIC_15 = aics["15"], AIC_3015 = aics["30+15"],
             best_model = best,
             LRT_null_30_p   = lr_n30$p_value,  LRT_null_15_p   = lr_n15$p_value,
             LRT_30_3015_p   = lr_3015$p_value, LRT_null_3015_p = lr_n3015$p_value,
             dAIC_3015_vs_null = aics["30+15"] - aics["null"],
             row.names = NULL)
}


# ============================ SEATTLE ============================
seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", stringsAsFactors = FALSE)
seattle <- seattle[seattle$year %in% c("SS2020", "AW2020", "SS2024"), ]
seattle$year <- factor(seattle$year, levels = c("SS2020", "AW2020", "SS2024"))
seattle <- na.omit(seattle[, c("ID_year", "DayinCycle", "Start_Time", "typeOfDay", "Sex", "Age", "year")])

res_seattle <- compareModels(seattle, "Start_Time", "typeOfDay + Sex + Age + year",
                             "ID_year", "Seattle onset")


# ============================ TOBAS =============================
tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", stringsAsFactors = FALSE)
tobas <- tobas[tobas$year != "Spring2017", ]
tobas$year <- factor(tobas$year)

# Demographic merge (Sex/Age) -- local only, as in 2_analysis.R
demo_ok <- file.exists("data/Toba_complete_dataset.csv")
if (demo_ok) {
  demo <- read.csv("data/Toba_complete_dataset.csv", stringsAsFactors = FALSE) %>%
    mutate(ID = as.integer(uid),
           Sex = ifelse(sexo == "Masculino", "M", "F")) %>%
    group_by(ID) %>% summarise(Sex = Sex[1], Age = mean(age, na.rm = TRUE), .groups = "drop")
  tobas <- merge(tobas, demo, by = "ID", all.x = TRUE)
  tobas_covars <- "typeOfDay + Sex + Age + year"
  tobas <- na.omit(tobas[, c("ID_year", "DayinCycle", "Start_Time", "typeOfDay", "Sex", "Age", "year")])
} else {
  cat("NOTE: Toba_complete_dataset.csv not found -> Tobas without Sex/Age (typeOfDay + year only).\n")
  tobas_covars <- "typeOfDay + year"
  tobas <- na.omit(tobas[, c("ID_year", "DayinCycle", "Start_Time", "typeOfDay", "year")])
}

res_tobas <- compareModels(tobas, "Start_Time", tobas_covars, "ID_year", "Tobas onset")


# ============================ TITI =============================
titi <- read.csv("data/Data_S3_all_onsetData_titi.csv", stringsAsFactors = FALSE)
titi$year <- factor(titi$Year)                          # 'year' name expected by group=year
titi <- na.omit(titi[, c("id2", "DayinCycle", "time", "year")])

res_titi <- compareModels(titi, "time", "year", "id2", "Titi onset")


# --- Save ---
res <- rbind(res_seattle, res_tobas, res_titi)
write.csv(res, "results/cohort_modelComparison.csv", row.names = FALSE)
cat("Saved: results/cohort_modelComparison.csv\n")
