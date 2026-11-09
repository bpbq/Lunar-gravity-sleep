library(dplyr)
library(TTR)
library(lubridate)
library(tidyverse)
library(mosaic)
library(zoo)
library(ggpubr)
library(here)

rm(list = ls())
setwd(here())


###############################################################################
# SPECTRAL ANALYSIS
#
# A direct spectral analysis demonstrating the existence of the ~30 and ~15 day
# periodicities. This script implements two complementary approaches:
#
# COHORT LEVEL:
#   FFT on the cohort mean profile averaged by DayinCycle (30 points). With 30
#   points, the discrete DFT frequencies correspond exactly to periods of 30,
#   15, 10, 7.5... days, which allows detecting both periodicities of interest.
#   Limitation: the signal already assumes a 30-day periodicity because of the
#   folding, so it is NOT independent evidence that a 30-day cycle exists; it is
#   a harmonic decomposition of the folded waveform.
#
# INDIVIDUAL LEVEL (averaged):
#   Periodogram of each individual's raw series (no folding), computed with
#   spec.pgram() with a 10% taper and linear detrending. The spectra are
#   interpolated onto a common period grid and averaged across individuals.
#   This does not assume a 30-day periodicity a priori and is more conservative.
#   Significance is assessed against an AR(1) red-noise null (see below).
#   Limitation: with 24-65 days of data the spectral resolution is low (~1/N
#   cycles/day), so the 30 and 15 day periods may be poorly resolved for
#   individuals with short series.
#
# Outputs saved to results/ and Figures/
###############################################################################


####### Functions #######

# Function: cohortSpectrum
# Computes the power spectrum (DFT) of a cohort mean profile.
# Input: vector of means per DayinCycle (length n_days)
# Output: data.frame with Period (days) and Power (spectral density)
cohortSpectrum <- function(mean_profile, n_days = 30) {
  x <- mean_profile - mean(mean_profile, na.rm = TRUE)  # Center (remove mesor)
  ft    <- fft(x)
  power <- (Mod(ft)^2) / n_days                         # Power density
  # Positive frequencies only (up to Nyquist)
  freqs  <- seq(0, n_days - 1) / n_days                 # cycles/day
  n_half <- floor(n_days / 2)
  data.frame(
    Period = n_days / seq_len(n_half),                   # days
    Power  = 2 * power[2:(n_half + 1)]                   # x2 for positive components
  )
}


# Function: specPipeline
# Common spectral pipeline: (optional SMA7) + linear detrend + 10% Hanning taper
# + periodogram, interpolated onto period_grid. Used by BOTH the observed series
# and the AR(1) surrogates so that the null has exactly the same filter bias.
#
# The 7-day SMA removes the weekly "free day" effect (greater variability on
# weekends) that would otherwise produce a spurious peak at 7 days and mask the
# lunar periodicities of interest. The SMA7 acts as a low-pass filter that
# attenuates frequencies > 1/7 day^-1 without affecting the 15 and 30 day
# periods. Because the SMA reddens the spectrum, the null must pass through the
# identical filter for the comparison to be fair.
specPipeline <- function(x_clean, period_grid, apply_sma7 = TRUE) {
  if (apply_sma7 && length(x_clean) >= 14) {
    x_clean <- as.numeric(TTR::SMA(x_clean, n = 7))
    x_clean <- x_clean[!is.na(x_clean)]
  }
  if (length(x_clean) < 10) return(NULL)
  pg <- tryCatch(
    spec.pgram(x_clean, detrend = TRUE, taper = 0.1, plot = FALSE, demean = TRUE),
    error = function(e) NULL
  )
  if (is.null(pg)) return(NULL)
  periods_pg <- 1 / pg$freq
  in_range   <- period_grid >= min(periods_pg) & period_grid <= max(periods_pg)
  out <- rep(NA_real_, length(period_grid))
  out[in_range] <- approx(periods_pg, pg$spec, xout = period_grid[in_range], rule = 1)$y
  out
}

# Function: individualSpectrum
# Observed spectrum of one individual (short NA gaps imputed; minimum 24 days).
individualSpectrum <- function(x_raw, period_grid, apply_sma7 = TRUE) {
  x_filled <- zoo::na.approx(x_raw, maxgap = 3, na.rm = FALSE)
  x_clean  <- x_filled[!is.na(x_filled)]
  if (length(x_clean) < 24) return(NULL)
  specPipeline(x_clean, period_grid, apply_sma7)
}

# Function: individualNullSpectra
# Generates n_surr AR(1) (red-noise) null spectra for one individual, passed
# through the SAME pipeline (specPipeline). The AR(1) is fit on the RAW detrended
# series (before the SMA) via its lag-1 autocorrelation (the Yule-Walker
# estimator for AR(1)). This way the null reproduces the real autocorrelation of
# the series and, after the same filter, its same spectral shape. A peak that
# exceeds this null is evidence of periodicity above the red-noise background --
# which the raw periodogram alone does NOT allow one to claim.
individualNullSpectra <- function(x_raw, period_grid, n_surr = 200, apply_sma7 = TRUE) {
  x_filled <- zoo::na.approx(x_raw, maxgap = 3, na.rm = FALSE)
  x_clean  <- x_filled[!is.na(x_filled)]
  n <- length(x_clean)
  if (n < 24) return(NULL)

  # AR(1) on the raw detrended series
  x_dt <- as.numeric(residuals(lm(x_clean ~ seq_len(n))))
  phi  <- tryCatch(stats::acf(x_dt, lag.max = 1, plot = FALSE)$acf[2], error = function(e) 0)
  if (is.na(phi)) phi <- 0
  phi   <- max(min(phi, 0.99), -0.99)
  sigma <- sqrt(max(var(x_dt) * (1 - phi^2), .Machine$double.eps))

  mat <- matrix(NA_real_, nrow = n_surr, ncol = length(period_grid))
  for (k in seq_len(n_surr)) {
    if (abs(phi) < 1e-6) {
      surr <- rnorm(n, 0, sigma)                       # AR(1) degenerates to white noise
    } else {
      surr <- as.numeric(arima.sim(list(ar = phi), n = n, sd = sigma))
    }
    sp <- specPipeline(surr, period_grid, apply_sma7)
    if (!is.null(sp)) mat[k, ] <- sp
  }
  mat
}


# Function: plotMeanSpectrum
# Plots the mean spectrum +/- standard error, marking the 30 and 15 day periods
# with vertical lines. If thresholds are supplied, overlays the AR(1) red-noise
# significance levels.
plotMeanSpectrum <- function(period_grid, mean_power, se_power, n_indiv,
                              title = "", color = "#D55E00",
                              thr95 = NULL, thr99 = NULL) {
  df <- data.frame(
    Period     = period_grid,
    Mean_Power = mean_power,
    SE         = se_power
  )
  if (!is.null(thr95)) df$Thr95 <- thr95
  if (!is.null(thr99)) df$Thr99 <- thr99
  # Keep a reasonable period range for the plot (5-45 days)
  df <- df[df$Period >= 5 & df$Period <= 45, ]

  p <- ggplot(df, aes(x = Period, y = Mean_Power)) +
    # Standard error band
    geom_ribbon(aes(ymin = Mean_Power - SE, ymax = Mean_Power + SE),
                fill = color, alpha = 0.25) +
    geom_line(color = color, linewidth = 1)

  # AR(1) red-noise null significance thresholds
  if (!is.null(thr95))
    p <- p + geom_line(aes(y = Thr95), color = "gray25", linetype = "dashed",
                       linewidth = 0.6)
  if (!is.null(thr99))
    p <- p + geom_line(aes(y = Thr99), color = "gray25", linetype = "dotted",
                       linewidth = 0.6)

  y_top <- max(df$Mean_Power + df$SE, df$Thr95, df$Thr99, na.rm = TRUE)
  p +
    # Reference lines for the periods of interest
    geom_vline(xintercept = 30, linetype = "dashed", color = "black",
               linewidth = 0.7) +
    geom_vline(xintercept = 15, linetype = "dotted", color = "black",
               linewidth = 0.7) +
    annotate("text", x = 30, y = y_top, label = "30d", hjust = -0.1, size = 3.2) +
    annotate("text", x = 15, y = y_top, label = "15d", hjust = -0.1, size = 3.2) +
    labs(
      title    = title,
      subtitle = paste0("n = ", n_indiv, " individuals | gray dashed line: AR(1) red-noise null 95%"),
      x        = "Period (days)",
      y        = "Spectral density (mean ± SE)"
    ) +
    scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 30, 35, 40, 45)) +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.8),
      panel.border     = element_blank(),
      panel.background = element_rect(fill = "white", colour = "lightgray"),
      panel.grid.major = element_line(linewidth = 0.4, linetype = "solid",
                                       colour = "lightgray"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 7)
    )
}


####### Common period grid for individual spectra #######
# Resolution of 0.5 days over the range 5-65 days
period_grid <- seq(5, 65, by = 0.5)

# AR(1) surrogates per individual for the red-noise null (significance).
# 200 is sufficient; lower to 50-100 for fast development (noisier in the tails).
N_SURR <- 200
set.seed(42)


################################################################################
####  COHORT-LEVEL SPECTRA (FFT on the averaged DayinCycle profile)  ###########
################################################################################

cat("=== Cohort-level spectra ===\n")

# --- Seattle ---
sleep_seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)
sleep_seattle <- sleep_seattle[sleep_seattle$year %in% c("SS2020", "AW2020", "SS2024"), ]
sleep_seattle$year <- factor(sleep_seattle$year, levels = c("SS2020", "AW2020", "SS2024"))

plotTitles_seattle <- c("Spring 2020", "Winter 2021", "Spring 2024")
colors_seattle     <- c("#33DDAA", "#33DDAA", "#FF7799")

spec_cohort_list_seattle <- list()

for (var in c("Start_Time", "Duration", "End_Time")) {
  plotlist_cohort <- list()

  for (i in seq_along(levels(sleep_seattle$year))) {
    cohort     <- levels(sleep_seattle$year)[i]
    df_cohort  <- sleep_seattle[sleep_seattle$year == cohort, ]

    # Mean profile per DayinCycle (the 30 days of the cycle)
    mean_profile <- df_cohort %>%
      group_by(DayinCycle) %>%
      summarise(val = mean(.data[[var]], na.rm = TRUE)) %>%
      arrange(DayinCycle)

    if (nrow(mean_profile) < 20) next  # Skip if too few days

    spec_df        <- cohortSpectrum(mean_profile$val, n_days = nrow(mean_profile))
    spec_df$cohort <- cohort
    spec_df$var    <- var
    spec_cohort_list_seattle[[paste(var, cohort, sep = "_")]] <- spec_df

    # Cohort spectrum plot
    p <- ggplot(spec_df, aes(x = Period, y = Power)) +
      geom_bar(stat = "identity", fill = colors_seattle[i], width = 0.5) +
      geom_vline(xintercept = 30, linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = 15, linetype = "dotted", linewidth = 0.7) +
      annotate("text", x = 30, y = max(spec_df$Power),
               label = "30d", hjust = -0.15, size = 3) +
      annotate("text", x = 15, y = max(spec_df$Power),
               label = "15d", hjust = -0.15, size = 3) +
      labs(title = plotTitles_seattle[i],
           x = "Period (days)", y = "Power") +
      scale_x_continuous(breaks = c(5, 10, 15, 20, 30)) +
      theme(
        axis.line        = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = "lightgray"),
        panel.grid.major = element_line(linewidth = 0.3, colour = "lightgray"),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(hjust = 0.5, size = 9)
      )
    plotlist_cohort[[i]] <- p
  }

  # Save cohort panel
  tiff(paste0("Figures/Seattle_CohortSpectrum_", var, ".tif"),
       width = 850, height = 2650, res = 300, compression = "lzw")
  print(ggarrange(plotlist = plotlist_cohort, ncol = 1, nrow = 3))
  dev.off()
  cat("  Saved: Seattle_CohortSpectrum_", var, "\n")
}


# --- Tobas ---
sleep_tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)
sleep_tobas <- sleep_tobas[sleep_tobas$year != "Spring2017", ]
sleep_tobas$year <- factor(sleep_tobas$year,
                            levels = c("Fall-Winter2016", "Fall-Winter2023",
                                       "Spring2018", "Spring2023",
                                       "Summer2024", "Fall-Winter2024"))

plotTitles_tobas <- c("Winter 2016", "Fall 2023", "Spring 2018",
                       "Spring 2023", "Summer 2024", "Fall 2024")
colors_tobas     <- c("#33DDAA", "#55BBFF", "#55BBFF", "#BB66FF", "#FF7799", "#FF7799")

spec_cohort_list_tobas <- list()

for (var in c("Start_Time", "Duration", "End_Time")) {
  plotlist_cohort <- list()

  for (i in seq_along(levels(sleep_tobas$year))) {
    cohort    <- levels(sleep_tobas$year)[i]
    df_cohort <- sleep_tobas[sleep_tobas$year == cohort, ]

    mean_profile <- df_cohort %>%
      group_by(DayinCycle) %>%
      summarise(val = mean(.data[[var]], na.rm = TRUE)) %>%
      arrange(DayinCycle)

    if (nrow(mean_profile) < 20) next

    spec_df        <- cohortSpectrum(mean_profile$val, n_days = nrow(mean_profile))
    spec_df$cohort <- cohort
    spec_df$var    <- var
    spec_cohort_list_tobas[[paste(var, cohort, sep = "_")]] <- spec_df

    p <- ggplot(spec_df, aes(x = Period, y = Power)) +
      geom_bar(stat = "identity", fill = colors_tobas[i], width = 0.5) +
      geom_vline(xintercept = 30, linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = 15, linetype = "dotted", linewidth = 0.7) +
      annotate("text", x = 30, y = max(spec_df$Power),
               label = "30d", hjust = -0.15, size = 3) +
      annotate("text", x = 15, y = max(spec_df$Power),
               label = "15d", hjust = -0.15, size = 3) +
      labs(title = plotTitles_tobas[i],
           x = "Period (days)", y = "Power") +
      scale_x_continuous(breaks = c(5, 10, 15, 20, 30)) +
      theme(
        axis.line        = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = "lightgray"),
        panel.grid.major = element_line(linewidth = 0.3, colour = "lightgray"),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(hjust = 0.5, size = 9)
      )
    plotlist_cohort[[i]] <- p
  }

  tiff(paste0("Figures/Tobas_CohortSpectrum_", var, ".tif"),
       width = 850, height = 5000, res = 300, compression = "lzw")
  print(ggarrange(plotlist = plotlist_cohort, ncol = 1, nrow = 6))
  dev.off()
  cat("  Saved: Tobas_CohortSpectrum_", var, "\n")
}

# Save all cohort spectra to CSV
all_cohort_spectra <- bind_rows(
  bind_rows(spec_cohort_list_seattle) %>% mutate(population = "Seattle"),
  bind_rows(spec_cohort_list_tobas)   %>% mutate(population = "Tobas")
)
write.csv(all_cohort_spectra, "results/cohort_spectra.csv", row.names = FALSE)


################################################################################
####  INDIVIDUAL-LEVEL SPECTRA (averaged across individuals)  ##################
################################################################################

cat("\n=== Individual-level spectra (averaged) ===\n")
cat("NOTE: with 24-65 days of data per individual, the spectral resolution is low.\n")
cat("The 30 and 15 day periods may be poorly resolved in short series.\n\n")

for (population in c("Seattle", "Tobas")) {

  if (population == "Seattle") {
    sleep_raw  <- sleep_seattle
    color_line <- "#D55E00"
  } else {
    sleep_raw  <- sleep_tobas
    color_line <- "#0072B2"
  }

  for (var in c("Start_Time", "Duration", "End_Time")) {
    cat(population, "-", var, "\n")

    spec_matrix <- list()
    null_list   <- list()   # AR(1) null spectra per individual [N_SURR x P]

    for (id in unique(sleep_raw$ID_year)) {
      df_id  <- sleep_raw[sleep_raw$ID_year == id, ]
      # Sort chronologically so the series is temporal
      df_id  <- df_id[order(df_id$night), ]
      x_raw  <- df_id[[var]]

      spec_i <- individualSpectrum(x_raw, period_grid)
      if (!is.null(spec_i)) {
        spec_matrix[[id]] <- spec_i
        null_list[[id]]   <- individualNullSpectra(x_raw, period_grid, n_surr = N_SURR)
      }
    }

    n_indiv <- length(spec_matrix)
    cat("  Individuals with a valid spectrum: ", n_indiv, "\n")

    if (n_indiv < 3) {
      cat("  Too few individuals. Skipping.\n")
      next
    }

    spec_mat   <- do.call(rbind, spec_matrix)
    mean_power <- colMeans(spec_mat, na.rm = TRUE)
    se_power   <- apply(spec_mat, 2, function(x) {
      n <- sum(!is.na(x))
      if (n < 2) return(NA_real_)
      sd(x, na.rm = TRUE) / sqrt(n)
    })
    n_per_period <- apply(spec_mat, 2, function(x) sum(!is.na(x)))

    # --- Red-noise null for the MEAN spectrum ---
    # At each iteration k we average the k-th surrogate of every individual,
    # building a null distribution of the mean spectrum. The 95/99% threshold is
    # the level the observed mean spectrum must exceed to be significant.
    ids_ok     <- names(null_list)[!vapply(null_list, is.null, logical(1))]
    null_means <- matrix(NA_real_, nrow = N_SURR, ncol = length(period_grid))
    for (k in seq_len(N_SURR)) {
      surr_k <- vapply(ids_ok, function(idn) null_list[[idn]][k, ],
                       numeric(length(period_grid)))          # [P x n_ids]
      null_means[k, ] <- rowMeans(surr_k, na.rm = TRUE)
    }
    thr95 <- apply(null_means, 2, quantile, 0.95, na.rm = TRUE)
    thr99 <- apply(null_means, 2, quantile, 0.99, na.rm = TRUE)

    # --- Per-individual significance in the 15 and 30 day bands ---
    # (each individual periodogram is very noisy; this is complementary to the
    # test on the mean spectrum, which is more powerful)
    sig_band <- function(target) {
      idx <- which.min(abs(period_grid - target))
      vals <- vapply(ids_ok, function(idn) {
        obs <- spec_matrix[[idn]][idx]
        nl  <- null_list[[idn]][, idx]
        if (is.na(obs) || all(is.na(nl))) return(NA_real_)
        as.numeric(obs > quantile(nl, 0.95, na.rm = TRUE))
      }, numeric(1))
      c(n = sum(!is.na(vals)), n_sig = sum(vals, na.rm = TRUE),
        frac = mean(vals, na.rm = TRUE))
    }
    s15 <- sig_band(15); s30 <- sig_band(30)
    cat(sprintf("  Per-individual significance (obs > null95): 15d %g/%g (%.1f%%), 30d %g/%g (%.1f%%); chance ~5%%\n",
                s15["n_sig"], s15["n"], 100 * s15["frac"],
                s30["n_sig"], s30["n"], 100 * s30["frac"]))
    # Exceedance of the MEAN spectrum at 15 and 30
    i15 <- which.min(abs(period_grid - 15)); i30 <- which.min(abs(period_grid - 30))
    cat(sprintf("  Mean spectrum vs null95: 15d %s, 30d %s\n",
                ifelse(mean_power[i15] > thr95[i15], "SIGNIFICANT", "n.s."),
                ifelse(mean_power[i30] > thr95[i30], "SIGNIFICANT", "n.s.")))

    # Save mean spectrum + thresholds to CSV
    spec_summary <- data.frame(
      Period        = period_grid,
      Mean_Power    = mean_power,
      SE            = se_power,
      N_individuals = n_per_period,
      Null95        = thr95,
      Null99        = thr99
    )
    write.csv(spec_summary,
              paste0("results/", population, "_meanSpectrum_", var, ".csv"),
              row.names = FALSE)

    # Mean spectrum plot with significance thresholds
    p <- plotMeanSpectrum(
      period_grid = period_grid,
      mean_power  = mean_power,
      se_power    = se_power,
      n_indiv     = n_indiv,
      title       = paste0(population, " - ", var, "\nMean individual spectrum"),
      color       = color_line,
      thr95       = thr95
    )

    tiff(paste0("Figures/", population, "_MeanSpectrum_", var, ".tif"),
         width = 1200, height = 900, res = 300, compression = "lzw")
    print(p)
    dev.off()
    cat("  Saved: ", population, "_MeanSpectrum_", var, "\n")
  }
}


################################################################################
####  COMPARATIVE PANEL: Cohort vs Individual (Start_Time)  ####################
################################################################################
# Summary figure showing both approaches side by side for the manuscript.

cat("\nBuilding comparative panel...\n")

# Cohort spectra: average power at each period across cohorts
spec_seattle_st <- all_cohort_spectra %>%
  filter(population == "Seattle", var == "Start_Time") %>%
  group_by(Period) %>%
  summarise(mean_power = mean(Power), se_power = sd(Power) / sqrt(n()), .groups = "drop")

spec_tobas_st <- all_cohort_spectra %>%
  filter(population == "Tobas", var == "Start_Time") %>%
  group_by(Period) %>%
  summarise(mean_power = mean(Power), se_power = sd(Power) / sqrt(n()), .groups = "drop")

# Load the saved individual spectra
ind_seattle <- tryCatch(
  read.csv("results/Seattle_meanSpectrum_Start_Time.csv"), error = function(e) NULL)
ind_tobas   <- tryCatch(
  read.csv("results/Tobas_meanSpectrum_Start_Time.csv"),   error = function(e) NULL)

# Build a panel of plots (one per population)
plot_list_panel <- list()

if (!is.null(ind_seattle)) {
  plot_list_panel[["Seattle_ind"]] <- plotMeanSpectrum(
    ind_seattle$Period, ind_seattle$Mean_Power, ind_seattle$SE,
    n_indiv = max(ind_seattle$N_individuals, na.rm = TRUE),
    title = "Seattle: Individual spectrum (averaged)",
    color = "#D55E00",
    thr95 = ind_seattle$Null95)
}
if (!is.null(ind_tobas)) {
  plot_list_panel[["Tobas_ind"]] <- plotMeanSpectrum(
    ind_tobas$Period, ind_tobas$Mean_Power, ind_tobas$SE,
    n_indiv = max(ind_tobas$N_individuals, na.rm = TRUE),
    title = "Tobas: Individual spectrum (averaged)",
    color = "#0072B2",
    thr95 = ind_tobas$Null95)
}

if (length(plot_list_panel) > 0) {
  tiff("Figures/ComparativeSpectra_StartTime.tif",
       width  = 2400, height = 900 * ceiling(length(plot_list_panel) / 2),
       res    = 300, compression = "lzw")
  print(ggarrange(plotlist = plot_list_panel,
                  ncol = 2, nrow = ceiling(length(plot_list_panel) / 2)))
  dev.off()
  cat("Comparative panel saved to Figures/ComparativeSpectra_StartTime.tif\n")
}

cat("\n=== Spectral analysis complete ===\n")
