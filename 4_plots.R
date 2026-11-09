library(dplyr)
library(TTR) # Moving Average function
library(lubridate) # Date/Time format
library(data.table) # for setDT function
library(tidyr)
library(cgwtools)
library(suncalc)
library(tidyverse) # to work with dataframes
library(mosaic) # to calculate zscores
library(nlstools) # to analyze NLS fits
library(data.table)
library(zoo)
library(stats)
library(AICcmodavg)
library(modelsummary)
library(GLMMcosinor)
library(tibble)
library(ggpubr) # For arranging multiple plots (ggarrange)
library(CircStats)
library(here)

# Clear environment
rm(list=ls())

# Set working directory
setwd(here())

##### Functions #####

# Function: transformPhase
# Converts acrophase (peak time) from radians to "Days in Cycle".
# Shifts by period/4 to align Sine reference (starts at 0) to Cosine reference (starts at peak).
transformPhase <- function(x, period){
  x = (((x*180/pi)*period)/360)-(period/4) 
  return(x)
}

# Function: sinusoid
# Generates the mathematical function for the fitted curve (Sine wave).
# Supports single component (n=1, e.g., 30-day) or dual component (n!=1, e.g., 30+15 day).
sinusoid <-function(data, results, n=1, period){
  if (n==1){
    fun = function (x){ results[1,"M"] + (results[1,"A"] * sin((2*pi*1/period)* (x - (results[1,"P"]))))}
  }else {
    fun = function (x){ results[1,"M"] + (results[1,"A"] * sin((2*pi*1/period)* (x - (results[1,"P"]))) + results[1,"A2"] * sin((2*pi*1/(period/2))* (x - (results[1,"P2"])))) 
    }
  }
  return(fun)
}

# Function: lunarPlot
# Standard plot for Sleep Data vs Lunar Cycle.
# - Plots raw data means (points).
# - Overlays the fitted sinusoidal model (line).
# - Marks Full/New moon phases with vertical lines.
lunarPlot <- function(data, variables, results, n=1, period, color){
  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle), y = variable)) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    stat_function(fun = sinusoid(data, results, n, period),
                  colour = color, size = 1.5) +
    stat_summary(fun="mean", colour = "black", pch=16, size=0.5, position = position_dodge(width = .1)) +
    labs(x = "Night in cycle") +
    scale_y_continuous(
      name   = if (grepl("hours from midnight", variables)) sub("hours from midnight", "clock time", variables) else variables,
      labels = if (grepl("hours from midnight", variables)) hours_to_clock else waiver()
    ) +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15),
                       labels = c("Full", "New","Full", "New", "Full"))+
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text(),
          axis.text = element_text(size=11))
  return(plot)
}

# Function: gravPlot
# Plot for Gravitational Variation only.
# Uses spline interpolation for a smooth curve instead of a fitted sine wave.
gravPlot <- function(data, variables, color){
  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle), y = variable)) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    geom_line(data=as.data.frame(spline(data$DayinCycle, data$variable)), aes(x=x,y=y),color=color,size=1.5)+
    stat_summary(fun="mean", colour = "black", pch=16, size=0.5, position = position_dodge(width = .1)) + 
    labs(x = "Night in cycle", y = variables) +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15), 
                       labels = c("Full", "New","Full", "New", "Full"))+
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text(),
          axis.text = element_text(size=11))
  return(plot)
}

# Function: scaling
# Helper function to scale a secondary variable (Gravity) to fit the primary y-axis (Sleep).
# Used for dual-axis plotting.
scaling = function(new, old) {y = (((max(new)-min(new)) * (old-min(old)))/(max(old)-min(old))) + min(new)
return(y)}

# Function: hours_to_clock
# Formats an axis of "hours from midnight" (continuous; may be negative or >24) as
# clock-time labels HH:MM, wrapping modulo 24 ONLY in the displayed label (the
# underlying scale stays continuous). Used so onset/offset axes read as e.g. 23:50
# / 00:10 / 02:50 while the models keep working in hours-from-midnight.
hours_to_clock <- function(x) {
  out <- rep(NA_character_, length(x))
  ok  <- !is.na(x)
  h   <- x[ok] %% 24
  hh  <- floor(h)
  mm  <- round((h - hh) * 60)
  hh  <- (hh + mm %/% 60) %% 24
  mm  <- mm %% 60
  out[ok] <- sprintf("%02d:%02d", as.integer(hh), as.integer(mm))
  out
}

# Function: lunarPlot_grav_nopoints
# Creates a dual-axis plot:
# - Primary Axis (Left): Sleep Variable (Fitted Sine Wave only, no raw points).
# - Secondary Axis (Right): Gravity (Scaled Line).
# This visualizes the phase relationship between Sleep and Gravity.
lunarPlot_grav_nopoints <- function(data, variables, results, n=1, period, color, var_grav, ylab_grav,color2){
  new_min=min(data$variable)
  new_max=max(data$variable)
  old_min=min(as.numeric(data[,var_grav]))
  old_max=max(as.numeric(data[,var_grav]))
  
  # Inverse scaling function for the secondary axis labels
  inv_scale_fun <- function(x) {
    ((x - new_min) * (old_max - old_min)) / (new_max - new_min) + old_min}
  
  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle))) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    # Sleep Model (Sine Wave)
    stat_function(fun = sinusoid(data, results, n, period),
                  colour = color, size = 1.5) +
    # Gravity Data (Scaled Line)
    geom_line(aes(y=scaling(data$variable, as.numeric(data[,var_grav]))), color=color2,size=1.5, linetype="81")+
    # Dual Axis Definition
    scale_y_continuous(
      name   = if (grepl("hours from midnight", variables)) sub("hours from midnight", "clock time", variables) else variables,
      labels = if (grepl("hours from midnight", variables)) hours_to_clock else waiver(),
      sec.axis = sec_axis(transform = ~ inv_scale_fun(.), name = ylab_grav)
    ) +
    labs(x = "Night in cycle") +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15), 
                       labels = c("Full", "New","Full", "New", "Full"))+
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text(),
          axis.text = element_text(size=11))
  return(plot)
}

# Function: lunarPlot_grav_nopoints2
# Similar to previous function, but includes:
# 1. Spline smoothing for the sleep data (instead of just fitted sine).
# 2. Raw mean points for sleep data.
# Used for Figure 4/5 where we show raw data + fit + gravity.
lunarPlot_grav_nopoints2 <- function(data, variables, results, n=1, period, color, var_grav, ylab_grav,color2){
  # [Same scaling logic as above]
  new_min=min(data$variable)
  new_max=max(data$variable)
  old_min=min(as.numeric(data[,var_grav]))
  old_max=max(as.numeric(data[,var_grav]))
  inv_scale_fun <- function(x) {
    ((x - new_min) * (old_max - old_min)) / (new_max - new_min) + old_min}
  
  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle))) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    # Spline smoothed sleep line
    geom_line(data=as.data.frame(spline(data$DayinCycle, data$variable)), aes(x=x,y=y),color=color,size=1.5)+
    # Scaled Gravity Line
    geom_line(aes(y=scaling(data$variable, as.numeric(data[,var_grav]))), color=color2,size=1.5, linetype="81")+
    # Raw Sleep Data Points
    stat_summary( aes(x=as.numeric(DayinCycle),y=data$variable), fun="mean", colour = "black", pch=16, size=0.45, position = position_dodge(width = .1)) + 
    
    scale_y_continuous(
      name   = if (grepl("hours from midnight", variables)) sub("hours from midnight", "clock time", variables) else variables,
      labels = if (grepl("hours from midnight", variables)) hours_to_clock else waiver(),
      sec.axis = sec_axis(transform = ~ inv_scale_fun(.), name = ylab_grav)
    ) +
    labs(x = "Night in cycle") +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15), 
                       labels = c("Full", "New","Full", "New", "Full"))+
    # [Theme settings same as above]
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text(),
          axis.text = element_text(size=11))
  return(plot)
}

# Function: lunarPlot_grav_pointsonly
# Dual-axis plot: gravity as a scaled line (right axis) + sleep raw mean points
# (left axis), with NO fitted/spline sleep line. Overlay format (b).
lunarPlot_grav_pointsonly <- function(data, variables, color, var_grav, ylab_grav, color2){
  new_min=min(data$variable); new_max=max(data$variable)
  old_min=min(as.numeric(data[,var_grav])); old_max=max(as.numeric(data[,var_grav]))
  inv_scale_fun <- function(x) {
    ((x - new_min) * (old_max - old_min)) / (new_max - new_min) + old_min}

  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle))) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    # Scaled Gravity Line
    geom_line(aes(y=scaling(data$variable, as.numeric(data[,var_grav]))), color=color2, size=1.5, linetype="81")+
    # Raw Sleep Data Points (no fitted line)
    stat_summary(aes(x=as.numeric(DayinCycle), y=data$variable), fun="mean", colour = "black", pch=16, size=0.45, position = position_dodge(width = .1)) +
    scale_y_continuous(
      name   = if (grepl("hours from midnight", variables)) sub("hours from midnight", "clock time", variables) else variables,
      labels = if (grepl("hours from midnight", variables)) hours_to_clock else waiver(),
      sec.axis = sec_axis(transform = ~ inv_scale_fun(.), name = ylab_grav)
    ) +
    labs(x = "Night in cycle") +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15),
                       labels = c("Full", "New","Full", "New", "Full"))+
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text(),
          axis.text = element_text(size=11))
  return(plot)
}

# Function: lunarPlot_grav_modelpoints
# Dual-axis plot: sleep COSINOR MODEL line (left axis) + raw sleep mean points +
# gravity as a scaled line (right axis). Overlay format (c): the model prediction
# overlaid on the points, together with gravity (NOT a spline through the points).
lunarPlot_grav_modelpoints <- function(data, variables, results, n=1, period, color, var_grav, ylab_grav, color2){
  new_min=min(data$variable); new_max=max(data$variable)
  old_min=min(as.numeric(data[,var_grav])); old_max=max(as.numeric(data[,var_grav]))
  inv_scale_fun <- function(x) {
    ((x - new_min) * (old_max - old_min)) / (new_max - new_min) + old_min}

  plot<-ggplot(data = data, aes(x = as.numeric(DayinCycle))) +
    geom_vline(aes(xintercept = 0))+
    geom_vline(aes(xintercept = 30))+
    geom_vline(aes(xintercept = 60))+
    # Sleep cosinor MODEL line
    stat_function(fun = sinusoid(data, results, n, period), colour = color, size = 1.5) +
    # Scaled Gravity Line
    geom_line(aes(y=scaling(data$variable, as.numeric(data[,var_grav]))), color=color2, size=1.5, linetype="81")+
    # Raw Sleep Data Points
    stat_summary(aes(x=as.numeric(DayinCycle), y=data$variable), fun="mean", colour = "black", pch=16, size=0.45, position = position_dodge(width = .1)) +
    scale_y_continuous(
      name   = if (grepl("hours from midnight", variables)) sub("hours from midnight", "clock time", variables) else variables,
      labels = if (grepl("hours from midnight", variables)) hours_to_clock else waiver(),
      sec.axis = sec_axis(transform = ~ inv_scale_fun(.), name = ylab_grav)
    ) +
    labs(x = "Night in cycle") +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15),
                       labels = c("Full", "New","Full", "New", "Full"))+
    theme(axis.line = element_line(colour = 'black', size = 1),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white",colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgray"),
          axis.title.x = element_text(),
          axis.text = element_text(size=11))
  return(plot)
}

# Function: extractResults_30
# Extracts amplitude and phase of the 30-day model from the GLMM results table.
# Handles the case where grepping returns an extra header row (length == 3).
extractResults_30 <- function(resultsAll, cohort, means_cohort) {
  results <- resultsAll[grepl(cohort, resultsAll$X), ]
  if (length(results[, 1]) == 3) results <- results[2:3, ]
  data.frame(M = mean(means_cohort$variable),
             A = results[1, 2],
             P = transformPhase(results[2, 2], 30))
}

# Function: extractResults_30_15
# Extracts amplitude and phase of BOTH components of the 30+15 model.
# GLMMcosinor reports the rows in order: amp1, amp2, acr1, acr2
# (all amplitudes first, then all acrophases).
# CHECK with names(summary(mixed_mod_30_15)$transformed.table) if it changes between versions.
extractResults_30_15 <- function(resultsAll, cohort, means_cohort) {
  results <- resultsAll[grepl(cohort, resultsAll$X), ]
  # Separate by parameter type using the row name
  amp_rows <- results[grepl("amp", results$X, ignore.case = TRUE), ]
  acr_rows <- results[grepl("acr", results$X, ignore.case = TRUE), ]

  if (nrow(amp_rows) < 2 || nrow(acr_rows) < 2) {
    warning(paste("extractResults_30_15: insufficient rows for", cohort,
                  "- amp:", nrow(amp_rows), "acr:", nrow(acr_rows),
                  ". Check the CSV format."))
    return(NULL)
  }
  # amp_rows[1] = amp1 (30d), amp_rows[2] = amp2 (15d)
  # acr_rows[1] = acr1 (30d), acr_rows[2] = acr2 (15d)
  data.frame(M  = mean(means_cohort$variable),
             A  = amp_rows[1, 2],
             P  = transformPhase(acr_rows[1, 2], 30),
             A2 = amp_rows[2, 2],
             P2 = transformPhase(acr_rows[2, 2], 15))
}

# Function: makeCohortFigs
# Produces the full, consistent cohort figure set for one population, with cohorts
# stacked in a single column, reusing the helpers above (methods unchanged):
#   - column cosinor, 30-day model           -> <pop>_<var>_30.tif
#   - column cosinor, 30+15 model            -> <pop>_<var>_30+15.tif
#   - overlay (a) lines only (sleep + grav)  -> <pop>_<var>_overlay_lines.tif
#   - overlay (b) gravity line + sleep points-> <pop>_<var>_overlay_points.tif
#   - overlay (c) gravity line + sleep line+pts -> <pop>_<var>_overlay_linespoints.tif
# for each variable in 'var_specs'. Results are read from
#   <results_dir>/GLMMCosinor_<csv>.csv  and  _30+15.csv  (the 30+15 file is optional).
# 'means' must have columns: year, DayinCycle, mean_Total_C and the variable columns.
makeCohortFigs <- function(means, cohorts, titles, colors, colors2,
                           results_dir, pop_label, outdir,
                           var_specs = list(
                             list(col = "Start_Time", tag = "onset",    csv = "StartTime", ylab = "Sleep onset (hours from midnight)"),
                             list(col = "End_Time",   tag = "offset",   csv = "EndTime",   ylab = "Sleep offset (hours from midnight)"),
                             list(col = "Duration",   tag = "duration", csv = "Duration",  ylab = "Sleep duration (minutes)")),
                           grav_lab = "Daily-mean gravitational deviation (\U00B5Gal)") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  nC <- length(cohorts)
  ttl <- function(p, i) p + ggtitle(titles[i]) + theme(plot.title = element_text(hjust = 0.5))

  for (vs in var_specs) {
    res30   <- read.csv(file.path(results_dir, paste0("GLMMCosinor_", vs$csv, ".csv")), header = TRUE)
    f3015   <- file.path(results_dir, paste0("GLMMCosinor_", vs$csv, "_30+15.csv"))
    res3015 <- if (file.exists(f3015)) read.csv(f3015, header = TRUE) else NULL

    has3015 <- !is.null(res3015)

    # Build each cohort's panels inside its OWN environment (lapply) so the cosinor
    # closures capture that cohort's parameters. A plain for-loop reuses one
    # environment and, because ggplot evaluates the closures lazily at print time,
    # every panel would otherwise be drawn with the LAST cohort's M/A/P.
    panels <- lapply(seq_along(cohorts), function(i) {
      co <- cohorts[i]
      mc <- as.data.frame(means[means$year == co, ])
      mc$variable <- mc[[vs$col]]
      rp30   <- extractResults_30(res30, co, mc)
      rp3015 <- if (has3015) extractResults_30_15(res3015, co, mc) else NULL

      # Common y-axis range for ALL plot types of this cohort+variable: covers the
      # data points AND both fitted curves, so the axis does not shift between the
      # plain, model and overlay versions. coord_cartesian zooms without clipping.
      xs   <- seq(0, 60, length.out = 200)
      vals <- c(mc$variable, sinusoid(mc, rp30, 1, 30)(xs))
      if (!is.null(rp3015)) vals <- c(vals, sinusoid(mc, rp3015, 2, 30)(xs))
      r    <- range(vals, na.rm = TRUE); pad <- diff(r) * 0.05
      yr   <- c(r[1] - pad, r[2] + pad)
      lim  <- function(p) p + coord_cartesian(ylim = yr)

      out <- list(
        p30 = lim(ttl(lunarPlot(mc, vs$ylab, rp30, 1, 30, colors[i]), i)),
        pa  = lim(ttl(lunarPlot_grav_nopoints(mc, vs$ylab, rp30, 1, 30, colors[i], "mean_Total_C", grav_lab, colors2[i]), i)),
        pb  = lim(ttl(lunarPlot_grav_pointsonly(mc, vs$ylab, colors[i], "mean_Total_C", grav_lab, colors2[i]), i)),
        pc  = lim(ttl(lunarPlot_grav_modelpoints(mc, vs$ylab, rp30, 1, 30, colors[i], "mean_Total_C", grav_lab, colors2[i]), i)),
        p3015 = NULL, pa3015 = NULL, pc3015 = NULL)

      if (!is.null(rp3015)) {
        out$p3015  <- lim(ttl(lunarPlot(mc, vs$ylab, rp3015, 2, 30, colors[i]), i))
        out$pa3015 <- lim(ttl(lunarPlot_grav_nopoints(mc, vs$ylab, rp3015, 2, 30, colors[i], "mean_Total_C", grav_lab, colors2[i]), i))
        out$pc3015 <- lim(ttl(lunarPlot_grav_modelpoints(mc, vs$ylab, rp3015, 2, 30, colors[i], "mean_Total_C", grav_lab, colors2[i]), i))
      } else if (has3015) {
        out$p3015 <- out$pa3015 <- out$pc3015 <- (ggplot() + theme_void())
      }
      # Realise the lazy cosinor closures now so each panel is frozen to its cohort.
      invisible(lapply(Filter(Negate(is.null), out), ggplot_build))
      out
    })
    pick <- function(nm) lapply(panels, `[[`, nm)

    H <- 850 * nC
    savecol <- function(L, suffix, w, vec = FALSE) {
      base <- file.path(outdir, paste0(pop_label, "_", vs$tag, "_", suffix))
      tiff(paste0(base, ".tif"),
           width = w, height = H, res = 300, compression = "lzw")
      print(ggarrange(plotlist = L, ncol = 1, nrow = nC)); dev.off()
      # vector twin at the SAME logical size (inches = pixels/300) so the layout,
      # axis ranges and tick marks are identical to the .tif, just resolution-free.
      if (vec) {
        svg(paste0(base, ".svg"), width = w/300, height = H/300)
        print(ggarrange(plotlist = L, ncol = 1, nrow = nC)); dev.off()
      }
    }
    savecol(pick("p30"), "30", 850)
    savecol(pick("pa"),  "overlay_lines",       1050)
    savecol(pick("pb"),  "overlay_points",      1050)
    savecol(pick("pc"),  "overlay_linespoints", 1050)
    if (has3015) {
      savecol(pick("p3015"),  "30+15", 850, vec = TRUE)
      savecol(pick("pa3015"), "overlay_lines_30+15",       1050)
      savecol(pick("pc3015"), "overlay_linespoints_30+15", 1050, vec = TRUE)
    }
  }
  cat(paste0("  ", pop_label, ": cohort cosinor + overlay figures saved in ", outdir, "/\n"))
}


# Function: circularShiftVec
# Circularly shifts a vector by 'shift' positions.
circularShiftVec <- function(x, shift) {
  n <- length(x)
  shift <- shift %% n
  if (shift == 0) return(x)
  c(x[(shift + 1):n], x[1:shift])
}


# Function: cohortCCFShiftPvals
# Permutation test for the cohort-level CCF using all possible circular rotations
# of the sleep series (exhaustive, not random).
# With N points there are N-1 unique shifts, giving a resolution of p = 1/(N-1).
# For N=60 (two cycles), p_min ~ 0.017.
#
# Parameters:
#   grav_vec   : gravity vector (length N)
#   sleep_vec  : sleep vector (length N)
#   max_lag    : maximum lag for the CCF (default 15)
#
# Returns: data.frame with Lag, obs_ACF, pval_perm
cohortCCFShiftPvals <- function(grav_vec, sleep_vec, max_lag = 15) {
  N <- length(grav_vec)
  shifts <- seq_len(N - 1)  # all non-trivial shifts

  # Observed CCF
  cc_obs <- ccf(grav_vec, sleep_vec, max_lag, plot = FALSE)
  obs_acf <- as.numeric(cc_obs$acf)
  lags    <- as.numeric(cc_obs$lag)

  # Matrix of null ACFs [n_shifts x n_lags]
  null_mat <- matrix(NA_real_, nrow = length(shifts), ncol = length(lags))
  for (k in seq_along(shifts)) {
    sleep_sh <- circularShiftVec(sleep_vec, shifts[k])
    cc_null  <- ccf(grav_vec, sleep_sh, max_lag, plot = FALSE)
    null_mat[k, ] <- as.numeric(cc_null$acf)
  }

  # P-value: proportion of |null ACF| >= |obs ACF|
  pvals <- sapply(seq_along(lags), function(j) {
    mean(abs(null_mat[, j]) >= abs(obs_acf[j]))
  })

  data.frame(Lag = lags, obs_ACF = obs_acf, pval_perm = pvals,
             n_shifts = length(shifts))
}


# Function: fitOLSCosinor
# Fits a simple (OLS) cosinor to a profile of n_days points and returns the
# fitted values at the same days.
# Equivalent to the 1-component cosinor but without mixed effects:
#   y ~ mesor + A*cos(2pi*t/T) + B*sin(2pi*t/T)
# Fast (just lm), suitable for running inside permutation loops.
fitOLSCosinor <- function(mean_profile, n_days = 30, n_comp = 1) {
  days <- seq_len(n_days) - 1  # 0:(n_days-1)
  df <- data.frame(
    y  = mean_profile,
    c1 = cos(2 * pi * days / n_days),
    s1 = sin(2 * pi * days / n_days)
  )
  if (n_comp >= 2) {
    # Second harmonic = 15-day component (n_days/2). Needed so the null is
    # comparable with the 30+15 curve overlaid on the plot.
    df$c2 <- cos(2 * pi * days / (n_days / 2))
    df$s2 <- sin(2 * pi * days / (n_days / 2))
    mod <- lm(y ~ c1 + s1 + c2 + s2, data = df)
  } else {
    mod <- lm(y ~ c1 + s1, data = df)
  }
  fitted(mod)
}


# Function: cohortShuffleBands
# Generates confidence bands (95% CI) for the cohort cosinor curve via
# individual-level phase shuffling.
#
# Permutation procedure:
#   1. Circularly shift each individual's DayinCycle by a random offset
#   2. Average across individuals to get the shuffled mean profile
#   3. Fit an OLS cosinor to the shuffled mean profile
#   4. Store the cosinor fitted values
# The bands are built on the distribution of fitted cosinor curves, not on the
# raw mean profiles. This is more consistent with the analysis: it tests whether
# the OBSERVED cosinor curve has an amplitude/phase that could not arise by
# chance given the heterogeneity of individual phases.
#
# Parameters:
#   sleep_df  : data.frame with columns ID_year, year, DayinCycle, and sleep_var
#   sleep_var : name of the column to use (e.g. "Start_Time")
#   cohort    : value of sleep_df$year to filter
#   n_perm    : number of permutations (default 500)
#   n_days    : cycle length (default 30)
#
# Returns: list with $lower, $upper (vectors of length n_days) and $days (0:n_days-1)
cohortShuffleBands <- function(sleep_df, sleep_var, cohort, n_perm = 500, n_days = 30,
                               n_comp = 1) {
  df_cohort <- sleep_df[sleep_df$year == cohort & !is.na(sleep_df[[sleep_var]]), ]

  # Matrix to store the fitted cosinor curve of each permutation [n_perm x n_days]
  perm_curves <- matrix(NA_real_, nrow = n_perm, ncol = n_days)

  for (perm in seq_len(n_perm)) {
    df_perm <- df_cohort

    # Step 1: circularly shift each individual's DayinCycle
    for (id in unique(df_cohort$ID_year)) {
      idx   <- df_perm$ID_year == id
      shift <- sample(0:(n_days - 1), 1)
      df_perm$DayinCycle[idx] <- (df_perm$DayinCycle[idx] + shift) %% n_days
    }

    # Step 2: cohort mean profile with shuffled data
    m_perm <- df_perm %>%
      group_by(DayinCycle) %>%
      summarise(val = mean(.data[[sleep_var]], na.rm = TRUE), .groups = "drop") %>%
      arrange(DayinCycle)

    if (nrow(m_perm) < n_days) next  # skip if days are missing

    # Step 3: fit an OLS cosinor to the shuffled mean profile
    fitted_curve <- tryCatch(
      fitOLSCosinor(m_perm$val, n_days, n_comp),
      error = function(e) rep(NA_real_, n_days)
    )
    perm_curves[perm, ] <- fitted_curve
  }

  list(
    lower = apply(perm_curves, 2, quantile, 0.025, na.rm = TRUE),
    upper = apply(perm_curves, 2, quantile, 0.975, na.rm = TRUE),
    days  = 0:(n_days - 1)
  )
}


# Function: lunarPlot_withShuffleBands
# Extends lunarPlot by adding the phase-shuffling confidence bands as a gray ribbon.
# Lets you visualize whether the observed curve exceeds the variation expected by chance.
lunarPlot_withShuffleBands <- function(data, variables, results, bands,
                                        n = 1, period, color) {
  # Duplicate the bands for the second turn of the cycle (DayinCycle 30-59)
  band_df <- data.frame(
    DayinCycle = c(bands$days, bands$days + 30),
    lower      = rep(bands$lower, 2),
    upper      = rep(bands$upper, 2)
  )

  lunarPlot(data, variables, results, n, period, color) +
    geom_ribbon(data = band_df,
                aes(x = DayinCycle, ymin = lower, ymax = upper),
                alpha = 0.15, fill = "gray50", inherit.aes = FALSE) +
    labs(caption = "Gray band: 95% CI individual phase shuffling")
}


###### Colors for cohorts #######
colors=c("#33DDAA","#33DDAA","#FF7799","#BB66FF","#FFBB67")
colors2=c("#117755","#117755","#CC2299","#6611AA","#EE8122")

# Open file with gravity data
df <- read.csv ("data/Data_S4_meanDaily_Gravitational_variation.csv", header = TRUE)

# Convert NightDate to a Date object (M/D/YYYY format in the gravity file).
# This is needed so the %in% with sleep$night works regardless of the string
# format of each sleep dataset (Seattle: YYYY-MM-DD, Tobas: M/D/YYYY).
df$NightDate <- as.Date(df$NightDate, format = "%m/%d/%Y")

# Subset data by city
tobas_full<-df[df$city=="Ingeniero Juarez",]
seattle_full<-df[df$city=="Seattle",]
titis_full<-df[df$city=="Davis",]


#### Plot: Gravity Variation Over Time (Figure S1) ####
# Visualizes how the amplitude of the semi-lunar tide changes over months.
df<-tobas_full[,c("NightDate","datenum", "DayinCycle", "mean_Total_C")]
data=df[120:300,] # Select a specific date range for visualization

# Add markers for Full (white circle) and New (black circle) moons
data$x_label <- ifelse(
  (data$NightDate %in% c("2016-04-22", "2016-05-21", "2016-06-20", "2016-07-19", "2016-08-18", "2016-09-16", "2016-10-16")),  "○",
  ifelse(data$NightDate %in% c("2016-04-07", "2016-05-06", "2016-06-05", "2016-07-04", "2016-08-02", "2016-09-01", "2016-10-01"), "●", "")
)

data$breaks<-data$x_label!="" 
data$Full<-data$NightDate %in% c("2016-04-22", "2016-05-21", "2016-06-20", "2016-07-19", "2016-08-18", "2016-09-16", "2016-10-16")
data$New<-data$NightDate %in% c("2016-04-07", "2016-05-06", "2016-06-05", "2016-07-04", "2016-08-02", "2016-09-01", "2016-10-01")

tiff(paste("Figures/ChangeInGravity.tif", sep=""), width = 2500, height = 900, res = 300,
     compression = "lzw")
ggplot(data=data)+geom_point(aes(x=datenum, y=mean_Total_C), data=data)+geom_line(aes(x=datenum, y=mean_Total_C), data=data)+
  geom_vline(xintercept = data$datenum[data$Full], color = "black")+
  geom_vline(xintercept = data$datenum[data$New], linetype = "dashed", color = "#888888")+
  scale_x_continuous(
    breaks = data$datenum[data$breaks],
    labels = data$x_label[data$breaks], 
  ) + ylab("Daily-mean gravitational deviation (\U00B5Gal)") +
  xlab("Julian date and night in cycle") +
  theme(
    axis.text.x = element_text(size = 20, hjust = 0.5, family = "Arial"),
    axis.line = element_line(colour = 'black', size = 1),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white",colour = "lightgray", size = 0.5, linetype = "solid"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  )
dev.off()


###### Seattle Analysis #####
# Plot Mean Gravity for Seattle cohorts
plotTitles<-c("Spring 2020","Winter 2021","Spring 2024")
sleep <- read.csv ("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)
sleep$night <- as.Date(sleep$night)  # Seattle: YYYY-MM-DD
sleep<-sleep[sleep$year!="SS2023" & sleep$year!="AW2024",] # Filter cohorts

# Set factor levels for correct ordering in plots
sleep$year <- factor(
  sleep$year,
  levels = c( "SS2020", "AW2020", "SS2024")
) 

sleep<-na.omit(sleep)
sleep$year<-as.factor(sleep$year)
plotlist<-list()

# Generate Gravity Plots for each cohort
for(i in (1:length(levels(sleep$year)))){
  cohort= levels(sleep$year)[i]
  color=colors2[i]
  seattle<-seattle_full[seattle_full$NightDate %in% sleep$night[sleep$year==cohort],]
  
  seattle$variable<-seattle$mean_Total_C
  
  # Pad data for double plotting (0-60 days)
  seattle2<-seattle
  seattle2$DayinCycle<-seattle$DayinCycle+30
  seattle<-rbind(seattle, seattle2)
  
  plot=gravPlot(seattle, "Daily-mean gravitational deviation (\U00B5Gal)", color) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist[[i]]<-plot
}

# Save composite gravity plot
tiff(paste("Figures/seattle_Group_MeanGravity.tif", sep=""), width = 850, height = 2650, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()


###### Cross-correlations at the Sample Level (Seattle) #####

# Load and clean sleep data
sleep <- read.csv ("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)
sleep$night <- as.Date(sleep$night)  # Seattle: YYYY-MM-DD
sleep<-sleep[sleep$year!="SS2023" & sleep$year!="AW2024",]
sleep$year<-as.factor(sleep$year)

# Calculate Daily Means for Sleep Variables
means<-sleep %>% group_by(DayinCycle, year) %>%
  summarise(Start_Time = mean(Start_Time), End_Time = mean(End_Time), Duration=mean(Duration))

# Pad data for circular smoothing
means2<-means
means2$DayinCycle<-means2$DayinCycle+30
means3<-means[means$DayinCycle<3,]
means3$DayinCycle<-means3$DayinCycle+60
means4<-means[means$DayinCycle>26,]
means4$DayinCycle<-means4$DayinCycle-30
means<-rbind(means2, means3, means4, means)

# Interpolate and Smooth (7-day MA)
filling=expand.grid(year=unique(means$year), DayinCycle=c(-3:62))
means<-merge(filling, means, by=c("year", "DayinCycle"), all.x=TRUE)  
means <- means[order(means$year, means$DayinCycle),]
means = means %>%
  group_by (year) %>%
  mutate_at(c("Start_Time","Duration","End_Time", "DayinCycle"), na.approx, na.rm=FALSE)%>%
  ungroup()
means<-na.omit(means)
means = means %>% group_by (year) %>%
  mutate(Start_Time = SMA(Start_Time, 7), End_Time = SMA(End_Time, 7), Duration = SMA(Duration, 7), DayinCycle = SMA(DayinCycle, 7))
means<-na.omit(means)

# Loop for Cross-Correlation Calculation
vars_sleep=c("Start_Time", "End_Time", "Duration")
crossCorrResults=data.frame()

for(i in (1:length(levels(sleep$year)))){
  cohort= levels(sleep$year)[i]
  color=colors[i]
  seattle<-seattle_full[seattle_full$NightDate %in% sleep$night[sleep$year==cohort],]
  means_sleep<-means[means$year==cohort,]
  
  # Prepare Gravity Means
  means_grav<-seattle %>% group_by(DayinCycle) %>%
    summarise(mean_Total_C = mean(mean_Total_C))
  # [Data Padding & Smoothing for Gravity - same logic as sleep]
  means_grav2<-means_grav
  means_grav2$DayinCycle<-means_grav2$DayinCycle+30
  means_grav3<-means_grav[means_grav$DayinCycle<3,]
  means_grav3$DayinCycle<-means_grav3$DayinCycle+60
  means_grav4<-means_grav[means_grav$DayinCycle>26,]
  means_grav4$DayinCycle<-means_grav4$DayinCycle-30
  means_grav<-rbind(means_grav2, means_grav3, means_grav4, means_grav)
  filling=expand.grid(DayinCycle=c(-3:62))
  means_grav<-merge(filling, means_grav, by=c("DayinCycle"), all.x=TRUE)  
  means_grav <- means_grav[order(means_grav$DayinCycle),]
  means_grav = means_grav %>%
    mutate_at(c("mean_Total_C", "DayinCycle"), na.approx, na.rm=FALSE)%>%
    ungroup()
  means_grav<-means_grav[means_grav$DayinCycle>=0 & means_grav$DayinCycle<=59,]
  
  # Merge Sleep and Gravity
  means_mix<-merge(means_sleep, means_grav, by="DayinCycle")
  
  # Calculate CCF + circular shift permutation test
  # With 60 points there are 59 unique shifts -> p resolution ~ 0.017
  for(vars in vars_sleep){
    tiff(paste("Figures/CrossCor/SampleLevel_Seattle_mean_Total_C_",vars,"_",cohort,".tif", sep=""), width = 900, height = 1200, res = 300,
         compression = "lzw") # file & size
    crossCorr<-(ccf(means_mix[,"mean_Total_C"], means_mix[,vars],15))
    dev.off()
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("group"=rep(cohort, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var2" = rep(vars, n), "var1" = rep("mean_Total_C", n)))

    # Permutation test (exhaustive circular shift over the sleep series)
    perm_res <- cohortCCFShiftPvals(means_mix[,"mean_Total_C"], means_mix[,vars], max_lag=15)
    perm_res$cohort <- cohort
    perm_res$var    <- vars
    write.csv(perm_res,
              paste0("results/seattle/sample/CCF_shiftPvals_", cohort, "_", vars, ".csv"),
              row.names = FALSE)
  }

  # Linear Regression on Lagged Data
  # Shifts the data based on the optimal lag observed for each cohort to calculate R-squared.
  means_mix$var<-NA
  if(cohort=="AW2020"){
    # Lag = +2 days
    means_mix$var[3:length(means_mix$Start_Time)]<-means_mix$Start_Time[1:(length(means_mix$Start_Time)-2)]
    means_mix$var[1:2]<-means_mix$Start_Time[59:60]
  }else if(cohort=="SS2020"){
    # Lag = 0 days
    means_mix$var<-means_mix$Start_Time
  }else{
    # Lag = -2 days
    means_mix$var[1:(length(means_mix$Start_Time)-2)]<-means_mix$Start_Time[3:length(means_mix$Start_Time)]
    means_mix$var[59:60]<-means_mix$Start_Time[1:2]
  }
  
  # Print Linear Model Stats
  print(paste(cohort,"and Start_Time"))
  mod<-lm(means_mix[1:30,]$var ~means_mix[1:30,]$mean_Total_C)
  summary(mod)
  print(paste("slope is", coef(mod)[2]))
  print(paste("r-square is",summary(mod)$r.squared))
  
  # Plot Linear Correlation
  tiff(paste("Figures/CrossCor/correlation_SampleLevel_Seattle_mean_Total_C_Start_Time_",cohort,".tif", sep=""), width = 1000, height = 850, res = 300,
       compression = "lzw") 
  print(ggplot(data= means_mix[1:30,], aes(x=mean_Total_C, y=var))+geom_point()+stat_smooth(method="lm", color="black")+
          labs(x = "Daily-mean gravitational deviation (\U00B5Gal)")+
          scale_y_continuous(name = "Sleep onset (clock time)", labels = hours_to_clock)+
          theme(axis.line = element_line(colour = 'black', size = 0.75),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()))
  dev.off()
}

r<-crossCorrResults %>% group_by(group, var1, var2) %>% filter(ACF==max(ACF))


##### Whole sample NO GROUPS #####

setwd(here())


##### Seattle #####
sleep <- read.csv ("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)
sleep$night <- as.Date(sleep$night)  # Seattle: YYYY-MM-DD

sleep$year <- factor(
  sleep$year,
  levels = c( "SS2020", "AW2020", "SS2024")
) 


sleep2<-sleep
sleep2$DayinCycle<-sleep2$DayinCycle+30
sleep3<-sleep[sleep$DayinCycle<3,]
sleep3$DayinCycle<-sleep3$DayinCycle+60
sleep4<-sleep[sleep$DayinCycle>26,]
sleep4$DayinCycle<-sleep4$DayinCycle-30

sleep<-rbind(sleep,sleep2, sleep3, sleep4)

sleep<-sleep[order(sleep$ID_year, sleep$DayinCycle),]

means<-sleep %>% group_by(DayinCycle) %>%
  summarise(Start_Time = mean(Start_Time), End_Time = mean(End_Time), Duration=mean(Duration))


means = means %>%
  mutate(Start_Time = SMA(Start_Time, 7), End_Time = SMA(End_Time, 7), Duration = SMA(Duration, 7), DayinCycle = SMA(DayinCycle, 7)) 

means<-na.omit(means)


####### Pooled (whole-sample) 30+15 model: onset / offset / duration in one column #######

color = "#D55E00"

# Sleep onset
means$variable<-means$Start_Time
results<-read.csv("results/seattle/sample/GLMMCosinor_StartTime_wholeSample_30+15.csv", header = TRUE)
results2<-data.frame(M=mean(means$variable),
                     A=results[5,2], A2=results[6,2],
                     P=transformPhase(results[7,2], 30),
                     P2=transformPhase(results[8,2], 15))
p_onset<-lunarPlot(means, "Sleep onset (hours from midnight)", results2, 2, 30, color); ggplot_build(p_onset)

# Sleep offset
means$variable<-means$End_Time
results<-read.csv("results/seattle/sample/GLMMCosinor_EndTime_wholeSample_30+15.csv", header = TRUE)
results2<-data.frame(M=mean(means$variable),
                     A=results[5,2], A2=results[6,2],
                     P=transformPhase(results[7,2], 30),
                     P2=transformPhase(results[8,2], 15))
p_offset<-lunarPlot(means, "Sleep offset (hours from midnight)", results2, 2, 30, color); ggplot_build(p_offset)

# Sleep duration
means$variable<-means$Duration
results<-read.csv("results/seattle/sample/GLMMCosinor_Duration_wholeSample_30+15.csv", header = TRUE)
results2<-data.frame(M=mean(means$variable),
                     A=results[5,2], A2=results[6,2],
                     P=transformPhase(results[7,2], 30),
                     P2=transformPhase(results[8,2], 15))
p_duration<-lunarPlot(means, "Sleep duration (minutes)", results2, 2, 30, color); ggplot_build(p_duration)

# Column figure (onset / offset / duration; no panel labels)
tiff("Figures/notGrouped/seattle_All_30+15.tif", width = 1050, height = 850 * 3, res = 300, compression = "lzw")
print(ggarrange(plotlist = list(p_onset, p_offset, p_duration), ncol = 1, nrow = 3))
dev.off()




##### SLEEP and GRAVITY by COHORT #####


###### BY GROUP ######
sleep<-sleep %>% group_by(year) %>%
  filter(length(unique(ID_year))>=10)


seattle_full<-seattle_full[seattle_full$NightDate %in% sleep$night,]

seattle_full$night<-seattle_full$NightDate
seattle_full<-as.data.frame(seattle_full)

seattle_full$DayinCycle<-NULL

sleep<-merge(x=sleep, y=seattle_full, by="night", x.all=TRUE)

means<-sleep %>% group_by(DayinCycle, year) %>%
  summarise(Start_Time = mean(Start_Time), End_Time = mean(End_Time), Duration=mean(Duration), mean_Total_C=mean(mean_Total_C))


filling=expand.grid(year=unique(means$year), DayinCycle=c(-3:62))
means<-merge(filling, means, by=c("year", "DayinCycle"), all.x=TRUE)  
means <- means[order(means$year, means$DayinCycle),]

means = means %>%
  group_by (year) %>%
  mutate_at(c("Start_Time","Duration","End_Time", "DayinCycle", "mean_Total_C"), na.approx, na.rm=FALSE)%>%
  ungroup()

means<-na.omit(means)



means = means %>% group_by(year) %>%
  mutate(Start_Time = SMA(Start_Time, 7), End_Time = SMA(End_Time, 7), Duration = SMA(Duration, 7), DayinCycle = SMA(DayinCycle, 7)) 

means$mean_Total_C[4:length(means$mean_Total_C)]<-means$mean_Total_C[1:(length(means$mean_Total_C)-3)]

means<-na.omit(means)


# Full consistent cohort figure set (onset/offset/duration x {30, 30+15} columns and
# the three gravity overlays: lines / points / lines+points) -> Figures/cohort/
makeCohortFigs(means, unique(means$year), c("Spring 2020", "Winter 2021", "Spring 2024"),
               colors, colors2, "results/seattle/sample", "Seattle", "Figures/cohort")


####### 30days plot by cohort ######
resultsAll<-read.csv("results/seattle/sample/GLMMCosinor_StartTime.csv", header = TRUE)


plotlist1<-list()
plotlist2<-list()
plotlist3<-list()
plotlist4<-list()
plotTitles<-c("Spring 2020","Winter 2021","Spring 2024")


for(i in c(1:length(unique(means$year)))){
means$variable<-means$Start_Time
cohort=unique(means$year)[i]
means_cohort<-means[means$year==cohort,]

results<-resultsAll[grepl(cohort, resultsAll$X),]
if(length(results[,1])==3){
  results<-results[2:3,]
}
means_cohort<-as.data.frame(means_cohort)

results_plot<-data.frame(M=mean(means_cohort$variable),
                     A=results[1,2],
                     P=transformPhase(results[2,2], 30))

plot<-lunarPlot(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
plotlist1[[i]]<-print(plot)


if(cohort=="SS2024"){
tiff(paste("Figures/notGrouped/Seattle_new_moon.tif", sep=""),       width = 1050, height = 850, res = 300,
     compression = "lzw") # file & size
print(lunarPlot(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, "#FF7799"))
dev.off()
}

plot<-lunarPlot_grav_nopoints(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
plotlist2[[i]]<-print(plot)

plot<-print(gravPlot(means_cohort, "Sleep onset (hours from midnight)", colors[i])) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
plotlist3[[i]]<-print(plot)

plot<-lunarPlot_grav_nopoints2(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
plotlist4[[i]]<-print(plot)
}


tiff(paste("Figures/seattle_Group_startTime.tif", sep=""), width = 850, height = 2650, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()


tiff(paste("Figures/seattle_Group_startTime_gravity_nopoint.tif", sep=""),width = 1050, height = 2650, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()



tiff(paste("Figures/seattle_Group_startTime_noModel.tif", sep=""), width = 850, height = 2650, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist3 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()

tiff(paste("Figures/seattle_Group_startTime_gravity_nopoint_noModel.tif", sep=""),width = 1050, height = 2650, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist4 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()

tiff(paste("Figures/seattle_Fullmoon_Group_startTime_gravity_nopoint.tif", sep=""),width = 1050, height = 850, res = 300,
     compression = "lzw") # file & size
print(plotlist4[1])
dev.off()

tiff(paste("Figures/seattle_NewMoon_Group_startTime_gravity_nopoint.tif", sep=""),width = 1050, height = 850, res = 300,
     compression = "lzw") # file & size
print(plotlist4[3])
dev.off()


####### 30+15 plots by cohort: Sleep Onset (Seattle) #######
# Equivalent to the 30-day plots above but using the combined 30+15 model.
# The overlaid sinusoidal curve is the sum of both components (lunar + semi-lunar),
# better reflecting the bimodal structure expected under the gravitational hypothesis.
# REQUIRES having run 2_analysis.R to generate GLMMCosinor_StartTime_30+15.csv

resultsAll_30_15 <- read.csv("results/seattle/sample/GLMMCosinor_StartTime_30+15.csv",
                              header = TRUE)

# Load raw sleep here to reuse in the 30+15 shuffle and in the next section
sleep_raw_seattle_tmp <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)
sleep_raw_seattle_tmp$year <- as.character(sleep_raw_seattle_tmp$year)

plotlist_30_15        <- list()
plotlist_30_15_shuffle <- list()
set.seed(42)

for (i in c(1:length(unique(means$year)))) {
  means$variable <- means$Start_Time
  cohort         <- unique(means$year)[i]
  means_cohort   <- as.data.frame(means[means$year == cohort, ])

  results_plot <- extractResults_30_15(resultsAll_30_15, cohort, means_cohort)
  if (is.null(results_plot)) next

  # 30+15 plot without bands
  plot <- lunarPlot(means_cohort, "Sleep onset (hours from midnight)",
                    results_plot, n = 2, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_30_15[[i]] <- print(plot)

  # 30+15 plot with shuffle bands (OLS cosinor over permutations)
  bands <- cohortShuffleBands(sleep_raw_seattle_tmp, "Start_Time", cohort,
                               n_perm = 500, n_days = 30, n_comp = 2)
  plot_sh <- lunarPlot_withShuffleBands(means_cohort,
                                         "Sleep onset (hours from midnight)",
                                         results_plot, bands,
                                         n = 2, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_30_15_shuffle[[i]] <- print(plot_sh)
}

tiff("Figures/seattle_Group_startTime_30+15.tif",
     width = 850, height = 2650, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_30_15, ncol = 1, nrow = 3, labels = c("", "", ""))
print(plot)
dev.off()

tiff("Figures/seattle_Group_startTime_30+15_shuffleBands.tif",
     width = 850, height = 2650, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_30_15_shuffle, ncol = 1, nrow = 3, labels = c("", "", ""))
print(plot)
dev.off()


####### Phase Shuffling: CI bands for the cohort mean (Seattle, Sleep Onset) #######
# Generates permutational confidence bands by circularly shifting each individual's
# DayinCycle before averaging. This creates a null distribution for the cohort mean
# curve, letting you visualize whether the observed pattern exceeds what is expected
# by chance.
#
# NOTE: sleep_raw (data without the double-cycle duplication) is used for the shuffling.
# The results are overlaid on the 30-day plots (without the 30+15 model so as not to
# overload the figure) and on the 30+15 ones.

sleep_raw_seattle <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE)
sleep_raw_seattle$night <- as.Date(sleep_raw_seattle$night)  # Seattle: YYYY-MM-DD
sleep_raw_seattle <- sleep_raw_seattle[sleep_raw_seattle$year %in%
                                         levels(sleep$year), ]
sleep_raw_seattle$year <- factor(sleep_raw_seattle$year, levels = levels(sleep$year))

set.seed(42)
plotlist_shuffle <- list()

for (i in c(1:length(levels(sleep_raw_seattle$year)))) {
  cohort       <- levels(sleep_raw_seattle$year)[i]
  means_cohort <- as.data.frame(means[means$year == cohort, ])
  means_cohort$variable <- means_cohort$Start_Time

  # Compute CI bands (may take ~1-2 min per cohort with n_perm=500)
  bands <- cohortShuffleBands(sleep_raw_seattle, "Start_Time", cohort,
                               n_perm = 500, n_days = 30)

  results_plot <- extractResults_30(
    read.csv("results/seattle/sample/GLMMCosinor_StartTime.csv", header = TRUE),
    cohort, means_cohort)

  plot <- lunarPlot_withShuffleBands(means_cohort,
                                     "Sleep onset (hours from midnight)",
                                     results_plot, bands,
                                     n = 1, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_shuffle[[i]] <- print(plot)
}

tiff("Figures/seattle_Group_startTime_shuffleBands.tif",
     width = 850, height = 2650, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_shuffle, ncol = 1, nrow = 3, labels = c("", "", ""))
print(plot)
dev.off()


###### Sleep offset ######
resultsAll<-read.csv("results/seattle/sample/GLMMCosinor_EndTime_30+15.csv", header = TRUE)



plotlist1<-list()
plotlist2<-list()

for(i in c(1:length(unique(means$year)))){
  means$variable<-means$End_Time
  cohort=unique(means$year)[i]
  means_cohort<-means[means$year==cohort,]
  
  means_cohort<-as.data.frame(means_cohort)

  results_plot <- extractResults_30_15(resultsAll, cohort, means_cohort)
  if (is.null(results_plot)) next

  plot<-lunarPlot(means_cohort, "Sleep offset (hours from midnight)", results_plot, 2, 30, colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist1[[i]]<-print(plot)

  plot<-lunarPlot_grav_modelpoints(means_cohort, "Sleep offset (hours from midnight)", results_plot, 2, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist2[[i]]<-print(plot)
}


svg(paste("Figures/seattle_Group_endTime.svg", sep=""), width = 2.8333, height = 8.8333)
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()


svg(paste("Figures/seattle_Group_endTime_gravity_nopoint.svg", sep=""),width = 3.5, height = 8.8333)
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()



###### Duration ######

resultsAll<-read.csv("results/seattle/sample/GLMMCosinor_Duration_30+15.csv", header = TRUE)



plotlist1<-list()
plotlist2<-list()

for(i in c(1:length(unique(means$year)))){
  means$variable<-means$Duration
  cohort=unique(means$year)[i]
  means_cohort<-means[means$year==cohort,]
  
  means_cohort<-as.data.frame(means_cohort)

  results_plot <- extractResults_30_15(resultsAll, cohort, means_cohort)
  if (is.null(results_plot)) next

  plot<-lunarPlot(means_cohort, "Sleep duration (minutes)", results_plot, 2, 30, colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist1[[i]]<-print(plot)

  plot<-lunarPlot_grav_modelpoints(means_cohort, "Sleep duration (minutes)", results_plot, 2, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist2[[i]]<-print(plot)
}


svg(paste("Figures/seattle_Group_Duration.svg", sep=""), width = 2.8333, height = 8.8333)
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()


svg(paste("Figures/seattle_Group_Duration_gravity_nopoint.svg", sep=""),width = 3.5, height = 8.8333)
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 3, labels = c("","",""))
print(plot)
dev.off()




###### Tobas #####
plotTitles<-c("Winter 2016","Fall 2023","Spring 2018","SS 2023/2024","Fall 2024")
# colors=c("#33DDAA","#55BBFF","#FF7799","#FFBB67","#BB66AA","#BB66FF")
# 
# colors2=c("#117755","#1133AA","#CC2299","#EE8122","#9922AA","#6611AA")


colors=c("#33DDAA","#55BBFF","#55BBFF","#BB66FF", "#FF7799","#FF7799")

colors2=c("#117755","#1133AA","#1133AA","#6611AA","#CC2299","#CC2299")

sleep <- read.csv ("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)
sleep$night <- as.Date(sleep$night, format = "%m/%d/%Y")  # Tobas: M/D/YYYY
sleep<-sleep[sleep$year!="Spring2017",] #n<10

sleep$year<-factor(sleep$year, levels = c("Fall-Winter2016", "Fall-Winter2023", "Spring2018", "Spring2023", "Summer2024", "Fall-Winter2024"))
plotTitles<-c("Winter 2016","Fall 2023","Spring 2018","Spring 2023", "Summer 2024","Fall 2024")


#sleep$year<-factor(sleep$year, levels = c("Fall-Winter2016", "Fall-Winter2023", "Spring2018", "Spring-Summer2023", "Fall-Winter2024"))



plotlist<-list()
for(i in c(1:length(levels(sleep$year)))){
  cohort= levels(sleep$year)[i]
  color=colors2[i]
  tobas<-tobas_full[tobas_full$NightDate %in% sleep$night[sleep$year==cohort],]
  
  tobas2<-tobas
  tobas2$DayinCycle<-tobas$DayinCycle+30
  tobas<-rbind(tobas, tobas2)
  
  tobas$variable<-tobas$mean_Total_C
  plot<-gravPlot(tobas, "Daily-mean gravitational deviation (\U00B5Gal)", color) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist[[i]]<-plot
  
}

tiff(paste("Figures/Tobas_Group_MeanGravity.tif", sep=""), width = 850, height = 5000, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()


###### Cross-correlations at the sample level by cohort #####

###### sleep means #####

means<-sleep %>% group_by(DayinCycle, year) %>%
  summarise(Start_Time = mean(Start_Time), End_Time = mean(End_Time), Duration=mean(Duration))

means2<-means
means2$DayinCycle<-means2$DayinCycle+30
means3<-means[means$DayinCycle<3,]
means3$DayinCycle<-means3$DayinCycle+60
means4<-means[means$DayinCycle>26,]
means4$DayinCycle<-means4$DayinCycle-30

means<-rbind(means2, means3, means4, means)

filling=expand.grid(year=unique(means$year), DayinCycle=c(-3:62))
means<-merge(filling, means, by=c("year", "DayinCycle"), all.x=TRUE)  
means <- means[order(means$year, means$DayinCycle),]

means = means %>%
  group_by (year) %>%
  mutate_at(c("Start_Time","Duration","End_Time", "DayinCycle"), na.approx, na.rm=FALSE)%>%
  ungroup()



means<-na.omit(means)

means = means %>% group_by (year) %>%
  mutate(Start_Time = SMA(Start_Time, 7), End_Time = SMA(End_Time, 7), Duration = SMA(Duration, 7), DayinCycle = SMA(DayinCycle, 7)) 

means<-na.omit(means)

missingRows<-means[means$year=="Summer2024" & (means$DayinCycle==30 | means$DayinCycle==31),]
missingRows$DayinCycle<-missingRows$DayinCycle-30
means<-rbind(means, missingRows)

vars_sleep=c("Start_Time", "End_Time", "Duration")

crossCorrResults=data.frame()

for(i in (1:length(levels(sleep$year)))){
  cohort= levels(sleep$year)[i]
  color=colors[i]
  tobas<-tobas_full[tobas_full$NightDate %in% sleep$night[sleep$year==cohort],]
  means_sleep<-means[means$year==cohort,]
  
  
  #### Gravity means ####
  means_grav<-tobas %>% group_by(DayinCycle) %>%
    summarise(mean_Total_C = mean(mean_Total_C))
  
  means_grav2<-means_grav
  means_grav2$DayinCycle<-means_grav2$DayinCycle+30
  means_grav3<-means_grav[means_grav$DayinCycle<3,]
  means_grav3$DayinCycle<-means_grav3$DayinCycle+60
  means_grav4<-means_grav[means_grav$DayinCycle>26,]
  means_grav4$DayinCycle<-means_grav4$DayinCycle-30
  
  means_grav<-rbind(means_grav2, means_grav3, means_grav4, means_grav)
  
  filling=expand.grid(DayinCycle=c(-3:62))
  means_grav<-merge(filling, means_grav, by=c("DayinCycle"), all.x=TRUE)  
  means_grav <- means_grav[order(means_grav$DayinCycle),]
  
  means_grav = means_grav %>%
    mutate_at(c("mean_Total_C", "DayinCycle"), na.approx, na.rm=FALSE)%>%
    ungroup()
  means_grav<-means_grav[means_grav$DayinCycle>=0 & means_grav$DayinCycle<=59,]
  
  means_mix<-merge(means_sleep, means_grav, by="DayinCycle")
  
  for(vars in vars_sleep){
    tiff(paste("Figures/CrossCor/SampleLevel_Tobas_MeanTotalC_",vars,"_",cohort,".tif", sep=""), width = 900, height = 1200, res = 300,
         compression = "lzw") # file & size
    crossCorr<-(ccf(means_mix[,"mean_Total_C"], means_mix[,vars],15))
    dev.off()
    n = length(crossCorr$acf)
    crossCorrResults<-rbind(crossCorrResults, data.frame("group"=rep(cohort, n), "ACF"= crossCorr$acf, "Lag"= crossCorr$lag, "var2" = rep(vars, n), "var1" = rep("mean_Total_C", n)))

    # Permutation test (exhaustive circular shift over the sleep series)
    perm_res <- cohortCCFShiftPvals(means_mix[,"mean_Total_C"], means_mix[,vars], max_lag=15)
    perm_res$cohort <- cohort
    perm_res$var    <- vars
    write.csv(perm_res,
              paste0("results/Tobas/sample/CCF_shiftPvals_", cohort, "_", vars, ".csv"),
              row.names = FALSE)
  }


  means_mix$var<-NA
  if(cohort=="Fall-Winter2016"){
    means_mix$var[6:length(means_mix$Start_Time)]<-means_mix$Start_Time[1:(length(means_mix$Start_Time)-5)]
    means_mix$var[1:5]<-means_mix$Start_Time[56:60]
  }else if(cohort=="Fall-Winter2023"){
    means_mix$var[12:length(means_mix$Start_Time)]<-means_mix$Start_Time[1:(length(means_mix$Start_Time)-11)]
    means_mix$var[1:11]<-means_mix$Start_Time[50:60]
  }else if(cohort=="Spring2018"){
    means_mix$var[10:length(means_mix$Start_Time)]<-means_mix$Start_Time[1:(length(means_mix$Start_Time)-9)]
    means_mix$var[1:9]<-means_mix$Start_Time[52:60]
  }else if(cohort=="Spring2023"){
    means_mix$var[14:length(means_mix$Start_Time)]<-means_mix$Start_Time[1:(length(means_mix$Start_Time)-13)]
    means_mix$var[1:13]<-means_mix$Start_Time[48:60]
  }else if(cohort=="Summer2024"){
    means_mix$var[1:(length(means_mix$Start_Time)-3)]<-means_mix$Start_Time[4:length(means_mix$Start_Time)]
    means_mix$var[58:60]<-means_mix$Start_Time[1:3]
  }else{
    means_mix$var[1:(length(means_mix$Start_Time)-7)]<-means_mix$Start_Time[8:length(means_mix$Start_Time)]
    means_mix$var[54:60]<-means_mix$Start_Time[1:7]
  }
  
  print(paste(cohort,"and Start_Time"))
  
  mod<-lm(means_mix[1:30,]$var ~means_mix[1:30,]$mean_Total_C)
  summary(mod)
  print(paste("slope is", coef(mod)[2]))        # correlation coefficient
  print(paste("95% CI is", confint(mod)["means_mix[1:30, ]$mean_Total_C", ]))
  print(paste("p value is", summary(mod)$coefficients[2, 4]))         # p-value
  print(paste("r-square is",summary(mod)$r.squared)) # R-squared
  
  tiff(paste("Figures/CrossCor/correlation_SampleLevel_Tobas_mean_Total_C_Onset_",cohort,".tif", sep=""), width = 1000, height = 850, res = 300,
       compression = "lzw") # file & size
  print(ggplot(data= means_mix[1:30,], aes(x=mean_Total_C, y=var))+geom_point()+stat_smooth(method="lm", color="black")+
          labs(x = "Daily-mean gravitational deviation (\U00B5Gal)")+
          scale_y_continuous(name = "Sleep onset (clock time)", labels = hours_to_clock)+
    theme(axis.line = element_line(colour = 'black', size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()))
  dev.off()

}
r<-crossCorrResults %>% group_by(group, var1, var2) %>% filter(ACF==max(ACF))


###### TOBAS #####

sleep2<-sleep
sleep2$DayinCycle<-sleep2$DayinCycle+30
sleep3<-sleep[sleep$DayinCycle<3,]
sleep3$DayinCycle<-sleep3$DayinCycle+60
sleep4<-sleep[sleep$DayinCycle>26,]
sleep4$DayinCycle<-sleep4$DayinCycle-30

sleep<-rbind(sleep,sleep2, sleep3, sleep4)

sleep<-sleep[order(sleep$ID_year, sleep$DayinCycle),]
#sleep$year<-factor(sleep$year, levels = c("Fall-Winter2016", "Fall-Winter2023", "Spring2018", "Spring-Summer2023", "Fall-Winter2024"))

### Added ###
plotTitles<-c("Winter 2016","Fall 2023","Spring 2018","Spring 2023", "Summer 2024","Fall 2024")


means<-sleep %>% group_by(DayinCycle) %>%
  summarise(Start_Time = mean(Start_Time), End_Time = mean(End_Time), Duration=mean(Duration))


means = means %>%
  mutate(Start_Time = SMA(Start_Time, 7), End_Time = SMA(End_Time, 7), Duration = SMA(Duration, 7), DayinCycle = SMA(DayinCycle, 7)) 

means<-na.omit(means)



##### Whole sample NO GROUPS ###

##### Pooled (whole-sample) 30+15 model: onset / offset / duration in one column #####
color = "#D55E00"

# Sleep onset
means$variable<-means$Start_Time
results<-read.csv("results/Tobas/sample/GLMMCosinor_StartTime_wholeSample_30+15.csv", header = TRUE)
results2<-data.frame(M=mean(means$variable),
                     A=results[5,2], A2=results[6,2],
                     P=transformPhase(results[7,2], 30),
                     P2=transformPhase(results[8,2], 15))
p_onset<-lunarPlot(means, "Sleep onset (hours from midnight)", results2, 2, 30, color); ggplot_build(p_onset)

# Sleep offset
means$variable<-means$End_Time
results<-read.csv("results/Tobas/sample/GLMMCosinor_EndTime_wholeSample_30+15.csv", header = TRUE)
results2<-data.frame(M=mean(means$variable),
                     A=results[5,2], A2=results[6,2],
                     P=transformPhase(results[7,2], 30),
                     P2=transformPhase(results[8,2], 15))
p_offset<-lunarPlot(means, "Sleep offset (hours from midnight)", results2, 2, 30, color); ggplot_build(p_offset)

# Sleep duration
means$variable<-means$Duration
results<-read.csv("results/Tobas/sample/GLMMCosinor_Duration_wholeSample_30+15.csv", header = TRUE)
results2<-data.frame(M=mean(means$variable),
                     A=results[5,2], A2=results[6,2],
                     P=transformPhase(results[7,2], 30),
                     P2=transformPhase(results[8,2], 15))
p_duration<-lunarPlot(means, "Sleep duration (minutes)", results2, 2, 30, color); ggplot_build(p_duration)

# Column figure (onset / offset / duration; no panel labels)
tiff("Figures/notGrouped/Tobas_All_30+15.tif", width = 1050, height = 850 * 3, res = 300, compression = "lzw")
print(ggarrange(plotlist = list(p_onset, p_offset, p_duration), ncol = 1, nrow = 3))
dev.off()


###### SLEEP and GRAVITY by COHORT ######
sleep<-sleep %>% group_by(year) %>%
  filter(length(unique(ID_year))>=10)


unique(sleep$year)

tobas_full<-tobas_full[tobas_full$NightDate %in% sleep$night,]

tobas_full$night<-tobas_full$NightDate
tobas_full<-as.data.frame(tobas_full)

tobas_full$DayinCycle<-NULL

sleep<-merge(x=sleep, y=tobas_full, by="night", x.all=TRUE)

means<-sleep %>% group_by(DayinCycle, year) %>%
  summarise(Start_Time = mean(Start_Time), End_Time = mean(End_Time), Duration=mean(Duration), mean_Total_C=mean(mean_Total_C))


filling=expand.grid(year=unique(means$year), DayinCycle=c(-3:62))
means<-merge(filling, means, by=c("year", "DayinCycle"), all.x=TRUE)  
means <- means[order(means$year, means$DayinCycle),]

means = means %>%
  group_by (year) %>%
  mutate_at(c("Start_Time","Duration","End_Time", "DayinCycle", "mean_Total_C"), na.approx, na.rm=FALSE)%>%
  ungroup()

means<-na.omit(means)



means = means %>% group_by(year) %>%
  mutate(Start_Time = SMA(Start_Time, 7), End_Time = SMA(End_Time, 7), Duration = SMA(Duration, 7), DayinCycle = SMA(DayinCycle, 7)) 

means$mean_Total_C[4:length(means$mean_Total_C)]<-means$mean_Total_C[1:(length(means$mean_Total_C)-3)]

means<-na.omit(means)


# Sleep onset: split into two groups of three cohorts
# grpA = Fall 2023 + Spring 2018 + Spring 2023 (cohorts 2,3,4)
# grpB = Winter 2016 + Summer 2024 + Fall 2024 (1,5,6)
# Sleep offset and Duration: all 6 cohorts together
tobas_cohorts6 <- c("Fall-Winter2016", "Fall-Winter2023", "Spring2018", "Spring2023", "Summer2024", "Fall-Winter2024")
tobas_titles6  <- c("Winter 2016", "Fall 2023", "Spring 2018", "Spring 2023", "Summer 2024", "Fall 2024")
tobas_gA <- c(2, 3, 4); tobas_gB <- c(1, 5, 6)

onset_spec <- list(list(col = "Start_Time", tag = "onset", csv = "StartTime", ylab = "Sleep onset (hours from midnight)"))
makeCohortFigs(means, tobas_cohorts6[tobas_gA], tobas_titles6[tobas_gA], colors[tobas_gA], colors2[tobas_gA],
               "results/Tobas/sample", "Tobas_grpA", "Figures/cohort", var_specs = onset_spec)
makeCohortFigs(means, tobas_cohorts6[tobas_gB], tobas_titles6[tobas_gB], colors[tobas_gB], colors2[tobas_gB],
               "results/Tobas/sample", "Tobas_grpB", "Figures/cohort", var_specs = onset_spec,
               grav_lab = "Mean daily gravity variation (\U00B5Gal)")  # Fig.3: keep current axis term unchanged

makeCohortFigs(means, tobas_cohorts6, tobas_titles6, colors, colors2,
               "results/Tobas/sample", "Tobas", "Figures/cohort",
               var_specs = list(
                 list(col = "End_Time",  tag = "offset",   csv = "EndTime",  ylab = "Sleep offset (hours from midnight)"),
                 list(col = "Duration",  tag = "duration", csv = "Duration", ylab = "Sleep duration (minutes)")))




##### Start Time #####
resultsAll<-read.csv("results/Tobas/sample/GLMMCosinor_StartTime.csv", header = TRUE)

plotlist1<-list()
plotlist2<-list()
plotlist3<-list()
plotlist4<-list()
#plotTitles<-c("Winter 2016","Fall 2023","Spring 2018","SS 2023/2024","Fall 2024")
plotTitles<-c("Winter 2016","Fall 2023","Spring 2018","Spring 2023", "Summer 2024","Fall 2024")

for(i in c(1:length(unique(means$year)))){
  means$variable<-means$Start_Time
  cohort=unique(means$year)[i]
  means_cohort<-means[means$year==cohort,]
  
  results<-resultsAll[grepl(cohort, resultsAll$X),]
  if(length(results[,1])==3){
    results<-results[2:3,]
  }
  means_cohort<-as.data.frame(means_cohort)
  
  results_plot<-data.frame(M=mean(means_cohort$variable),
                           A=results[1,2],
                           P=transformPhase(results[2,2], 30))
  
  plot<-lunarPlot(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist1[[i]]<-print(plot)
  
  if(cohort=="Fall-Winter2024"){
    tiff(paste("Figures/notGrouped/Tobas_new_moon.tif", sep=""), width = 1050, height = 850, res = 300,
         compression = "lzw") # file & size
    print(lunarPlot(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, "#FF7799"))
    dev.off()
  }
  
  plot<-lunarPlot_grav_nopoints(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist2[[i]]<-print(plot)
  
  plot<-gravPlot(means_cohort, "Sleep onset (hours from midnight)", colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist3[[i]]<-print(plot)
  
  plot<-lunarPlot_grav_nopoints2(means_cohort, "Sleep onset (hours from midnight)", results_plot, 1, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist4[[i]]<-print(plot)
}


tiff(paste("Figures/Tobas_Group_startTime.tif", sep=""), width = 850, height = 5000, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()


tiff(paste("Figures/Tobas_Group_startTime_gravity_nopoint.tif", sep=""),width = 1050, height = 5000, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()




tiff(paste("Figures/Tobas_Group_startTime_noModel.tif", sep=""), width = 850, height = 5000, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist3 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()

tiff(paste("Figures/Tobas_Group_startTime_gravity_nopoint_noModel.tif", sep=""),width = 1050, height = 5000, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist4 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()

tiff(paste("Figures/Tobas_Fullmoon_Group_startTime_gravity_nopoint.tif", sep=""),width = 1050, height = 850, res = 300,
     compression = "lzw") # file & size
print(plotlist4[1])
dev.off()

tiff(paste("Figures/Tobas_NewMoon_Group_startTime_gravity_nopoint.tif", sep=""),width = 1050, height = 850, res = 300,
     compression = "lzw") # file & size
print(plotlist4[6])
dev.off()


####### 30+15 plots by cohort: Sleep Onset (Tobas) #######
resultsAll_30_15 <- read.csv("results/Tobas/sample/GLMMCosinor_StartTime_30+15.csv",
                              header = TRUE)

sleep_raw_tobas_tmp <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)
sleep_raw_tobas_tmp <- sleep_raw_tobas_tmp[sleep_raw_tobas_tmp$year != "Spring2017", ]
sleep_raw_tobas_tmp$year <- as.character(sleep_raw_tobas_tmp$year)

plotlist_30_15        <- list()
plotlist_30_15_shuffle <- list()
set.seed(42)

for (i in c(1:length(unique(means$year)))) {
  means$variable <- means$Start_Time
  cohort         <- unique(means$year)[i]
  means_cohort   <- as.data.frame(means[means$year == cohort, ])

  results_plot <- extractResults_30_15(resultsAll_30_15, cohort, means_cohort)
  if (is.null(results_plot)) next

  # Without bands
  plot <- lunarPlot(means_cohort, "Sleep onset (hours from midnight)",
                    results_plot, n = 2, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_30_15[[i]] <- print(plot)

  # With shuffle bands (OLS cosinor)
  bands <- cohortShuffleBands(sleep_raw_tobas_tmp, "Start_Time", cohort,
                               n_perm = 500, n_days = 30)
  plot_sh <- lunarPlot_withShuffleBands(means_cohort,
                                         "Sleep onset (hours from midnight)",
                                         results_plot, bands,
                                         n = 2, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_30_15_shuffle[[i]] <- print(plot_sh)
}

tiff("Figures/Tobas_Group_startTime_30+15.tif",
     width = 850, height = 5000, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_30_15, ncol = 1, nrow = 6, labels = rep("", 6))
print(plot)
dev.off()

tiff("Figures/Tobas_Group_startTime_30+15_shuffleBands.tif",
     width = 850, height = 5000, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_30_15_shuffle, ncol = 1, nrow = 6, labels = rep("", 6))
print(plot)
dev.off()


####### Phase Shuffling: CI bands for the cohort mean (Tobas, Sleep Onset) #######
sleep_raw_tobas <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE)
sleep_raw_tobas$night <- as.Date(sleep_raw_tobas$night, format = "%m/%d/%Y")  # Tobas: M/D/YYYY
sleep_raw_tobas <- sleep_raw_tobas[sleep_raw_tobas$year != "Spring2017", ]
sleep_raw_tobas$year <- factor(sleep_raw_tobas$year,
                                levels = c("Fall-Winter2016", "Fall-Winter2023",
                                           "Spring2018", "Spring2023",
                                           "Summer2024", "Fall-Winter2024"))

set.seed(42)
plotlist_shuffle <- list()
resultsAll_30_tobas <- read.csv("results/Tobas/sample/GLMMCosinor_StartTime.csv",
                                 header = TRUE)

for (i in c(1:length(levels(sleep_raw_tobas$year)))) {
  cohort       <- levels(sleep_raw_tobas$year)[i]
  means_cohort <- as.data.frame(means[means$year == cohort, ])
  if (nrow(means_cohort) == 0) next
  means_cohort$variable <- means_cohort$Start_Time

  bands <- cohortShuffleBands(sleep_raw_tobas, "Start_Time", cohort,
                               n_perm = 500, n_days = 30)

  results_plot <- extractResults_30(resultsAll_30_tobas, cohort, means_cohort)

  plot <- lunarPlot_withShuffleBands(means_cohort,
                                     "Sleep onset (hours from midnight)",
                                     results_plot, bands,
                                     n = 1, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_shuffle[[i]] <- print(plot)
}

tiff("Figures/Tobas_Group_startTime_shuffleBands.tif",
     width = 850, height = 5000, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_shuffle, ncol = 1, nrow = 6, labels = rep("", 6))
print(plot)
dev.off()


##### End Time #####
resultsAll<-read.csv("results/Tobas/sample/GLMMCosinor_EndTime_30+15.csv", header = TRUE)


plotlist1<-list()
plotlist2<-list()

for(i in c(1:length(unique(means$year)))){
  means$variable<-means$End_Time
  cohort=unique(means$year)[i]
  means_cohort<-means[means$year==cohort,]
  
  means_cohort<-as.data.frame(means_cohort)

  results_plot <- extractResults_30_15(resultsAll, cohort, means_cohort)
  if (is.null(results_plot)) next

  plot<-lunarPlot(means_cohort, "Sleep offset (hours from midnight)", results_plot, 2, 30, colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist1[[i]]<-print(plot)

  plot<-lunarPlot_grav_modelpoints(means_cohort, "Sleep offset (hours from midnight)", results_plot, 2, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist2[[i]]<-print(plot)
}


svg(paste("Figures/Tobas_Group_endTime.svg", sep=""), width = 2.8333, height = 16.6667)
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()


svg(paste("Figures/tobas_Group_endTime_gravity_nopoint.svg", sep=""),width = 3.5, height = 16.6667)
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()


##### Duration #####
resultsAll<-read.csv("results/Tobas/sample/GLMMCosinor_Duration_30+15.csv", header = TRUE)


plotlist1<-list()
plotlist2<-list()

for(i in c(1:length(unique(means$year)))){
  means$variable<-means$Duration
  cohort=unique(means$year)[i]
  means_cohort<-means[means$year==cohort,]
  
  means_cohort<-as.data.frame(means_cohort)

  results_plot <- extractResults_30_15(resultsAll, cohort, means_cohort)
  if (is.null(results_plot)) next

  plot<-lunarPlot(means_cohort, "Sleep duration (minutes)", results_plot, 2, 30, colors[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist1[[i]]<-print(plot)

  plot<-lunarPlot_grav_modelpoints(means_cohort, "Sleep duration (minutes)", results_plot, 2, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist2[[i]]<-print(plot)
}


svg(paste("Figures/Tobas_Group_Duration.svg", sep=""), width = 2.8333, height = 16.6667)
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()


svg(paste("Figures/Tobas_Group_Duration_gravity_nopoint.svg", sep=""),width = 3.5, height = 16.6667)
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 6, labels = c("","","","",""))
print(plot)
dev.off()


######## Monkeys ######
colors=c("#33DDAA", "#FF7799")
colors2=c("#117755","#CC2299")

data<-read.csv("data/Data_S3_all_onsetData_titi.csv", header = TRUE)
data<-data[,c(2:(length(colnames(data))))]


###### Gravity and Sleep #####

resultsAll<-read.csv("results/Titis/GLMMCosinor_ActivityOnset.csv", header = TRUE)

means<-data%>% group_by(DayinCycle, Year) %>%
  summarise(variable=mean(time), mean_Total_C=mean(mean_Total_C))


means1<-means
means1$DayinCycle<-means1$DayinCycle+30
means2<-means[means$DayinCycle<=2,]
means2$DayinCycle<-means2$DayinCycle+60
means3<-means[means$DayinCycle>=27,]
means3$DayinCycle<-means3$DayinCycle-30

means<-rbind(means3,means,means1,means2)

filling=expand.grid(DayinCycle=c(-3:62), Year=c(2023, 2024))
means<-merge(filling, means, by=c("DayinCycle", "Year"), all.x=TRUE)  
means <- means[order(means$Year, means$DayinCycle),]

means = means %>% group_by(Year) %>%
  mutate_at(c("mean_Total_C", "variable"), na.approx, na.rm=FALSE)%>%
  ungroup()


means = means %>%  group_by(Year) %>%
  mutate(variable = SMA(variable, 7), DayinCycle=SMA(DayinCycle,7))
means$mean_Total_C[4:length(means$mean_Total_C)]<-means$mean_Total_C[1:(length(means$mean_Total_C)-3)]
means<-na.omit(means)

means<-as.data.frame(means)
means$DayinCycle<-as.numeric(means$DayinCycle)
means$variable<-as.numeric(means$variable)


# Full consistent cohort figure set (titi: activity onset only) -> Figures/cohort/
means_titi <- means
means_titi$year  <- means_titi$Year
means_titi$onset <- means_titi$variable
makeCohortFigs(means_titi, unique(means_titi$year), c("Summer 2023", "Spring 2024"),
               colors, colors2, "results/Titis", "Titi", "Figures/cohort",
               var_specs = list(list(col = "onset", tag = "onset", csv = "ActivityOnset",
                                      ylab = "Activity onset (hours from midnight)")))


plotlist1<-list()
plotlist2<-list()
plotlist3<-list()
plotTitles<-c("Summer 2023","Spring 2024")

crossCorrResults=data.frame()

for(i in c(1:length(unique(means$Year)))){
  cohort=unique(means$Year)[i]
  means_cohort<-means[means$Year==cohort,]
  
  results<-resultsAll[grepl(cohort, resultsAll$X),]
  if(length(results[,1])==3){
    results<-results[2:3,]
  }
  means_cohort<-as.data.frame(means_cohort)
  
  results_plot<-data.frame(M=mean(means_cohort$variable),
                           A=results[1,2],
                           P=transformPhase(results[2,2], 30))
  
  plot<-lunarPlot(means_cohort, "Activity onset (hours from midnight)", results_plot, 1, 30, colors[i])+ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist1[[i]]<-print(plot)
  
  plot<-lunarPlot_grav_nopoints(means_cohort, "Activity onset (hours from midnight)", results_plot, 1, 30, colors[i], "mean_Total_C", "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist2[[i]]<-print(plot)
  
  means_cohort$onset<-means_cohort$variable
  means_cohort$variable<-means_cohort$mean_Total_C
  plot= gravPlot(means_cohort, "Daily-mean gravitational deviation (\U00B5Gal)", colors2[i]) +ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist3[[i]]<-print(plot)
  
  
}


tiff(paste("Figures/Titis_ActivityOnset.tif", sep=""), width = 850, height = 1800, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist1 , ncol = 1, nrow= 2, labels = c("",""))
print(plot)
dev.off()


tiff(paste("Figures/Titis_ActivityOnset_gravity_nopoint.tif", sep=""),width = 1050, height = 1800, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist2 , ncol = 1, nrow= 2, labels = c("",""))
print(plot)
dev.off()


tiff(paste("Figures/Titis_Group_MeanGravity.tif", sep=""), width = 850, height = 1800, res = 300,
     compression = "lzw") # file & size
plot<-ggarrange( plotlist = plotlist3 , ncol = 1, nrow= 2, labels = c("",""))
print(plot)
dev.off()


####### Plots 30+15: Activity Onset (Titis) #######
resultsAll_30_15 <- read.csv("results/Titis/GLMMCosinor_ActivityOnset_30+15.csv",
                              header = TRUE)

plotlist_30_15 <- list()

for (i in c(1:length(unique(means$Year)))) {
  cohort       <- unique(means$Year)[i]
  means_cohort <- as.data.frame(means[means$Year == cohort, ])

  results_plot <- extractResults_30_15(resultsAll_30_15, as.character(cohort), means_cohort)
  if (is.null(results_plot)) next

  plot <- lunarPlot(means_cohort, "Activity onset (hours from midnight)",
                    results_plot, n = 2, period = 30, colors[i]) +
    ggtitle(plotTitles[i]) + theme(plot.title = element_text(hjust = 0.5))
  plotlist_30_15[[i]] <- print(plot)
}

tiff("Figures/Titis_ActivityOnset_30+15.tif",
     width = 850, height = 1800, res = 300, compression = "lzw")
plot <- ggarrange(plotlist = plotlist_30_15, ncol = 1, nrow = 2, labels = c("", ""))
print(plot)
dev.off()


###############################################################################
######  CIRCULAR PHASE PLOTS (Option A: one plot per component)  ##############
###############################################################################
#
# For each sleep variable (Start_Time) and each population (Seattle, Tobas),
# the individual phases are plotted as points on a circle where:
#   - 0 degrees / "Full Moon" = DayinCycle 0 (full moon)
#   - 180 degrees / "New Moon" = DayinCycle 15 (new moon)
#
# TWO plots are generated per population:
#   Plot 1 (30-day component): individual phases from the 30-day model
#   Plot 2 (15-day component): individual phases from the 30+15-day model
#
# Points are weighted by amplitude (point size proportional to amplitude).
# The mean vector (Rayleigh) and the p-value of the amplitude-weighted Rayleigh
# test (implemented via bootstrapping) are included.
#
# NOTE: the phase of the 15-day component is expressed on the 15-day cycle
# (0 = maximum, 7.5 = minimum). This is more informative than expressing it on
# the 30-day cycle, where two equidistant peaks may look like one.

library(circular)

# --- Helper function: weighted Rayleigh test via permutations ---
# Function: unweightedRayleigh
# Standard Rayleigh test (unweighted). Permutational p-value for consistency.
unweightedRayleigh <- function(phases_rad, n_perm = 5000) {
  n     <- length(phases_rad)
  C     <- mean(cos(phases_rad))
  S     <- mean(sin(phases_rad))
  R     <- sqrt(C^2 + S^2)
  theta <- atan2(S, C)
  R_null <- replicate(n_perm, {
    phi <- runif(n, 0, 2 * pi)
    sqrt(mean(cos(phi))^2 + mean(sin(phi))^2)
  })
  list(R = R, theta = theta, p_value = mean(R_null >= R),
       theta_days_30 = (theta / (2 * pi)) * 30)
}


# Computes the amplitude-weighted mean resultant length (R) and its permutational
# p-value (under the hypothesis of circular uniformity).
weightedRayleigh <- function(phases_rad, weights, n_perm = 5000) {
  # Weighted mean resultant length
  C <- sum(weights * cos(phases_rad)) / sum(weights)
  S <- sum(weights * sin(phases_rad)) / sum(weights)
  R_obs <- sqrt(C^2 + S^2)
  theta_obs <- atan2(S, C)

  # Null distribution: random uniform phases with the same weights
  R_null <- replicate(n_perm, {
    phi <- runif(length(phases_rad), 0, 2 * pi)
    Cn  <- sum(weights * cos(phi)) / sum(weights)
    Sn  <- sum(weights * sin(phi)) / sum(weights)
    sqrt(Cn^2 + Sn^2)
  })

  p_val <- mean(R_null >= R_obs)

  list(R       = R_obs,
       theta   = theta_obs,       # angle of the mean vector in radians
       p_value = p_val,
       theta_days_30 = (theta_obs / (2 * pi)) * 30)  # converted to days of the cycle
}


# --- Function to build the circular plot ---
# phases_days : phases in days of the cycle (0 to period)
# amplitudes  : weights (individual amplitudes; rescaled to [0.5, 3] for point size)
# period      : cycle length in days (30 or 15)
# rayleigh    : object returned by weightedRayleigh()
# title       : plot title
# color_pts   : color of the individual points
# Function: circularPhasePlot
# weighted = TRUE  : point size proportional to amplitude (default)
# weighted = FALSE : all points the same size (standard Rayleigh)
circularPhasePlot <- function(phases_days, amplitudes, period, rayleigh,
                               title = "", color_pts = "#D55E00",
                               colors_vec = NULL,   # per-point colors (overrides color_pts)
                               weighted = TRUE, method_lbl = "weighted Rayleigh",
                               pole_top = "Full Moon", pole_bottom = "New Moon") {
  # Rotation: Full Moon (day=0) at the TOP, clockwise toward New Moon.
  phases_rad <- pi/2 - (phases_days / period) * 2 * pi

  # Point size: proportional to amplitude (weighted) or uniform (unweighted)
  if (weighted) {
    amp_norm <- 0.5 + 2.5 * (amplitudes - min(amplitudes)) /
      (max(amplitudes) - min(amplitudes) + 1e-9)
  } else {
    amp_norm <- rep(1.5, length(amplitudes))
  }

  use_per_point <- !is.null(colors_vec) && length(colors_vec) == length(phases_days)

  df_pts <- data.frame(
    x    = cos(phases_rad),
    y    = sin(phases_rad),
    size = amp_norm,
    col  = if (use_per_point) colors_vec else rep(color_pts, length(phases_days))
  )

  # Mean vector in the rotated system (Rayleigh's theta is in the original system)
  mean_x <- rayleigh$R * sin(rayleigh$theta)
  mean_y <- rayleigh$R * cos(rayleigh$theta)
  df_arrow <- data.frame(x = 0, y = 0, xend = mean_x, yend = mean_y)

  # Labels for the circle poles (top / bottom). Full/new moon by default;
  # for the 15-day component "Syzygy"/"Quadrature" are passed.
  label_df <- data.frame(
    x     = c(0, 0),
    y     = c(1.15, -1.15),
    label = c(pole_top, pole_bottom)
  )

  p_label <- ifelse(rayleigh$p_value < 0.001, "p < 0.001",
                    sprintf("p = %.3f", rayleigh$p_value))
  subtitle_txt <- sprintf("R = %.3f,  %s", rayleigh$R, p_label)

  ggplot() +
    # Reference circle
    annotate("path",
             x = cos(seq(0, 2 * pi, length.out = 200)),
             y = sin(seq(0, 2 * pi, length.out = 200)),
             color = "gray60", linewidth = 0.5) +
    # Reference lines (full/new moon)
    geom_segment(aes(x = 0, y = -1.05, xend = 0, yend = 1.05),
                 color = "gray80", linetype = "dashed", linewidth = 0.4) +
    # Individual points (size = amplitude, color per-point or uniform)
    geom_point(data = df_pts, aes(x = x, y = y, size = size, color = col), alpha = 0.6) +
    scale_color_identity() +
    scale_size_identity() +
    # Mean vector (Rayleigh)
    geom_segment(data = df_arrow,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 color = "black", linewidth = 1.2) +
    # Full/new moon labels
    geom_text(data = label_df, aes(x = x, y = y, label = label),
              size = 3.5, fontface = "bold") +
    coord_fixed(xlim = c(-1.3, 1.3), ylim = c(-1.4, 1.4)) +
    labs(title = title, subtitle = subtitle_txt) +
    theme_void() +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
      plot.margin   = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}


# Gravity acrophase (peak day) of the 30- or 15-day COMPONENT, estimated from the
# COMBINED 30+15 cosinor fit so each component's phase is isolated from the other
# (the same reason the sleep analysis uses the 30+15 model: a single-component fit
# to bimodal gravity gives a biased peak). 'period' selects which component (30 or
# 15) peak to return. Same convention as the sleep acrophase (peak = phi*T/2pi).
gravPeakDay <- function(t, g, period) {
  ok <- !is.na(t) & !is.na(g)
  if (sum(ok) < 6) return(NA_real_)              # need >= 6 points for the 5-parameter fit
  t <- t[ok]; g <- g[ok]
  df <- data.frame(g = g,
                   c30 = cos(2*pi*t/30), s30 = sin(2*pi*t/30),
                   c15 = cos(2*pi*t/15), s15 = sin(2*pi*t/15))
  b <- coef(lm(g ~ c30 + s30 + c15 + s15, data = df))
  if (period == 30) { phi <- atan2(b["s30"], b["c30"]); Tp <- 30 }
  else              { phi <- atan2(b["s15"], b["c15"]); Tp <- 15 }
  as.numeric((phi * Tp / (2 * pi)) %% Tp)
}

# DIAGNOSTIC: per-cohort gravity profile (points) with the fitted cosinors and the
# model-defined 30-day / 15-day maxima vs the empirical maximum. Lets you check
# whether the acrophases that define the gravity markers/references match the data.
# Remove if not needed.
gravFitPanel <- function(grav_sub, sleep_night, cohorts, cohort_titles, pop_label) {
  pl <- list()
  for (ci in seq_along(cohorts)) {
    co  <- cohorts[ci]
    nts <- unique(sleep_night$night[sleep_night$year == co])
    gg  <- grav_sub[grav_sub$NightDate %in% nts, ]
    if (nrow(gg) < 6) {
      pl[[ci]] <- ggplot() + labs(title = paste0(cohort_titles[ci], " — no data")) + theme_void(); next
    }
    d <- aggregate(mean_Total_C ~ DayinCycle, data = gg, FUN = mean)
    names(d) <- c("t", "g")
    d$c30 <- cos(2*pi*d$t/30); d$s30 <- sin(2*pi*d$t/30)
    d$c15 <- cos(2*pi*d$t/15); d$s15 <- sin(2*pi*d$t/15)
    m30   <- lm(g ~ c30 + s30, data = d)
    m15   <- lm(g ~ c15 + s15, data = d)
    m3015 <- lm(g ~ c30 + s30 + c15 + s15, data = d)

    xx <- seq(0, 30, by = 0.2)
    nd <- data.frame(c30 = cos(2*pi*xx/30), s30 = sin(2*pi*xx/30),
                     c15 = cos(2*pi*xx/15), s15 = sin(2*pi*xx/15))
    curves <- rbind(
      data.frame(t = xx, g = predict(m30,   nd), model = "30-day"),
      data.frame(t = xx, g = predict(m15,   nd), model = "15-day"),
      data.frame(t = xx, g = predict(m3015, nd), model = "30+15")
    )
    peak30   <- gravPeakDay(d$t, d$g, 30)
    peak15   <- gravPeakDay(d$t, d$g, 15)
    data_max <- d$t[which.max(d$g)]

    pl[[ci]] <- ggplot() +
      geom_vline(xintercept = data_max, color = "gray55", linetype = "dotted", linewidth = 0.7) +
      geom_vline(xintercept = peak30, color = "#D55E00", linewidth = 0.9) +
      geom_vline(xintercept = c(peak15, peak15 + 15), color = "#1B9E77",
                 linetype = "dashed", linewidth = 0.9) +
      geom_point(data = d, aes(t, g), color = "black", size = 1.4) +
      geom_line(data = curves, aes(t, g, color = model), linewidth = 0.9) +
      scale_color_manual(values = c("30-day" = "#D55E00", "15-day" = "#1B9E77", "30+15" = "black")) +
      scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
      labs(title = paste0(pop_label, " ", cohort_titles[ci]),
           subtitle = "vlines: 30d & 15d COMPONENT peaks of the 30+15 fit (orange / green dashed); data max (gray dotted)",
           x = "Day in cycle (0 = Full, 15 = New)", y = "Mean gravity (\U00B5Gal)", color = NULL) +
      theme(axis.line = element_line(colour = "black", linewidth = 0.7),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "white", colour = "lightgray"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 6, color = "gray40"))
  }
  nC <- length(pl)
  tiff(paste0("Figures/GravFitDiagnostic_", pop_label, ".tif"),
       width = 1500, height = 700 * nC, res = 300, compression = "lzw")
  print(ggarrange(plotlist = pl, ncol = 1, nrow = nC, common.legend = TRUE, legend = "bottom"))
  dev.off()
  cat(paste0("  Saved: Figures/GravFitDiagnostic_", pop_label, ".tif\n"))
}


# NOTE: the POOLED circular plots referenced to the lunar phase (all cohorts
# combined vs the moon) were removed. Pooling cohorts against the moon smears the
# heterogeneous per-cohort phases and is uninformative. The by-cohort plots
# (below) and the gravity-referenced pooled plots (further down) are kept instead.


###############################################################################
######  CIRCULAR PLOTS BY COHORT: weighted + unweighted, 30 and 15 days     ##
###############################################################################
#
# One circle per cohort (rows = cohorts, same order and colors as the other
# figures), in two versions: amplitude-weighted Rayleigh and standard unweighted
# Rayleigh. A green marker shows the cohort's gravity-peak day on each circle.
#
# 15-day component (semilunar): the two syzygies (full moon AND new moon) fall at
# the SAME point of the 15-day cycle, so the top pole is labeled "Syzygy" and the
# bottom "Quadrature" (quarters), instead of Full/New Moon.
#
# All circular figures are written to Figures/circular/.

cat("Building circular plots by cohort (30 and 15 days)...\n")
set.seed(42)

# --- BY-COHORT driver: one row per cohort, weighted + unweighted columns ---
makeByCohortCirc <- function(data_indiv, sleep_map, cohorts, cohort_titles,
                             cohort_colors, period, pop_label, file_tag, title_tag,
                             min_n = 3, pole_top = "Full Moon",
                             pole_bottom = "New Moon",
                             amp_pat = "^amp1", acr_pat = "^acr1",
                             disp_label = pop_label) {
  # pop_label  -> used for the FILE NAME (must be path-safe, no "/")
  # disp_label -> used for the plot TITLE (can be "Toba/Qom")
  # amp_pat / acr_pat default to "^amp1" / "^acr1" to extract only the 30-day
  # component (component 1) from the 30+15-day combined cosinor model.
  data_indiv$id     <- as.character(data_indiv$id)
  sleep_map$ID_year <- as.character(sleep_map$ID_year)
  sleep_map$year    <- as.character(sleep_map$year)
  d <- merge(data_indiv, sleep_map, by.x = "id", by.y = "ID_year", all.x = TRUE)

  L_w <- list(); L_uw <- list()
  for (ci in seq_along(cohorts)) {
    cohort  <- cohorts[ci]
    dc      <- d[d$year == cohort, ]
    amp     <- dc$estimate[grepl(amp_pat, dc$rowname)]
    acr     <- dc$estimate[grepl(acr_pat, dc$rowname)]
    ph      <- (acr * period / (2 * pi)) %% period
    valid   <- !is.na(amp) & !is.na(ph)
    title_c <- paste0(disp_label, " - ", cohort_titles[ci])

    if (sum(valid) >= min_n) {
      ph_rad <- (ph[valid] / period) * 2 * pi
      p_w  <- circularPhasePlot(ph[valid], amp[valid], period,
                                weightedRayleigh(ph_rad, amp[valid]),
                                title = title_c, color_pts = cohort_colors[ci],
                                weighted = TRUE, method_lbl = "weighted Rayleigh",
                                pole_top = pole_top, pole_bottom = pole_bottom)
      p_uw <- circularPhasePlot(ph[valid], amp[valid], period,
                                unweightedRayleigh(ph_rad),
                                title = title_c,
                                color_pts = cohort_colors[ci], weighted = FALSE,
                                method_lbl = "unweighted Rayleigh",
                                pole_top = pole_top, pole_bottom = pole_bottom)
    } else {
      empty <- ggplot() + labs(title = paste0(title_c, " - n < ", min_n)) + theme_void()
      p_w <- empty; p_uw <- empty
    }
    L_w[[ci]] <- p_w; L_uw[[ci]] <- p_uw
  }

  nC   <- length(cohorts)
  base <- paste0("Figures/circular/", pop_label, "_CircularPhase_ByCohort_", file_tag, "_StartTime")
  tiff(paste0(base, "_weighted.tif"),   width = 1200, height = 1200 * nC, res = 300, compression = "lzw")
  print(ggarrange(plotlist = L_w,  ncol = 1, nrow = nC)); dev.off()
  tiff(paste0(base, "_unweighted.tif"), width = 1200, height = 1200 * nC, res = 300, compression = "lzw")
  print(ggarrange(plotlist = L_uw, ncol = 1, nrow = nC)); dev.off()
  invisible(NULL)
}

# --- Load individual results (30 and 15 days) ---
data_30_seattle_v2   <- read.csv("results/seattle/individuals/modelResults_individual_30+15days_StartTime.csv", stringsAsFactors = FALSE)
data_15_seattle_v2   <- read.csv("results/seattle/individuals/modelResults_individual_15days_StartTime.csv",    stringsAsFactors = FALSE)
data_30_tobas_v2     <- read.csv("results/Tobas/individuals/modelResults_individual_30+15days_StartTime.csv",   stringsAsFactors = FALSE)
data_15_tobas_v2     <- read.csv("results/Tobas/individuals/modelResults_individual_15days_StartTime.csv",      stringsAsFactors = FALSE)

# (POOLED-vs-moon plots removed; see note above. By-cohort plots kept below.)

# --- BY COHORT: 30 and 15 days (cohorts in the same order as the rest) ---
sleep_s_map_v2 <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", header = TRUE, stringsAsFactors = FALSE)
sleep_s_map_v2 <- sleep_s_map_v2[!duplicated(sleep_s_map_v2$ID_year), c("ID_year", "year")]
sleep_t_map_v2 <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", header = TRUE, stringsAsFactors = FALSE)
sleep_t_map_v2 <- sleep_t_map_v2[!duplicated(sleep_t_map_v2$ID_year), c("ID_year", "year")]

seattle_cohorts_v2       <- c("SS2020", "AW2020", "SS2024")
seattle_cohort_titles_v2 <- c("Spring 2020", "Winter 2021", "Spring 2024")
seattle_cohort_colors_v2 <- c("#33DDAA", "#33DDAA", "#FF7799")

tobas_cohorts_v2       <- c("Fall-Winter2016", "Fall-Winter2023", "Spring2018",
                            "Spring2023",      "Summer2024",      "Fall-Winter2024")
tobas_cohort_titles_v2 <- c("Winter 2016", "Fall 2023", "Spring 2018",
                            "Spring 2023", "Summer 2024", "Fall 2024")
tobas_cohort_colors_v2 <- c("#33DDAA", "#55BBFF", "#55BBFF", "#BB66FF", "#FF7799", "#FF7799")

dir.create("Figures/circular", showWarnings = FALSE, recursive = TRUE)

# Per-cohort gravity peak day (for the green gravity-max marker on each circle).
# Loaded here because the time-series section that defines grav_fresh runs later.
grav_circ <- read.csv("data/Data_S4_meanDaily_Gravitational_variation.csv", stringsAsFactors = FALSE)
grav_circ$NightDate <- as.Date(grav_circ$NightDate, format = "%m/%d/%Y")
seattle_grav_circ <- grav_circ[grav_circ$city == "Seattle", ]
tobas_grav_circ   <- grav_circ[grav_circ$city == "Ingeniero Juarez", ]
sleep_s_night <- read.csv("data/Data_S1_all_sleepData_Seattle.csv", stringsAsFactors = FALSE)
sleep_s_night$night <- as.Date(sleep_s_night$night)
sleep_t_night <- read.csv("data/Data_S2_all_sleepData_Tobas.csv", stringsAsFactors = FALSE)
sleep_t_night$night <- as.Date(sleep_t_night$night, format = "%m/%d/%Y")

# Diagnostic: gravity fits + model-defined maxima per cohort (remove if not needed)
gravFitPanel(seattle_grav_circ, sleep_s_night, seattle_cohorts_v2, seattle_cohort_titles_v2, "Seattle")
gravFitPanel(tobas_grav_circ,   sleep_t_night, tobas_cohorts_v2,   tobas_cohort_titles_v2,   "Tobas")

makeByCohortCirc(data_30_seattle_v2, sleep_s_map_v2, seattle_cohorts_v2,
                 seattle_cohort_titles_v2, seattle_cohort_colors_v2,
                 30, "Seattle", "30day", "30-day")
# Tobas is split into two groups: grpA = Winter 2016 + Summer 2024 + Fall 2024
# (cohorts 1,5,6); grpB = Fall 2023 + Spring 2018 + Spring 2023 (cohorts 2,3,4).
# Order within each group follows the cohort order used in the other figures.
tobas_gA <- c(1, 5, 6)
tobas_gB <- c(2, 3, 4)
makeByCohortCirc(data_30_tobas_v2, sleep_t_map_v2, tobas_cohorts_v2[tobas_gA],
                 tobas_cohort_titles_v2[tobas_gA], tobas_cohort_colors_v2[tobas_gA],
                 30, "Tobas", "30day_grpA", "30-day", disp_label = "Toba/Qom")
makeByCohortCirc(data_30_tobas_v2, sleep_t_map_v2, tobas_cohorts_v2[tobas_gB],
                 tobas_cohort_titles_v2[tobas_gB], tobas_cohort_colors_v2[tobas_gB],
                 30, "Tobas", "30day_grpB", "30-day", disp_label = "Toba/Qom")
makeByCohortCirc(data_15_seattle_v2, sleep_s_map_v2, seattle_cohorts_v2,
                 seattle_cohort_titles_v2, seattle_cohort_colors_v2,
                 15, "Seattle", "15day", "15-day",
                 pole_top = "Syzygy", pole_bottom = "Quadrature")
makeByCohortCirc(data_15_tobas_v2, sleep_t_map_v2, tobas_cohorts_v2[tobas_gA],
                 tobas_cohort_titles_v2[tobas_gA], tobas_cohort_colors_v2[tobas_gA],
                 15, "Tobas", "15day_grpA", "15-day",
                 pole_top = "Syzygy", pole_bottom = "Quadrature", disp_label = "Toba/Qom")
makeByCohortCirc(data_15_tobas_v2, sleep_t_map_v2, tobas_cohorts_v2[tobas_gB],
                 tobas_cohort_titles_v2[tobas_gB], tobas_cohort_colors_v2[tobas_gB],
                 15, "Tobas", "15day_grpB", "15-day",
                 pole_top = "Syzygy", pole_bottom = "Quadrature", disp_label = "Toba/Qom")

cat("Circular plots by cohort complete.\n")


###############################################################################
######  SUPPLEMENTARY FIGURE: GRAVITY-REFERENCED vs LUNAR CIRCULAR PLOTS    ##
###############################################################################
#
# For each individual, the 30-day sleep phase is expressed either relative to
# the full moon (lunar reference) or relative to the cohort's gravity peak
# (gravity reference). The gravity peak is the 30-day component peak from the
# combined 30+15 cosinor fit on mean daily gravity for that cohort's dates.
#
# The supplementary figure contains:
#   - By-cohort panels: 2 columns (lunar | gravity reference) × N cohort rows
#   - Pooled panels: lunar reference (expected non-significant) vs
#                    gravity reference (expected significant)
#
# Only the 30-day component is used (15-day is too variable across cohorts).

# Helper: compute gravity peak day for one cohort from the grav_circ data
gravPeakDayCohort <- function(grav_sub, sleep_night_df, cohort, period = 30) {
  nts <- unique(sleep_night_df$night[sleep_night_df$year == cohort])
  gg  <- grav_sub[grav_sub$NightDate %in% nts, ]
  if (nrow(gg) < 6) return(NA_real_)
  d   <- aggregate(mean_Total_C ~ DayinCycle, data = gg, FUN = mean)
  gravPeakDay(d$DayinCycle, d$mean_Total_C, period)
}

# Helper: extract individual phases and amplitudes for one cohort
extractIndivPhases <- function(data_indiv, sleep_map, cohort, period = 30,
                               amp_pat = "^amp1", acr_pat = "^acr1") {
  # amp_pat / acr_pat target the 30-day component (component 1) of the 30+15-day
  # cosinor model. Patterns match amp1, amp11, amp12, ... (one row per participant).
  data_indiv$id     <- as.character(data_indiv$id)
  sleep_map$ID_year <- as.character(sleep_map$ID_year)
  sleep_map$year    <- as.character(sleep_map$year)
  d     <- merge(data_indiv, sleep_map, by.x = "id", by.y = "ID_year", all.x = TRUE)
  dc    <- d[d$year == cohort, ]
  amp   <- dc$estimate[grepl(amp_pat, dc$rowname)]
  acr   <- dc$estimate[grepl(acr_pat, dc$rowname)]
  ph    <- (acr * period / (2 * pi)) %% period
  valid <- !is.na(amp) & !is.na(ph)   # include all participants with valid estimates
  cat(sprintf("    %s: n = %d\n", cohort, sum(valid)))
  list(phase = ph[valid], amplitude = amp[valid])
}

# Compute gravity peak days for all cohorts (Seattle and Tobas)
seattle_grav_peaks_30 <- setNames(
  sapply(seattle_cohorts_v2, gravPeakDayCohort,
         grav_sub = seattle_grav_circ, sleep_night_df = sleep_s_night, period = 30),
  seattle_cohorts_v2)

tobas_grav_peaks_30 <- setNames(
  sapply(tobas_cohorts_v2, gravPeakDayCohort,
         grav_sub = tobas_grav_circ, sleep_night_df = sleep_t_night, period = 30),
  tobas_cohorts_v2)

cat("Gravity peak days (30-day component):\n")
print(seattle_grav_peaks_30)
print(tobas_grav_peaks_30)

# Helper: builds panels (lunar + gravity) for one population and accumulates
# pooled data. Returns list with L_moon, L_grav, ph_moon, ph_grav, amp vectors.
buildCircPanels <- function(data_30, sleep_map, cohorts, titles, colors,
                             grav_peaks, period = 30, pop_label = "") {
  L_moon <- list(); L_grav <- list()
  all_ph_moon <- c(); all_ph_grav <- c(); all_amp <- c(); all_cols <- c()

  for (ci in seq_along(cohorts)) {
    co       <- cohorts[ci]
    res      <- extractIndivPhases(data_30, sleep_map, co, period)
    ph       <- res$phase;  amp <- res$amplitude
    if (length(ph) < 3) next

    gpeak <- grav_peaks[co]
    if (is.na(gpeak)) next

    ph_grav <- (ph - gpeak) %% period

    all_ph_moon <- c(all_ph_moon, ph)
    all_ph_grav <- c(all_ph_grav, ph_grav)
    all_amp     <- c(all_amp, amp)
    all_cols    <- c(all_cols, rep(colors[ci], length(ph)))

    ph_rad      <- (ph      / period) * 2 * pi
    ph_grav_rad <- (ph_grav / period) * 2 * pi

    panel_title <- if (nchar(pop_label) > 0) paste0(pop_label, " - ", titles[ci]) else titles[ci]

    L_moon[[length(L_moon) + 1]] <- circularPhasePlot(
      ph, amp, period, weightedRayleigh(ph_rad, amp),
      title = panel_title, color_pts = colors[ci],
      weighted = TRUE, method_lbl = "weighted Rayleigh",
      pole_top = "Full Moon", pole_bottom = "New Moon"
    )
    L_grav[[length(L_grav) + 1]] <- circularPhasePlot(
      ph_grav, amp, period, weightedRayleigh(ph_grav_rad, amp),
      title = panel_title, color_pts = colors[ci],
      weighted = TRUE, method_lbl = "weighted Rayleigh",
      pole_top = "Gravity Peak", pole_bottom = "Gravity Trough"
    )
  }
  list(L_moon = L_moon, L_grav = L_grav,
       ph_moon = all_ph_moon, ph_grav = all_ph_grav,
       amp = all_amp, cols = all_cols)
}

# Helper: saves a 1-column (lunar phase only) tiff for one population
saveCircPopFig <- function(L_moon, L_grav, pop_label, outdir, period = 30) {
  nC <- length(L_moon)
  tiff(file.path(outdir,
       paste0("SupplFig_Circular_", pop_label, "_LunarRef_30day.tif")),
       width = 1200, height = 1200 * nC, res = 300, compression = "lzw")
  print(ggarrange(plotlist = L_moon, ncol = 1, nrow = nC))
  dev.off()
  cat(sprintf("  Saved %s (%d cohorts).\n", pop_label, nC))
}

# Helper: builds one pair of pooled circular plots (lunar | gravity) from
# pre-extracted phase/amplitude vectors and returns them as a list.
# cols: per-point color vector (uses cohort colors when provided).
makeCircPoolPair <- function(ph_moon, ph_grav, amp, title, cols = NULL,
                              period = 30) {
  ph_moon_rad <- (ph_moon / period) * 2 * pi
  ph_grav_rad <- (ph_grav / period) * 2 * pi
  list(
    moon = circularPhasePlot(
      ph_moon, amp, period, weightedRayleigh(ph_moon_rad, amp),
      title = title, colors_vec = cols,
      weighted = TRUE, method_lbl = "weighted Rayleigh",
      pole_top = "Full Moon", pole_bottom = "New Moon"),
    grav = circularPhasePlot(
      ph_grav, amp, period, weightedRayleigh(ph_grav_rad, amp),
      title = title, colors_vec = cols,
      weighted = TRUE, method_lbl = "weighted Rayleigh",
      pole_top = "Gravity Peak", pole_bottom = "Gravity Trough")
  )
}

# Main function: builds the supplementary circular figure, Seattle and Tobas separate
makeCircSupplFig <- function(data_30_seattle, data_30_tobas,
                              sleep_s_map, sleep_t_map,
                              seattle_cohorts, seattle_titles, seattle_colors,
                              tobas_cohorts, tobas_titles, tobas_colors,
                              seattle_grav_peaks, tobas_grav_peaks,
                              period = 30, outdir = "Figures/circular") {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # --- By-cohort figures (lunar | gravity), one per population ---
  s <- buildCircPanels(data_30_seattle, sleep_s_map,
                        seattle_cohorts, seattle_titles, seattle_colors,
                        seattle_grav_peaks, period, pop_label = "Seattle")
  saveCircPopFig(s$L_moon, s$L_grav, "Seattle", outdir, period)

  t <- buildCircPanels(data_30_tobas, sleep_t_map,
                        tobas_cohorts, tobas_titles, tobas_colors,
                        tobas_grav_peaks, period, pop_label = "Toba/Qom")
  saveCircPopFig(t$L_moon, t$L_grav, "Tobas", outdir, period)

  # --- Pooled figures ---
  # Tobas subsets:
  #   grp1  = cohort 1 only  (Fall-Winter2016)
  #   grp24 = cohorts 2-4    (Fall-Winter2023, Spring2018, Spring2023)
  #   grp56 = cohorts 5-6    (Summer2024, Fall-Winter2024)
  idx1  <- 1
  idx24 <- 2:4
  idx56 <- 5:6

  t1  <- buildCircPanels(data_30_tobas, sleep_t_map,
                          tobas_cohorts[idx1],  tobas_titles[idx1],
                          tobas_colors[idx1],   tobas_grav_peaks, period, pop_label = "Toba/Qom")
  t24 <- buildCircPanels(data_30_tobas, sleep_t_map,
                          tobas_cohorts[idx24], tobas_titles[idx24],
                          tobas_colors[idx24],  tobas_grav_peaks, period, pop_label = "Toba/Qom")
  t56 <- buildCircPanels(data_30_tobas, sleep_t_map,
                          tobas_cohorts[idx56], tobas_titles[idx56],
                          tobas_colors[idx56],  tobas_grav_peaks, period, pop_label = "Toba/Qom")

  # Each point is colored by its own cohort color (via cols vector from buildCircPanels)
  pair_seattle <- makeCircPoolPair(s$ph_moon, s$ph_grav, s$amp,
                                    "Seattle - All cohorts", cols = s$cols, period)
  pair_t1      <- makeCircPoolPair(t1$ph_moon,  t1$ph_grav,  t1$amp,
                                    "Toba/Qom - Winter 2016", cols = t1$cols, period)
  pair_t24     <- makeCircPoolPair(t24$ph_moon, t24$ph_grav, t24$amp,
                                    "Toba/Qom - Fall 2023 /\nSpring 2018 / Spring 2023",
                                    cols = t24$cols, period)
  pair_t56     <- makeCircPoolPair(t56$ph_moon, t56$ph_grav, t56$amp,
                                    "Toba/Qom - Summer 2024 / Fall 2024",
                                    cols = t56$cols, period)
  pair_tall    <- makeCircPoolPair(t$ph_moon,   t$ph_grav,   t$amp,
                                    "Toba/Qom - All cohorts", cols = t$cols, period)

  # Seattle pooled: 1 row × 2 cols (lunar | gravity)
  tiff(file.path(outdir, "SupplFig_CircularPooled_Seattle_LunarVsGrav_30day.tif"),
       width = 2400, height = 1200, res = 300, compression = "lzw")
  print(ggarrange(plotlist = list(pair_seattle$moon, pair_seattle$grav),
                  ncol = 2, nrow = 1))
  dev.off()

  # Tobas pooled: 4 rows × 2 cols (grp1, grp24, grp56, Tobas all) × (lunar | gravity)
  tobas_pool_list <- list(
    pair_t1$moon,   pair_t1$grav,
    pair_t24$moon,  pair_t24$grav,
    pair_t56$moon,  pair_t56$grav,
    pair_tall$moon, pair_tall$grav
  )
  tiff(file.path(outdir, "SupplFig_CircularPooled_Tobas_LunarVsGrav_30day.tif"),
       width = 2400, height = 1200 * 4, res = 300, compression = "lzw")
  print(ggarrange(plotlist = tobas_pool_list, ncol = 2, nrow = 4))
  dev.off()

  # Summary lunar-only: Seattle (1 row) and Tobas (1 row) saved as separate files
  tiff(file.path(outdir, "SupplFig_CircularLunarRef_Seattle_30day.tif"),
       width = 1200, height = 1200, res = 300, compression = "lzw")
  print(pair_seattle$moon)
  dev.off()

  tiff(file.path(outdir, "SupplFig_CircularLunarRef_Tobas_30day.tif"),
       width = 1200, height = 1200, res = 300, compression = "lzw")
  print(pair_tall$moon)
  dev.off()

  n_total <- length(c(s$ph_moon, t$ph_moon))
  cat(sprintf("Supplementary circular figures complete (%d individuals total).\n", n_total))
}

set.seed(42)
makeCircSupplFig(
  data_30_seattle_v2, data_30_tobas_v2,
  sleep_s_map_v2, sleep_t_map_v2,
  seattle_cohorts_v2, seattle_cohort_titles_v2, seattle_cohort_colors_v2,
  tobas_cohorts_v2,   tobas_cohort_titles_v2,   tobas_cohort_colors_v2,
  seattle_grav_peaks_30, tobas_grav_peaks_30
)


###############################################################################
######  SHARED DATA: gravity + raw sleep series  ############################
###############################################################################
# Loaded here (gravity file + raw sleep, full/new moon dates, cohort labels and
# colors) and reused by the per-person and gravity-referenced plots below.

# Function: smooth7_calendar
# 7-day moving average by CALENDAR WINDOW (+-3 days around each night).
# Unlike SMA() over the ordered vector, this version is robust to gaps: it
# averages only the nights that fall within the real 7-day window, so it does
# not "jump" over missing days by counting them as consecutive.
#   dates  : Date vector (for a single participant, ordered or not)
#   values : numeric vector of the same length
smooth7_calendar <- function(dates, values) {
  dd <- as.numeric(dates)
  vapply(seq_along(dd), function(i) {
    w <- abs(dd - dd[i]) <= 3      # centered 7-day window
    mean(values[w], na.rm = TRUE)
  }, numeric(1))
}

# --- Seattle: gravity + raw sleep ---

# Reload raw data -- seattle_full and tobas_full were modified earlier in the
# script (DayinCycle was removed), so we re-read the gravity file.
grav_fresh <- read.csv("data/Data_S4_meanDaily_Gravitational_variation.csv",
                        stringsAsFactors = FALSE)
grav_fresh$NightDate <- as.Date(grav_fresh$NightDate, format = "%m/%d/%Y")
seattle_grav_fresh   <- grav_fresh[grav_fresh$city == "Seattle", ]
tobas_grav_fresh     <- grav_fresh[grav_fresh$city == "Ingeniero Juarez", ]

sleep_raw_ts <- read.csv("data/Data_S1_all_sleepData_Seattle.csv",
                          stringsAsFactors = FALSE)
sleep_raw_ts$night <- as.Date(sleep_raw_ts$night)

seattle_cohorts_ts <- c("SS2020", "AW2020", "SS2024")
seattle_titles_ts  <- c("Spring 2020", "Winter 2021", "Spring 2024")
seattle_colors_ts  <- c("#117755", "#117755", "#CC2299")

# Full and new moon dates (DayinCycle == 0 -> full moon, == 15 -> new moon)
full_moon_s <- seattle_grav_fresh$NightDate[seattle_grav_fresh$DayinCycle == 0]
new_moon_s  <- seattle_grav_fresh$NightDate[seattle_grav_fresh$DayinCycle == 15]

# --- Tobas: gravity + raw sleep ---

sleep_raw_ts_t <- read.csv("data/Data_S2_all_sleepData_Tobas.csv",
                             stringsAsFactors = FALSE)
sleep_raw_ts_t$night <- as.Date(sleep_raw_ts_t$night, format = "%m/%d/%Y")
sleep_raw_ts_t <- sleep_raw_ts_t[sleep_raw_ts_t$year != "Spring2017", ]

tobas_cohorts_ts <- c("Fall-Winter2016", "Fall-Winter2023", "Spring2018",
                       "Spring2023",      "Summer2024",      "Fall-Winter2024")
tobas_titles_ts  <- c("Winter 2016", "Fall 2023", "Spring 2018",
                       "Spring 2023", "Summer 2024", "Fall 2024")
tobas_colors_ts  <- c("#117755", "#1133AA", "#1133AA",
                       "#6611AA", "#CC2299", "#CC2299")

# Use tobas_grav_fresh: tobas_full already lost DayinCycle (line 1308).
# NightDate is already a Date, no need for as.Date.
full_moon_t <- tobas_grav_fresh$NightDate[tobas_grav_fresh$DayinCycle == 0]
new_moon_t  <- tobas_grav_fresh$NightDate[tobas_grav_fresh$DayinCycle == 15]

# ============================================================
# ONE PLOT PER PERSON: smoothed onset (7d MA) vs gravity
# Folder: Figures/IndivxGravity/
# File name: <Seattle|Tobas>_<cohort>_<ID_year>.tif
# ============================================================

cat("Building individual plots (IndivxGravity)...\n")
dir.create("Figures/IndivxGravity", showWarnings = FALSE, recursive = TRUE)

# Function: indivGravityPlot
# Plot of ONE person: their smoothed onset (sleep_ma) + the cohort gravity
# (scaled to the secondary axis, drawn on top). Returns NULL if there is not
# enough range to scale (flat series or a single point).
indivGravityPlot <- function(person_df, grav_df, plot_title,
                             color_sleep = "#CC2299", color_grav = "black",
                             full_moon_dates = NULL, new_moon_dates = NULL,
                             y_name = "Sleep onset (clock time)") {
  s_min <- min(person_df$sleep_ma, na.rm = TRUE)
  s_max <- max(person_df$sleep_ma, na.rm = TRUE)
  g_min <- min(grav_df$mean_grav, na.rm = TRUE)
  g_max <- max(grav_df$mean_grav, na.rm = TRUE)
  if (!all(is.finite(c(s_min, s_max, g_min, g_max))) ||
      s_max == s_min || g_max == g_min) return(NULL)

  inv_scale_fun <- function(x) {
    ((x - s_min) * (g_max - g_min)) / (s_max - s_min) + g_min
  }
  grav_df$grav_scaled <- s_min + (s_max - s_min) *
    (grav_df$mean_grav - g_min) / (g_max - g_min + 1e-9)

  ggplot() +
    {if (!is.null(full_moon_dates))
      geom_vline(xintercept = as.numeric(full_moon_dates),
                 color = "black", linewidth = 0.4)} +
    {if (!is.null(new_moon_dates))
      geom_vline(xintercept = as.numeric(new_moon_dates),
                 linetype = "dotted", color = "gray50", linewidth = 0.4)} +
    geom_line(data = person_df, aes(x = night, y = sleep_ma),
              color = color_sleep, linewidth = 0.8) +
    geom_point(data = person_df, aes(x = night, y = sleep_ma),
               color = color_sleep, size = 1) +
    # Gravity on top of everything
    geom_line(data = grav_df, aes(x = night, y = grav_scaled),
              color = color_grav, linewidth = 1.1, linetype = "81") +
    scale_y_continuous(
      name     = y_name,
      labels   = hours_to_clock,
      sec.axis = sec_axis(transform = ~ inv_scale_fun(.),
                          name = "Daily-mean gravitational deviation (\U00B5Gal)")
    ) +
    scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
    labs(title = plot_title, x = NULL) +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.8),
      panel.border     = element_blank(),
      panel.background = element_rect(fill = "white", colour = "lightgray"),
      panel.grid.major = element_line(linewidth = 0.3, colour = "lightgray"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 30, hjust = 1, size = 7),
      axis.text.y      = element_text(size = 8),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 9)
    )
}

# Function: makeIndivGravityPlots
# Loops over the cohorts and participants of a population and saves one .tif per person.
makeIndivGravityPlots <- function(pop_label, sleep_raw, cohorts, titles, colors,
                                  grav_sub, full_moon, new_moon) {
  n_done <- 0
  for (ci in seq_along(cohorts)) {
    cohort <- cohorts[ci]
    dfc <- sleep_raw[sleep_raw$year == cohort & !is.na(sleep_raw$Start_Time), ]
    if (nrow(dfc) == 0) next

    # Onset smoothed with a 7-day MA per participant (calendar window)
    sm <- dfc %>%
      arrange(ID_year, night) %>%
      group_by(ID_year) %>%
      mutate(sleep_ma = smooth7_calendar(night, Start_Time)) %>%
      ungroup()

    for (id in unique(sm$ID_year)) {
      person_df <- sm[sm$ID_year == id, ]
      person_df <- person_df[order(person_df$night), ]
      if (nrow(person_df) < 2) next

      # Cohort gravity restricted to this person's date range
      rng <- range(person_df$night)
      grav_df <- grav_sub[grav_sub$NightDate >= rng[1] &
                            grav_sub$NightDate <= rng[2], ] %>%
        group_by(NightDate) %>%
        summarise(mean_grav = mean(mean_Total_C, na.rm = TRUE), .groups = "drop") %>%
        mutate(night = NightDate) %>%
        arrange(night)
      if (nrow(grav_df) < 2) next

      fm <- full_moon[full_moon >= rng[1] & full_moon <= rng[2]]
      nm <- new_moon[new_moon  >= rng[1] & new_moon  <= rng[2]]

      p <- indivGravityPlot(person_df, grav_df,
                            paste0(titles[ci], " — ", id),
                            colors[ci], "black", fm, nm)
      if (is.null(p)) next

      id_clean <- gsub("[^A-Za-z0-9]+", "_", as.character(id))
      tiff(paste0("Figures/IndivxGravity/", pop_label, "_", cohort, "_",
                  id_clean, ".tif"),
           width = 1500, height = 850, res = 300, compression = "lzw")
      print(p)
      dev.off()
      n_done <- n_done + 1
    }
  }
  cat(paste0("  ", pop_label, ": ", n_done, " plots saved in Figures/IndivxGravity/\n"))
}

makeIndivGravityPlots("Seattle", sleep_raw_ts, seattle_cohorts_ts,
                      seattle_titles_ts, seattle_colors_ts,
                      seattle_grav_fresh, full_moon_s, new_moon_s)

makeIndivGravityPlots("Tobas", sleep_raw_ts_t, tobas_cohorts_ts,
                      tobas_titles_ts, tobas_colors_ts,
                      tobas_grav_fresh, full_moon_t, new_moon_t)


# ============================================================
# TITI MONKEYS: individual raw activity onset vs gravity
# One plot per monkey (id2): smoothed activity onset + daily gravity (Davis).
# The titi file already contains the Davis gravity (mean_Total_C) and the date.
# NOTE: titi recordings are short (~21 days max, battery-limited), so these are
# qualitative; no individual or pooled CCF is computed for the monkeys.
# Output: Figures/IndivxGravity/Titi_<id2>.tif
# ============================================================

cat("Building individual titi plots (activity onset vs gravity)...\n")
dir.create("Figures/IndivxGravity", showWarnings = FALSE, recursive = TRUE)

titi <- read.csv("data/Data_S3_all_onsetData_titi.csv", stringsAsFactors = FALSE)
titi$night <- as.Date(titi$night)              # titi file: YYYY-MM-DD

# Davis full/new moon dates from the gravity file
davis_grav  <- grav_fresh[grav_fresh$city == "Davis", ]
full_moon_d <- davis_grav$NightDate[davis_grav$DayinCycle == 0]
new_moon_d  <- davis_grav$NightDate[davis_grav$DayinCycle == 15]

n_titi <- 0
for (id in unique(titi$id2)) {
  pf <- titi[titi$id2 == id & !is.na(titi$time), ]
  pf <- pf[order(pf$night), ]
  if (nrow(pf) < 2) next
  pf$sleep_ma <- smooth7_calendar(pf$night, pf$time)   # smoothed activity onset

  # Daily gravity for this monkey (mean_Total_C is already in the titi file)
  gd <- pf[, c("night", "mean_Total_C")]
  gd <- gd[!duplicated(gd$night), ]
  names(gd)[names(gd) == "mean_Total_C"] <- "mean_grav"
  gd <- gd[order(gd$night), ]
  if (nrow(gd) < 2) next

  rng <- range(pf$night)
  fm <- full_moon_d[full_moon_d >= rng[1] & full_moon_d <= rng[2]]
  nm <- new_moon_d[new_moon_d  >= rng[1] & new_moon_d  <= rng[2]]

  p <- indivGravityPlot(pf, gd, paste0("Titi ", id, " (", pf$Year[1], ")"),
                        color_sleep = "#117733", color_grav = "black",
                        full_moon_dates = fm, new_moon_dates = nm,
                        y_name = "Activity onset (clock time)")
  if (is.null(p)) next
  id_clean <- gsub("[^A-Za-z0-9]+", "_", as.character(id))
  tiff(paste0("Figures/IndivxGravity/Titi_", id_clean, ".tif"),
       width = 1500, height = 850, res = 300, compression = "lzw")
  print(p); dev.off()
  n_titi <- n_titi + 1
}
cat(paste0("  Titi: ", n_titi, " plots saved in Figures/IndivxGravity/\n"))


# ============================================================
# CIRCULAR PHASE RELATIVE TO GRAVITY (pooled across cohorts)
# Each individual's sleep acrophase MINUS its cohort's gravity acrophase, per
# component. By referencing every cohort to its OWN gravity (instead of to the
# absolute lunar phase), heterogeneous cohort phases should collapse if sleep is
# organized by gravity: the gravity-referenced phases then cluster (Rayleigh).
# This is the phase-domain analog of the pooled cross-correlation and is a clean
# single-figure test of "sleep tracks gravity, not absolute lunar phase".
# Weighted and unweighted versions of each circle. Seattle (and Titi) pooled across
# cohorts; Tobas split into 3 groups stacked in one column. All -> Figures/circular/.
# ============================================================

cat("Building gravity-referenced circular plots...\n")
set.seed(42)

# (gravPeakDay is defined earlier, near circularPhasePlot.)

# Per-individual sleep phase referenced to its cohort's gravity peak (component
# acrophase), pooled across the requested cohorts. Returns a data.frame of phases
# (days from gravity max) and amplitudes for all individuals with a valid fit.
relGravPhases <- function(data_indiv, sleep_map, sleep_night, grav_sub, cohorts, period) {
  data_indiv$id     <- as.character(data_indiv$id)
  sleep_map$ID_year <- as.character(sleep_map$ID_year)
  sleep_map$year    <- as.character(sleep_map$year)
  d <- merge(data_indiv, sleep_map, by.x = "id", by.y = "ID_year", all.x = TRUE)

  amp <- d$estimate[grepl("^amp1", d$rowname)]
  acr <- d$estimate[grepl("^acr1", d$rowname)]
  yr  <- as.character(d$year[grepl("^amp1", d$rowname)])
  sleep_peak <- (acr * period / (2 * pi)) %% period

  grav_peak <- vapply(cohorts, function(co) {
    nts <- unique(sleep_night$night[sleep_night$year == co])
    gg  <- grav_sub[grav_sub$NightDate %in% nts, ]
    if (nrow(gg) < 6) return(NA_real_)
    prof <- aggregate(mean_Total_C ~ DayinCycle, data = gg, FUN = mean)
    gravPeakDay(prof$DayinCycle, prof$mean_Total_C, period)
  }, numeric(1))
  names(grav_peak) <- cohorts

  gp   <- grav_peak[yr]
  dday <- (sleep_peak - gp) %% period            # sleep phase relative to gravity max
  keep <- !is.na(amp) & amp > 0 & !is.na(dday)
  data.frame(ph = as.numeric(dday[keep]), amp = amp[keep],
             year = yr[keep], stringsAsFactors = FALSE)
}

# Weighted + unweighted circular plots from precomputed gravity-referenced phases.
makeRelGravCircPlots <- function(ph_days, amp, period, title, color,
                                 year = NULL, highlight = NULL, color_map = NULL) {
  if (length(ph_days) < 3) {
    empty <- ggplot() + labs(title = paste0(title, " - n < 3")) + theme_void()
    return(list(weighted = empty, unweighted = empty))
  }
  # Per-point colors inside a pooled circle: a per-cohort color map takes
  # precedence; otherwise a single cohort can be highlighted; otherwise one color.
  pt_col <- color
  if (!is.null(color_map) && !is.null(year)) {
    mapped <- unname(color_map[year])
    pt_col <- ifelse(is.na(mapped), color, mapped)
  } else if (!is.null(highlight) && !is.null(year)) {
    pt_col <- ifelse(year == highlight$cohort, highlight$color, color)
  }
  ph_rad <- (ph_days / period) * 2 * pi
  p_w  <- circularPhasePlot(ph_days, amp, period, weightedRayleigh(ph_rad, amp),
                            title = title, color_pts = pt_col, weighted = TRUE,
                            method_lbl = "weighted Rayleigh",
                            pole_top = "Gravity max", pole_bottom = "Gravity min")
  p_uw <- circularPhasePlot(ph_days, amp, period, unweightedRayleigh(ph_rad),
                            title = title, color_pts = pt_col,
                            weighted = FALSE, method_lbl = "unweighted Rayleigh",
                            pole_top = "Gravity max", pole_bottom = "Gravity min")
  list(weighted = p_w, unweighted = p_uw)
}

# Single pooled population (Seattle, Titi): one panel, weighted + unweighted files.
saveRelGravPooled <- function(data_indiv, sleep_map, sleep_night, grav_sub, cohorts,
                              period, pop_label, file_tag, title_tag, color, color_map = NULL) {
  ph <- relGravPhases(data_indiv, sleep_map, sleep_night, grav_sub, cohorts, period)
  pl <- makeRelGravCircPlots(ph$ph, ph$amp, period,
                             paste0(pop_label, " - All cohorts"), color,
                             year = ph$year, color_map = color_map)
  base <- paste0("Figures/circular/", pop_label, "_CircularPhase_relGravity_", file_tag, "_StartTime")
  tiff(paste0(base, "_weighted.tif"),   width = 1200, height = 1200, res = 300, compression = "lzw")
  print(pl$weighted);   dev.off()
  tiff(paste0(base, "_unweighted.tif"), width = 1200, height = 1200, res = 300, compression = "lzw")
  print(pl$unweighted); dev.off()
  cat(paste0("  Saved: ", pop_label, "_CircularPhase_relGravity_", file_tag,
             "_StartTime (weighted / unweighted)\n"))
}

# --- Seattle: pooled across cohorts, 30 and 15 days. Points colored by cohort
#     (SS2020 + AW2020 green, SS2024 pink) so their origin stays visible. ---
seattle_cmap <- c("SS2020" = "#33DDAA", "AW2020" = "#33DDAA", "SS2024" = "#FF7799")
saveRelGravPooled(data_30_seattle_v2, sleep_s_map_v2, sleep_raw_ts, seattle_grav_fresh,
                  seattle_cohorts_v2, 30, "Seattle", "30day", "30-day", "#D55E00", color_map = seattle_cmap)
saveRelGravPooled(data_15_seattle_v2, sleep_s_map_v2, sleep_raw_ts, seattle_grav_fresh,
                  seattle_cohorts_v2, 15, "Seattle", "15day", "15-day", "#D55E00", color_map = seattle_cmap)

# --- Titi: pooled across cohorts (years 2023, 2024), 30 and 15 days ---
# Requires the individual titi cosinor fits produced by 1_individual_plotsAndModels.R.
titi_f30 <- "results/Titis/individuals/modelResults_individual_30days_ActivityOnset.csv"
titi_f15 <- "results/Titis/individuals/modelResults_individual_15days_ActivityOnset.csv"
if (file.exists(titi_f30) && file.exists(titi_f15)) {
  data_30_titi <- read.csv(titi_f30, stringsAsFactors = FALSE)
  data_15_titi <- read.csv(titi_f15, stringsAsFactors = FALSE)

  titi_onset <- read.csv("data/Data_S3_all_onsetData_titi.csv", stringsAsFactors = FALSE)
  titi_onset$night <- as.Date(titi_onset$night)
  titi_sleep_map   <- data.frame(ID_year = as.character(titi_onset$id2),
                                 year    = as.character(titi_onset$Year), stringsAsFactors = FALSE)
  titi_sleep_map   <- titi_sleep_map[!duplicated(titi_sleep_map$ID_year), ]
  titi_sleep_night <- data.frame(night = titi_onset$night,
                                 year  = as.character(titi_onset$Year), stringsAsFactors = FALSE)
  titi_grav_fresh  <- grav_circ[grav_circ$city == "Davis", ]   # grav_circ loaded above (NightDate as Date)

  saveRelGravPooled(data_30_titi, titi_sleep_map, titi_sleep_night, titi_grav_fresh,
                    c("2023", "2024"), 30, "Titi", "30day", "30-day", "#117733")
  saveRelGravPooled(data_15_titi, titi_sleep_map, titi_sleep_night, titi_grav_fresh,
                    c("2023", "2024"), 15, "Titi", "15day", "15-day", "#117733")
} else {
  cat("  Titi gravity-referenced circles skipped: run 1_individual_plotsAndModels.R first.\n")
}

# --- Tobas: 3-group split in a single column, based on the gravity 30/15 structure.
#     Order top -> bottom: green (Winter 2016), celeste (other three), pink (Summer + Fall 2024).
tobas_groups <- list(
  list(cohorts = c("Fall-Winter2016"),                             color = "#33DDAA", lab = "Winter 2016"),
  list(cohorts = c("Fall-Winter2023", "Spring2018", "Spring2023"), color = "#55BBFF",
       lab = "Fall 2023 +\nSpring 2018 + Spring 2023",
       highlight = list(cohort = "Spring2023", color = "#BB66FF")),
  list(cohorts = c("Summer2024", "Fall-Winter2024"),               color = "#FF7799", lab = "Summer + Fall 2024")
)
for (per in c(30, 15)) {
  dat <- if (per == 30) data_30_tobas_v2 else data_15_tobas_v2
  Lw <- list(); Luw <- list()
  for (gi in seq_along(tobas_groups)) {
    grp <- tobas_groups[[gi]]
    ph  <- relGravPhases(dat, sleep_t_map_v2, sleep_raw_ts_t, tobas_grav_fresh, grp$cohorts, per)
    pl  <- makeRelGravCircPlots(ph$ph, ph$amp, per,
                                paste0("Toba/Qom - ", grp$lab), grp$color,
                                year = ph$year, highlight = grp$highlight)
    Lw[[gi]] <- pl$weighted; Luw[[gi]] <- pl$unweighted
  }
  base <- paste0("Figures/circular/Tobas_CircularPhase_relGravity_", per, "day_3group_StartTime")
  tiff(paste0(base, "_weighted.tif"),   width = 1200, height = 1200 * 3, res = 300, compression = "lzw")
  print(ggarrange(plotlist = Lw,  ncol = 1, nrow = 3)); dev.off()
  tiff(paste0(base, "_unweighted.tif"), width = 1200, height = 1200 * 3, res = 300, compression = "lzw")
  print(ggarrange(plotlist = Luw, ncol = 1, nrow = 3)); dev.off()
  cat(paste0("  Saved: Tobas_CircularPhase_relGravity_", per, "day_3group_StartTime (weighted / unweighted)\n"))
}

cat("Gravity-referenced circular plots complete.\n")

