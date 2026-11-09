# Lunar Gravitational Cycles and Sleep Timing

This repository contains the data processing and analysis code used in the study:

**“Lunar gravity predicts sleep timing”**
(preprint on bioRxiv; journal submission in progress)

The goal of this repository is to ensure transparency and reproducibility of the analyses reported in the manuscript.

---

📁 Repository Structure
├── data
|   ├── all_onsetData_titi.csv
|   ├── all_sleepData_Seattle.csv
|   ├── all_sleepData_Tobas.csv
|   ├── Combined_Gravity_SunMoon_JuarezAndSeattle.csv
|
├── Figures
|   ├── CrossCor
|   ├── individual_Plots
|       ├──crossCorr
|   ├── notGrouped
|
├── results
|   ├── seattle
|       ├── individuals
|       ├── sample
|   ├── Titis
|   ├── Tobas
|       ├── individuals
|       ├── sample
|
├── 1_individual_plotsAndModels.r
├── 2_analysis.r
├── 3_crossCor_grav_sleep.r
├── 4_plots.r
├── UsefulFunctions.r



---

## 📊 Data Availability

Processed and aggregated datasets sufficient to reproduce the main analyses and figures are provided in `data/`.

Access to raw data may be granted upon reasonable request and in accordance with IRB and data use agreements.

---

## 🔁 Reproducibility

Analyses were conducted in **R**.

To reproduce the main results and figures:
1. Clone this repository
2. Install required R packages
3. Run scripts in the following order:

```r
1_individual_plotsAndModels.r
2_analysis.r
3_crossCor_grav_sleep.r
4_plots.r
