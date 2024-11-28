# Modeling Wage Determinants: A Multivariate Analysis of CPS Data

## Overview
This project analyzes wage determinants using multivariate statistical techniques on the Current Population Survey (CPS) Outgoing Rotation Group (ORG) data. It is part of the Advanced Multivariate Statistics (AMS) course project.

## Objectives
- Build statistical models to predict wages.
- Analyze wage disparities by studying:
  - Individual characteristics (e.g., education, age, gender).
  - Occupational and regional factors.
  - Household-level influences.

## Repository Contents
- `wage_analysis.R`: Main R script for data preprocessing, analysis, and visualization.
- `.gitignore`: Ensures sensitive and large files are excluded.
- `README.md`: Project documentation.

## Dependencies
Install required R packages:

1. Install `epiextractr` from a custom repository:
   ```R
   install.packages("epiextractr", repos = c("https://economic.r-universe.dev", "https://cloud.r-project.org"))
   ```
2. Install other required packages:
   ```R
   install.packages(c(
      "dplyr", "ggplot2", "lme4", "robustbase", "car", "boot", 
      "lmtest", "sandwich", "MASS", "DescTools", "MuMIn", 
      "performance", "VIM", "psych", "GGally", "lmerTest", 
      "corrplot", "ggcorrplot", "leaps", "rrcov", "RobStatTM", 
      "robustlmm", "stevemisc", "lattice", "tidyverse"
   ))
   ```
## Data Source

The data for this project comes from the CPS ORG, January 2022. It provides detailed information on earnings and demographic characteristics, collected by the U.S. Census Bureau and Bureau of Labor Statistics.

### Download Instructions

To replicate the results, follow these steps:
1. Install `epiextractr` as described in Dependencies.
2. Download the CPS ORG data:
   ```R
   download_cps("org", "/path/to/your/directory")
   ```
   Replace /path/to/your/directory with your desired location for saving the data.
3. Load the data for analysis:
   ```R
   cps_data <- load_org(
      2022,
      year, female, wage, orgwgt, age, education, race, industry, job_class, metro,
      .extracts_dir = "/path/to/your/directory"
   )
   ```

## Methodology

The project uses a combination of statistical models and robust validation techniques:

- **Models:**
  - **OLS Regression:** Baseline model.
  - **Robust Regression:** Addressed outliers and heteroscedasticity (Huber’s M-estimator, MM-estimator, LTS).
  - **Linear Mixed Effects Models (LME):** Incorporated household-level random effects (~25% of total unexplained variance).

- **Validation:**
  - 70-30 train-test split for performance evaluation.
  - .632 bootstrap for unbiased test Mean Squared Error (MSE) estimates.

## Contact

For questions or suggestions, contact Aleksandr Dudakov via [GitHub Profile](https://github.com/aleksandr-dudakov).
