# install.packages("epiextractr", repos = c("https://economic.r-universe.dev", "https://cloud.r-project.org"))
library(epiextractr)  # For importing CPS data using epiextractr
library(VIM)          # For performing hot deck imputation
library(psych)        # For descriptive statistics
library(GGally)       # For pairwise plots and correlations
library(car)          # For VIF calculation
library(corrplot)     # For correlation plots
library(ggcorrplot)   # For enhanced correlation plots
library(lmtest)       # For heteroskedasticity tests
library(sandwich)     # For robust standard errors
library(leaps)        # For best subset selection
library(DescTools)    # For Winsorized mean and trimmed mean
library(robustbase)   # For robust regression methods
library(MASS)         # For robust regression (rlm)
library(rrcov)        # For robust covariance estimation
library(RobStatTM)    # For MM-estimator with Peña–Yohai initial estimator
library(lme4)         # For linear mixed-effects models
library(lmerTest)     # For p-values in lmer models
library(MuMIn)        # For model selection in mixed models
library(robustlmm)    # For robust mixed-effects models
library(performance)  # For model diagnostics
library(stevemisc)    # For caterpillar plot of random effects
library(lattice)      # For caterpillar plot of random effects (robust)
library(boot)         # For bootstrap methods
library(parallel)     # For parallel computing in bootstrap
library(tidyverse)    # Includes dplyr, ggplot2, etc., for data manipulation and visualization

#---------------------------------------------------------------
# Data Loading and Preprocessing
#---------------------------------------------------------------
# download_cps("org", "/advanced multivariate statistics/project") uncomment for initial data loading
# Load CPS data (2021-2023) and select relevant variables
cps_data <- load_org(
  2021:2023,
  year, 
  month, 
  hrhhid, 
  hrhhid2,
  pulineno,
  age, 
  female,
  citizen,
  wbhao,
  married,
  metstat,
  division,
  gradeatn,
  cow1,
  emp,
  mind03,
  wageotc_noadj,
  .extracts_dir = "data"
) %>%
  as_factor() %>%
  filter(
    emp == "Employed",             # Include only employed individuals
    !cow1 %in% c("Without pay")     # Exclude individuals without pay
  )

# Data Transformation
cps_data <- cps_data %>%
  mutate(
    hhid = str_c(hrhhid, hrhhid2, sep = ""),
    age = as.integer(ifelse(age == '80+', 80, as.character(age))),
    education_level = gradeatn,
    gradeatn = case_when(
      gradeatn == 'Less than 1st grade' ~ 0,
      gradeatn == '1st-4th grade' ~ 2,
      gradeatn == '5th-6th grade' ~ 5,
      gradeatn == '7th-8th grade' ~ 7,
      gradeatn == '9th grade' ~ 9,
      gradeatn == '10th grade' ~ 10,
      gradeatn == '11th grade' ~ 11,
      gradeatn == '12th grade-no diploma' ~ 12,
      gradeatn == 'HS graduate, GED' ~ 12,
      gradeatn == 'Some college but no degree' ~ 13,
      gradeatn == 'Associate degree-occupational/vocational' ~ 14,
      gradeatn == 'Associate degree-academic program' ~ 14,
      gradeatn == 'Bachelor\'s degree' ~ 16,
      gradeatn == 'Master\'s degree' ~ 18,
      gradeatn == 'Professional school' ~ 19,
      gradeatn == 'Doctorate' ~ 21,
      TRUE ~ NA_integer_
    ),
    report_date = make_date(year, month, 1)
  ) %>%
  rename(
    personid = pulineno,
    sex = female,
    race = wbhao,
    metropolitan = metstat,
    education = gradeatn,
    job_class = cow1,
    industry = mind03,
    wage = wageotc_noadj
  )

# Filter Data for January 2022
cps_jan_2022 <- cps_data %>%
  filter(year == 2022, month == 1)

# Initial Missing Value Assessment
missing_value_rate <- cps_jan_2022 %>%
  summarise_all(~ mean(is.na(.)))
print(missing_value_rate, width = Inf)

# Impute Missing Values Using Other Months of the Same Individuals
missing_records <- cps_jan_2022 %>%
  filter(is.na(wage) | is.na(metropolitan)) %>%
  dplyr::select(hhid, personid, report_date)

# Impute "wage" and "metropolitan" variables
variable_list <- c("wage", "metropolitan")
for (var in variable_list) {
  imputation_data <- cps_data %>%
    inner_join(missing_records, by = c("hhid", "personid"), suffix = c(".full", ".missing")) %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(date_difference = abs(difftime(report_date.full, report_date.missing, units = "days"))) %>%
    group_by(hhid, personid) %>%
    slice_min(date_difference, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(hhid, personid, all_of(var)) %>%
    rename(!!paste0("imputed_", var) := !!sym(var))
  
  cps_jan_2022 <- cps_jan_2022 %>%
    left_join(imputation_data, by = c("hhid", "personid")) %>%
    mutate(
      !!var := ifelse(is.na(.data[[var]]), .data[[paste0("imputed_", var)]], .data[[var]])
    ) %>%
    dplyr::select(-all_of(paste0("imputed_", var)))
}

# Post-Imputation Missing Value Assessment
missing_rate <- cps_jan_2022 %>%
  summarise_all(~ mean(is.na(.)))
print(missing_rate, width = Inf)

# Hot Deck Imputation for Remaining Missing Values
set.seed(123)
cps_jan_2022 <- hotdeck(
  cps_jan_2022,
  variable = c("wage", "metropolitan"),
  domain_var = c("age", "race", "sex", "education")
) %>%
  as_tibble()

# Final Missing Value Assessment
missing_rate <- cps_jan_2022 %>%
  summarise_all(~ mean(is.na(.)))
print(missing_rate, width = Inf)

# Remove Records with Missing Values and Zero Wages (145 records)
cps <- cps_jan_2022 %>%
  filter(!if_any(c(wage, metropolitan), is.na), wage != 0) %>%
  mutate(
    metropolitan = factor(metropolitan, levels = c(1, 2), labels = c('Nonmetropolitan', 'Metropolitan')),
    job_class = fct_drop(job_class),
    industry = fct_drop(industry)
  ) %>%
  dplyr::select(hhid, age, sex, citizen, race, married, metropolitan, 
                division, education, job_class, industry, wage, education_level)

# Final Data Structure
str(cps)

#---------------------------------------------------------------
# Exploratory Data Analysis (EDA)
#---------------------------------------------------------------

## 1. Descriptive Statistics

# Create log-transformed wage variable
cps <- cps %>% mutate(log_wage = log(wage))

### 1.1 Summary Statistics

#### Continuous Variables
continuous_vars <- cps %>% select(age, education, wage, log_wage)

# Summary Statistics for Continuous Variables
summary(continuous_vars)
psych::describe(continuous_vars)

#### Categorical Variables
categorical_vars <- cps %>% select(sex, citizen, race, married, metropolitan, division, job_class, industry)

# Frequency and Percentage Tables for Categorical Variables
for (var in names(categorical_vars)) {
  cat("Frequency Table for", var, ":\n")
  print(table(categorical_vars[[var]]))
  cat("\nPercentage Table for", var, ":\n")
  print(prop.table(table(categorical_vars[[var]])) * 100)
  cat("\n-------------------------------------\n")
}

### 1.2 Visualization

#### Histograms and Density Plots for Wage and Log(Wage)

# Define a reusable theme for ggplot
plot_theme <- theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# Histogram and Density Plot for Wage
ggplot(cps, aes(x = wage)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Wage", x = "Wage", y = "Frequency") +
  plot_theme

ggplot(cps, aes(x = wage)) +
  geom_density(fill = "lightblue") +
  labs(title = "Density Plot of Wage", x = "Wage", y = "Density") +
  plot_theme

# Histogram and Density Plot for Log(Wage)
ggplot(cps, aes(x = log_wage)) +
  geom_histogram(binwidth = 0.2, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Log(Wage)", x = "Log(Wage)", y = "Frequency") +
  plot_theme

ggplot(cps, aes(x = log_wage)) +
  geom_density(fill = "lightblue") +
  labs(title = "Density Plot of Log(Wage)", x = "Log(Wage)", y = "Density") +
  plot_theme

#### Boxplots of Log(Wage) by Categorical Variables

# Define a function for creating boxplots
plot_boxplot <- function(var) {
  ggplot(cps, aes(x = .data[[var]], y = log_wage)) +
    geom_boxplot(fill = "lightblue") +
    labs(title = paste("Boxplot of Log(Wage) by", var), x = var, y = "Log(Wage)") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Apply boxplot function for selected categorical variables
lapply(c("sex", "citizen", "race", "married", "metropolitan", "division", "education_level", "job_class", "industry"), plot_boxplot)

#### Bar Plots for Categorical Variables

# Define a function for creating bar plots
plot_bar <- function(var) {
  ggplot(cps, aes(x = fct_infreq(.data[[var]]))) +
    geom_bar(fill = "lightblue") +
    labs(title = paste("Bar Plot of", var), x = var, y = "Count") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Apply bar plot function for selected categorical variables
lapply(names(categorical_vars), plot_bar)

#### Scatter Plots of Wage Against Continuous Predictors

# Wage vs Age
ggplot(cps, aes(x = age, y = log(wage))) +
  geom_jitter(alpha = 0.2) +
  geom_smooth(method = "loess", color = "red") +
  ylim(1, 5) +
  labs(title = "Scatter Plot of Wage vs Age", x = "Age", y = "Wage") +
  plot_theme

# Wage vs Education (Years)
ggplot(cps, aes(x = education, y = log(wage)))+
  geom_jitter(alpha = 0.2) +
  geom_smooth(method = "loess", color = "red") +
  ylim(1, 5) +
  labs(title = "Scatter Plot of Wage vs Years of Education", x = "Years of Education", y = "Wage") +
  plot_theme

# Pairwise Scatter Plots and Correlations for Continuous Variables
continuous_vars_with_logwage <- cps %>% select(age, education, wage, log_wage)
GGally::ggpairs(continuous_vars_with_logwage)

#---------------------------------------------------------------
# Advanced EDA Techniques
#---------------------------------------------------------------

### 2.1 Correlation Analysis and Multicollinearity

# Prepare Numeric Data for Correlation Analysis
numeric_vars <- cps %>%
  mutate(
    citizen_num = ifelse(citizen == "US citizen", 1, 0),
    married_num = ifelse(married == "Married", 1, 0),
    female = ifelse(sex == "Female", 1, 0),
    metro_num = ifelse(metropolitan == "Metropolitan", 1, 0)
  ) %>%
  select(age, education, log_wage, female, citizen_num, married_num, metro_num)

# Correlation Matrix and Heatmap
cor_matrix <- cor(numeric_vars, use = "complete.obs")
ggcorrplot::ggcorrplot(cor_matrix, lab = TRUE, title = "Correlation Heatmap")

# Variance Inflation Factor (VIF) Calculation
vif_model <- lm(log_wage ~ age + sex + citizen + race + married + metropolitan +
                  division + education + job_class + industry, data = cps)
vif_values <- car::vif(vif_model)
print(vif_values)

### 2.2 Normality Testing of Wage and Log(Wage)

# Q-Q Plot for Wage
ggplot(cps, aes(sample = wage)) +
  stat_qq(shape = 16, size = 2) +
  stat_qq_line(color = "red", size = 1) +
  labs(title = "Q-Q Plot of Wage", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Q-Q Plot for Log(Wage)
ggplot(cps, aes(sample = log_wage)) +
  stat_qq(shape = 16, size = 2) +
  stat_qq_line(color = "red", size = 1) +
  labs(title = "Q-Q Plot of Log(Wage)", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Shapiro-Wilk Normality Test on a Sample of Wage and Log(Wage)
wage_sample <- sample(cps$wage, 5000)
log_wage_sample <- sample(cps$log_wage, 5000)

shapiro_test_wage <- shapiro.test(wage_sample)
shapiro_test_log_wage <- shapiro.test(log_wage_sample)

cat("Shapiro-Wilk Test for Wage: p-value =", shapiro_test_wage$p.value, "\n")
cat("Shapiro-Wilk Test for Log(Wage): p-value =", shapiro_test_log_wage$p.value, "\n")

#---------------------------------------------------------------
# Splitting Data into Training and Test Sets
#---------------------------------------------------------------
# Set train-test split proportion and create indices
train_proportion <- 0.7
train_indices <- sample(seq_len(nrow(cps)), size = floor(train_proportion * nrow(cps)))

# Split Data into Training and Test Sets
train_data <- cps[train_indices, ]
test_data <- cps[-train_indices, ]

# Print Split Sizes
cat("Training Set Size:", nrow(train_data), "\n")
cat("Test Set Size:", nrow(test_data), "\n")

#---------------------------------------------------------------
# Initialize Dataframe to Store MSE for Models
#---------------------------------------------------------------
# Dataframe to store model names and MSE results
mse_results <- data.frame(
  Model = character(),
  Test_MSE = numeric(),
  stringsAsFactors = FALSE
)

#---------------------------------------------------------------
# Function to Evaluate Model and Store Test MSE
#---------------------------------------------------------------
evaluate_model <- function(model_name, actuals, predictions, mse_df) {
  mse <- mean((actuals - predictions)^2)
  mse_df <- rbind(mse_df, data.frame(Model = model_name, Test_MSE = mse, stringsAsFactors = FALSE))
  print(mse_df[nrow(mse_df), ])
  return(mse_df)
}

#---------------------------------------------------------------
# Baseline OLS Model
#---------------------------------------------------------------
# Define the formula excluding non-predictor variables
ols_formula <- log_wage ~ age + I(age^2) + sex + citizen + race + married + metropolitan +
  division + education + job_class + industry

# Fit the model on training data
ols_model <- lm(ols_formula, data = train_data)

# Model summary
summary(ols_model)

# Predictions on test data
test_predictions <- predict(ols_model, newdata = test_data)

# Evaluate and store Test MSE
mse_results <- evaluate_model("Ordinary Least Squares (OLS) Regression", test_data$log_wage, test_predictions, mse_results)

#---------------------------------------------------------------
# Model Diagnostics
#---------------------------------------------------------------
## Enhanced Residual Analysis Plots
# Residuals vs Fitted Values
ggplot(data.frame(Fitted = ols_model$fitted.values, Residuals = ols_model$residuals), aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# Scale-Location Plot
ggplot(data.frame(Fitted = ols_model$fitted.values, Scale = sqrt(abs(ols_model$residuals))), aes(x = Fitted, y = Scale)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Square Root of Absolute Residuals") +
  theme_minimal()

# Q-Q Plot of Residuals
ggplot(data.frame(Residuals = ols_model$residuals), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "Q-Q Plot of Residuals",    
    x = "Theoretical Quantiles",
    y = "Sample Quantiles") +
  theme_minimal()

# Residuals vs Leverage Plot
# Outlier Maps and Studentized Residuals
# Compute studentized residuals for the OLS model
studentized_residuals <- rstudent(ols_model)
leverage_values <- hatvalues(ols_model)
cutoff_leverage <- 2 * mean(leverage_values)

# Plot Studentized Residuals vs Leverage (Outlier Map) using ggplot2
outlier_data <- data.frame(
  Leverage = leverage_values,
  Studentized_Residuals = studentized_residuals
)

ggplot(outlier_data, aes(x = Leverage, y = Studentized_Residuals)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(cutoff_leverage, 1.5 * cutoff_leverage), color = "blue", linetype = "dashed") +
  labs(title = "Outlier Map for OLS Model", x = "Leverage", y = "Studentized Residuals") +
  theme_minimal()

# Histogram of Residuals
ggplot(data.frame(Residuals = ols_model$residuals), aes(x = Residuals)) +
  geom_histogram(bins = 80, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Residuals", x = "Residuals") +
  theme_minimal()

# Shapiro-Wilk Test for Normality of Residuals (using a sample due to large N)
shapiro_test <- shapiro.test(sample(ols_model$residuals, min(5000, length(ols_model$residuals))))
print(shapiro_test)

## Heteroskedasticity Test (Breusch-Pagan Test)
bptest_result <- lmtest::bptest(ols_model)
print(bptest_result)

#---------------------------------------------------------------
# Robust Standard Errors
#---------------------------------------------------------------
# Calculate robust standard errors
robust_se <- sqrt(diag(sandwich::vcovHC(ols_model, type = "HC1")))

# Robust summary with robust standard errors
robust_summary <- lmtest::coeftest(ols_model, vcov = sandwich::vcovHC(ols_model, type = "HC1"))
print(robust_summary)

#---------------------------------------------------------------
# Best Subset Selection (BSS)
#---------------------------------------------------------------
# Perform Best Subset Selection with `regsubsets`
nvmax <- 37  # Set maximum number of predictors to consider
bss_fit <- leaps::regsubsets(
  log_wage ~ age + I(age^2) + sex + citizen + race + married + metropolitan +
    division + education + job_class + industry,
  data = train_data,
  nvmax = nvmax,
  method = "exhaustive"
)

# Custom predict function for regsubsets to handle new data
predict.regsubsets <- function(object, newdata, id) {
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id = id)
  predictions <- mat[, names(coefi)] %*% coefi
  return(predictions)
}

# Calculate and store test MSE for each model size
test_mse_bss <- sapply(1:nvmax, function(i) {
  predictions <- predict.regsubsets(bss_fit, newdata = test_data, id = i)
  mean((test_data$log_wage - predictions)^2)
})

# Identify the optimal model size with the lowest test MSE
best_model_size <- which.min(test_mse_bss)
cat("Optimal Model Size (Test MSE):", best_model_size, "\n")
cat("Lowest Test MSE for BSS Model:", test_mse_bss[best_model_size], "\n")

# Visualization: Test MSE vs. Model Size using ggplot2
# Create a data frame for ggplot
bss_mse_df <- data.frame(
  Model_Size = 1:nvmax,
  Test_MSE = test_mse_bss
)

# Plot Test MSE vs Model Size
ggplot(bss_mse_df, aes(x = Model_Size, y = Test_MSE)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  annotate("point", x = best_model_size, y = test_mse_bss[best_model_size], 
           color = "red", size = 3) +
  labs(
    title = "Best Subset Selection (Test MSE by Model Size)",
    x = "Number of Predictors",
    y = "Test MSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Extract and Evaluate the Best BSS Model
# Retrieve coefficients of the best BSS model
best_bss_coefs <- coef(bss_fit, id = best_model_size)
print(best_bss_coefs)

# Generate predictions on test data using the best BSS model
best_bss_predictions <- predict.regsubsets(bss_fit, newdata = test_data, id = best_model_size)

# Evaluate and store Test MSE for the best BSS model
mse_results <- evaluate_model("Best Subset Selection (BSS) Regression", test_data$log_wage, best_bss_predictions, mse_results)

#---------------------------------------------------------------
# Robust Methods
#---------------------------------------------------------------
## 1. Univariate Robust Estimators on Continuous Variables
# Select continuous variables
continuous_vars <- cps %>% dplyr::select(age, education, wage, log_wage)
continuous_vars_no_wage <- continuous_vars %>% dplyr::select(-wage)

# Function to compute robust statistics
compute_univariate_statistics <- function(x, alpha = 0.1) {
  x <- na.omit(x)
  list(
    Mean = mean(x),
    Median = median(x),
    Trimmed_Mean = mean(x, trim = alpha),
    Winsorized_Mean = mean(DescTools::Winsorize(x, val = quantile(x, probs = c(alpha, 1 - alpha)))),
    Variance = var(x),
    Standard_Deviation = sd(x),
    Mean_Absolute_Deviation = mean(abs(x - mean(x))),
    Median_Absolute_Deviation = mad(x, constant = 1),
    Normalized_MAD = mad(x),
    IQR = IQR(x),
    Huber_Mean = MASS::huber(x)$mu,
    Huber_Scale = MASS::huber(x)$s,
    Tukey_Bisquare_Mean = coef(MASS::rlm(x ~ 1, psi = psi.bisquare))[[1]],
    Tukey_Scale = MASS::rlm(x ~ 1, psi = psi.bisquare)$s
  )
}

# Compute and format robust statistics
univariate_stats <- lapply(continuous_vars, compute_univariate_statistics)
univariate_stats_df <- as.data.frame(do.call(rbind, lapply(univariate_stats, unlist)))
rownames(univariate_stats_df) <- names(continuous_vars)
print(round(univariate_stats_df, 4))

## 2. Multivariate Robust Estimators of Location and Scatter
### 2.1 Classical and Robust Mahalanobis Distance
# Classical Mahalanobis distances
mu_classical <- colMeans(continuous_vars_no_wage)
cov_classical <- cov(continuous_vars_no_wage)
md_classical <- mahalanobis(continuous_vars_no_wage, center = mu_classical, cov = cov_classical)

# Robust Mahalanobis Distance using Minimum Covariance Determinant (MCD)
mcd_result <- rrcov::CovMcd(continuous_vars_no_wage)
mu_rmcd <- mcd_result@center
cov_rmcd <- mcd_result@cov
md_rmcd <- mahalanobis(continuous_vars_no_wage, center = mu_rmcd, cov = cov_rmcd)

robust_stats <- data.frame(
  Variable = names(continuous_vars_no_wage),
  Robust_Mean = mu_rmcd,
  Robust_Variance = diag(cov_rmcd)
)

## Combine univariate and multivariate statistics
robust_stats_df <- univariate_stats_df %>%
  rownames_to_column("Variable") %>%
  left_join(robust_stats, by = "Variable")

# Display the combined statistics
robust_stats_df

## Distance-Distance Plot
distance_data <- data.frame(
  Classical_Distance = md_classical,
  Robust_Distance = md_rmcd
)

ggplot(distance_data, aes(x = Classical_Distance, y = Robust_Distance)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  xlim(0, 50) +
  ylim(0, 50) +
  labs(
    title = "Distance-Distance Plot: Classical vs Robust Mahalanobis Distances",
    x = "Classical Mahalanobis Distance",
    y = "Robust Mahalanobis Distance"
  ) +
  theme_minimal()

### 2.3 Tolerance Ellipsoids
# Pairwise scatter plots with classical and robust tolerance ellipsoids
pairs(continuous_vars, 
      panel = function(x, y) {
        points(x, y, pch = 20, cex = 0.6)
        # Add classical ellipse
        dataEllipse(x, y, levels = 0.95, add = TRUE, plot.points = FALSE, col = "blue", lty = 2)
        # Add robust ellipse (using RMCD estimates)
        dataEllipse(x, y, levels = 0.95, add = TRUE, plot.points = FALSE, robust = TRUE, col = "red")
      },
      main = "Pairwise Scatter Plots with Tolerance Ellipsoids")

## 3. Robust Regression Models
## 3.1 Robust Regression using M-Estimators (Huber's Method)
# Fit robust regression model using Huber's M-estimator
huber_model <- MASS::rlm(ols_formula, data = train_data, psi = psi.huber)
summary(huber_model)

# Predictions on test data
huber_predictions <- predict(huber_model, newdata = test_data)

# Evaluate and store Test MSE for Huber's M-estimator
mse_results <- evaluate_model("Huber's M-estimator Robust Regression", test_data$log_wage, huber_predictions, mse_results)

# Diagnostic Plots for Huber Model
par(mfrow = c(2, 2))
plot(huber_model, main = "Diagnostic Plots for Huber's M-estimator")

## 3.2 Robust Regression using MM-Estimator with Peña–Yohai Initial Estimator
# Control parameters for MM-estimator
control_lmrobdetMM <- lmrobdet.control(
  bb = 0.5,            # Breakdown point for initial S-estimator
  efficiency = 0.95,   # Desired efficiency
  family = "mopt",     # Optimal psi function
  max.it = 100,        # Max iterations
  rel.tol = 1e-07      # Convergence tolerance
)

# Fit MM-estimator with Peña–Yohai initial estimator
lmrobdetMM_model <- RobStatTM::lmrobdetMM(ols_formula, data = train_data, control = control_lmrobdetMM)
summary(lmrobdetMM_model)

# Predictions on test data
lmrobdetMM_predictions <- predict(lmrobdetMM_model, newdata = test_data)

# Evaluate and store Test MSE for MM-estimator
mse_results <- evaluate_model("MM-estimator Robust Regression with Peña–Yohai Initialization", test_data$log_wage, lmrobdetMM_predictions, mse_results)

# Diagnostic Plots for MM-Estimator
par(mfrow = c(2, 2))
plot(lmrobdetMM_model, which = c(1, 2, 4, 5), main = "Diagnostic Plots for MM-Estimator with PY Init")

## 3.3 Least Trimmed Squares (LTS) Regression
# Fit LTS regression model (alpha = 0.5), using nsamp = 5000
lts_model <- robustbase::ltsReg(ols_formula, data = train_data, alpha = 0.5, nsamp = 5000)
summary(lts_model)

# Predictions on test data (LTS doesn't have a predict method)
lts_predictions <- model.matrix(ols_formula, data = test_data) %*% coef(lts_model)

# Evaluate and store Test MSE for LTS model
mse_results <- evaluate_model("Least Trimmed Squares (LTS) Robust Regression", test_data$log_wage, lts_predictions, mse_results)

# Diagnostic Plots for LTS Model
# Plot 1: Response vs Fitted Values
ggplot(data.frame(Fitted = lts_model$fitted.values, Observed = train_data$log_wage), aes(x = Fitted, y = Observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Response vs Fitted Values (LTS Regression)", x = "Fitted Values", y = "Observed Log(Wage)") +
  theme_minimal()

# Plot 2: Residuals vs Fitted Values
ggplot(data.frame(Fitted = lts_model$fitted.values, Residuals = residuals(lts_model)), aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values (LTS Regression)", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# Plot 3: Scale-Location Plot
sqrt_abs_resid <- sqrt(abs(residuals(lts_model)))
mean_sqrt_resid <- mean(sqrt_abs_resid)  
ggplot(data.frame(Fitted = lts_model$fitted.values, Sqrt_Abs_Residuals = sqrt_abs_resid), 
       aes(x = Fitted, y = Sqrt_Abs_Residuals)) +
  geom_point() +
  geom_hline(yintercept = mean_sqrt_resid, color = "red", linetype = "dashed") +
  labs(title = "Scale-Location Plot (LTS Regression)", x = "Fitted Values", y = "Sqrt(|Residuals|)") +
  theme_minimal()

# Plot 4: Residuals vs Robust Distances
# Compute robust Mahalanobis distances on continuous predictors
continuous_vars <- train_data %>% dplyr::select(age, education)
cov_robust <- covMcd(continuous_vars)
center <- cov_robust$center
covariance <- cov_robust$cov

# Check for singularity and compute robust distances if non-singular
if (det(covariance) != 0) {
  robust_distances <- mahalanobis(continuous_vars, center = center, cov = covariance)
  std_resid_lts <- residuals(lts_model) / lts_model$scale
  
  ggplot(data.frame(Robust_Distance = robust_distances, Std_Residuals = std_resid_lts), aes(x = Robust_Distance, y = Std_Residuals)) +
    geom_point() +
    geom_hline(yintercept = c(-2, 0, 2), color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(2 * mean(robust_distances)), color = "red", linetype = "dashed") +
    labs(title = "Residuals vs Robust Distances (LTS Regression)", x = "Robust Distances", y = "Standardized Residuals") +
    theme_minimal()
} else {
  cat("Covariance matrix is singular. Cannot compute robust distances.\n")
}

# 3.5 Outlier Maps for Robust Models
## 3.5.1 Outlier Map for Huber's M-estimator Model
std_resid_huber <- residuals(huber_model) / huber_model$s
huber_hatvalues <- hatvalues(lm(ols_formula, data = train_data)) 
huber_weight <- huber_model$weights

outlier_data_huber <- data.frame(
  Leverage = huber_hatvalues,
  Std_Residuals = std_resid_huber,
  Weight = huber_weight
)

ggplot(outlier_data_huber, aes(x = Leverage, y = Std_Residuals, alpha = Weight)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(2, 3) * mean(huber_hatvalues), color = "blue", linetype = "dashed") +
  labs(title = "Outlier Map for Huber's M-estimator", x = "Leverage", y = "Standardized Residuals") +
  theme_minimal()

## 3.5.2 Outlier Map for MM-Estimator with PY init
std_resid_mm <- residuals(lmrobdetMM_model) / lmrobdetMM_model$scale
mm_distances <- lmrobdetMM_model$MD  # Robust Mahalanobis distances from MM model
mm_weight <- lmrobdetMM_model$rweights

outlier_data_mm <- data.frame(
  Mahalanobis_Distances = mm_distances,
  Std_Residuals = std_resid_mm,
  Weight = mm_weight
)

ggplot(outlier_data_mm, aes(x = Mahalanobis_Distances, y = Std_Residuals, alpha = Weight)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(2, 3) * mean(mm_distances), color = "blue", linetype = "dashed") +
  labs(title = "Outlier Map for MM-Estimator with PY init", x = "Mahalanobis Distances", y = "Standardized Residuals") +
  theme_minimal()

## 3.5.3 Outlier Map for LTS Regression Model
std_resid_lts <- residuals(lts_model) / lts_model$scale
lts_hatvalues <- hatvalues(lm(ols_formula, data = train_data))  # Approximate leverage
lts_weight <- lts_model$raw.weights

outlier_data_lts <- data.frame(
  Leverage = lts_hatvalues,
  Std_Residuals = std_resid_lts,
  Weight = lts_weight
)

ggplot(outlier_data_lts, aes(x = Leverage, y = Std_Residuals, alpha = Weight)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(2, 3) * mean(lts_hatvalues), color = "blue", linetype = "dashed") +
  labs(title = "Outlier Map for LTS Regression", x = "Leverage", y = "Standardized Residuals") +
  theme_minimal()

# 3.6 QQ Plots of Residuals for Robust Models
## 3.6.1 QQ Plot for Huber's M-estimator Residuals
# Extract residuals and weights
residuals_huber <- residuals(huber_model)
weights_huber <- huber_model$weights

# Create QQ plot
ggplot(data.frame(Residuals = residuals_huber, Weight = weights_huber), aes(sample = Residuals, alpha = Weight)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  scale_alpha_continuous(name = "Weight", range = c(0.1, 1), guide = guide_colourbar(title.position = "top")) +
  labs(
    title = "QQ Plot of Residuals (Huber's M-estimator)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

## 3.6.2 QQ Plot for MM-Estimator Residuals
# Extract residuals and weights
residuals_mm <- residuals(lmrobdetMM_model)
weights_mm <- lmrobdetMM_model$rweights

# Compute theoretical quantiles and sort residuals
n <- length(residuals_mm)
theoretical_quantiles <- qnorm(ppoints(n))
sorted_residuals <- sort(residuals_mm)
sorted_weights <- weights_mm[order(residuals_mm)]

# Create data frame
qq_data_mm <- data.frame(
  Theoretical = theoretical_quantiles,
  Residuals = sorted_residuals,
  Weight = sorted_weights
)

# Calculate slope and intercept for the reference line
slope <- sd(qq_data_mm$Residuals) / sd(qq_data_mm$Theoretical)
intercept <- mean(qq_data_mm$Residuals) - slope * mean(qq_data_mm$Theoretical)

# Create QQ plot
ggplot(qq_data_mm, aes(x = Theoretical, y = Residuals, alpha = Weight)) +
  geom_point() +
  geom_abline(slope = slope, intercept = intercept, color = "red", linewidth = 1) +
  scale_alpha_continuous(name = "Weight", range = c(0.1, 1), guide = guide_colourbar(title.position = "top")) +
  labs(
    title = "QQ Plot of Residuals (MM-Estimator with PY init)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

## 3.6.3 QQ Plot for LTS Residuals
# Extract residuals
# Extract residuals and weights
residuals_lts <- residuals(lts_model)
weights_lts <- lts_model$raw.weights  # Use raw weights from the LTS model

# Compute theoretical quantiles and sort residuals
n <- length(residuals_lts)
theoretical_quantiles <- qnorm(ppoints(n))
sorted_residuals <- sort(residuals_lts)
sorted_weights <- weights_lts[order(residuals_lts)]  # Match weights to sorted residuals

# Create data frame
qq_data_lts <- data.frame(
  Theoretical = theoretical_quantiles,
  Residuals = sorted_residuals,
  Weight = sorted_weights
)

# Calculate slope and intercept for the reference line
slope <- sd(qq_data_lts$Residuals) / sd(qq_data_lts$Theoretical)
intercept <- mean(qq_data_lts$Residuals) - slope * mean(qq_data_lts$Theoretical)

# Create QQ plot
ggplot(qq_data_lts, aes(x = Theoretical, y = Residuals, alpha = Weight)) +
  geom_point() +
  geom_abline(slope = slope, intercept = intercept, color = "red", linewidth = 1) +
  scale_alpha_continuous(name = "Weight", range = c(0.1, 1), guide = guide_colourbar(title.position = "top")) +
  labs(
    title = "QQ Plot of Residuals (LTS Regression)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
#---------------------------------------------------------------
# Linear Mixed Effects Models
#---------------------------------------------------------------
# Calculate the distribution of household sizes
cps %>%
  group_by(hhid) %>%
  summarise(household_size = n()) %>%
  count(household_size, name = "number_of_households") %>%
  arrange(household_size) |>
  mutate(
    household_size = household_size, 
    number_of_households = number_of_households,
    number_of_people = household_size * number_of_households,
    proportion_of_people = number_of_people / 12148
  )

# Define Candidate Mixed Effects Models
lme_formulas <- list(
  # Model1: Random intercept for households
  Model1 = log_wage ~ age + I(age^2) + sex + race + citizen + married + 
    metropolitan + division + education + job_class + industry + 
    (1 | hhid),
  
  # Model2: Nested random effects for division and households
  Model2 = log_wage ~ age + I(age^2) + sex + race + citizen + married + 
    metropolitan + education + job_class + industry + 
    (1 | division/hhid)
)

# Fit each model and store results
lme_models <- lapply(names(lme_formulas), function(name) {
  cat("Fitting", name, "\n")
  lmer(lme_formulas[[name]], data = train_data)
})
names(lme_models) <- names(lme_formulas)

# Model Comparison and Selection
# Collect model selection metrics (AIC, BIC, Log-Likelihood, Test MSE)
model_selection_metrics <- data.frame(Model = names(lme_models), stringsAsFactors = FALSE)

for (i in seq_along(lme_models)) {
  model <- lme_models[[i]]
  test_predictions <- predict(model, newdata = test_data, allow.new.levels = TRUE)
  model_selection_metrics[i, "AIC"] <- AIC(model)
  model_selection_metrics[i, "BIC"] <- BIC(model)
  model_selection_metrics[i, "LogLik"] <- as.numeric(logLik(model))
  model_selection_metrics[i, "Test_MSE"] <- mean((test_data$log_wage - test_predictions)^2)
}

# Display Model Selection Metrics
print(model_selection_metrics)

# Select the Best Model Based on AIC
# Although 2nd model is (very-very) slightly better in AIC, 
# the 1st is simplier and better in MSE
best_model_name <- model_selection_metrics$Model[1]
cat("Best Model Based on AIC:", best_model_name, "\n")
best_lme_model <- lme_models[[best_model_name]]

summary(best_lme_model)

# Evaluate Best Model and Update MSE Results
best_lme_predictions <- predict(best_lme_model, newdata = test_data, allow.new.levels = TRUE)
mse_results <- evaluate_model("Linear Mixed Effects (LME) Regression", test_data$log_wage, best_lme_predictions, mse_results)

# Model Diagnostics for the Best Mixed Effects Model
# Extract residuals and fitted values for diagnostics
residuals_best <- residuals(best_lme_model, scaled = TRUE)
fitted_best <- fitted(best_lme_model)

# Plot 1: Standardized Residuals vs Fitted Values
ggplot(data.frame(Fitted = fitted_best, Residuals = residuals_best), aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Standardized Residuals vs Fitted Values for LME model",
    x = "Fitted Values",
    y = "Standardized Residuals"
  ) +
  theme_minimal()

# Plot 2: Q-Q Plot of Standardized Residuals
ggplot(data.frame(Residuals = residuals_best), aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "Q-Q Plot of Standardized Residuals for LME model",
    x = "Theoretical Quantiles",
    y = "Standardized Residuals"
  ) +
  theme_minimal()

# Plot 3.1: Q-Q Plot of Random Effects 
ranef_hh <- ranef(best_lme_model)$hhid[, "(Intercept)"]
ggplot(data.frame(RandomEffects = ranef_hh), aes(sample = RandomEffects)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "Q-Q Plot of Random Effects (Household)",
    x = "Theoretical Quantiles",
    y = "Random Effects (Intercept)"
  ) +
  theme_minimal()

# Plot 4: Caterpillar Plot of Random Effects
# Extract random effects and conditional variances
# Extract random effects and conditional variances
ranef_data <- ranef(best_lme_model, condVar = TRUE)$hhid
post_var <- attr(ranef(best_lme_model, condVar = TRUE)$hhid, "postVar")

# Prepare data frame with random effects and confidence intervals
qq <- qnorm(0.975)
ranef_df <- as.data.frame(ranef_data) %>%
  mutate(
    hhid = rownames(ranef_data),
    variance = sapply(1:dim(post_var)[3], function(i) post_var[1, 1, i]),
    CI_low = `(Intercept)` - qq * sqrt(variance),
    CI_high = `(Intercept)` + qq * sqrt(variance)
  ) %>%
  arrange(desc(abs(`(Intercept)`)))

# Plot top 50 random effects
ranef_df %>% 
  slice_head(n = 50) %>%
  ggplot(aes(x = reorder(hhid, `(Intercept)`), y = `(Intercept)`)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  coord_flip() +
  labs(
    title = "Top 50 Random Effects for hhid",
    x = "Household ID",
    y = "Random Effect (Intercept)"
  ) +
  theme_minimal()

# Percentage of Variation due to Random Effects (PVRE)
# PVRE for the Standard LME Model
var_components <- as.data.frame(VarCorr(best_lme_model))
var_random <- sum(var_components$vcov[var_components$grp != "Residual"])
var_residual <- var_components$vcov[var_components$grp == "Residual"]
pvre <- var_random / (var_random + var_residual)
cat("PVRE for Standard LME Model:", round(pvre * 100, 2), "%\n")

#---------------------------------------------------------------
# Bootstrap Methods
#---------------------------------------------------------------
## 1. Bootstrap Confidence Intervals (BCa) on Continuous Variables

### 1.1 Bootstrapping the Mean of `wage` and `log_wage`

# Define the statistic function for bootstrapping the mean
boot_mean <- function(data, indices) {
  mean(data[indices])  # Resample data and calculate mean
}

# Perform bootstrap for `wage` and `log_wage`
boot_wage <- boot(data = cps$wage, statistic = boot_mean, R = 1000)
boot_logwage <- boot(data = cps$log_wage, statistic = boot_mean, R = 1000)

# Bootstrap results and BCa confidence intervals
wage_ci_bca <- boot.ci(boot_wage, type = "bca", L = empinf(boot_wage, index = 1L, type = "jack"))
logwage_ci_bca <- boot.ci(boot_logwage, type = "bca", L = empinf(boot_logwage, index = 1L, type = "jack"))

# Print results
cat("Bootstrap Mean and 95% BCa CI for `wage`:\n")
print(boot_wage)
print(wage_ci_bca)

cat("\nBootstrap Mean and 95% BCa CI for `log_wage`:\n")
print(boot_logwage)
print(logwage_ci_bca)

# Compare with Normal Theory and t-based Confidence Intervals for Mean

# Normal Theory and t-based CI for `wage`
wage_mean <- mean(cps$wage)
wage_sd <- sd(cps$wage)
n_wage <- length(cps$wage)
alpha <- 0.05
z <- qnorm(1 - alpha / 2)
t_value_wage <- qt(1 - alpha / 2, df = n_wage - 1)

wage_ci_normal <- wage_mean + c(-1, 1) * z * wage_sd / sqrt(n_wage)
wage_ci_t <- wage_mean + c(-1, 1) * t_value_wage * wage_sd / sqrt(n_wage)
cat("Normal Theory 95% CI for mean of `wage`:", wage_ci_normal, "\n")
cat("t-based 95% CI for mean of `wage`:", wage_ci_t, "\n")

# Normal Theory and t-based CI for `log_wage`
logwage_mean <- mean(cps$log_wage)
logwage_sd <- sd(cps$log_wage)
n_logwage <- length(cps$log_wage)

logwage_ci_normal <- logwage_mean + c(-1, 1) * z * logwage_sd / sqrt(n_logwage)
logwage_ci_t <- logwage_mean + c(-1, 1) * t_value_wage * logwage_sd / sqrt(n_logwage)
cat("Normal Theory 95% CI for mean of `log_wage`:", logwage_ci_normal, "\n")
cat("t-based 95% CI for mean of `log_wage`:", logwage_ci_t, "\n")

### 1.2 Bootstrapping the Median of `wage` and `log_wage`

# Define the statistic function for bootstrapping the median
boot_median <- function(data, indices) {
  median(data[indices])
}

# Perform bootstrap for median of `wage` and `log_wage`
boot_median_wage <- boot(data = cps$wage, statistic = boot_median, R = 1000)
boot_median_logwage <- boot(data = cps$log_wage, statistic = boot_median, R = 1000)

# Compute BCa confidence intervals for medians
wage_median_ci_bca <- boot.ci(boot_median_wage, type = "bca", L = empinf(boot_wage, index = 1L, type = "jack"))
logwage_median_ci_bca <- boot.ci(boot_median_logwage, type = "bca", L = empinf(boot_logwage, index = 1L, type = "jack"))

# Print median bootstrap results
cat("\nBootstrap Median and 95% BCa CI for `wage`:\n")
print(boot_median_wage)
print(wage_median_ci_bca)

cat("\nBootstrap Median and 95% BCa CI for `log_wage`:\n")
print(boot_median_logwage)
print(logwage_median_ci_bca)

# Note: No parametric confidence intervals for the median under normal theory.

### 1.3 Visualize Bootstrap Distributions

hist_bins <- 50
hist_fill <- "lightblue"
hist_color <- "black"
vline_color <- "red"

# Plot 1: Bootstrap Distribution of Mean Wage
ggplot(data.frame(Value = boot_wage$t), aes(x = Value)) +
  geom_histogram(bins = hist_bins, fill = hist_fill, color = hist_color, alpha = 0.7) +
  geom_vline(xintercept = wage_mean, color = vline_color, linetype = "dashed") +
  labs(
    title = "Bootstrap Distribution of Mean Wage",
    x = "Bootstrap Mean Wage",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10)
  )

# Plot 2: Bootstrap Distribution of Mean Log Wage
ggplot(data.frame(Value = boot_logwage$t), aes(x = Value)) +
  geom_histogram(bins = hist_bins, fill = hist_fill, color = hist_color, alpha = 0.7) +
  geom_vline(xintercept = logwage_mean, color = vline_color, linetype = "dashed") +
  labs(
    title = "Bootstrap Distribution of Mean Log Wage",
    x = "Bootstrap Mean Log Wage",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10)
  )

# Plot 3: Bootstrap Distribution of Median Wage
ggplot(data.frame(Value = boot_median_wage$t), aes(x = Value)) +
  geom_histogram(bins = hist_bins, fill = hist_fill, color = hist_color, alpha = 0.7) +
  geom_vline(xintercept = median(cps$wage), color = vline_color, linetype = "dashed") +
  labs(
    title = "Bootstrap Distribution of Median Wage",
    x = "Bootstrap Median Wage",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10)
  )

# Plot 4: Bootstrap Distribution of Median Log Wage
ggplot(data.frame(Value = boot_median_logwage$t), aes(x = Value)) +
  geom_histogram(bins = hist_bins, fill = hist_fill, color = hist_color, alpha = 0.7) +
  geom_vline(xintercept = median(cps$log_wage), color = vline_color, linetype = "dashed") +
  labs(
    title = "Bootstrap Distribution of Median Log Wage",
    x = "Bootstrap Median Log Wage",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10)
  )

#---------------------------------------------------------------
# Bootstrap OLS Regression using the .632 Bootstrap Estimator
#---------------------------------------------------------------
## 2.1 Define Bootstrap Function for OLS Regression
# Bootstrap function to obtain OLS coefficients
boot_ols_full <- function(data, indices) {
  d <- data[indices, ]  # Resample data (bootstrapping pairs)
  model <- lm(ols_formula, data = d)
  return(coef(model))  # Return coefficients
}

# Set number of bootstrap samples
R = 300

# Perform bootstrap
boot_ols_results_full <- boot(data = train_data, statistic = boot_ols_full, R = R)
print(boot_ols_results_full)

## 2.3 Compute .632 Bootstrap Estimates of Test MSE
# Extract bootstrap sample indices and initialize storage vectors
bootstrap_indices <- boot.array(boot_ols_results_full)
bootstrap_mse_train <- numeric(boot_ols_results_full$R)
bootstrap_mse_test <- numeric(boot_ols_results_full$R)
bootstrap_mse_632 <- numeric(boot_ols_results_full$R)
n_total <- nrow(cps)

# Calculate in-bag and out-of-bag MSEs
for (i in 1:boot_ols_results_full$R) {
  freq <- bootstrap_indices[i, ]
  in_bag <- which(freq > 0)
  out_of_bag <- which(freq == 0)
  coef_i <- boot_ols_results_full$t[i, ]
  
  # MSE for in-bag
  pred_inbag <- model.matrix(ols_formula, data = cps[in_bag, ]) %*% coef_i
  mse_inbag <- mean((cps$log_wage[in_bag] - pred_inbag)^2)
  bootstrap_mse_train[i] <- mse_inbag
  
  # MSE for out-of-bag (if any)
  if (length(out_of_bag) > 0) {
    pred_oob <- model.matrix(ols_formula, data = cps[out_of_bag, ]) %*% coef_i
    mse_oob <- mean((cps$log_wage[out_of_bag] - pred_oob)^2)
    bootstrap_mse_test[i] <- mse_oob
  } else {
    bootstrap_mse_test[i] <- NA
  }
  bootstrap_mse_632[i] <- 0.368 * bootstrap_mse_train[i] + 0.632 * bootstrap_mse_test[i]
}

# Compute .632 bootstrap estimator of MSE
mse_632 <- mean(bootstrap_mse_632)
cat(".632 Bootstrap Estimate of Test MSE:", mse_632, "\n")

# Compare with original MSE on train full data
original_mse <- mean(residuals(ols_model)^2)
cat("Original train OLS MSE:", original_mse, "\n")

# Append result to `mse_results`
mse_results <- rbind(mse_results, data.frame(Model = ".632 Bootstrap OLS Regression", Test_MSE = mse_632))

## 2.4 Diagnostics and Visualization
### 2.4.1 Plot Bootstrap Distributions of Coefficients
# Extract coefficient names and bootstrap values for plotting
coef_names <- names(coef(ols_model))
coef_values <- as.data.frame(boot_ols_results_full$t)
colnames(coef_values) <- coef_names

# Generate individual histograms for each coefficient
for (coef_name in coef_names) {
  p <- ggplot(coef_values, aes(x = .data[[coef_name]])) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = coef(ols_model)[coef_name], color = "red", linetype = "dashed") +
    labs(
      title = paste("Bootstrap Distribution of", coef_name),
      x = coef_name,
      y = "Frequency"
    ) +
    theme_minimal()
  
  # Print each plot to ensure it renders in a loop
  print(p)
}

### 2.4.2 Residual Analysis
# Compute average residuals from bootstrapped models
residuals_matrix <- sapply(1:boot_ols_results_full$R, function(i) {
  coef_i <- boot_ols_results_full$t[i, ]
  pred_i <- model.matrix(ols_formula, data = cps) %*% coef_i
  cps$log_wage - pred_i
})
avg_residuals <- rowMeans(residuals_matrix)

# Plot histogram of average residuals from bootstrapped models
ggplot(data.frame(AverageResiduals = avg_residuals), aes(x = AverageResiduals)) +
  geom_histogram(bins = 100, fill = "lightblue", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Average Residuals from Bootstrapped Models",
    x = "Residuals",
    y = "Frequency"
  ) +
  #xlim(-5, 5) +
  theme_minimal()

#---------------------------------------------------------------
# Bootstrap MM-Estimator Regression using the .632 Bootstrap Estimator
#---------------------------------------------------------------
## 1. Define Bootstrap Function for MM-Estimator Regression
# Bootstrap function to obtain MM coefficients
boot_mm_full <- function(data, indices) {
  d <- data[indices, ]  # Resample data (bootstrapping pairs)
  model <- RobStatTM::lmrobdetMM(ols_formula, data = d, control = control_lmrobdetMM)
  return(coef(model))  # Return coefficients
}

# Perform bootstrap
# Detect number of available cores
n_cores <- detectCores() - 1  

# Set up a cluster for parallel computing
cl <- makeCluster(n_cores)
clusterExport(cl, c("cps", "ols_formula", "control_lmrobdetMM", "boot_mm_full"))  # Export required variables and functions

# Enable parallel processing for the 'boot' function
boot_mm_results_full <- boot(
  data = cps,
  statistic = boot_mm_full,
  R = R,
  parallel = "snow",
  ncpus = n_cores,
  cl = cl
)

# Stop the cluster
stopCluster(cl)
print(boot_mm_results_full)

## 2. Compute .632 Bootstrap Estimates of Test MSE
# Extract bootstrap sample indices and initialize storage vectors
bootstrap_mm_indices <- boot.array(boot_mm_results_full)
bootstrap_mse_mm_train <- numeric(boot_mm_results_full$R)
bootstrap_mse_mm_test <- numeric(boot_mm_results_full$R)
bootstrap_mse_mm_632 <- numeric(boot_mm_results_full$R)
n_total <- nrow(cps)

# Calculate in-bag and out-of-bag MSEs
for (i in 1:boot_mm_results_full$R) {
  freq <- bootstrap_mm_indices[i, ]
  in_bag <- which(freq > 0)
  out_of_bag <- which(freq == 0)
  coef_i <- boot_mm_results_full$t[i, ]
  
  # MSE for in-bag
  pred_inbag <- model.matrix(ols_formula, data = cps[in_bag, ]) %*% coef_i
  mse_inbag <- mean((cps$log_wage[in_bag] - pred_inbag)^2)
  bootstrap_mse_mm_train[i] <- mse_inbag
  
  # MSE for out-of-bag (if any)
  if (length(out_of_bag) > 0) {
    pred_oob <- model.matrix(ols_formula, data = cps[out_of_bag, ]) %*% coef_i
    mse_oob <- mean((cps$log_wage[out_of_bag] - pred_oob)^2)
    bootstrap_mse_mm_test[i] <- mse_oob
  } else {
    bootstrap_mse_mm_test[i] <- NA
  }
  bootstrap_mse_mm_632[i] <- 0.368 * bootstrap_mse_mm_train[i] + 0.632 * bootstrap_mse_mm_test[i]
}

# Compute .632 bootstrap estimator of MSE
mse_mm_632 <- mean(bootstrap_mse_mm_632, na.rm = TRUE)
cat(".632 Bootstrap Estimate of Test MSE for MM-Estimator:", mse_mm_632, "\n")

# Append result to `mse_results`
mse_results <- rbind(mse_results, data.frame(Model = ".632 Bootstrap MM-Estimator Regression", Test_MSE = mse_mm_632))

## 3.4 Compare OLS and Robust Test MSEs and Test Hypothesis of Equal MSEs
# Ensure that the number of bootstrap iterations is the same and align the MSEs
# Compute the MSE difference between Robust and OLS models
bootstrap_mse_diff <- bootstrap_mse_mm_632 - bootstrap_mse_632 
mse_difference <- mse_mm_632 - mse_632 

# Compute 5% Bootstrap Confidence Interval for the difference in MSEs
alpha <- 0.05
ci_lower <- quantile(bootstrap_mse_diff, probs = alpha / 2)
ci_upper <- quantile(bootstrap_mse_diff, probs = 1 - alpha / 2)

cat("95% Bootstrap Confidence Interval for Difference in MSEs:", ci_lower, "to", ci_upper, "\n")

# Plot histogram of bootstrap MSE differences
ggplot(data.frame(mse_diff = bootstrap_mse_diff), aes(x = mse_diff)) +
  geom_histogram(binwidth = 0.005, color = "black", fill = "skyblue") +
  labs(
    title = "Histogram of Bootstrap MSE Differences",
    x = "Bootstrap MSE Difference",
    y = "Frequency"
  ) +
  theme_minimal()

# Interpret the results
if ((ci_lower > 0) || (ci_upper < 0)) {
  cat("There is a significant difference between the MSEs of the LME and OLS models at the 5% significance level.\n")
} else {
  cat("There is no significant difference between the MSEs of the LME and OLS models at the 5% significance level.\n")
}

#---------------------------------------------------------------
# Comparing MSE Results of Different Models
#---------------------------------------------------------------
# View the MSE results
print(mse_results)

# Plot the Test MSEs for each model
ggplot(mse_results[order(mse_results$Test_MSE), ], aes(x = reorder(Model, Test_MSE), y = Test_MSE)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(Test_MSE, 6)), vjust = -0.5, size = 3) +
  labs(
    title = "Comparison of Test MSE Across Models",
    x = "Model",
    y = "Test MSE"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
