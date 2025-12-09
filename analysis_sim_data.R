rm(list = ls())

setwd("Desktop/KR-project/")
data <- readRDS("cs_cohort_simulated_data.rds")
#source("kennedy_erf_binary.R") # uncomment to run simple complete-case analysis
source("kennedy_erf_binary_censoring.R")
source("plotting.R")

# Rework of censoring variables to ensure once censored, always censored

# Create a function to fix censoring
fix_censoring <- function(data) {
  # Work with a copy
  df_fixed <- data
  
  # For each row, once censored (0), remain censored
  for(i in 1:nrow(df_fixed)) {
    if(df_fixed$cens_oud_period_1[i] == 0) {
      df_fixed$cens_oud_period_2[i] <- 0
      df_fixed$cens_oud_period_3[i] <- 0
      df_fixed$cens_oud_period_4[i] <- 0
    } else if(df_fixed$cens_oud_period_2[i] == 0) {
      df_fixed$cens_oud_period_3[i] <- 0
      df_fixed$cens_oud_period_4[i] <- 0
    } else if(df_fixed$cens_oud_period_3[i] == 0) {
      df_fixed$cens_oud_period_4[i] <- 0
    }
  }
  
  return(df_fixed)
}

# Apply the fix
df <- fix_censoring(data)

# Verify the fix
# Check that censoring is monotonic (once 0, always 0)
check_censoring <- function(data) {
  issues <- 0
  for(i in 1:nrow(data)) {
    cens_vals <- c(data$cens_oud_period_1[i], 
                   data$cens_oud_period_2[i],
                   data$cens_oud_period_3[i], 
                   data$cens_oud_period_4[i])
    
    # Check if there's a 1 after a 0
    for(j in 2:length(cens_vals)) {
      if(cens_vals[j-1] == 0 && cens_vals[j] == 1) {
        issues <- issues + 1
        cat("Issue in row", i, "\n")
      }
    }
  }
  cat("Total issues found:", issues, "\n")
}

check_censoring(df)


# Analysis function incorporating censoring weights
analyze_continuous_treatment_effect <- function(df, period, 
                                                bw.seq = seq(0.1, 2, length.out = 50),
                                                use_censoring_weights = TRUE) {
  
  # Extract variables for this period
  cens_var <- paste0("cens_oud_period_", period)
  oud_var <- paste0("oud_period_", period)
  
  # Covariates
  W_vars <- c("W_01", paste0("W_", 2:14))
  
  # Apply censoring
  cens <- df[[cens_var]]
  y <- df[[oud_var]]
  a <- df$A 
  x <- as.matrix(df[, W_vars])
  
  # Run modified Kennedy method (now with censoring weights by default)
  results <- ctseff_binary(
    y = y,
    a = a,
    x = x,
    cens = cens,
    bw.seq = bw.seq,
    n.pts = 100,
    use_censoring_weights = use_censoring_weights,  # Can toggle on/off
    stabilize_weights = TRUE,
    truncate_weights = TRUE,
    sl.lib = c("SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
    verbose = (period == 1)  # Only verbose for first period
  )
  
  return(results)
}

# Run analysis for all periods (with censoring correction)
all_results <- list()
for (period in 1:4) {
  cat("\nAnalyzing Period", period, "\n")
  all_results[[period]] <- analyze_continuous_treatment_effect(df, period)
  
  # Check bandwidth selection
  cat("Optimal bandwidth:", 
      all_results[[period]]$bw.risk$bw[which.min(all_results[[period]]$bw.risk$risk)], "\n")
}


# Customize density appearance
plot_dose_response(all_results, df, 
                            xlim = c(50, 350), 
                            ylim = c(0, 0.1),
                            density_height = 0.2,
                            density_color = "#2E86AB",
                            density_alpha = 0.2)

# OPTIONAL: Uncomment to Run without censoring weights for comparison

# # If you want to compare with uncorrected version
# all_results_uncorrected <- list()
# for (period in 1:4) {
#   cat("\nAnalyzing Period", period, " (no censoring weights)\n")
#   all_results_uncorrected[[period]] <- analyze_continuous_treatment_effect(
#     df, period, 
#     use_censoring_weights = FALSE  # Turn off censoring weights
#   )
# }

# Below I played around with possible binnings of the exposure to see if results were sensitive to binning (TLDR: results did not change much)
# Binning Exposure
bin_width = 20

df$A_bin <- floor(df$A / bin_width) * bin_width + bin_width/2

plot(density(df$A))
plot(density(df$A_bin))

analyze_continuous_treatment_effect <- function(df, period, 
                                                bw.seq = seq(0.1, 2, length.out = 50),
                                                use_censoring_weights = TRUE) {
  
  # Extract variables for this period
  cens_var <- paste0("cens_oud_period_", period)
  oud_var <- paste0("oud_period_", period)
  
  # Covariates
  W_vars <- c("W_01", paste0("W_", 2:14))
  
  # Apply censoring
  cens <- df[[cens_var]]
  y <- df[[oud_var]]
  a <- df$A_bin 
  x <- as.matrix(df[, W_vars])
  
  # Run modified Kennedy method (now with censoring weights by default)
  results <- ctseff_binary(
    y = y,
    a = a,
    x = x,
    cens = cens,
    bw.seq = bw.seq,
    n.pts = 100,
    use_censoring_weights = use_censoring_weights,  # Can toggle on/off
    stabilize_weights = TRUE,
    truncate_weights = TRUE,
    sl.lib = c("SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
    verbose = (period == 1)  # Only verbose for first period
  )
  
  return(results)
}

# Run analysis for all periods
all_results <- list()
for (period in 1:4) {
  cat("\nAnalyzing Period", period, "\n")
  all_results[[period]] <- analyze_continuous_treatment_effect(df, period)
  
  # Check bandwidth selection
  cat("Optimal bandwidth:", 
      all_results[[period]]$bw.risk$bw[which.min(all_results[[period]]$bw.risk$risk)], "\n")
}


# Customize density appearance
plot_dose_response(all_results, df, 
                            xlim = c(50, 350), 
                            ylim = c(0, 0.1),
                            density_height = 0.2,
                            density_color = "#2E86AB",
                            density_alpha = 0.2)
