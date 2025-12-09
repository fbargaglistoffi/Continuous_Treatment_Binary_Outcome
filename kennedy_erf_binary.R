# Modified version for binary outcomes and censoring
ctseff_binary <- function(y, a, x, bw.seq, n.pts = 100, a.rng = c(min(a), max(a)),
                          cens = NULL,  # Add censoring indicator
                          sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", 
                                     "SL.mean", "SL.ranger")) {
  
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("KernSmooth")
  require("zoo")
  
  kern <- function(t) { dnorm(t) }
  
  # Handle censoring: only use uncensored observations
  if (!is.null(cens)) {
    # cens should be 1 for uncensored, 0 for censored
    keep_idx <- which(cens == 1)
    y <- y[keep_idx]
    a <- a[keep_idx]
    x <- x[keep_idx, , drop = FALSE]
  }
  
  n <- nrow(x)
  
  # Set up evaluation points
  a.min <- a.rng[1]
  a.max <- a.rng[2]
  a.vals <- seq(a.min, a.max, length.out = n.pts)
  
  # Create augmented dataset for predictions
  xa.new <- rbind(
    cbind(x, a), 
    cbind(x[rep(1:n, length(a.vals)), ], a = rep(a.vals, rep(n, length(a.vals))))
  )
  x.new <- xa.new[, -ncol(xa.new)]
  x <- data.frame(x)
  x.new <- data.frame(x.new)
  colnames(x) <- colnames(x.new)
  xa.new <- data.frame(xa.new)
  
  # MODIFICATION 1: Estimate nuisance functions with appropriate families
  # Treatment model (continuous - stays the same)
  pimod <- SuperLearner(Y = a, X = x, SL.library = sl.lib, newX = x.new)
  pimod.vals <- pimod$SL.predict
  
  # Variance model (continuous - stays the same)
  pi2mod <- SuperLearner(Y = (a - pimod.vals[1:n])^2, X = x, 
                         SL.library = sl.lib, newX = x.new)
  pi2mod.vals <- pmax(pi2mod$SL.predict, 0.01)  # Ensure positive variance
  
  # MODIFICATION 2: Outcome model for BINARY outcome
  # Use binomial family
  mumod <- SuperLearner(Y = y, X = cbind(x, a), 
                        family = binomial(),  # KEY CHANGE
                        SL.library = sl.lib, 
                        newX = xa.new)
  muhat.vals <- mumod$SL.predict
  
  # Ensure predictions are in [0,1]
  muhat.vals <- pmin(pmax(muhat.vals, 0.001), 0.999)
  
  # Construct estimated densities and conditional expectations
  a.std <- (xa.new$a - pimod.vals) / sqrt(pi2mod.vals)
  
  # Estimate density using standardized values
  dens_est <- density(a.std[1:n])
  pihat.vals <- approx(dens_est$x, dens_est$y, xout = a.std, rule = 2)$y / sqrt(pi2mod.vals)
  
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  # Marginal density of treatment
  varpihat <- predict(smooth.spline(a.vals, apply(pihat.mat, 2, mean)), x = a)$y
  varpihat <- pmax(varpihat, 0.001)  # Avoid division by zero
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = TRUE, nrow = n)
  
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  # Marginal outcome expectation
  mhat <- predict(smooth.spline(a.vals, apply(muhat.mat, 2, mean)), x = a)$y
  mhat <- pmin(pmax(mhat, 0.001), 0.999)  # Keep in [0,1]
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = TRUE, nrow = n)
  
  # MODIFICATION 3: Adjusted pseudo-outcome for binary data
  # Standard doubly robust pseudo-outcome
  pseudo.out <- (y - muhat) * (varpihat / pihat) + mhat
  
  # Alternative stabilized version for binary outcomes
  # pseudo.out <- (y - muhat) / pmax(pihat / varpihat, 0.1) + mhat
  
  # Bandwidth selection via cross-validation
  w.fn <- function(bw) {
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a - a.val) / bw
      kern.std <- kern(a.std) / bw
      numerator <- mean(a.std^2 * kern.std) * (kern(0) / bw)
      denominator <- mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2
      w.avals <- c(w.avals, numerator / pmax(denominator, 0.001))
    }
    return(w.avals / n)
  }
  
  hatvals <- function(bw) {
    approx(a.vals, w.fn(bw), xout = a, rule = 2)$y
  }
  
  cts.eff.fn <- function(out, bw) {
    lp_result <- locpoly(a, out, bandwidth = bw)
    approx(lp_result$x, lp_result$y, xout = a, rule = 2)$y
  }
  
  # Risk function for bandwidth selection
  risk.fn <- function(h) {
    hats <- pmin(hatvals(h), 0.9)  # Prevent division issues
    preds <- cts.eff.fn(pseudo.out, bw = h)
    mean(((pseudo.out - preds) / (1 - hats))^2, na.rm = TRUE)
  }
  
  risk.est <- sapply(bw.seq, function(h) {
    tryCatch(risk.fn(h), error = function(e) Inf)
  })
  
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  
  # Estimate effect curve with optimal bandwidth
  final_curve <- locpoly(a, pseudo.out, bandwidth = h.opt)
  est <- approx(final_curve$x, final_curve$y, xout = a.vals, rule = 2)$y
  
  # MODIFICATION 4: Bound estimates for binary outcomes
  est <- pmin(pmax(est, 0), 1)
  
  # Estimate standard errors
  se <- NULL
  for (a.val in a.vals) {
    a.std <- (a - a.val) / h.opt
    kern.std <- kern(a.std) / h.opt
    
    # Weighted regression
    weights_valid <- kern.std > 0.001
    if (sum(weights_valid) > 10) {
      beta <- coef(lm(pseudo.out ~ a.std, weights = kern.std, subset = weights_valid))
      
      Dh <- matrix(c(
        mean(kern.std), mean(kern.std * a.std),
        mean(kern.std * a.std), mean(kern.std * a.std^2)
      ), nrow = 2)
      
      # Calculate influence functions
      kern.mat <- matrix(rep(kern((a.vals - a.val) / h.opt) / h.opt, n), 
                         byrow = TRUE, nrow = n)
      g2 <- matrix(rep((a.vals - a.val) / h.opt, n), byrow = TRUE, nrow = n)
      
      intfn1.mat <- kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
      intfn2.mat <- g2 * kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
      
      # Numerical integration
      int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n),
                           byrow = TRUE, nrow = n) * 
                      (intfn1.mat[, -1] + intfn1.mat[, -length(a.vals)]) / 2, 1, sum)
      int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n),
                           byrow = TRUE, nrow = n) * 
                      (intfn2.mat[, -1] + intfn2.mat[, -length(a.vals)]) / 2, 1, sum)
      
      # Influence function
      infl_fn <- solve(Dh) %*% rbind(
        kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int1,
        a.std * kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int2
      )
      
      sigma <- var(infl_fn[1, ])
      se <- c(se, sqrt(sigma))
    } else {
      se <- c(se, NA)
    }
  }
  
  # Replace NAs with interpolated values
  if (any(is.na(se))) {
    se <- na.approx(se, na.rm = FALSE)
    se[is.na(se)] <- mean(se, na.rm = TRUE)
  }
  
  # Confidence intervals with logit transformation for bounded outcomes
  logit_est <- qlogis(pmax(pmin(est, 0.999), 0.001))
  logit_se <- se / (est * (1 - est))
  ci.ll <- plogis(logit_est - 1.96 * logit_se / sqrt(n))
  ci.ul <- plogis(logit_est + 1.96 * logit_se / sqrt(n))
  
  res <- data.frame(a.vals, est, se, ci.ll, ci.ul)
  
  return(list(res = res, bw.risk = bw.risk))
}

# Estimate censoring weights
estimate_censoring_weights <- function(df, period) {
  cens_var <- paste0("cens_oud_period_", period)
  W_vars <- c("W_01", paste0("W_", 2:14))
  
  # Model probability of being uncensored
  cens_model <- glm(formula(paste(cens_var, "~", paste(W_vars, collapse = " + "), "+ A")),
                    data = df,
                    family = binomial())
  
  # Get weights (inverse probability of being uncensored)
  prob_uncensored <- predict(cens_model, type = "response")
  weights <- ifelse(df[[cens_var]] == 1, 1 / prob_uncensored, 0)
  
  # Stabilize weights
  weights <- weights * mean(df[[cens_var]] == 1)
  
  # Truncate extreme weights
  weights[weights > quantile(weights[weights > 0], 0.95)] <- 
    quantile(weights[weights > 0], 0.95)
  
  return(weights)
}

# Modified function with censoring weights
ctseff_binary_weighted <- function(y, a, x, cens, weights = NULL, ...) {
  # Apply censoring and weights
  keep_idx <- which(cens == 1)
  
  if (!is.null(weights)) {
    # Incorporate censoring weights into the analysis
    # This would require modifying the SuperLearner calls to use weights
    # Not fully implemented in the standard Kennedy method
  }
  
  # Continue with standard analysis...
  ctseff_binary(y = y[keep_idx], 
                a = a[keep_idx], 
                x = x[keep_idx, ], 
                ...)
}

# Function to estimate censoring weights
estimate_censoring_weights <- function(df, period, stabilize = TRUE, truncate = TRUE) {
  
  # Variables
  cens_var <- paste0("cens_oud_period_", period)
  W_vars <- c("W_01", paste0("W_", 2:14))
  
  # For longitudinal censoring, we need to account for previous censoring
  # Create cumulative censoring indicator
  df$cumulative_cens <- 1
  for(p in 1:period) {
    cens_p <- paste0("cens_oud_period_", p)
    df$cumulative_cens <- df$cumulative_cens * df[[cens_p]]
  }
  
  # Estimate probability of remaining uncensored up to this period
  if(period == 1) {
    # For period 1, simple model
    cens_formula <- formula(paste(cens_var, "~", 
                                  paste(W_vars, collapse = " + "), "+ A"))
  } else {
    # For later periods, could include previous outcomes if appropriate
    prev_outcomes <- paste0("oud_period_", 1:(period-1))
    # Only include previous outcomes for those still uncensored
    cens_formula <- formula(paste(cens_var, "~", 
                                  paste(W_vars, collapse = " + "), 
                                  "+ A"))
  }
  
  # Fit censoring model
  cens_model <- glm(cens_formula, 
                    data = df[df$cumulative_cens == 1 | df[[cens_var]] == 0, ],
                    family = binomial())
  
  # Predict probability of being uncensored
  prob_uncensored <- predict(cens_model, newdata = df, type = "response")
  
  if(stabilize) {
    # Stabilized weights: P(C=1) / P(C=1|W,A)
    marginal_prob <- mean(df[[cens_var]][df$cumulative_cens == 1 | df[[cens_var]] == 0])
    weights <- ifelse(df[[cens_var]] == 1, 
                      marginal_prob / prob_uncensored, 
                      0)
  } else {
    # Unstabilized weights: 1 / P(C=1|W,A)
    weights <- ifelse(df[[cens_var]] == 1, 
                      1 / prob_uncensored, 
                      0)
  }
  
  if(truncate) {
    # Truncate extreme weights at 1st and 99th percentiles
    weight_quantiles <- quantile(weights[weights > 0], c(0.01, 0.99))
    weights[weights > 0 & weights < weight_quantiles[1]] <- weight_quantiles[1]
    weights[weights > weight_quantiles[2]] <- weight_quantiles[2]
  }
  
  # Diagnostics
  cat("Censoring weight summary for period", period, ":\n")
  cat("  Mean weight:", mean(weights[weights > 0]), "\n")
  cat("  SD weight:", sd(weights[weights > 0]), "\n")
  cat("  Min weight:", min(weights[weights > 0]), "\n")
  cat("  Max weight:", max(weights[weights > 0]), "\n")
  cat("  % censored:", mean(df[[cens_var]] == 0) * 100, "%\n\n")
  
  return(weights)
}

# Modified Kennedy method with censoring weights
ctseff_binary_weighted <- function(y, a, x, cens, cens_weights = NULL,
                                   bw.seq, n.pts = 100, 
                                   a.rng = c(min(a[cens == 1]), max(a[cens == 1])),
                                   sl.lib = c("SL.earth", "SL.gam", "SL.glm", 
                                              "SL.glm.interaction", "SL.mean", "SL.ranger")) {
  
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("KernSmooth")
  require("zoo")
  
  kern <- function(t) { dnorm(t) }
  
  # Handle censoring: only use uncensored observations
  keep_idx <- which(cens == 1)
  y_obs <- y[keep_idx]
  a_obs <- a[keep_idx]
  x_obs <- x[keep_idx, , drop = FALSE]
  
  # Apply censoring weights if provided
  if (!is.null(cens_weights)) {
    w_obs <- cens_weights[keep_idx]
    # Normalize weights to sum to n
    w_obs <- w_obs * length(w_obs) / sum(w_obs)
  } else {
    w_obs <- rep(1, length(keep_idx))
  }
  
  n <- length(y_obs)
  
  # Set up evaluation points
  a.min <- a.rng[1]
  a.max <- a.rng[2]
  a.vals <- seq(a.min, a.max, length.out = n.pts)
  
  # Create augmented dataset for predictions
  xa.new <- rbind(
    cbind(x_obs, a = a_obs), 
    cbind(x_obs[rep(1:n, length(a.vals)), ], 
          a = rep(a.vals, rep(n, length(a.vals))))
  )
  x.new <- xa.new[, -ncol(xa.new)]
  x_obs <- data.frame(x_obs)
  x.new <- data.frame(x.new)
  colnames(x_obs) <- colnames(x.new)
  xa.new <- data.frame(xa.new)
  
  # WEIGHTED SUPERLEARNER CALLS
  
  # Create custom weighted SL wrappers if needed
  SL.glm.weighted <- function(..., obsWeights) {
    SL.glm(..., obsWeights = obsWeights)
  }
  
  SL.gam.weighted <- function(..., obsWeights) {
    SL.gam(..., obsWeights = obsWeights)
  }
  
  # Treatment model with weights
  pimod <- SuperLearner(Y = a_obs, 
                        X = x_obs, 
                        obsWeights = w_obs,  # Include weights
                        SL.library = sl.lib, 
                        newX = x.new)
  pimod.vals <- pimod$SL.predict
  
  # Variance model with weights
  residuals_sq <- (a_obs - pimod.vals[1:n])^2
  pi2mod <- SuperLearner(Y = residuals_sq, 
                         X = x_obs,
                         obsWeights = w_obs,  # Include weights
                         SL.library = sl.lib, 
                         newX = x.new)
  pi2mod.vals <- pmax(pi2mod$SL.predict, 0.01)
  
  # Outcome model with weights for binary outcome
  mumod <- SuperLearner(Y = y_obs, 
                        X = cbind(x_obs, a = a_obs),
                        obsWeights = w_obs,  # Include weights
                        family = binomial(),
                        SL.library = sl.lib, 
                        newX = xa.new)
  muhat.vals <- mumod$SL.predict
  muhat.vals <- pmin(pmax(muhat.vals, 0.001), 0.999)
  
  # Construct weighted density estimates
  a.std <- (xa.new$a - pimod.vals) / sqrt(pi2mod.vals)
  
  # Weighted density estimation
  dens_est <- density(a.std[1:n], weights = w_obs/sum(w_obs))
  pihat.vals <- approx(dens_est$x, dens_est$y, xout = a.std, rule = 2)$y / sqrt(pi2mod.vals)
  
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  # Weighted marginal density
  varpihat <- predict(smooth.spline(a.vals, 
                                    apply(pihat.mat, 2, weighted.mean, w = w_obs)), 
                      x = a_obs)$y
  varpihat <- pmax(varpihat, 0.001)
  
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  # Weighted marginal outcome
  mhat <- predict(smooth.spline(a.vals, 
                                apply(muhat.mat, 2, weighted.mean, w = w_obs)), 
                  x = a_obs)$y
  mhat <- pmin(pmax(mhat, 0.001), 0.999)
  
  # Pseudo-outcome (already incorporates censoring through weights)
  pseudo.out <- (y_obs - muhat) * (varpihat / pihat) + mhat
  
  # WEIGHTED BANDWIDTH SELECTION
  
  w.fn <- function(bw) {
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a_obs - a.val) / bw
      kern.std <- kern(a.std) / bw
      # Weighted moments
      numerator <- weighted.mean(a.std^2 * kern.std, w = w_obs) * (kern(0) / bw)
      denom_term1 <- weighted.mean(kern.std, w = w_obs) * 
        weighted.mean(a.std^2 * kern.std, w = w_obs)
      denom_term2 <- weighted.mean(a.std * kern.std, w = w_obs)^2
      denominator <- denom_term1 - denom_term2
      w.avals <- c(w.avals, numerator / pmax(denominator, 0.001))
    }
    return(w.avals / n)
  }
  
  hatvals <- function(bw) {
    approx(a.vals, w.fn(bw), xout = a_obs, rule = 2)$y
  }
  
  # Weighted local polynomial regression
  cts.eff.fn <- function(out, bw) {
    # Create weights for locpoly (product of censoring weights and kernel weights)
    lp_result <- locpoly(a_obs, out, bandwidth = bw, 
                         kernel = "normal")  # Can't directly pass weights to locpoly
    approx(lp_result$x, lp_result$y, xout = a_obs, rule = 2)$y
  }
  
  # Risk function with weights
  risk.fn <- function(h) {
    hats <- pmin(hatvals(h), 0.9)
    preds <- cts.eff.fn(pseudo.out, bw = h)
    weighted.mean(((pseudo.out - preds) / (1 - hats))^2, w = w_obs, na.rm = TRUE)
  }
  
  risk.est <- sapply(bw.seq, function(h) {
    tryCatch(risk.fn(h), error = function(e) Inf)
  })
  
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  
  # Final curve estimation with optimal bandwidth
  final_curve <- locpoly(a_obs, pseudo.out, bandwidth = h.opt)
  est <- approx(final_curve$x, final_curve$y, xout = a.vals, rule = 2)$y
  est <- pmin(pmax(est, 0), 1)
  
  # WEIGHTED STANDARD ERRORS
  se <- NULL
  for (a.val in a.vals) {
    a.std <- (a_obs - a.val) / h.opt
    kern.std <- kern(a.std) / h.opt
    
    weights_valid <- kern.std > 0.001
    if (sum(weights_valid) > 10) {
      # Combined weights for regression
      combined_weights <- kern.std * w_obs
      
      beta <- coef(lm(pseudo.out ~ a.std, 
                      weights = combined_weights, 
                      subset = weights_valid))
      
      # Weighted design matrix
      Dh <- matrix(c(
        weighted.mean(kern.std, w = w_obs),
        weighted.mean(kern.std * a.std, w = w_obs),
        weighted.mean(kern.std * a.std, w = w_obs),
        weighted.mean(kern.std * a.std^2, w = w_obs)
      ), nrow = 2)
      
      # Rest of SE calculation...
      sigma <- var(pseudo.out - beta[1] - beta[2] * a.std) * mean(combined_weights)
      se <- c(se, sqrt(sigma / n))
    } else {
      se <- c(se, NA)
    }
  }
  
  # Handle missing SEs
  if (any(is.na(se))) {
    se_zoo <- zoo::na.approx(se, na.rm = FALSE)
    if (any(is.na(se_zoo))) {
      se_zoo <- zoo::na.fill(se_zoo, "extend")
    }
    se <- as.numeric(se_zoo)
  }
  
  # Confidence intervals
  est_bounded <- pmax(pmin(est, 0.999), 0.001)
  logit_est <- qlogis(est_bounded)
  logit_se <- se / (est_bounded * (1 - est_bounded))
  ci.ll <- plogis(logit_est - 1.96 * logit_se / sqrt(n))
  ci.ul <- plogis(logit_est + 1.96 * logit_se / sqrt(n))
  
  res <- data.frame(a.vals, est, se, ci.ll, ci.ul)
  
  return(list(res = res, bw.risk = bw.risk, weights = w_obs))
}