# Kennedy Doubly Robust Method with Incorporated Censoring Weigths

ctseff_binary <- function(y, a, x, bw.seq, 
                          n.pts = 100, 
                          a.rng = NULL,  # Changed to NULL default
                          cens = NULL,
                          use_censoring_weights = TRUE,
                          stabilize_weights = TRUE,
                          truncate_weights = TRUE,
                          weight_quantiles = c(0.01, 0.99),
                          sl.lib = c("SL.earth", "SL.gam", "SL.glm", 
                                     "SL.glm.interaction", "SL.mean", "SL.ranger"),
                          verbose = FALSE) {
  
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("KernSmooth")
  require("zoo")
  
  kern <- function(t) { dnorm(t) }
  
  # Store full data for censoring weight estimation
  full_y <- y
  full_a <- a
  full_x <- x
  full_cens <- cens
  
  # STEP 1: Handle censoring and determine range
  
  # First, identify uncensored observations
  if (!is.null(cens)) {
    keep_idx <- which(cens == 1)
  } else {
    keep_idx <- 1:length(y)
  }
  
  # Set range based on uncensored data if not provided
  if (is.null(a.rng)) {
    if (!is.null(cens)) {
      a.rng <- c(min(a[keep_idx]), max(a[keep_idx]))
    } else {
      a.rng <- c(min(a), max(a))
    }
  }
  
  # STEP 2: Estimate Censoring Weights (if requested)
  
  w_cens <- rep(1, length(y))  # Default: no weighting
  
  if (!is.null(cens) && use_censoring_weights) {
    if(verbose) cat("Estimating censoring weights...\n")
    
    # Fit censoring model on full data
    cens_data <- data.frame(cens = cens, x, a = a)
    cens_model <- glm(cens ~ ., 
                      data = cens_data,
                      family = binomial())
    
    # Predict censoring probabilities
    prob_uncensored <- predict(cens_model, type = "response")
    
    # Calculate weights
    if(stabilize_weights) {
      marginal_prob <- mean(cens)
      w_cens <- ifelse(cens == 1, marginal_prob / prob_uncensored, 0)
    } else {
      w_cens <- ifelse(cens == 1, 1 / prob_uncensored, 0)
    }
    
    # Ensure no extreme weights
    w_cens[!is.finite(w_cens)] <- 0
    w_cens <- pmax(w_cens, 0)
    
    # Truncate extreme weights
    if(truncate_weights && sum(w_cens > 0) > 0) {
      wq <- quantile(w_cens[w_cens > 0], weight_quantiles)
      w_cens[w_cens > 0 & w_cens < wq[1]] <- wq[1]
      w_cens[w_cens > wq[2]] <- wq[2]
    }
    
    if(verbose) {
      cat("Censoring rate:", round(mean(cens == 0) * 100, 1), "%\n")
      cat("Weight summary: mean =", round(mean(w_cens[w_cens > 0]), 2),
          ", sd =", round(sd(w_cens[w_cens > 0]), 2), "\n")
    }
  }
  
  # Filter to uncensored observations
  if (!is.null(cens)) {
    y <- y[keep_idx]
    a <- a[keep_idx]
    x <- x[keep_idx, , drop = FALSE]
    w_cens <- w_cens[keep_idx]
  }
  
  # Normalize weights
  if(sum(w_cens) > 0) {
    w_cens <- w_cens * length(w_cens) / sum(w_cens)
  } else {
    w_cens <- rep(1, length(y))
  }
  
  n <- length(y)
  
  # Set up evaluation points
  a.min <- a.rng[1]
  a.max <- a.rng[2]
  
  # Check for valid range
  if(!is.finite(a.min) || !is.finite(a.max) || a.min >= a.max) {
    stop("Invalid treatment range. Check that there are uncensored observations.")
  }
  
  a.vals <- seq(a.min, a.max, length.out = n.pts)
  
  # Create augmented dataset
  xa.new <- rbind(
    cbind(x, a), 
    cbind(x[rep(1:n, length(a.vals)), ], a = rep(a.vals, rep(n, length(a.vals))))
  )
  x.new <- xa.new[, -ncol(xa.new)]
  x <- data.frame(x)
  x.new <- data.frame(x.new)
  colnames(x) <- colnames(x.new)
  xa.new <- data.frame(xa.new)
  
  # STEP 3: Estimate Nuisance Functions with Weights
  
  # Treatment model
  pimod <- SuperLearner(Y = a, X = x, 
                        obsWeights = w_cens,
                        SL.library = sl.lib, 
                        newX = x.new)
  pimod.vals <- pimod$SL.predict
  
  # Variance model
  pi2mod <- SuperLearner(Y = (a - pimod.vals[1:n])^2, X = x,
                         obsWeights = w_cens,
                         SL.library = sl.lib, 
                         newX = x.new)
  pi2mod.vals <- pmax(pi2mod$SL.predict, 0.01)
  
  # Outcome model
  mumod <- SuperLearner(Y = y, X = cbind(x, a),
                        obsWeights = w_cens,
                        family = binomial(),
                        SL.library = sl.lib, 
                        newX = xa.new)
  muhat.vals <- mumod$SL.predict
  muhat.vals <- pmin(pmax(muhat.vals, 0.001), 0.999)
  
  # STEP 4: Construct Densities and Pseudo-outcome
  
  a.std <- (xa.new$a - pimod.vals) / sqrt(pi2mod.vals)
  
  # Weighted density
  dens_weights <- w_cens / sum(w_cens)
  dens_est <- density(a.std[1:n], weights = dens_weights)
  pihat.vals <- approx(dens_est$x, dens_est$y, xout = a.std, rule = 2)$y / sqrt(pi2mod.vals)
  
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  # Weighted marginals
  varpihat <- predict(smooth.spline(a.vals, 
                                    apply(pihat.mat, 2, weighted.mean, w = w_cens)), 
                      x = a)$y
  varpihat <- pmax(varpihat, 0.001)
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, weighted.mean, w = w_cens), n), 
                         byrow = TRUE, nrow = n)
  
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  mhat <- predict(smooth.spline(a.vals, 
                                apply(muhat.mat, 2, weighted.mean, w = w_cens)), 
                  x = a)$y
  mhat <- pmin(pmax(mhat, 0.001), 0.999)
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, weighted.mean, w = w_cens), n), 
                     byrow = TRUE, nrow = n)
  
  # Pseudo-outcome
  pseudo.out <- (y - muhat) * (varpihat / pihat) + mhat
  
  # STEP 5: Bandwidth Selection
  
  w.fn <- function(bw) {
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a - a.val) / bw
      kern.std <- kern(a.std) / bw
      numerator <- weighted.mean(a.std^2 * kern.std, w = w_cens) * (kern(0) / bw)
      denom1 <- weighted.mean(kern.std, w = w_cens)
      denom2 <- weighted.mean(a.std^2 * kern.std, w = w_cens)
      denom3 <- weighted.mean(a.std * kern.std, w = w_cens)^2
      denominator <- denom1 * denom2 - denom3
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
  
  risk.fn <- function(h) {
    hats <- pmin(hatvals(h), 0.9)
    preds <- cts.eff.fn(pseudo.out, bw = h)
    weighted.mean(((pseudo.out - preds) / (1 - hats))^2, w = w_cens, na.rm = TRUE)
  }
  
  risk.est <- sapply(bw.seq, function(h) {
    tryCatch(risk.fn(h), error = function(e) Inf)
  })
  
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  

  # STEP 6: Final Estimation
  
  final_curve <- locpoly(a, pseudo.out, bandwidth = h.opt)
  est <- approx(final_curve$x, final_curve$y, xout = a.vals, rule = 2)$y
  est <- pmin(pmax(est, 0), 1)
  

  # STEP 7: Standard Errors

  se <- NULL
  for (a.val in a.vals) {
    a.std <- (a - a.val) / h.opt
    kern.std <- kern(a.std) / h.opt
    
    weights_valid <- kern.std > 0.001
    if (sum(weights_valid) > 10) {
      combined_weights <- kern.std * w_cens
      
      tryCatch({
        beta <- coef(lm(pseudo.out ~ a.std, 
                        weights = combined_weights, 
                        subset = weights_valid))
        
        Dh <- matrix(c(
          weighted.mean(kern.std, w = w_cens),
          weighted.mean(kern.std * a.std, w = w_cens),
          weighted.mean(kern.std * a.std, w = w_cens),
          weighted.mean(kern.std * a.std^2, w = w_cens)
        ), nrow = 2)
        
        if(det(Dh) > 0.0001) {
          residuals <- pseudo.out - beta[1] - beta[2] * a.std
          sigma <- weighted.mean(residuals^2, w = combined_weights)
          se <- c(se, sqrt(sigma / n))
        } else {
          se <- c(se, NA)
        }
      }, error = function(e) {
        se <- c(se, NA)
      })
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
  
  return(list(res = res, bw.risk = bw.risk))
}