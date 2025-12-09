# Plotting function with density plot
plot_dose_response <- function(all_results, df, 
                                        xlim = c(0, 500), 
                                        ylim = c(0, 0.2),
                                        x_tick_interval = NULL,  # Auto-calculate if NULL
                                        y_tick_interval = NULL,  # Auto-calculate if NULL
                                        grid_x_interval = NULL,  # Auto-calculate if NULL
                                        grid_y_interval = NULL,  # Auto-calculate if NULL
                                        show_percent = TRUE,     # Show % on y-axis
                                        x_label = "Cumulative Treatment (A)",
                                        y_label = "Probability of OUD",
                                        main_title = "Period",
                                        show_density = TRUE,     # Show density plot
                                        density_height = 0.15,   # Proportion of y-range for density
                                        density_color = "#666666",
                                        density_alpha = 0.3,
                                        colors = list(main = "#2E86AB",
                                                      ci = "#2E86AB",
                                                      grid = "#E0E0E0",
                                                      rug = "#666666")) {
  
  # Auto-calculate intervals if not provided
  if(is.null(x_tick_interval)) {
    x_range <- xlim[2] - xlim[1]
    x_tick_interval <- round(x_range / 6, -floor(log10(x_range / 6)))
  }
  
  if(is.null(y_tick_interval)) {
    y_range <- ylim[2] - ylim[1]
    y_tick_interval <- round(y_range / 5, 2)
  }
  
  if(is.null(grid_x_interval)) {
    grid_x_interval <- x_tick_interval
  }
  
  if(is.null(grid_y_interval)) {
    grid_y_interval <- y_tick_interval
  }
  
  # Set up better margins and parameters
  par(mfrow = c(2, 2), 
      mar = c(4.5, 4.5, 3, 2),  
      mgp = c(3, 0.7, 0),       
      las = 1,                   
      family = "sans")           
  
  for(period in 1:4) {
    res <- all_results[[period]]$res
    
    # Subset to specified range
    plot_idx <- res$a.vals >= xlim[1] & res$a.vals <= xlim[2]
    
    # Count observations
    cens_var <- paste0("cens_oud_period_", period)
    a_values <- df$A[df[[cens_var]] == 1]
    a_values_in_range <- a_values[a_values >= xlim[1] & a_values <= xlim[2]]
    n_in_range <- length(a_values_in_range)
    n_total <- sum(df[[cens_var]] == 1)
    pct_shown <- round(100 * n_in_range / n_total)
    
    # Create base plot
    plot(res$a.vals[plot_idx], res$est[plot_idx], 
         type = "n",  
         ylim = ylim,
         xlim = xlim,
         xlab = "",
         ylab = "",
         main = "",
         axes = FALSE,
         frame.plot = FALSE)
    
    # Add custom grid
    y_grid_seq <- seq(ylim[1], ylim[2], by = grid_y_interval)
    abline(h = y_grid_seq, col = colors$grid, lty = 1, lwd = 0.5)
    
    x_grid_seq <- seq(xlim[1], xlim[2], by = grid_x_interval)
    abline(v = x_grid_seq, col = colors$grid, lty = 1, lwd = 0.5)
    
    # Add density plot if requested
    if(show_density && length(a_values_in_range) > 10) {
      # Calculate density
      dens <- density(a_values_in_range, from = xlim[1], to = xlim[2], n = 512)
      
      # Scale density to fit in bottom portion of plot
      dens_scaled <- dens$y / max(dens$y) * (ylim[2] - ylim[1]) * density_height
      
      # Create polygon for density
      polygon(c(dens$x, rev(dens$x)), 
              c(rep(ylim[1], length(dens$x)), rev(dens_scaled + ylim[1])),
              col = adjustcolor(density_color, alpha.f = density_alpha),
              border = NA)
      
      # Add a subtle line at the top of density
      lines(dens$x, dens_scaled + ylim[1], 
            col = adjustcolor(density_color, alpha.f = 0.5), 
            lwd = 0.5)
    }
    
    # Add confidence band
    xx <- c(res$a.vals[plot_idx], rev(res$a.vals[plot_idx]))
    yy <- c(pmax(res$ci.ll[plot_idx], ylim[1]), 
            rev(pmin(res$ci.ul[plot_idx], ylim[2])))
    polygon(xx, yy, col = adjustcolor(colors$ci, alpha.f = 0.2), border = NA)
    
    # Add main line
    lines(res$a.vals[plot_idx], res$est[plot_idx], 
          lwd = 2.5, col = colors$main)
    
    # Add horizontal reference line at 0
    abline(h = 0, lty = 1, col = "black", lwd = 0.5)
    
    # Custom axes
    x_tick_seq <- seq(xlim[1], xlim[2], by = x_tick_interval)
    axis(1, at = x_tick_seq, 
         labels = x_tick_seq,
         cex.axis = 0.9, col = "black", col.axis = "black", lwd = 0.5)
    
    # Y-axis with optional percentage labels
    y_tick_seq <- seq(ylim[1], ylim[2], by = y_tick_interval)
    if(show_percent) {
      y_labels <- paste0(round(y_tick_seq * 100, 1), "%")
    } else {
      y_labels <- round(y_tick_seq, 3)
    }
    axis(2, at = y_tick_seq, 
         labels = y_labels,
         cex.axis = 0.9, col = "black", col.axis = "black", lwd = 0.5)
    
    # Add axis labels
    mtext(x_label, side = 1, line = 3, cex = 0.95)
    mtext(y_label, side = 2, line = 3, cex = 0.95, las = 0)
    
    # Add title with sample size info
    mtext(paste0(main_title, " ", period), side = 3, line = 1.5, cex = 1.1, font = 2)
    mtext(paste0("n = ", format(n_in_range, big.mark = ","), 
                 " (", pct_shown, "% of uncensored)"), 
          side = 3, line = 0.3, cex = 0.8, col = "gray40")
    
    # Add subtle rug plot for additional detail
    rug(a_values_in_range, col = adjustcolor(colors$rug, alpha.f = 0.2), 
        lwd = 0.3, ticksize = 0.01)
    
    # Add box around plot
    box(lwd = 0.5)
  }
}
