################################################################################
#
# Helper functions used by other scripts
#
################################################################################

# Takes a UTM sf object in m and converts the units to km
convert_to_km <- function(dat) {
  # Do conversion
  dat_km <- st_geometry(dat)/1000
  st_crs(dat_km) <- st_crs(dat)
  st_geometry(dat) <- dat_km
  
  # Save coordinates in a column
  dat[c("x", "y")] <- st_coordinates(dat)
  
  return(dat)
}

# Generates environmental covariate conditional effects plots with 95% confidence 
# bounds for a model fit... fit is the model object, xvars is a vector containing 
# the names of covariates to plot effects for, xlabs is a vector of labels 
# corresponding to those covariates, ylimits is a list containing a ylim vector 
# for each covariate, and ylab is a label for all y-axes to share... nrow and ncol 
# define how the plot facets are arranged... no_ylab is a vector of + integers 
# (e.g., c(1,2,3)) corresponding to the covariates to drop y-axes labels for... 
# yr is the year to set f_year/year to... exp = T will convert log-transformed 
# covariates back to their natural scale... function will always reverse mean-centering 
# and scaling
plot_effects <- function(fit, xvars, yr = 2023, xlabs, ylab, no_ylab = NULL, ylimits,
                         nrow = NULL, ncol = NULL,
                         text_size = 10, linewidth = 1, ptsize = 1.6, 
                         color_est = "deepskyblue", color_se = "grey",
                         alpha_est = 1, alpha_se = 0.5, exp = T)
{
  # Create a list to store the plots
  plot_list <- list()
  
  # Generate plots and store them in plot_list
  for(i in 1:length(xvars))
  {
    # Set up initial dataframe to expand
    pred_dat <- data.frame(year = yr,
                           snowcrab = TRUE,
                           BtmTempScaled = 0, 
                           BtmTempRangeScaled = 0,
                           ShearVelocityScaled = 0,
                           LogDepthScaled = 0,
                           PhiScaled = 0)
    
    # Expand dataframe using covariate of interest
    xvar_data <- seq(min(fit$data[[xvars[i]]], na.rm = TRUE), 
                     max(fit$data[[xvars[i]]], na.rm = TRUE), by = 0.1)
    pred_dat <- pred_dat[rep(1, length(xvar_data)), ]
    pred_dat[[xvars[i]]] <- xvar_data
    pred_dat$unscaled_covar <- pred_dat[[xvars[i]]]*sd(fit$data[[gsub("Scaled", "", xvars[i])]]) + 
      mean(fit$data[[gsub("Scaled", "", xvars[i])]]) # Reverse mean-centering and scaling
    if(grepl("Log", xvars[i]) && exp == T) # Reverse log if required
    {
      pred_dat$unscaled_covar <- exp(pred_dat$unscaled_covar)
    }
    
    # Make prediction
    prediction <- predict(fit, pred_dat, re_form = NA, se_fit = TRUE)
    
    # Generate plot
    rug_col_name <- if (grepl("Log", xvars[i]) && exp == TRUE) { # Reverse log if required
      gsub("Log|Scaled", "", xvars[i])
    } else {
      gsub("Scaled", "", xvars[i])
    }
    plot_list[[i]] <- ggplot(prediction, aes(x = unscaled_covar, exp(est),
                                             ymin = exp(est - 1.96 * est_se),
                                             ymax = exp(est + 1.96 * est_se))) +
      geom_rug(data = fit$data, inherit.aes = F, 
               aes(x = .data[[rug_col_name]]), sides = "b", alpha = 0.5) +
      geom_ribbon(fill = color_se, col = NA, alpha = alpha_se) + 
      geom_line(lwd = linewidth, col = color_est, alpha = alpha_est) +
      coord_cartesian(expand = FALSE, ylim = ylimits[[i]]) + 
      theme_classic() +
      theme(text = element_text(size = text_size), 
            axis.text = element_text(size = text_size)) +
      labs(x = xlabs[i], y = ylab)
    # If plot doesn't have a ylabel, remove
    if(i %in% c(no_ylab))
    {
      plot_list[[i]] <- plot_list[[i]] + theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank())
    }
  }
  
  # Return the final combined plot
  wrap_plots(plot_list, ncol = ncol, nrow = nrow) + 
    plot_layout(axis_titles = "collect")
}