#' plot.tvem:  Produce plots for a tvem model.
#' 
#' Produces plots from a tvem object produced by 
#' the tvem function.  These plots will be shown on the default
#' output device (likely the screen);  they can of course be 
#' written to a file instead, by preceding the call to plot.tvem 
#' with a call to png(), pdf(), or other R graphic file output functions.
#' 
#' @param x The TVEM object to be plotted.
#' @param use_panes Whether to plot multiple coefficient
#' functions in a single image.
#' @param which_plot The coefficient number to plot,
#' if only one plot is desired at a time. 
#' @param diagnostics If this is set to TRUE, then 
#' instead of plotting coefficient functions,
#' the function will show a histogram of residuals 
#' and a plot of fitted values versus
#' residuals.  These may be useful in checking for 
#' outliers or skew in TVEM with a numeric
#' outcome.  They are not likely to be as useful 
#' in TVEM with a binary or other discrete outcome.
#' @param exponentiate If this is set to TRUE and if
#' the TVEM had a binary outcome, then the exponentiated 
#' coefficient functions (representing odds and odds
#' ratios) will be plotted rather than the usual
#' coefficient functions (representing log odds and log
#' odds ratios).
#' @param ... Further arguments currently not supported
#' 
#' @export
#' @method plot tvem

plot.tvem <- function(x,
                      use_panes=TRUE,
                      which_plot=NULL,
                      diagnostics=FALSE, 
                      exponentiate=FALSE, ...) {
  old_par <- par(no.readonly = TRUE);
  on.exit(par(old_par));
  num_tv_coefs <- length(x$grid_fitted_coefficients);
  if (use_panes) {
    if (diagnostics) {
      panel_dims <- c(1,2);
    } else {
      if (num_tv_coefs==1) {panel_dims <- c(1,1);}
      if (num_tv_coefs==2) {panel_dims <- c(1,2);}
      if (num_tv_coefs==3 | num_tv_coefs==4) {panel_dims <- c(2,2);}
      if (num_tv_coefs==5 | num_tv_coefs==6) {panel_dims <- c(3,2);}
      if (num_tv_coefs==7 | num_tv_coefs==8) {panel_dims <- c(4,2);}
      if (num_tv_coefs==9) {panel_dims <- c(3,3);}
      if(num_tv_coefs>10) {stop("Too many functions to plot in panes.");}  
    }
    par(mfrow=panel_dims);
  };
  if (diagnostics) {
    if (exponentiate) {stop("Error:  This function cannot currently provide diagnostic plots on the odds ratio scale.");}
    hist(x$back_end_model$fitted,
         main="Residuals",
         xlab="Residual");
    plot(x$back_end_model$fitted,
         x$back_end_model$residuals,
         main="Fitted versus residuals",
         xlab="Fitted",
         ylab="Residuals");
  } else {
    the_grid <- x$time_grid;
    temp_plot_function <- function(the_var_name,
                                   the_grid,
                                   the_coef,
                                   ymin,
                                   ymax,
                                   exponentiate) {
      if (ymax - ymin < 1e-3) {
         warning(paste("At least one of the coefficient functions has been",
                       "estimated as essentially equal to zero across the,",
                       "interval ('shrunken to zero').  It should be treated",
                       "as removed from the model, and its confidence intervals",
                       "should be ignored."));
      }
      if (exponentiate) {
        plot(x=the_grid,
             y=exp(the_coef$estimate),
             main=paste("Exponentiated TVEM coefficient:\n",the_var_name),
             xlab=expression(t),
             ylab=expression(exp(beta(t))),
             ylim=c(exp(ymin),exp(ymax)),
             type="l",
             lty="solid",
             lwd=2, 
             mgp=c(2,1,0));
        abline(h=1,lty="dotted");
        lines(x=the_grid,
              y=exp(the_coef$lower));
        lines(x=the_grid,
              y=exp(the_coef$upper));
        
      } else {
        plot(x=the_grid,
             y=the_coef$estimate,
             main=paste("TVEM coefficient:\n",the_var_name),
             xlab=expression(t),
             ylab=expression(beta(t)),
             ylim=c(ymin,ymax),
             type="l",
             lty="solid",
             lwd=2, 
             mgp=c(2,1,0));
        abline(h=0,lty="dotted");
        lines(x=the_grid,
              y=the_coef$lower);
        lines(x=the_grid,
              y=the_coef$upper);
        
      }
      return(0);
    }
    if (is.null(which_plot)) {
      for (which_plot in 1:num_tv_coefs) {
        the_coef <- x$grid_fitted_coefficients[[which_plot]];
        the_var_name <- names(x$grid_fitted_coefficients)[which_plot];
        ymin <- min(0,min(the_coef$lower));
        ymax <- max(0,max(the_coef$upper));
        temp_plot_function(the_var_name,
                           the_grid,
                           the_coef,
                           ymin,
                           ymax,
                           exponentiate);
      }
    } else {
      the_coef <- x$grid_fitted_coefficients[[which_plot]];
      the_var_name <- names(x$grid_fitted_coefficients)[which_plot];
      ymin <- min(0,min(the_coef$lower));
      ymax <- max(0,max(the_coef$upper));
      temp_plot_function(the_var_name,
                         the_grid,
                         the_coef,
                         ymin,
                         ymax,
                         exponentiate);
    }
  }
}