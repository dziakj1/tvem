select_tvem <- function(max_knots=5,
                        keep_going_if_too_few=FALSE,
                        use_bic=FALSE,
                        penalize=FALSE,
                        print_output=TRUE,
                        ...) { 
  #' select_tvem:  Select number of interior knots for an unpenalized TVEM.
  #' @param max_knots  The maximum number of interior knots to try (0 through max_knots)
  #' @param keep_going_if_too_few  Whether to continue in a stepwise fashion if the max_knots 
  #' does not seem to be high enough
  #' @param use_bic Whether to use BIC (TRUE) instead of AIC (FALSE) when selecting the best 
  #' model. Note that both of these IC's are calculated from the working-independence 
  #' pseudolikelihood rather than the unknown true likelihood. However, for BIC, the
  #' sample size is taken to be the number of subjects, not the number of observations.
  #' @param penalize Whether to include a penalty function in estimation
  #' @param print_output Whether to print the pseudolikelihoods obtained for each candidate number of interior knots.
  #' @param  ... Other inputs to be sent along to each call to the tvem function. 
  #' 
  #' @return A TVEM object for the fitted model, with an additional component containing
  #' a table of information criteria.
  #' @export
  
  args1 <- match.call();  
  num_knots_values <- 0:max_knots;  # will try at least each of these values for num_knots;
  IC_values <- rep(NA,length(num_knots_values)); 
  done_looking <- FALSE;
  num_knots_values_index <- 0;
  while (done_looking==FALSE) {
    num_knots_values_index <- num_knots_values_index + 1;
    this_num_knots <- num_knots_values[num_knots_values_index]; 
    more_args <- as.list(args1)[-1];
    more_args$num_knots <- this_num_knots;
    more_args$use_naive_se <- TRUE;
    more_args$print_gam_formula <- FALSE; 
    more_args <- more_args[ (names(more_args)!="max_knots")  &
                            (names(more_args)!="keep_going_if_too_few") & 
                            (names(more_args)!="use_bic")];
    # I got this trick from https://statisticsglobe.com/remove-element-from-list-in-r ;
    #    ans1 <- try(suppressWarnings(do.call(tvem, more_args)));
    ans1 <- do.call(tvem, more_args); 
    if (class(ans1)=="try-error") {
      IC_values[num_knots_values_index] <- Inf;
    } 
    if (use_bic==FALSE) {
      # use AIC;
      IC_values[num_knots_values_index] <- ans1$model_information$pseudo_aic;
    }
    if (use_bic==TRUE) {
      # use BIC;
      IC_values[num_knots_values_index] <- ans1$model_information$pseudo_bic;
    }
    if (num_knots_values_index==length(num_knots_values)) {
      # If you have come to end of list of values to try for number of knots
      if ((which.min(IC_values)==num_knots_values_index) &
          (keep_going_if_too_few==TRUE)) {
        # If there still might be more knots needed then try more
        num_knots_values <- c(num_knots_values, max(num_knots_values)+1);
        IC_values <- c(IC_values, NA);
      }  else {
        done_looking <- TRUE;
      }
    } 
  } 
  best_num_knots <- num_knots_values[which.min(IC_values)];  
  more_args <- as.list(args1)[-1];
  more_args$num_knots <- best_num_knots;
  more_args$use_naive_se <- FALSE; 
  more_args <- more_args[ (names(more_args)!="max_knots")  &
                          (names(more_args)!="keep_going_if_too_few") & 
                          (names(more_args)!="use_bic")];
  ans1 <- do.call(tvem, more_args);
  ans1$ICs_table <- cbind(knots=num_knots_values, ic=IC_values);
  if (print_output) {
    print(ans1$ICs_table);
    print(paste("Selected",best_num_knots,"interior knots."));
  }
  return(ans1);
}