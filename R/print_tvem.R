#' print.tvem:  Print output from a model that was fit by the tvem function.
#' 
#' @param x The tvem object (output of the tvem or select_tvem function)
#' @param ornate Whether to print lines between different sections of the output for easier reading.
#' @param ... Further arguments currently not supported
#' 
#' @export
#' @method print tvem




print.tvem <- function(x, ornate=TRUE, ...) {
  if (ornate) {
    divider <- "======================================================= \n";
  } else {
    divider <- "\n";
  }
  if (ornate) {
    cat(divider);
    cat("Time-Varying Effects Modeling (TVEM) Function Output \n");
    cat(divider);
  }
  cat(paste("Response variable:  ",
            x$model_information$response_name,
            "\n"));
  cat(paste("Response outcome distribution type:",
            x$model_information$outcome_family,
            "\n"));
  cat(paste("Time interval:  ",
            round(min(x$time_grid),4),
            "to",
            round(max(x$time_grid),4),"\n"));
  cat("Number of subjects:  ");
  cat(paste(x$model_information$n_subjects));
  cat("\nEffects specified as time-varying:  ");
  cat(paste(names(x$grid_fitted_coefficients),sep=" ",collapse=", ")); 
  cat("\nYou can use the plot function to view the plots.\n");
  if (!is.null(x$invar_effects_estimates)) {
    cat(divider)
    cat("Effects specified as non-time-varying: \n");
    print(x$invar_effects_estimates);
  }
  cat(divider);
  cat("Back-end model fitted in mgcv::bam function: \n")
  cat(paste("Method",x$back_end_model$method))
  cat("\nFormula:\n");
  print(x$back_end_model$formula); 
  cat(paste("Pseudolikelihood AIC:",round(x$model_information$pseudo_aic,2)));
  cat(paste("\nPseudolikelihood BIC:",round(x$model_information$pseudo_bic,2),"\n"));
  if (x$model_information$used_listwise_deletion) {
    cat("Note: Used listwise deletion for missing data.\n");
  }
  if (!is.null(x$ICs_table)) {
    cat("Model selection table for number of interior knots:\n")
    print(x$ICs_table);
  }
  cat(divider);
}
