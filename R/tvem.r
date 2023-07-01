#' tvem:  Fit a time-varying effect model.
#' 
#' Fits a time-varying effect model (Tan et al., 2012);  that is,
#' a varying-coefficients model (Hastie & Tibshirani, 1993) for 
#' longitudinal data.
#' 
#' @note The interface is based somewhat on the TVEM 3.1.1 SAS macro by 
#' the Methodology Center (Li et al., 2017).  However, that macro uses 
#' either "P-splines" (penalized truncated power splines) or "B-splines" 
#' (unpenalized B[asic]-splines, like those of Eilers and Marx, 1996, 
#' but without the smoothing penalty).  The current function uses 
#' penalized B-splines, much more like those of Eilers and Marx (1996).
#' However, their use is more like the "P-spline" method than the "B-spline" method   
#' in the TVEM 3.1.1 SAS macro, in that the precise choice of knots 
#' is not critical, the tuning is done automatically, and the fitted model
#' is intended to be interpreted in a population-averaged (i.e., marginal) 
#' way.  Thus, random effects are not allowed, but sandwich standard 
#' errors are used in attempt to account for within-subject correlation,
#' similar to working-independence GEE (Liang and Zeger, 1986).
#' 
#' @note Note that as in ordinary parametric regression, if the range
#' of the covariate does not include values near zero, then the 
#' interpretation of the intercept coefficient may be somewhat 
#' difficult and its standard errors may be large (i.e., due to extrapolation).
#' 
#' @note The bam ("Big Additive Models") function in the 
#' mgcv package ("Mixed GAM Computation Vehicle with GCV/AIC/REML smoothness
#' estimation and GAMMs by REML/PQL") by Simon Wood is used for back-end
#' calculations (see Wood, Goude, & Shaw, 2015).
#' 
#' @references 
#' Eilers, P. H. C., & Marx, B.  D. (1996). Flexible smoothing with B-splines
#' and penalties. Statistical Science, 11: 89-121. <doi:10.1214/ss/1038425655>
#' @references
#' Hastie, T, Tibshirani, R. (1993). Varying-coefficient models. Journal 
#' of the Royal Statistical Socety, B, 55:757-796. <doi:10.1057/9780230280830_39>
#' @references
#' Li, R., Dziak, J. J., Tan, X., Huang, L., Wagner, A. T., & Yang, J. (2017).
#' TVEM (time-varying effect model) SAS macro users' guide (Version 3.1.1).
#' University Park: The Methodology Center, Penn State. Retrieved from
#' <http://methodology.psu.edu>. Available online at
#' <https://aimlab.psu.edu/tvem/tvem-sas-macro/> and archived at  
#' <https://github.com/dziakj1/MethodologyCenterTVEMmacros>
#' and
#' <https://scholarsphere.psu.edu/collections/v41687m23q>.
#' @references
#' Liang, K. Y., Zeger, S. L. Longitudinal data analysis using generalized linear
#' models. Biometrika. 1986; 73:13-22. <doi:10.1093/biomet/73.1.13>
#' @references
#' Tan, X., Shiyko, M. P., Li, R., Li, Y., & Dierker, L. (2012). A time-varying
#' effect model for intensive longitudinal data. Psychological Methods, 17: 
#' 61-77. <doi:10.1037/a0025814>
#' @references
#' Wood, S. N., Goude, Y., & Shaw, S. (2015). Generalized additive models 
#' for large data sets. Applied Statistics, 64: 139-155. ISBN 10 1498728332, 
#' ISBN 13 978-1498728331.
#' 
#' @param data The dataset containing the observations, assumed
#' to be in long form (i.e., one row per observation, potentially
#' multiple rows per subject).
#' @param formula A formula listing the outcome on the left side, 
#' and the time-varying effects on the right-side.  For a time-
#' varying intercept only, use y~1, where y is the name of the 
#' outcome.  For a single time-varying-effects covariate, use
#' y~x, where x is the name of the covariate.  For multiple 
#' covariates, use syntax like y~x1+x2.  Do not include the non-
#' time-varying-effects covariates here.  Note that the values 
#' of these covariates themselves may either be time-varying 
#' or time-invariant. For example, time-invariant biological sex may have
#' a time-varying effect on time-varying height during childhood.
#' @param id The name of the variable in the dataset which represents 
#' subject (participant) identity.  Observations are considered 
#' to be correlated within subject (although the correlation 
#' structure is not explicitly modeled) but are assumed independent
#' between subjects.
#' @param time The name of the variable in the dataset which 
#' represents time. The regression coefficient functions representing 
#' the time-varying effects are assumed to be smooth functions 
#' of this variable.  
#' @param invar_effects  Optionally, the names of one or more 
#' variables in the dataset assumed to have a non-time-varying (i.e., time-invariant)
#' regression effect on the outcome.  The values of these covariates 
#' themselves may either be time-varying or time-invariant.  The
#' covariates should be specified as the right side of a formula, e.g.,
#' ~x1 or ~x1+x2.
#' @param family  The outcome family, as specified in functions like 
#' glm.  For a numerical outcome you can use the default of gaussian(). 
#' For a binary outcome, use binomial().  For a count outcome, you can 
#' use poisson().  The parentheses after the family name are there because it is 
#' actually a built-in R object.
#' @param weights An optional sampling weight variable.
#' @param num_knots The number of interior knots assumed per spline function,
#' not counting exterior knots. This is assumed to be the same for each function.
#' If penalized=TRUE is used, it is probably okay to leave num_knots at its default.
#' @param spline_order The shape of the function between knots, with a
#' default of 3 representing cubic spline.
#' @param penalty_function_order The order of the penalty function (see 
#' Eilers and Marx, 1996), with a default of 1 for first-order 
#' difference penalty.  Eilers and Marx (1996) used second-order difference
#' but we found first-order seemed to perform parsimoniously in this setting.
#' Please feel free to consider setting this to 2 to explore other possible 
#' results. The penalty function is something analogous to a prior distribution
#' describing how smooth or flat the estimated coefficient functions should be,  
#' with 1 being smoothest.
#' @param grid The number of points at which the spline coefficients
#' will be estimated, for the purposes of the pointwise estimates and 
#' pointwise standard errors to be included in the output object. The
#' grid points will be generated as equally spaced over the observed 
#' interval. Alternatively, grid can be specified as a vector instead, in which 
#' each number in the vector is interpreted as a time point for the grid
#' itself.
#' @param penalize Whether to add a complexity penalty; TRUE or FALSE
#' @param alpha  One minus the nominal coverage for 
#' the pointwise confidence intervals to be constructed.  Note that a 
#' multiple comparisons correction is not applied.  Also, in some cases
#' the nominal coverage may not be exactly achieved even pointwise, 
#' because of uncertainty in the tuning parameter and risk of overfitting.
#' These problems are not unique to TVEM but are found in many curve-
#' fitting situations.  
#' @param basis Form of function basis (an optional argument about computational 
#' details passed on to the mgcv::s function as bs=).  We strongly recommend
#' leaving it at the default value.
#' @param method Fitting method (an optional argument about computational 
#' details passed on to the mgcv::bam function as method).  We strongly recommend
#' leaving it at the default value.
#' @param use_naive_se  Whether to save time by using a simpler, less valid formula for 
#' standard errors. Only do this if you are doing TVEM inside a loop for
#' bootstrapping or model selection and plan to ignore these standard errors.
#' @param print_gam_formula  whether to print the formula used to do the back-end calculations
#' in the bam (large data gam) function in the mgcv package.
#' @param normalize_weights Whether to rescale (standardize) the weights variable to have a mean 
#' of 1 for the dataset used in the analysis. Setting this to FALSE might lead to invalid 
#' standard errors caused by misrepresentation of the true sample size. This 
#' option is irrelevant and ignored if a weight variable is not specified, because
#' in that case all the weights are effectively 1 anyway.  An error will result if the function is asked to rescale weights and any of the weights are negative; however, it is very rare for
#' sampling weights to be negative.
#'
#' @return An object of type tvem. The components of an object of 
#' type tvem are as follows:
#' \describe{
#' \item{time_grid}{A vector containing many evenly spaced time 
#' points along the interval between the lowest and highest observed 
#' time value. The exact number of points is determined by the 
#' input parameter 'grid'.}
#' \item{grid_fitted_coefficients}{A list of data frames, one for 
#' each smooth function which was fit (including the intercept). Each 
#' data frame contains the fitted estimates of the function
#' at each point of time_grid, along with pointwise standard 
#' errors and pointwise confidence intervals.}
#' \item{invar_effects_estimates}{If any variables are specified in
#' invar_effects, their estimated regression coefficients and 
#' standard errors are shown here.}
#' \item{model_information}{A list summarizing the options
#' specified in the call to the function, as well as fit statistics
#' based on the log-pseudo-likelihood function. The term pseudo
#' here means that the likelihood function is evaluated as though
#' the correct knot locations were known, as though the 
#' observations were independent and, if applicable, as though sampling 
#' weights were multiples of a participant rather than 
#' inverse probabilities. This allows tvem to be used without
#' specifying a fully parametric probability model.}
#' \item{back_end_model}{The full output from the bam() 
#' function from the mgcv package, which was used to fit the 
#' penalized spline regression model underlying the TVEM.}
#' }  
#' 
#' @keywords Statistics|smooth
#' @keywords Statistics|models|regression
#'
#'@import mgcv
#'@importFrom graphics abline hist lines par plot text
#'@importFrom stats AIC as.formula binomial coef
#'           gaussian glm plogis poisson qnorm rbinom
#'           rnorm sd terms update var
#' 
#'@examples
#' set.seed(123)
#' the_data <- simulate_tvem_example()
#' tvem_model <- tvem(data=the_data,
#'               formula=y~x1,
#'               invar_effects=~x2,
#'               id=subject_id,
#'               time=time)
#' print(tvem_model)
#' plot(tvem_model)
#' 
#'@export

tvem <- function(data,
                 formula,
                 id,
                 time,
                 invar_effects=NULL,
                 family=gaussian(), # use binomial() for binary;
                 weights=NULL,
                 num_knots=20,  
                 # number of interior knots per function (not counting exterior knots);
                 spline_order=3,
                 # shape of function between knots, with default of 3 for cubic;
                 penalty_function_order=1,
                 grid=100,
                 penalize=TRUE,
                 alpha=.05,
                 basis="ps",
                 method="fREML",
                 use_naive_se=FALSE,
                 print_gam_formula=FALSE,
                 normalize_weights=TRUE) 
{   
  ##################################
  # Process the input;
  ################################## 
  m <- match.call(expand.dots = FALSE);
  weights_argument_number <- match("weights",names(m));
  if (is.na(weights_argument_number)) {
    use_weights <- FALSE;
  } else {
    use_weights <- TRUE;
    weights_variable_name <- as.character(m[[weights_argument_number]]);
  } 
  m$formula <- NULL;
  m$invar_effects <- NULL;
  m$family <- NULL;
  m$grid <- NULL;
  m$num_knots <- NULL;
  m$spline_order <- NULL;
  m$penalty_function_order <- NULL;
  m$alpha <- NULL;
  m$penalize <- NULL;
  m$basis <- NULL;
  m$method <- NULL;
  m$use_naive_se <- NULL;
  m$print_gam_formula <- NULL;
  m$normalize_weights <- NULL;
  if (is.matrix(eval.parent(m$data))) {
    m$data <- as.data.frame(data);
  }
  m[[1]] <- quote(stats::model.frame);
  m <- eval.parent(m); 
  id_variable_name <- as.character(substitute(id));
  time_variable_name <- as.character(substitute(time)); 
  if (is.character(family)) {
    # Handle the possibility that the user specified family
    # as a string instead of an object of class family.
    family <- tolower(family);
    if (family=="normal" | 
        family=="gaussian" | 
        family=="linear" | 
        family=="numerical") {
      family <- gaussian();
    }
    if (family=="binary" | 
        family=="bernoulli" | 
        family=="logistic" | 
        family=="binomial") {
      family <- binomial();
    }
    if (family=="count" | 
        family=="poisson"|
        family=="loglinear") {
      family <- poisson();
    }
  }
  if ((family$family=="gaussian") &
      (family$link!="identity")) {
    stop("Currently the tvem function only supports the identity link for Gaussian outcomes.")
  }
  if ((family$family=="binomial") &
      (family$link!="logit")) {
    stop("Currently the tvem function only supports the logit link for binomial outcomes.")
  }
  if ((family$family=="poisson") &
      (family$link!="log")) {
    stop("Currently the tvem function only supports the log link for Poisson outcomes.")
  }
  if ((family$family!="gaussian")&
      (family$family!="binomial")&
      (family$family!="poisson")) {
    stop("This version of tvem only handles the gaussian(), binomial() and poisson() families.")
  }
  orig_formula <- formula; 
  the_terms <- terms(formula);
  whether_intercept <- attr(the_terms,"intercept");
  if (whether_intercept != 1) {
    stop("An intercept function is required in the current version of this function")
  }
  formula_variable_names <- all.vars(formula); 
  if (is.null(invar_effects)) {
    num_invar_effects <- 0;
    invar_effects_names <- NA;
  } else {
    invar_effects_names <- attr(terms(invar_effects),"term.labels");
    num_invar_effects <- length(invar_effects_names);
  }
  data_for_analysis <- as.data.frame(eval.parent(data)); 
  used_listwise_deletion <- FALSE;
  names_to_check <- c(id_variable_name,
                      time_variable_name,
                      formula_variable_names);
  if (num_invar_effects>0) {
    names_to_check <- c(names_to_check,
                        invar_effects_names);
  }
  if (use_weights) {
    names_to_check <- c(names_to_check, weights_variable_name);
  }
  for (variable_name in names_to_check) {
    if (sum(is.na(data_for_analysis[,variable_name]))>0) {
      data_for_analysis <- data_for_analysis[which(!is.na(data_for_analysis[,variable_name])),]; 
      used_listwise_deletion <- TRUE;
    }
  } 
  id_variable <- data_for_analysis[,id_variable_name];
  time_variable <- data_for_analysis[,time_variable_name];
  response_name <- formula_variable_names[[1]];
  if (length(formula_variable_names)<=1) {
    num_varying_effects <- 0;
  } else {
    num_varying_effects <- length(formula_variable_names)-1; # not including intercept;
  }
  if (num_varying_effects>0) {
    varying_effects_names <- c(formula_variable_names[2:(1+num_varying_effects)]);
  } else {
    varying_effects_names <- NA;
  }
  if (length(num_knots)>1) {
    stop(paste("Please provide a single number for num_knots."));
  };
  crit_value <- qnorm(1-alpha/2);
  if (use_weights) {
    if (max(abs(data_for_analysis))<1e-8) {
      stop("The weights must not all be zero.");
    }
    if (normalize_weights) {
      if (min(data_for_analysis[,weights_variable_name])<0) {
        stop("Weights cannot be rescaled if any of them are negative."); 
      } 
      data_for_analysis$weights_for_analysis <- data_for_analysis[,weights_variable_name] / mean(data_for_analysis[,weights_variable_name]);
    } else {
      data_for_analysis$weights_for_analysis <- data_for_analysis[,weights_variable_name];
    }
  };
  if (num_knots + spline_order + 1 > length(unique(time_variable))) {
    stop(paste("Because of the limited number of unique time points, \n",
               "please provide a smaller value for num_knots."));
  }
  ####################################################################
  # Construct regular grid for plotting fitted coefficient functions;
  ####################################################################
  if (is.null(grid)) {
    grid <- 100;
  }
  if (length(grid)==1) {
    grid_size <- as.integer(grid);
    min_time <- min(time_variable, na.rm = TRUE);
    max_time <- max(time_variable, na.rm = TRUE);
    time_grid <- seq(min_time, max_time, length=grid_size);
  }
  if (length(grid)>1) {
    time_grid <- grid;
  }
  ##########################################################################
  # Construct formula to send in to the back-end computation function.
  # This function is bam (big generalized additive models) in the  
  # mgcv (Mixed GAM Computation Vehicle) package by Simon Wood of R Project.
  ##########################################################################
  if (num_invar_effects>0 | num_varying_effects>0) {
    bam_formula <- as.formula(paste(response_name," ~ ", "1"));
  } else {
    bam_formula <- orig_formula;
  }   
  if (num_invar_effects>0 & num_varying_effects>0) {
    if (length(intersect(invar_effects_names,
                         varying_effects_names))>0) {
      stop(paste("Please do not specify the same variable",
                 "as having time-varying and time-invariant effects."));
    }
  }
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      new_text <- paste("~ . +",varying_effects_names[i]); 
      bam_formula <- update(bam_formula,as.formula(new_text));
      # add linear effect of each time-varying-effect to the formula
    }
  }
  if (num_invar_effects>0) {
    for (i in 1:num_invar_effects) {
      new_text <- paste("~ . +",invar_effects_names[i]);
      bam_formula <- update(bam_formula,as.formula(new_text));
      # add effect of each non-time-varying-effect to the formula
    }
  }
  new_text <- paste("~ . + s(",time_variable_name,",bs='",basis,"',by=NA,pc=0,",
                    "k=",
                    num_knots+spline_order+1,",fx=",
                    ifelse(penalize,"FALSE","TRUE"),")",sep=""); 
  # for time-varying intercept; 
  bam_formula <- update(bam_formula,as.formula(new_text)); 
  # for time-varying intercept;
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      this_covariate_name <- varying_effects_names[i];
      new_text <- paste("~ . + s(",time_variable_name,",bs='",basis,"',by=",
                        this_covariate_name,
                        ", pc=0,",
                        "m=c(",
                        spline_order-1,
                        ",",
                        penalty_function_order,"),",
                        "k=",
                        num_knots+spline_order+1,",fx=",ifelse(penalize,"FALSE","TRUE"),
                        ")",sep=""); 
      bam_formula <- update(bam_formula,as.formula(new_text));
    }
  } 
  if (family$family=="binomial") {
    if (max(data_for_analysis[,response_name],na.rm=TRUE)>1) {
      stop(paste("In this version of tvem, binomial data must be specified",
                 "as 0s and 1s; multinomial data is not allowed."));
    }
  }
  if (family$family=="poisson") {
    if (max(abs(data_for_analysis[,response_name]-round(data_for_analysis[,response_name])),na.rm=TRUE)>1e-10) {
      stop("Please provide only integer (whole number) count data for the outcome.");
    } else {
      data_for_analysis[,response_name] <- round(data_for_analysis[,response_name]);
    }
    if (min(data_for_analysis[,response_name],na.rm=TRUE)<0) {
      stop(paste("Count data must be nonnegative."));
    }
  }
  if (print_gam_formula) {print(bam_formula);}
  if (use_weights) {
    weights_for_analysis <- NULL; # this is a clumsy workaround
    # to tell the R syntax checker that weights_for_analysis 
    # is not an undeclared object (but is actually a column of 
    # data_for_analysis);
    model1 <- mgcv::bam(bam_formula,
                        data=data_for_analysis,
                        family=family,
                        weights=weights_for_analysis,
                        method=method);
  } else {
    model1 <- mgcv::bam(bam_formula,
                        data=data_for_analysis,
                        family=family,
                        method=method);
  }
  within_subject_variance <- sapply(X=unique(data_for_analysis[,id_variable_name]),
                                    function(X){var(data_for_analysis[which(data_for_analysis[,id_variable_name]==X),response_name],na.rm=TRUE)});
  if (length(within_subject_variance)>0) {
    if (sum(!is.na(within_subject_variance))>0) {
      highest_within_subject_variance_not_counting_singletons <- 
        max(within_subject_variance,na.rm=TRUE)
      if (highest_within_subject_variance_not_counting_singletons<1e-10) {
        warning(paste("The variable specified as the output seems to be",
                      "time-invariant within subject \n despite multiple measurements",
                      "per subject.  Results may not be interpretable."));
      }
    } 
  }
  ##################################
  # Extract coefficient estimates;
  ##################################
  grand_intercept <- model1$coefficients["(Intercept)"];
  bam_coefs <- predict.bam(model1,type="terms"); 
  estimated_b0 <- grand_intercept + 
    bam_coefs[,paste("s(",time_variable_name,")",sep="")];
  estimated_b <- NULL;
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      this_covariate_name <- varying_effects_names[i];
      this_covariate_values <- model1$model[,this_covariate_name];
      this_spline_term <- paste("s(",time_variable_name,"):",
                                this_covariate_name,sep=""); 
      this_covariate_function_values <-  (bam_coefs[,this_covariate_name]+
                                            bam_coefs[,this_spline_term])/this_covariate_values;
      estimated_b <- cbind(estimated_b,
                           this_covariate_function_values);
    }
    colnames(estimated_b) <- varying_effects_names;
  }
  ###############################################################
  # Calculate confidence intervals using sandwich formula
  # to better take into account the existence of within-subject
  # correlation;
  ###############################################################
  # Consult design matrix of spline bases on grid of observed times:
  design <- predict.bam(model1,type="lpmatrix");
  # Construct design matrix of spline bases on regular grid of times:
  temp_data <- data.frame(as.matrix(0*model1$model[1,])%x%rep(1,length(time_grid)));
  colnames(temp_data) <- colnames(as.data.frame(model1$model));
  temp_data[,time_variable_name] <- time_grid;
  grid_design <- predict.bam(model1,type="lpmatrix",newdata = temp_data); 
  # Working independence variance estimates (will make sandwich):
  if (penalize) {
    bread <- model1$Vc;
  } else {
    bread <- model1$Vp;
  }
  npar <- length(model1$coefficients);
  if (use_naive_se) {
    sandwich <- bread;
  } else {
    meat <- matrix(0,npar,npar);  # apologies to vegetarians -- 
    # peanut butter or vegetables are also fine!
    working_sigsqd <- var(model1$residuals);
    for (i in unique(id_variable[which(!is.na(id_variable))])) {
      these <- which(id_variable==i);
      residuals_these <- model1$y[these] - model1$fitted.values[these];
      if (use_weights) {
        stopifnot(length(residuals_these) == length(model1$weights[these]));
        residuals_these <- residuals_these * model1$weights[these];
      }
      if (family$family=="gaussian") {
        multiplier <- 1/working_sigsqd;
      }
      if (family$family=="binomial" | family$family=="poisson") {
        multiplier <- 1;
      }
      if (length(these)>0) {
        meat <- meat + t(design[these,,drop=FALSE])%*%
          (outer(multiplier*residuals_these,
                 multiplier*residuals_these))%*%
          design[these,,drop=FALSE];
      }
    }  
    sandwich <- bread %*% meat %*% bread;  
  }
  general_term_name <- sub(" .*","", 
                           gsub(x=names(model1$coefficients),
                                pattern="[.]",
                                replacement=" "));
  # removes the .1, .2, .3, etc., from after the term name in the names
  # of coefficients.
  indices_for_b0 <- which((general_term_name=="(Intercept)") | 
                            (general_term_name == paste("s(",
                                                        time_variable_name,
                                                        ")",
                                                        sep="")));
  basis_for_b0 <-  design[,indices_for_b0];
  estimated_b0_hard_way <-  basis_for_b0%*%model1$coefficients[indices_for_b0];  
  stopifnot(max(abs(estimated_b0_hard_way-estimated_b0), na.rm =TRUE)<1e-10); #just for double checking;
  grid_basis_for_b0 <- grid_design[,indices_for_b0]; 
  grid_estimated_b0 <- grid_basis_for_b0 %*% model1$coefficients[indices_for_b0];
  cov_mat_for_b0 <- sandwich[indices_for_b0,indices_for_b0];
  # Fitted b0 on observed times:
  temp_function <- function(i){return(t(basis_for_b0[i,]) %*%
                                        cov_mat_for_b0 %*%
                                        basis_for_b0[i,])};
  standard_error_b0 <- drop(sqrt(sapply(X=1:nrow(basis_for_b0),
                                        FUN=temp_function)));
  upper_b0 <- estimated_b0 + crit_value*standard_error_b0;
  lower_b0 <- estimated_b0 - crit_value*standard_error_b0;
  fitted_coefficients <- list("(Intercept)"=data.frame(estimate = estimated_b0,
                                                       standard_error = standard_error_b0,
                                                       upper = upper_b0,
                                                       lower = lower_b0));
  # Fitted b0 on regular grid:
  temp_function <- function(i){return(t(grid_basis_for_b0[i,]) %*%
                                        cov_mat_for_b0 %*%
                                        grid_basis_for_b0[i,])};
  grid_standard_error_b0 <- drop(sqrt(sapply(X=1:nrow(grid_basis_for_b0),
                                             FUN=temp_function)));
  grid_upper_b0 <- grid_estimated_b0 + crit_value*grid_standard_error_b0;
  grid_lower_b0 <- grid_estimated_b0 - crit_value*grid_standard_error_b0;
  grid_fitted_coefficients <- list("(Intercept)"=data.frame(estimate = grid_estimated_b0,
                                                            standard_error = grid_standard_error_b0,
                                                            upper = grid_upper_b0,
                                                            lower = grid_lower_b0));
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      this_covariate_name <- varying_effects_names[i];
      this_covariate_values <- model1$model[,this_covariate_name];
      this_spline_term <- paste("s(",time_variable_name,"):",
                                this_covariate_name,sep="");
      indices_for_this_b <- which((general_term_name == this_covariate_name) | 
                                    (general_term_name == this_spline_term)); 
      estimated_this_b <- (basis_for_b0%*%model1$coefficients[indices_for_this_b]);
      grid_estimated_this_b <- (grid_basis_for_b0%*%model1$coefficients[indices_for_this_b]);
      # Using the b0 basis for b1 makes code simpler but assumes that number of knots, 
      # spline order, and coefficient order are the same between 
      # intercept and substantive effect; perhaps remove this assumption
      # in some later version;   
      estimated_this_b <- estimated_b[,i]; 
      estimated_this_b_hard_way <- as.vector(basis_for_b0%*%model1$coefficients[indices_for_this_b]);
      stopifnot(max(abs(estimated_this_b_hard_way-estimated_b[,i]), na.rm=TRUE)<1e-10); #just for double checking;
      cov_mat_for_this_b <- sandwich[indices_for_this_b,indices_for_this_b];
      temp_function <- function(i){return(t(basis_for_b0[i,]) %*%
                                            cov_mat_for_this_b %*% 
                                            basis_for_b0[i,])};
      standard_error_this_b <- drop(sqrt(sapply(X=1:nrow(basis_for_b0),
                                                FUN=temp_function)));
      upper_this_b <- estimated_this_b + crit_value*standard_error_this_b;
      lower_this_b <- estimated_this_b - crit_value*standard_error_this_b;
      temp_list <- list(data.frame(estimate = estimated_this_b,
                                   standard_error = standard_error_this_b,
                                   upper = upper_this_b,
                                   lower = lower_this_b));
      names(temp_list) <- this_covariate_name;
      fitted_coefficients <- c(fitted_coefficients,
                               temp_list);
      temp_function <- function(i){return(t(grid_basis_for_b0[i,]) %*%
                                            cov_mat_for_this_b %*%
                                            grid_basis_for_b0[i,])};
      grid_standard_error_this_b <- drop(sqrt(sapply(X=1:nrow(grid_basis_for_b0),
                                                     FUN=temp_function)));
      grid_upper_this_b <- grid_estimated_this_b + crit_value*grid_standard_error_this_b;
      grid_lower_this_b <- grid_estimated_this_b - crit_value*grid_standard_error_this_b;
      temp_list_grid <- list(data.frame(estimate = grid_estimated_this_b,
                                        standard_error = grid_standard_error_this_b,
                                        upper = grid_upper_this_b,
                                        lower = grid_lower_this_b));
      names(temp_list_grid) <- this_covariate_name;
      grid_fitted_coefficients <- c(grid_fitted_coefficients,
                                    temp_list_grid);
    }
    names(fitted_coefficients[2:length(fitted_coefficients)]) <- varying_effects_names;
    names(grid_fitted_coefficients[2:length(grid_fitted_coefficients)]) <- varying_effects_names;
  }
  invar_effects_estimates <- NULL;
  if (num_invar_effects>0) {
    invar_effects_indices <- which(general_term_name %in% invar_effects_names);
    invar_effects_est <- coef(model1)[invar_effects_names];
    invar_effects_se <- sqrt(diag(sandwich[invar_effects_indices,
                                           invar_effects_indices,
                                           drop=FALSE]));
    invar_effects_estimates <- data.frame(estimate=invar_effects_est,
                                          standard_error=invar_effects_se);
  }
  ##################################
  # Return answers;
  ##################################
  nsub <- length(unique(id_variable[which(!is.na(id_variable))]));
  model_information <- list(outcome_family=family$family,
                            response_name=response_name,
                            num_varying_effects=num_varying_effects,
                            varying_effects_names=varying_effects_names,
                            num_invar_effects=num_invar_effects,
                            invar_effects_names=invar_effects_names,
                            n_subjects=nsub,
                            pseudo_aic=AIC(model1),
                            pseudo_bic=AIC(model1,k=log(nsub)),
                            used_listwise_deletion=used_listwise_deletion);
  if (use_weights) {
    model_information <- c(model_information,
                           weights_variable_name = weights_variable_name);
  }
  answer <- list(time_grid=time_grid,
                 grid_fitted_coefficients=grid_fitted_coefficients,
                 invar_effects_estimates=invar_effects_estimates,
                 model_information=model_information,
                 back_end_model=model1);
  class(answer) <- "tvem";
  return(answer);
}  