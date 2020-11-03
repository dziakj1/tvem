#' simulate_tvem_example: Simulate a dataset for demonstrating the tvem function. 
#' 
#' By default, the data-generating model has a time-varying intercept,
#' and two time-varying covariates named x1 and x2. 
#' x1 has a time-varying effect and x2 has a time-invariant effect.
#' 
#' @param n_subjects Number of subjects in dataset 
#' @param max_time The time point at the end of the simulated time interval 
#' @param simulate_binary Whether the simulated data should be binary
#' @param n_obs_possible Total number of possible measurement times per subject
#' @param prop_obs_observed Proportion of these that are actually observed
#' @param sigma_x1 Standard deviation of covariate 1, assumed homoskedastic over time
#' @param sigma_x2 Standard deviation of covariate 2, assumed homoskedastic over time
#' @param truncate_for_realism Whether to prevent simulated values from going below 0 or above 10, 
#' in order to imitate survey data; used only for normally distributed outcomes.
#' Set this to FALSE if you are running a simulation, in order to avoid losing coverage
#' due to departing from the parametric model.
#' @param round_digits Number of digits at which to round the generated data; used only for normally distributed outcomes 
#' @param x1_short_term_rho Correlation between adjacent measurements of covariate 1
#' @param x2_short_term_rho Correlation between adjacent measurements of covariate 2
#' @param sigma_y  Error standard deviation of y, only used if the outcomes are to be normal rather than binary 
#' @param mu_x1_function Mean of covariate 1 as function of time
#' @param mu_x2_function Mean of covariate 2 as function of time 
#' @param beta0_y_function TVEM intercept as function of time
#' @param beta1_y_function TVEM coefficient of covariate 1 as function of time  
#' @param beta2_y_function TVEM coefficient of covariate 2 as function of time 
#' @param y_short_term_rho Correlation between adjacent measurements of y, only used if it is normal and not binary 
#' @return A simulated dataset with the following variables:
#' \describe{
#' \item{subject_id}{Subject ID}
#' \item{time}{Observation time}
#' \item{x1}{First covariate}
#' \item{x2}{Second covariate}
#' \item{y}{Outcome variable}
#' }   
#'
#' @keywords Statistics|datagen
#'
#' @importFrom stats rbinom rnorm
#'
#' @export

simulate_tvem_example <- function(
  n_subjects = 300, 
  max_time = 7,
  simulate_binary = FALSE,
  n_obs_possible = 141,
  prop_obs_observed = .3, 
  sigma_x1 = 2,
  sigma_x2 = 2,
  truncate_for_realism = TRUE,
  round_digits = 3,
  x1_short_term_rho = .7,
  x2_short_term_rho = .7,
  sigma_y = 1.5,
  mu_x1_function = function(t) {6 - 2*sqrt(pmax(0,(t/7)-.2))},
  mu_x2_function = function(t) {3 + sqrt(pmax(0,(t/7)-.5))},
  beta0_y_function = function(t) {1 - .3*sqrt(t/7)},
  beta1_y_function = function(t) {.5*(t/7)^2}, 
  beta2_y_function = function(t) {rep(.2,length(t))},
  y_short_term_rho = 0.8 
  ) {
  possible_observation_times <- seq(0,max_time,length=n_obs_possible);
  mu_x1 <- mu_x1_function(possible_observation_times);
  mu_x2 <- mu_x2_function(possible_observation_times); 
  beta0_y <- beta0_y_function(possible_observation_times);
  beta1_y <- beta1_y_function(possible_observation_times);
  beta2_y <- beta2_y_function(possible_observation_times);
  n_obs_per_subject <- rbinom(n_subjects,n_obs_possible,prop_obs_observed);
  # Generate X1;
  x1_error_term_AR1 <- matrix(NA, n_subjects, n_obs_possible);
  x1_error_term_AR1[,1] <- rnorm(n_subjects);
  for (j in 2:n_obs_possible) {
    x1_error_term_AR1[,j] <- x1_short_term_rho*x1_error_term_AR1[,j-1] + 
      sqrt(1-x1_short_term_rho^2)*rnorm(n_subjects);
  }
  x1 <- t(apply(sigma_x1*x1_error_term_AR1,1,"+",mu_x1));
  x1 <- round(x1,round_digits);
  # Generate X2;
  x2_error_term_AR1 <- matrix(NA, n_subjects, n_obs_possible);
  x2_error_term_AR1[,1] <- rnorm(n_subjects);
  for (j in 2:n_obs_possible) {
    x2_error_term_AR1[,j] <- x2_short_term_rho*x2_error_term_AR1[,j-1] + 
      sqrt(1-x2_short_term_rho^2)*rnorm(n_subjects);
  }
  x2 <- t(apply(sigma_x2*x2_error_term_AR1,1,"+",mu_x2));
  x2 <- round(x2,round_digits);
  # Generate Y;
  matrix_beta0_y <- matrix(rep(beta0_y,n_subjects),byrow=TRUE,nrow=n_subjects);
  matrix_beta1_y <- matrix(rep(beta1_y,n_subjects),byrow=TRUE,nrow=n_subjects);
  matrix_beta2_y <- matrix(rep(beta2_y,n_subjects),byrow=TRUE,nrow=n_subjects); 
  if (simulate_binary) {
    eta_y <- matrix_beta0_y + matrix_beta1_y*x1 + matrix_beta2_y*x2;
    mu_y <- plogis(eta_y); 
    y <- apply(mu_y,MARGIN=c(1,2),FUN=rbinom,n=1,size=1);
  } else {
      mu_y <- matrix_beta0_y + matrix_beta1_y*x1 + matrix_beta2_y*x2;
      y_error_term_AR1 <- matrix(NA, n_subjects, n_obs_possible);
      y_error_term_AR1[,1] <- rnorm(n_subjects);
      for (j in 2:n_obs_possible) {
        y_error_term_AR1[,j] <- y_short_term_rho*y_error_term_AR1[,j-1] + 
          sqrt(1-y_short_term_rho^2)*rnorm(n_subjects);
      }
      y <- mu_y + sigma_y * y_error_term_AR1;  
      y <- round(y,round_digits);
    } 
  for (this_subject in 1:n_subjects) {
    unobserved_for_this_subject <- sample(1:ncol(y),size=n_obs_possible - n_obs_per_subject[this_subject]);
    y[this_subject,unobserved_for_this_subject] <- NA;
  }
  entire_long_dataset <- data.frame(subject_id=rep(1:n_subjects,each=n_obs_possible),
                                    time=rep(possible_observation_times,times=n_subjects),
                                    x1=as.vector(t(x1)),
                                    x2=as.vector(t(x2)),
                                    y=as.vector(t(y)));
  long_dataset <- entire_long_dataset[which(is.na(entire_long_dataset$y)==FALSE),]
  if (truncate_for_realism) {
    entire_long_dataset$x1 <- round(pmin(pmax(entire_long_dataset$x1,0),10),1);
    entire_long_dataset$x2 <- round(pmin(pmax(entire_long_dataset$x2,0),10),1);
    entire_long_dataset$y <- round(pmin(pmax(entire_long_dataset$y,0),10),1);
  }
  return(entire_long_dataset); 
}