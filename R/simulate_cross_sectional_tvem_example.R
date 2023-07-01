#' simulate_cross_sectional_tvem_example: Simulate a cross-sectional dataset for 
#' demonstrating tvem. 
#' 
#' This differs from what simulate_tvem_example does, in that each
#' participant (subject) is only measured at a single measurement time,
#' but the measurement time of each subject is random.  This might
#' simulate a sample of people of many different ages, where age 
#' is treated as a continuous number.
#' 
#' By default, the data-generating model has a time-varying intercept,
#' and two time-varying covariates named x1 and x2. 
#' x1 has a time-varying effect and x2 has a time-invariant effect.
#' 
#' @param n_subjects Number of subjects in dataset 
#' @param min_time The time point at the end of the simulated time interval 
#' @param max_time The time point at the end of the simulated time interval 
#' @param simulate_binary Whether the simulated data should be binary
#' @param sigma_x1 Standard deviation of covariate 1, assumed homoskedastic over time
#' @param sigma_x2 Standard deviation of covariate 2, assumed homoskedastic over time
#' @param truncate_for_realism Whether to prevent simulated values from going below 0 or above 10, 
#' in order to imitate survey data; used only for normally distributed outcomes.
#' Set this to FALSE if you are running a simulation, in order to avoid losing coverage
#' due to departing from the parametric model.
#' @param round_digits Number of digits at which to round the generated data; used only for normally distributed outcomes 
#' @param sigma_y  Error standard deviation of y, only used if the outcomes are to be normal rather than binary 
#' @param mu_x1_function Mean of covariate 1 as function of time
#' @param mu_x2_function Mean of covariate 2 as function of time 
#' @param beta0_y_function TVEM intercept as function of time
#' @param beta1_y_function TVEM coefficient of covariate 1 as function of time  
#' @param beta2_y_function TVEM coefficient of covariate 2 as function of time 
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
#' @importFrom stats rbinom rnorm runif
#'
#' @examples 
#' set.seed(16802)
#' the_data <- simulate_tvem_example(simulate_binary=TRUE)
#'
#' @export

simulate_cross_sectional_tvem_example <- function(
    n_subjects = 300, 
    min_time = 7,
    max_time = 7,
    simulate_binary = FALSE,
    sigma_x1 = 2,
    sigma_x2 = 2,
    truncate_for_realism = TRUE,
    round_digits = 3,
    sigma_y = 1.5,
    mu_x1_function = function(t) {6 - 2*sqrt(pmax(0,((t-min(t))/7)-.2))},
    mu_x2_function = function(t) {3 + sqrt(pmax(0,((t-min(t))/7)-.5))},
    beta0_y_function = function(t) {1 - .3*sqrt((t-min(t))/7)},
    beta1_y_function = function(t) {.5*((t-min(t))/7)^2}, 
    beta2_y_function = function(t) {rep(.2,length(t))}
) {
  observation_times <- round(runif(n_subjects, min=min_time, max=max_time),round_digits);
  mu_x1 <- mu_x1_function(observation_times);
  mu_x2 <- mu_x2_function(observation_times); 
  beta0_y <- beta0_y_function(observation_times);
  beta1_y <- beta1_y_function(observation_times);
  beta2_y <- beta2_y_function(observation_times);
  # Generate X1;
  x1_error_term <- rnorm(n_subjects);
  x1 <- sigma_x1*x1_error_term + mu_x1;
  x1 <- round(x1,round_digits);
  # Generate X2;
  x2_error_term <- rnorm(n_subjects);
  x2 <- sigma_x2*x2_error_term + mu_x2;
  x2 <- round(x2,round_digits);
  # Generate Y;
  if (simulate_binary) {
    eta_y <- beta0_y + beta1_y*x1 +beta2_y*x2;
    mu_y <- plogis(eta_y); 
    y <- rbinom(n=n_subjects, size=1, prob=mu_y);
  } else {
    mu_y <- beta0_y + beta1_y*x1 + beta2_y*x2;
    y_error_term <- rnorm(n_subjects);
    y <- mu_y + sigma_y * y_error_term;  
    y <- round(y,round_digits);
  } 
  entire_dataset <- data.frame(subject_id=1:n_subjects,
                                    time=observation_times,
                                    x1=x1,
                                    x2=x2,
                                    y=y);
  if (truncate_for_realism) {
    entire_dataset$x1 <- round(pmin(pmax(entire_dataset$x1,0),10),1);
    entire_dataset$x2 <- round(pmin(pmax(entire_dataset$x2,0),10),1);
    entire_dataset$y <- round(pmin(pmax(entire_dataset$y,0),10),1);
  }
  return(entire_dataset); 
}