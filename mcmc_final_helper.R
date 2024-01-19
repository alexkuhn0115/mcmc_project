#############################################################################
# Task: Helper functions for the script `mcmc_final_script.R`.
# Author: Alex Kuhn
#############################################################################

##################################################################
# Data generation helper functions
##################################################################

# Setup bins for given range and number of bins
bin_func <- function(L, U, K) {
  delta <- (U - L)/K                   # bin widths
  bin_lbs <- L + 0:(K-1) * delta       # collection of bin lower bds
  bin_ubs <- L + 1:K * delta           # collection of bin upper bds
  lambdas <- (bin_lbs + bin_ubs)/2     # collection of bin centers "lambda_k"'s
  return(list(delta=delta, bin_lbs=bin_lbs, bin_ubs=bin_ubs, lambdas=lambdas))
}


# Double check this: 
# https://stackoverflow.com/questions/63165956/r-select-n-evenly-spaced-out-elements-in-vector-including-first-and-last
split_func <- function(x, by) {
  r <- diff(range(x))
  out <- seq(0, r - by - 1, by = by)
  c(round(min(x) + c(0, out - 0.51 + (max(x) - max(out)) / 2), 0), max(x))
}

# Background function, takes vector of lambdas and gamma parameter
b_func <- function(lambdas, gamma) {
  rep(gamma, K) # constant for now
}

# Signal function, takes vector of lambdas, mean and sd for normal density
s_func <- function(lambdas, mu, sigma) {
  dnorm(lambdas, mean = mu, sd = sigma) # maybe change to truncated normal
                                        # since values can only be in [L,U]
}

# Function to construct mean function directly 
m_func <- function(lambdas, gamma, eta, mu, sigma){
  b_func(lambdas, gamma) + eta * s_func(lambdas, mu, sigma)
}

# Function to construct variance function for vector of lambdas
sd_func <- function(lambdas, sd) {
  rep(sd, K) # constant for now
}

##################################################################
# M2 Model: Log posterior densities and derivatives 
##################################################################

### Define log-unnormalized density function for m2 ("alternative model")
### TODO: ADD INPUTS; currently based on R session
lupost_alt <- function(theta) {
  
  ## Following relies on variables in environment
  # Unnormalized log density ; return -Inf if state out of param space
  if (theta[1] > 0 && theta[2] > 0 && L <= theta[3] && theta[3] <= U) {
    
    # Model defines mean structure; all else the same
    m <- theta[1] + theta[2] * dnorm(lambdas, mean = theta[3], sd = c)
    
    sum(dnorm(y, mean = m, sd = d, log = TRUE)) + 
      dgamma(theta[1], shape = a, rate = b, log = TRUE) + 
      dgamma(theta[2], shape = e, rate = f, log = TRUE) - 
      log(U - L)
  } else {
    -Inf
  }
}

### Functions below all partial derivatives for the gradient
phi_k <- function(theta) {
  dnorm(lambdas, mean = theta[3], sd = c)
}

dgam_alt <- function(theta) {
  -mean(y)/d^2 + theta[1]/d^2 + theta[2]/(K * d^2) * sum(phi_k(theta)) + 
    (1-a)/(K * theta[1]) + b/K
}

deta_alt <- function(theta) {
  -1/(K * d^2) * sum(y * phi_k(theta)) + theta[1]/(K * d^2) * sum(phi_k(theta)) + 
    theta[2]/(K * d^2) * sum(phi_k(theta)^2) + (1 - e)/(K * theta[2]) + f/K
}

dmu_alt <- function(theta) {
  1/(K * c^2 * d^2) * 
    (
      -theta[2] * sum(y * phi_k(theta) * (lambdas - theta[3])) + 
        theta[1] * theta[2] * sum(phi_k(theta) * (lambdas - theta[3])) + 
        theta[2]^2 * sum(phi_k(theta)^2 * (lambdas - theta[3])) 
    )
}

gradient_alt <- function(theta) {
  c(dgam_alt(theta), deta_alt(theta), dmu_alt(theta))
}

##################################################################
# M2 Model: Log posterior densities and derivatives 
##################################################################

lupost_null <- function(theta) {
  
  ## Following relies on variables in environment
  # Unnormalized log density ; return -Inf if state out of param space
  if (theta > 0) {

    sum(dnorm(y, mean = theta, sd = d, log = TRUE)) + 
      dgamma(theta[1], shape = a, rate = b, log = TRUE)
    
  } else {
    -Inf
  }
}

##################################################################
# Define function to use with `temper()`
##################################################################

# Define log-unnormalized density function(s)
ludfun <- function(state, log.pseudo.prior) {
  stopifnot(is.numeric(state))
  stopifnot(length(state) == 1 + 3) # model indicator + 3 params (gamma, eta, mu)
  icomp <- state[1]
  stopifnot(icomp == as.integer(icomp))
  stopifnot(1 <= icomp && icomp <= 2) # only have 2 models
  stopifnot(is.numeric(log.pseudo.prior))
  stopifnot(length(log.pseudo.prior) == 2) # 2 models; need pseudoprior for each
  
  ## Following relies on variables in environment
  
  # Unnormalized log density h_m(theta); return -Inf if state out of param space
  if (state[2] > 0 && state[3] > 0 && L <= state[4] && state[4] <= U) {
    
    # Model defines mean structure; all else the same
    if (icomp == 1) {
      m <- state[2]
    } else {
      m <- state[2] + state[3] * dnorm(lambdas, mean = state[4], sd = c)
    }
    
    sum(dnorm(y, mean = m, sd = d, log = TRUE)) + 
      dgamma(state[2], shape = a, rate = b, log = TRUE) + 
      dgamma(state[3], shape = e, rate = f, log = TRUE) - 
      log(U - L) + 
      log.pseudo.prior[icomp]
    
  } else {
    -Inf
  }
}

##################################################################
# Helper functions for MLE calculation in frequentist setting 
##################################################################

# TODO: Add sd as an argument for the function
alt_ll <- function(theta, mu_r, y, lambdas) {
  gamma <- theta[1]
  eta <- theta[2]
  m <- gamma + eta * dnorm(lambdas, mean = mu_r, sd = a)
  sum(dnorm(y, mean = m, sd = b, log = TRUE))
}

# TODO: - Add sd as an argument for the function
#       - Check that these calculations don't have numeric issues
get_alt_MLE <- function(y, lambdas, mu_r){
  
  phi_vec <- dnorm(lambdas, mean = mu_r, sd = a)
  
  ybar <- mean(y)
  phi_bar <- mean(phi_vec)
  yphi_bar <- mean(y * phi_vec)
  phi2_bar <- mean(phi_vec^2)
  
  eta_hat <- (ybar * phi_bar - yphi_bar)/((phi_bar)^2 - phi2_bar)
  eta_hat <- ifelse(eta_hat < 0, 0, eta_hat) # enforce nonnegative
  
  gamma_hat_1 <- ybar - eta_hat * phi_bar
  gamma_hat_1 <- ifelse(gamma_hat_1 < 0, 0, gamma_hat_1) # enforce nonnegative
  
  return(c(gamma_hat_1 = gamma_hat_1, eta_hat = eta_hat))
}






































