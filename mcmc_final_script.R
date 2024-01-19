#############################################################################
# Task: Code to run results shared in the mcmc final report. 
# Author: Alex Kuhn
#############################################################################

source('mcmc_final_helper.R')

library(tidyverse)
library(mcmc)
library(SimTools)
library(mcmcplots)
library(mcmcse)
library(R2jags)

##################################################################
# Simulate data for use below
##################################################################
set.seed(0202)

# True parameter values
L <- 12000                       # (known) wavelength lower bound 
U <- 18000                       # (known) wavelength upper bound
K <- 500                         # (known) number of bins (sample size) 
gamma <- 100                     # (unknown) true background rate
mu <- 17000                      # (unknown) true location of signal 
c <- 20                          # (known) sd for density in mean function 
eta <- 7000                      # (unknown) signal intensity
d <- 30                          # (known) sd for Gaussian process

# Set up bin structure for data
bin_obj <- bin_func(L=L, U=U, K=K)  
delta <- bin_obj$delta
bin_lbs <- bin_obj$bin_lbs
bin_ubs <- bin_obj$bin_ubs
lambdas <- bin_obj$lambdas       # (known) vector of centers of each bin

# Vectors of mean and sd for each observation in the GP
m <- m_func(lambdas, gamma, eta, mu, c)  # mean of GP
sd <- sd_func(lambdas, sd = d)         # sd of GP

# Generate GP as *independent* normal draws according to the above
y <- rnorm(n = K, mean = m, sd = sd)
df <- data.frame(y, lambdas)

# Plot simulated data
df %>% ggplot(aes(x = lambdas, y = y)) + 
  geom_line() + 
  geom_vline(xintercept = mu, color = "red") + 
  ylab("Flux (Y)") + 
  xlab("wavelength (lambda)") + 
  theme_minimal()

##################################################################
# Define Priors
##################################################################

# For gamma prior, want mean of 100 and sd about 10 (smaller than realistic)
#a <- 10
#b <- 1/10
a <- 100
b <- 1

# For eta prior, want mean of 7000 and sd of 100 (smaller than realistic)
#e <- 10
#f <- 1/700
e <- 4900
f <- 7/10

# Prior on mu is always Uniform([L,U])

##################################################################
# Laplace Approximation: Model m1 (reduced model)
##################################################################

# Max gamma is a root of quadratic
a_prime <- 1/(d^2)
b_prime <- (b/K - mean(y)/(d^2))
c_prime <- (1 - a)/K

gamma_hat <- (-b_prime + sqrt(b_prime^2 - 4 * a_prime * c_prime))/(2 * a_prime)
gamma_hat # 101.1754

hess_null_hat <- 1/(d^2) + (a - 1)/(K * gamma_hat^2)
log_det_hess_null <- log(abs(hess_null_hat)) # -6.785136

# Compute LA based on expansion
log_LA_null <- lupost_null(theta = gamma_hat) + (1/2)*log(2*pi) + 
  (1/2)*log_det_hess_alt - (1/2) * K # -2711.441

##################################################################
# Laplace Approximation: Model m2 (full model)
##################################################################

# (Contstrained) optimization to obtain max of posterior density
bounds <- matrix(c(0,Inf, 0,Inf, L, U), nc=2, byrow=TRUE)
colnames(bounds) <- c("lower", "upper")
D <- nrow(bounds)
ui <- rbind( diag(D), -diag(D) )
ci <- c( bounds[,1], - bounds[,2] )
finite <- as.vector(is.finite(bounds))
ui <- ui[finite,]
ci <- ci[finite]

# Initial values for optimization
params_init <- c(median(y),             # central y value, ignoring outliers
                 qgamma(0.5, shape = e, rate = f), # median of prior
                 lambdas[which.max(y)]) # Where it *looks* like there is a peak

# Run BFGS 
opt_res <- constrOptim(theta = params_init, 
                       f = lupost_alt,
                       grad = gradient_alt,
                       #method = "Nelder-Mead",
                       hessian = TRUE,
                       ui = ui, ci = ci)
theta_alt_hat <- opt_res$par     # 99.06569  6999.52382 17010.00000
hess_alt_hat <- opt_res$hessian  # Seems reasonable
#               [,1]          [,2]          [,3]
# [1,]  1.131816e-03  1.851852e-07 -5.421011e-17
# [2,]  1.851852e-07  2.168874e-07 -1.320222e-09
# [3,] -5.421011e-17 -1.320222e-09  7.306703e-05
log_det_hess_alt <- determinant(hess_alt_hat, logarithm = TRUE)$modulus[1]

# Compute LA based on expansion 
log_LA_alt <- lupost_alt(theta = theta_alt_hat) + (3/2)*log(2*pi) + 
            (1/2)*log_det_hess_alt - (3/2) * K  # -3201.247

##################################################################
# Bayes Factor estimate based on two LA's
##################################################################

log_BF_12 <- log_LA_null - log_LA_alt # 489.8063
log_10_BF_12 <- log(exp(log_LA_null - log_LA_alt), base = 10) # may not be great
log_10_BF_12 # 212.7202; Decisive evidence in favor of Model 2 (correct)

##################################################################
# Serial Tempering; now that lupost_alt() is corrected
##################################################################
set.seed(0202)

### Set up requirements for `temper()` function following Charlie's writeup
# at https://cran.r-project.org/web/packages/mcmc/vignettes/bfst.pdf

# Each model is a neighbor to the other, but not itself
neighbors <- matrix(c(FALSE, TRUE, TRUE, FALSE), 
                    nrow = 2, ncol = 2, byrow = TRUE)

# Initial values; use output from optimization in LA
state.initial <- c(1, theta_alt_hat)
qux <- rep(0, 2) # initial pseudopriors

# Run initial mcmc chain
out <- temper(ludfun, initial = state.initial, neighbors = neighbors, 
              nbatch = 1000, blen = 1,
              log.pseudo.prior = qux)

ibar <- colMeans(out$ibatch)
ibar # Poor mixing

# Initialize log pseudoprior until it "converges" roughly
qux.save <- qux
time.save <- out$time
repeat{
  out <- temper(out, log.pseudo.prior = qux)
  ibar <- colMeans(out$ibatch)
  qux <- qux + pmin(log(max(ibar)/ibar), 10)
  qux <- qux - min(qux)
  qux.save <- rbind(qux.save, qux, deparse.level = 0)
  time.save <- rbind(time.save, out$time, deparse.level = 0)
  if (max(ibar) / min(ibar) < 2) break
}

qux # Look at what we got

# Look at acceptance rates; a bit high
out$accepti
#        NA 0.7947368
# 0.4733542        NA

out$acceptx
# 0.7752809 0.7316294; high as well 

# Adjust acceptance rates for within model moves by adjusting scaling
out <- temper(out, scale = 10, log.pseudo.prior = qux)
out$acceptx 
# 0.1776316 0.2083333; looks alright

### Do a longer run now
out <- temper(out, nbatch = 5e3,
              scale = 10, log.pseudo.prior = qux)

### Some initial diagnostics
mix_post_sample <- as.mcmc(out$batch)
head(mix_post_sample)

plot(ts(out$batch)) # appears to be working!
plot(ts(out$ibatch)) # not a super helpful plot; can see we stay in single 
                     # model for quite a while

acf(mix_post_sample[,1])
acf(mix_post_sample[,2]) # pretty bad
acf(mix_post_sample[,3]) # pretty bad

apply(mix_post_sample[,1:3], 2, ess)
# var1        var2       var3 
# 225.486784  6.350786   6.935909

S_post<-Smcmc(mix_post_sample) # Which model is this associated with?
plot(S_post, Q=c(0.05, 0.95),  xlab="", ylab="")

minESS(p=1 + 3 * 3)
# 8831 

### Do a much longer run 
out2 <- temper(out, nbatch = 5e5,
               scale = 10, log.pseudo.prior = qux)

full_post_sample <- rbind(out$batch, out2$batch)
full_post_isample <- rbind(out$ibatch, out2$ibatch)

full_post_sample <- as.mcmc(full_post_sample)
head(full_post_sample)

apply(full_post_sample[,1:3], 2, ess)
# var1         var2      var3 
# 23124.44049  88.33316  15.23771

S_post<-Smcmc(full_post_sample)
plot(S_post, Q=c(0.05, 0.95),  xlab="", ylab="")

#save(full_post_sample, file = "serial_temper_result.rda")
#save(full_post_isample, file = "serial_temper_iresult.rda")
load("serial_temper_result.rda") ## LOAD RESULTS IN HERE
load("serial_temper_iresult.rda")

# Note: Still large MCSE for mu, but let's proceed for 
#       the sake of demonstration. 
#       See figure /Dropbox/line_detect/serial_temper_plot.pdf.

##################################################################
# Bayes factor estimated by Serial Tempering
##################################################################

# Calculate log 10 Bayes factor of Model 1 vs. Model 2
log_10_BF_terms_ST <- (qux - log(colMeans(full_post_isample))) / log(10)

log_10_BF_12_ST <- log_10_BF_terms_ST[1] - log_10_BF_terms_ST[2]
log_10_BF_12_ST # 8.183456; way smaller, but same decision as LA

# What is the MCSE for the model state parameter i? 

##################################################################
# MCSE for BF w/ Serial Tempering
##################################################################

# Calculations still following Charlie's write up here:
# https://cran.r-project.org/web/packages/mcmc/vignettes/bfst.pdf

### Delta method approach (relies on naive sample variance)
fred <- var(full_post_isample)/nrow(full_post_isample)
sally <- colMeans(full_post_isample)
mcse_log_10_BF <- (1 / log(10)) * sqrt(diag(fred) / sally^2 - 
                                         2 * fred[, 2] / (sally * sally[2]) +
                                         fred[2, 2] / sally[2]^2)
mcse_log_10_BF
# 0.004960854 0.000000000

ST_output_summary <- c(log_10_BF_12_ST, mcse_log_10_BF[1])
ST_output_summary
# 8.183455507 0.004960854

### Batch means approach w/ best linear approximation to BF
ibar <- colMeans(full_post_isample)
herman <- sweep(full_post_isample, 2, ibar, "/")
herman <- sweep(herman, 1, herman[ , 2], "-")
mcse_log_10_BF_batch <- (1/log(10)) * 
  apply(herman, 2, sd) / sqrt(nrow(full_post_isample))
# 0.004960854 0.000000000; same values 





