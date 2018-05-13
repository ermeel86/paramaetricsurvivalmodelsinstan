# Parametric proportional hazard survival model that uses M-splines for the baseline hazard or
# equivalently I-splines for the cumulative baseline hazard.
# In the Parametric and penalized generalized survival models lingua of Liu et al. we use the -log link function.
# The splines are functions of t (not log(t)) as done by Royston-Parmar (equivalently when using log -log link)
library(rstan)
library(flexsurv)
library(bayesplot)
library(splines2)
library(tidyverse)
library(survival)
library(simsurv)

scale_X <- FALSE # whether or not to scale the feature (covariate) matrix X
order <- 3 #equals degree+1 of spline basis functions
ninterior_knots <- 2 # needs to be >= 1

sm <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/pgsm_log_link_2_qr.stan")

# get the breast cancer survival data set used by Royston & Parmar
# see https://cran.r-project.org/web/packages/flexsurv/flexsurv.pdf#Rfn.bc.1
df <- as.tibble(flexsurv::bc)
N <- nrow(df)
msk_censored <- df$censrec == 0
is_censored <- as.integer(msk_censored) # censrec==1 means dead, 0 means censored in the original data
X <- as.matrix(cbind(as.integer(df$group == "Medium"), as.integer(df$group=="Poor")))
if(scale_X) X <- scale(X)
X_qr <- qr(X)
Q <- qr.Q(X_qr)*sqrt(N-1)
Q_censored <- as.matrix(Q[msk_censored,])
Q_uncensored <- as.matrix(Q[!msk_censored,])
R <- qr.R(X_qr)/sqrt(N-1)

times <- as.vector(df$rectime)

N_censored <- sum(is_censored)
N_uncensored <- N-N_censored
times_uncensored <- times[!msk_censored]

knots <- quantile(times_uncensored,head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
nknots <- length(knots)
isOut <- iSpline(times, knots = knots, degree = order-1, intercept = TRUE)
deriv_isOut <- deriv(isOut)
isOut_censored <- isOut[msk_censored,]
isOut_uncensored <- isOut[!msk_censored,]
deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]
nbasis <- dim(isOut)[2]
stan_data <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                   m=nbasis, Q_censored=Q_censored, Q_uncensored=Q_uncensored,R=R, 
                   NC=ncol(X),
                   basis_evals_censored=t(isOut_censored), 
                   basis_evals_uncensored=t(isOut_uncensored),
                   deriv_basis_evals_uncensored=t(deriv_isOut_uncensored)
)
fit <- sampling(sm, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
post <- as.array(fit)

mcmc_dens_chains(post, regex_pars = "betas")+
  vline_at(c(.847, 1.672))
mcmc_intervals(post, regex_pars = "betas")+
  vline_at(c(.847, 1.672))
mcmc_dens_chains(post, regex_pars = "gammas")

# =====================================================================================================================

# c.f. https://cran.r-project.org/web/packages/simsurv/vignettes/simsurv_usage.html


sim_run <- function() {
  # Create a data frame with the subject IDs and treatment covariate
  cov <- data.frame(id = 1:200,
                    trt = rbinom(200, 1, 0.5))
  
  # Simulate the event times
  dat <- simsurv(lambdas = 0.1, 
                 gammas = 1.5, 
                 betas = c(trt = -0.5), 
                 x = cov, 
                 maxt = 5)
  
  # Merge the simulated event times onto covariate data frame
  dat <- merge(cov, dat)
  
  df <- as.tibble(dat)
  is_censored <- df$status==0
  times <- as.vector(df$eventtime)
  msk_censored <- is_censored == 1
  
  N <- nrow(df)
  X <- as.matrix(df[c("trt")])
  X_qr <- qr(X)
  Q <- qr.Q(X_qr)*sqrt(N-1)
  Q_censored <- as.matrix(Q[msk_censored,])
  Q_uncensored <- as.matrix(Q[!msk_censored,])
  R <- qr.R(X_qr)/sqrt(N-1)
  
  N_uncensored <- N-sum(msk_censored)
  N_censored <- sum(msk_censored)
  times_censored <- times[msk_censored]
  times_uncensored <- times[!msk_censored]
  
  
  ninterior_knots <- 3 # needs to be > 1
  knots <- quantile(times_uncensored,head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
  nknots <- length(knots)
  order<- 3
  isOut <- iSpline(times, knots = knots, degree = order-1, intercept = TRUE)
  deriv_isOut <- deriv(isOut)
  #isOut <- ibs(log_times, knots = knots, degree = order-1, intercept = TRUE)
  #deriv_isOut <- bSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
  isOut_censored <- isOut[msk_censored,]
  isOut_uncensored <- isOut[!msk_censored,]
  deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]
  nbasis <- dim(isOut)[2]
  
  stan_data  <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                     m=nbasis, Q_censored=Q_censored, Q_uncensored=Q_uncensored,R=R,
                     NC=ncol(X),
                     basis_evals_censored=t(isOut_censored), 
                     basis_evals_uncensored=t(isOut_uncensored),
                     deriv_basis_evals_uncensored=t(deriv_isOut_uncensored)
  )
  
  
  fit <- sampling(sm, data=stan_data, seed=42, chains=2, cores=2, iter=1000)
  
  
  # Obtain estimates, standard errors and 95% CI limits
  est <- summary(fit)$summary["betas[1]","mean"]
  ses <- summary(fit)$summary["betas[1]","sd"]
  cil <- summary(fit)$summary["betas[1]","2.5%"]
  ciu <- summary(fit)$summary["betas[1]","97.5%"]
  
  # Return bias and coverage indicator for treatment effect
  c(bias = est - (-0.5), 
    coverage = ((-0.5 > cil) && (-0.5 < ciu)))
}

# Set seed for simulations
set.seed(908070)
 
# Perform 100 replicates in simulation study (this might take some time)
rowMeans(replicate(100, sim_run()))