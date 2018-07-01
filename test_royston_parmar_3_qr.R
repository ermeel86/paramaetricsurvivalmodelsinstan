library(rstan)
library(bayesplot)
library(splines2)
library(tidyverse)
library(survival)
sm_qr0 <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar_2_qr.stan")

sm_qr <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar_3_qr.stan")
sm <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar_2.stan")

#fifth example Northeren California Oncology Group  data
# c.f. https://web.stanford.edu/%7Ehastie/CASI_files/DATA/ncog.html

df <- as.tibble(read_delim("./data/ncog.txt", delim = " "))

fit_coxph <- coxph(Surv(t, d) ~ arm, data=df)


N <- nrow(df)
X <- as.matrix(as.integer(df$arm=="B"))
X_qr <- qr(X)
Q <- qr.Q(X_qr)*sqrt(N-1)
R <- qr.R(X_qr)/sqrt(N-1)
is_censored <- df$d==0
times <- as.vector(df$t)
log_times <- log(times)
msk_censored <- is_censored == 1

N_uncensored <- N-sum(msk_censored)
N_censored <- sum(msk_censored)
X_censored =  as.matrix(X[msk_censored,])
X_uncensored = as.matrix(X[!msk_censored,])
Q_censored =  as.matrix(Q[msk_censored,])
Q_uncensored = as.matrix(Q[!msk_censored,])
log_times_censored <- log_times[msk_censored]
log_times_uncensored <- log_times[!msk_censored]

ninterior_knots <- 3 # needs to be > 1
knots <- quantile(log_times,head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
nknots <- length(knots)
order<- 3
isOut <- iSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
deriv_isOut <- deriv(isOut)
#isOut <- ibs(log_times, knots = knots, degree = order-1, intercept = TRUE)
#deriv_isOut <- bSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
isOut_censored <- isOut[msk_censored,]
isOut_uncensored <- isOut[!msk_censored,]
deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]
nbasis <- dim(isOut)[2]

stan_data_qr <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                   m=nbasis, Q_censored=Q_censored,
                   X_censored=X_censored, X_uncensored=X_uncensored,
                   Q_uncensored=Q_uncensored, R=R,
                   log_times_censored=log_times_censored,
                   log_times_uncensored = log_times_uncensored,
                   NC=ncol(X),
                   basis_evals_censored=isOut_censored, 
                   basis_evals_uncensored=isOut_uncensored,
                   deriv_basis_evals_uncensored=deriv_isOut_uncensored
)

stan_data <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                     m=nbasis, Q_censored=Q_censored,
                     X_censored=X_censored, X_uncensored=X_uncensored,
                     Q_uncensored=Q_uncensored, R=R,
                     log_times_censored=log_times_censored,
                     log_times_uncensored = log_times_uncensored,
                     NC=ncol(X),
                     basis_evals_censored=t(isOut_censored), 
                     basis_evals_uncensored=t(isOut_uncensored),
                     deriv_basis_evals_uncensored=t(deriv_isOut_uncensored)
)
fit <- sampling(sm, data=stan_data, seed=42, chains=4, cores=2, iter=4000,control=list(adapt_delta=.99))
post <- as.array(fit)
mcmc_intervals(post, regex_pars = "betas")+vline_at(fit_coxph$coefficients)
mcmc_intervals(post,pars=c(sprintf("gammas[%d]", 1:nbasis)) )
#gammas_median <- summary(fit)$summary[sprintf("gammas[%d]", 1:nbasis), "50%"]


fitqr <- sampling(sm_qr, data=stan_data_qr, seed=42, chains=4, cores=2, iter=4000,control=list(adapt_delta=.99))
postqr <- as.array(fitqr)
mcmc_intervals(postqr, pars = c("betas[1]"))+vline_at(fit_coxph$coefficients)
mcmc_intervals(postqr,pars=c(sprintf("gammas[%d]", 1:nbasis)) )

gammas_medianqr <- summary(fitqr)$summary[sprintf("gammas[%d]", 1:nbasis), "50%"]


fitqr0 <- sampling(sm_qr0, data=stan_data, seed=42, chains=4, cores=2, iter=4000,control=list(adapt_delta=.99))
postqr0 <- as.array(fitqr0)
mcmc_intervals(postqr0, pars = c("betas[1]"))+vline_at(fit_coxph$coefficients)
mcmc_intervals(postqr0,pars=c(sprintf("gammas[%d]", 1:nbasis)) )

gammas_medianqr0 <- summary(fitqr0)$summary[sprintf("gammas[%d]", 1:nbasis), "50%"]


