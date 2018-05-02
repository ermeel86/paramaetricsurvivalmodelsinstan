library(rstan)
library(flexsurv)
library(bayesplot)
library(splines2)
library(tidyverse)
sm <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar.stan")
sm2 <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar_2.stan")

# get the breast cancer survival data set used by Royston & Parmar
# see https://cran.r-project.org/web/packages/flexsurv/flexsurv.pdf#Rfn.bc.1
df <- flexsurv::bc
df_censored <- filter(df, censrec==0)
df_uncensored <- filter(df, censrec==1)
N <- nrow(df)
X <- as.matrix(cbind(as.integer(df$group == "Medium"), as.integer(df$group=="Poor")))
is_censored <- as.vector(as.integer(df$censrec == 0)) # censrec==1 means dead, 0 means censored in the original data
times <- as.vector(df$rectime)
log_times <- log(times)

N_uncensored <- nrow(df_uncensored)
N_censored <- nrow(df_censored)
X_censored = as.matrix(cbind(as.integer(df_censored$group == "Medium"), as.integer(df_censored$group=="Poor")))
X_uncensored = as.matrix(cbind(as.integer(df_uncensored$group == "Medium"), as.integer(df_uncensored$group=="Poor")))
log_times_censored <- log(as.vector(df_censored$rectime))
log_times_uncensored <- log(as.vector(df_uncensored$rectime))

nknots <- 2 # needs to be > 1
knots <- quantile(log_times,head(tail(seq(0,1, length.out = nknots+2),-1),-1))
nknots <- length(knots)
order<- 3
isOut <- iSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
deriv_isOut <- deriv(isOut)
msk_censored <- is_censored == 1
isOut_censored <- isOut[msk_censored,]
isOut_uncensored <- isOut[!msk_censored,]
deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]

############
ggplot(mapping=aes(x=x, y=y))+
  geom_line(data=tibble(x=log_times, y=isOut[,1]))+
  geom_line(data=tibble(x=log_times, y=isOut[,2]))+
  geom_line(data=tibble(x=log_times, y=isOut[,3]))+
  geom_line(data=tibble(x=log_times, y=isOut[,4]))+
  geom_line(data=tibble(x=log_times, y=isOut[,5]))
isOut2 <- iSpline(log_times, knots = knots, degree = order-1, intercept = FALSE)
ggplot(mapping=aes(x=x, y=y))+
  geom_line(data=tibble(x=log_times, y=isOut2[,1]))+
  geom_line(data=tibble(x=log_times, y=isOut2[,2]))+
  geom_line(data=tibble(x=log_times, y=isOut2[,3]))+
  geom_line(data=tibble(x=log_times, y=isOut2[,4]))
############
nbasis <- dim(isOut)[2]
stan_data <- list(N=N, m=nbasis, X=X, is_censored=is_censored, log_times=log_times, knots=knots, NC=ncol(X),
                  ninterior_knots=nknots,
                  basis_evals=t(isOut), deriv_basis_evals=t(deriv_isOut)
                  )
stan_data2 <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                   m=nbasis, X_censored=X_censored, X_uncensored=X_uncensored, 
                   log_times_censored=log_times_censored,
                   log_times_uncensored = log_times_uncensored,
                   knots=knots, NC=ncol(X),
                   ninterior_knots=nknots,
                   basis_evals_censored=t(isOut_censored), 
                   basis_evals_uncensored=t(isOut_uncensored),
                   deriv_basis_evals_uncensored=t(deriv_isOut_uncensored)
)
fit2 <- sampling(sm2, data=stan_data2, seed=42, chains=4, cores=2, iter=4000,control=list(adapt_delta=.9))
post2 <- as.array(fit2)
fit <- sampling(sm, data=stan_data, seed=42, chains=4, cores=1, iter=4000, control=list(adapt_delta=.99))
post <- as.array(fit)
mcmc_dens_chains(post2, regex_pars = "betas")+
  vline_at(c(.847, 1.672))
mcmc_intervals(post2, regex_pars = "betas")+
  vline_at(c(.847, 1.672))
mcmc_dens_chains(post2, regex_pars = "gammas")
mcmc_dens_chains(post, regex_pars = "gamma_intercept")
#########################
gammas <- summary(fit2)$summary[sprintf("gammas[%d]", 1:nbasis), "50%"]
gamma_icpt <- summary(fit2)$summary["gamma_intercept","50%"]
betas  <-  summary(fit2)$summary[sprintf("betas[%d]", 1:2), "50%"]
ss <- as.vector(isOut %*% gammas + X %*% betas +gamma_icpt)
Ss <- exp(-exp(ss))
# survival curves
ggplot(data=tibble(x=exp(log_times)/365, y=Ss), aes(x=x, y=y))+
  geom_point()
#Hz <- as.vector(exp(ss))
#ggplot(data=tibble(x=exp(log_times)/365, y=Hz), aes(x=x, y=y))+
#  geom_point()





