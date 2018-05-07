library(rstan)
library(flexsurv)
library(bayesplot)
library(splines2)
library(tidyverse)
# =====================================================================================================================
sm <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar.stan")
sm2 <- stan_model("~/Desktop/Stan/Parametric_Survival_Models/royston_parmar_2.stan")
# =====================================================================================================================
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

ninterior_knots <- 2 # needs to be > 1
knots <- quantile(log_times,head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
nknots <- length(knots)
order<- 3
isOut <- iSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
deriv_isOut <- deriv(isOut)
msk_censored <- is_censored == 1
isOut_censored <- isOut[msk_censored,]
isOut_uncensored <- isOut[!msk_censored,]
deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]
# =====================================================================================================================
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
# =====================================================================================================================
nbasis <- dim(isOut)[2]
stan_data <- list(N=N, m=nbasis, X=X, is_censored=is_censored, log_times=log_times, NC=ncol(X),
                  basis_evals=t(isOut), deriv_basis_evals=t(deriv_isOut)
                  )
stan_data2 <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                   m=nbasis, X_censored=X_censored, X_uncensored=X_uncensored, 
                   log_times_censored=log_times_censored,
                   log_times_uncensored = log_times_uncensored,
                   NC=ncol(X),
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
# =====================================================================================================================
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

# =====================================================================================================================

# second example https://web.stanford.edu/~hastie/CASI_files/DATA/pediatric.html

df_pediatric <- read_delim("~/Desktop/Stan/Parametric_Survival_Models/data/pediatric.txt", delim = " ")


df_censored <- filter(df_pediatric, d==0)
df_uncensored <- filter(df_pediatric, d==1)
N <- nrow(df_pediatric)
X <- as.matrix(df_pediatric[c("sex","race","age","entry","far")])
X[,3:5] <- scale(X[,3:5])
is_censored <- as.vector(as.integer(df_pediatric$d == 0)) # censrec==1 means dead, 0 means censored in the original data
times <- as.vector(df_pediatric$t)
log_times <- log(times)
msk_censored <- is_censored == 1

N_uncensored <- nrow(df_uncensored)
N_censored <- nrow(df_censored)
X_censored =  X[msk_censored,]
X_uncensored = X[!msk_censored,]
log_times_censored <- log(as.vector(df_censored$t))
log_times_uncensored <- log(as.vector(df_uncensored$t))

ninterior_knots <- 2 # needs to be > 1
knots <- quantile(log_times,head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
nknots <- length(knots)
order<- 3
isOut <- iSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
deriv_isOut <- deriv(isOut)
isOut_censored <- isOut[msk_censored,]
isOut_uncensored <- isOut[!msk_censored,]
deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]
nbasis <- dim(isOut)[2]

stan_data2 <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                   m=nbasis, X_censored=X_censored, X_uncensored=X_uncensored, 
                   log_times_censored=log_times_censored,
                   log_times_uncensored = log_times_uncensored,
                   NC=ncol(X),
                   basis_evals_censored=t(isOut_censored), 
                   basis_evals_uncensored=t(isOut_uncensored),
                   deriv_basis_evals_uncensored=t(deriv_isOut_uncensored)
)
fit2 <- sampling(sm2, data=stan_data2, seed=42, chains=4, cores=2, iter=4000,control=list(adapt_delta=.9))
post2 <- as.array(fit2)
mcmc_dens_chains(post2, regex_pars = "betas")
# values below copied from table 9.7 in https://web.stanford.edu/~hastie/CASI/index.html
mcmc_intervals(post2, regex_pars = "betas")+vline_at(c(-.023,.282, -.235, -.460,.296 ))
mcmc_dens_chains(post2, regex_pars = "gammas")
mcmc_dens_chains(post2, regex_pars = "gamma_intercept")

# =====================================================================================================================
# third example http://www.openbugs.net/Examples/Leuk.html

leuk = read.table("http://www.karlin.mff.cuni.cz/~pesta/prednasky/NMFM404/Data/leuk.dat", sep=" ", head=T, skip=7) # read data
head(leuk)


N <- nrow(leuk)
X <- leuk$trtmt
X <- as.matrix(scale(leuk$trtmt,scale=FALSE))
is_censored <- leuk$remiss==0
times <- as.vector(leuk$weeks)
log_times <- log(times)
msk_censored <- is_censored == 1

N_uncensored <- N-sum(msk_censored)
N_censored <- sum(msk_censored)
X_censored =  as.matrix(X[msk_censored,])
X_uncensored = as.matrix(X[!msk_censored,])
log_times_censored <- log_times[msk_censored]
log_times_uncensored <- log_times[!msk_censored]

ninterior_knots <- 3 # needs to be > 1
knots <- quantile(log_times,head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
nknots <- length(knots)
order<- 3
isOut <- iSpline(log_times, knots = knots, degree = order-1, intercept = TRUE)
deriv_isOut <- deriv(isOut)
isOut_censored <- isOut[msk_censored,]
isOut_uncensored <- isOut[!msk_censored,]
deriv_isOut_uncensored <- deriv_isOut[!msk_censored,]
nbasis <- dim(isOut)[2]

stan_data2 <- list(N_uncensored=N_uncensored, N_censored=N_censored, 
                   m=nbasis, X_censored=X_censored, X_uncensored=X_uncensored, 
                   log_times_censored=log_times_censored,
                   log_times_uncensored = log_times_uncensored,
                    NC=ncol(X),
                   basis_evals_censored=t(isOut_censored), 
                   basis_evals_uncensored=t(isOut_uncensored),
                   deriv_basis_evals_uncensored=t(deriv_isOut_uncensored)
)
fit2 <- sampling(sm2, data=stan_data2, seed=42, chains=4, cores=2, iter=4000,control=list(adapt_delta=.99))
post2 <- as.array(fit2)
mcmc_dens_chains(post2, regex_pars = "betas")
# reference values from http://www.karlin.mff.cuni.cz/~pesta/NMFM404/ph.html#Assessing_the_PH_assumption
# they come from coxph with different methods to resolve ties
mcmc_intervals(post2, regex_pars = "betas") + 
  vline_at(c(-1.5721,-1.5092,-1.6282))
mcmc_dens_chains(post2, regex_pars = "gammas")
mcmc_dens_chains(post2, regex_pars = "gamma_intercept")

isOut_test <- iSpline(log(leuk_data$t), knots = knots, degree = order-1, intercept = TRUE)
gammas_median <- summary(fit2)$summary[sprintf("gammas[%d]",1:nbasis), "50%"]
gamma_icpt <- summary(fit2)$summary["gamma_intercept","50%"]
betas  <-  summary(fit2)$summary["betas[1]", "50%"]

ss_treat <- as.vector(isOut_test %*% gammas_median + 0.5*betas + gamma_icpt)#treatment
ss_placebo <- as.vector(isOut_test %*% gammas_median- 0.5*betas+ gamma_icpt)#placebo

Ss_treat <- exp(-exp(ss_treat))
Ss_placebo <- exp(-exp(ss_placebo))


# survival curves
ggplot(data=tibble(x=leuk_data$t, y1=Ss_treat, y2=Ss_placebo))+
  geom_line(aes(x=x, y=y1), color='green')+ 
  geom_point(aes(x=x, y=y1), color='green')+
  geom_line(aes(x=x, y=y2), color='blue')+
  geom_point(aes(x=x, y=y2), color='blue')

  
