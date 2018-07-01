/************************************************************************************************************************/
data {
    int<lower=0> N_uncensored;                                      // number of uncensored data points
    int<lower=0> N_censored;                                        // number of censored data points
    int<lower=1> m;                                                 // number of basis splines
    int<lower=1> NC;                                                // number of covariates
    matrix[N_censored,NC] Q_censored;                               // Q-transf. of design matrix (censored)
    matrix[N_uncensored,NC] Q_uncensored;                           // Q-transf. of design matrix (uncensored)
    matrix[NC, NC] R;
    vector[N_censored] log_times_censored;                          // x=log(t) in the paper (censored)
    vector[N_uncensored] log_times_uncensored;                      // x=log(t) in the paper (uncensored)
    matrix[N_censored,m] basis_evals_censored;                      // ispline basis matrix (censored)
    matrix[N_uncensored,m] basis_evals_uncensored;                  // ispline basis matrix (uncensored)
    matrix[N_uncensored,m] deriv_basis_evals_uncensored;            // derivatives of isplines matrix (uncensored)
}
/************************************************************************************************************************/
transformed data {
    matrix[NC,NC] R_inv = inverse(R);
}
/************************************************************************************************************************/
parameters {
    vector<lower=0>[m] gammas;                                  // regression coefficients for splines
    vector[NC] betas_tr;                                               // regression coefficients for covariates
    real gamma_intercept;                                           // \gamma_0 in the paper
}
/************************************************************************************************************************/
transformed parameters {
    vector[NC] betas = R_inv * betas_tr;
}
/************************************************************************************************************************/
model {
    vector[N_censored] etas_censored;
    vector[N_uncensored] etas_uncensored;
    gammas ~ normal(0, 1);
    betas ~ normal(0,1);
    gamma_intercept   ~ normal(0,1);
    
    etas_censored = Q_censored*betas_tr + basis_evals_censored*gammas  + gamma_intercept;
    etas_uncensored = Q_uncensored*betas_tr + basis_evals_uncensored*gammas  + gamma_intercept;
    
    //etas_censored = Q_censored*betas_tr + (gammas*basis_evals_censored)' + gamma_intercept;
    //etas_uncensored = Q_uncensored*betas_tr + (gammas*basis_evals_uncensored)' + gamma_intercept;
    target += -exp(etas_censored);
    target += etas_uncensored - exp(etas_uncensored) - log_times_uncensored + log(deriv_basis_evals_uncensored*gammas);
}
/************************************************************************************************************************/
