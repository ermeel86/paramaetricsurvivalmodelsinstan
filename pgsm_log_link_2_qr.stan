/************************************************************************************************************************/
data {
    int<lower=0> N_uncensored;                                      // number of uncensored data points
    int<lower=0> N_censored;                                        // number of censored data points
    int<lower=1> m;                                                 // number of basis splines
    int<lower=1> NC;                                                // number of covariates
    matrix[N_censored,NC] Q_censored;                               // Q-transf. of design matrix (censored)
    matrix[N_uncensored,NC] Q_uncensored;                           // Q-transf. of design matrix (uncensored)
    matrix[NC, NC] R;
    matrix[m,N_censored] basis_evals_censored;                      // ispline basis matrix (censored)
    matrix[m,N_uncensored] basis_evals_uncensored;                  // ispline basis matrix (uncensored)
    matrix[m,N_uncensored] deriv_basis_evals_uncensored;            // derivatives of isplines matrix (uncensored)
}
/************************************************************************************************************************/
transformed data {
    matrix[NC,NC] R_inv = inverse(R);
}
/************************************************************************************************************************/
parameters {
    row_vector<lower=0>[m] gammas;                                  // regression coefficients for splines
    vector[NC] betas_tr;                                            // regression coefficients for covariates
    //real gamma_intercept;                                         // \gamma_0 in the paper
}
/************************************************************************************************************************/
transformed parameters {
    vector[NC] betas = R_inv * betas_tr;
}
/************************************************************************************************************************/
model {
    gammas ~ normal(0, 2);
    betas ~ normal(0,1);
    //gamma_intercept   ~ normal(0,1);                              // we use no intercept on the I-splines, since this would lead to identifiability problems
    
    target += -exp(Q_censored*betas_tr).*(gammas*basis_evals_censored)';
    target += -exp(Q_uncensored*betas_tr).*(gammas*basis_evals_uncensored)';
    target +=  Q_uncensored*betas_tr + log(gammas*deriv_basis_evals_uncensored)';
}
/************************************************************************************************************************/
