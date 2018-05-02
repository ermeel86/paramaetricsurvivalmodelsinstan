/************************************************************************************************************************/
data {
    int<lower=1> N;                                                 // number of data points
    int<lower=1> m;                                                 // number of basis splines
    int<lower=1> ninterior_knots;                                   // number of interior knots
    ordered[ninterior_knots] knots;                                 // location of knots
    int<lower=1> NC;                                                // number of covariates
    matrix[N,NC] X;                                                 // design matrix
    int<lower=0, upper=1> is_censored[N];                           // delta in the paper
    vector[N] log_times;                                            // x=log(t) in the paper
    matrix[m,N] basis_evals;
    matrix[m,N] deriv_basis_evals;
}
/************************************************************************************************************************/
parameters {
    row_vector<lower=0>[m] gammas;                                  // regression coefficients for splines
    vector[NC] betas;                                               // regression coefficients for covariates
    real gamma_intercept;                                           // \gamma_0 in the paper
}
/************************************************************************************************************************/
model {
    vector[N] etas;
    gammas ~ normal(0, 2);
    betas ~ normal(0,1);
    gamma_intercept   ~ normal(0,1);

    etas = X*betas + (gammas*basis_evals)' + gamma_intercept;

    for(i in 1:N) {
        if(is_censored[i] == 0) {
          target +=  etas[i] - exp(etas[i]) - log_times[i] + log(gammas*deriv_basis_evals[,i]);
        }
        else {
           target += -exp(etas[i]);                                    // we work with log-likelihood
        }
   }
}
// TODO: vectorize target+= statements, which can be done. requires that we split data into censored and non-censored
/************************************************************************************************************************/
