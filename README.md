# Parametric survival models in Stan

## I-spline based implementation of Royston-Parmar (proportional hazard model)

- https://www.ncbi.nlm.nih.gov/pubmed/12210632
- https://projecteuclid.org/euclid.ss/1177012761

Uses a QR decomposition of the covariates.

See also http://discourse.mc-stan.org/t/survival-models-in-rstanarm

## Parametric hazard model that uses M-splines for the baseline hazard function

We use M-splines (in time) to parametrise the baseline hazard. This is equivalent to parametrising the baseline cumulative hazard function by
I-splines (in time). In the language of the paper below, we use a link function that equals -log.

The I-splines have the nice property that I(0) = 0, which suffices to guarantee that S(0)=exp(I(0))^exp(beta*X) =1

- See [Parametric and penalized generalized survival models](https://www.ncbi.nlm.nih.gov/pubmed/27587596)

Uses a QR decomposition of the covariates.


See also http://discourse.mc-stan.org/t/survival-models-in-rstanarm

## Todo
- Think whether it would be useful to add an intercept to M(t), equivalently replacing I(t) --> M(t) + intercept*t
- Clean-up tests and make them more systematic
- Work with NAs
- Implement the proportional odds model
- Work with  random-walk priors for $\gamma$ in order to allow for more knots, c.f. http://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
- Implement time-varying hazard and odds ratios (requires to solve the problem of constraining the spline coefficients to be non-negative)
