# Parametric survival models in Stan

## I-spline based implementation of Royston-Parmar (proportional hazard model)

- https://www.ncbi.nlm.nih.gov/pubmed/12210632
- https://projecteuclid.org/euclid.ss/1177012761


See also http://discourse.mc-stan.org/t/survival-models-in-rstanarm

## Todo

- Try to reproduce further examples from the literature
    - https://web.stanford.edu/~hastie/CASI_files/DATA/ncog.html
    - https://web.stanford.edu/~hastie/CASI_files/DATA/pediatric.html
- Work with NAs
- Implement the proportional odds model
- Work with (smooth) random-walk priors for $\gamma$ in order to allow for more knots, c.f. http://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
- Implement time-varying hazard and odds ratios
