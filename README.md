# longpower: Sample size calculations for longitudinal data

The longpower package contains functions for computing power and sample size for linear models of longitudinal data based on the formula due to Liu and Liang (1997) and Diggle et al (2002). Either formula is expressed in terms of marginal model or Generalized Estimating Equations (GEE) parameters. This package contains functions which translate pilot mixed effect model parameters (e.g. random intercept and/or slope) into marginal model parameters so that the formulas of Diggle et al or Liu and Liang formula can be applied to produce sample size calculations for two sample longitudinal designs assuming known variance.

## To install via CRAN:

```r
install.packages("longpower")
```

## To install via bitbucket:

```r
# install.packages('devtools')
devtools::install_bitbucket("longpower", "mdonohue")
```
