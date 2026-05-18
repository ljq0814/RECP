Single change point test by multivariate ranks-based energy statistics
(RECP)
================

## Overview

`RECP` provides nonparametric methods for detecting and locating a
single change point in multivariate data sequences. This package
implements a rank energy statistic constructed from multivariate ranks
defined via optimal transport theory. The asymptotic null distribution
is exactly distribution-free. This property enables us to develop a
computationally feasible algorithm for computing universal rejection
thresholds across different sample sizes, which is applicable for long
data sequences. Two testing procedures are provided:

- `rank_energy_perm()`: Rank Energy Change Point Test via permutation
  (REp)
- `rank_energy_eig()`: Rank Energy Change Point Test via limiting
  spectral distribution (REs)

## Installation

You can install the development version of `RECP` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ljq0814/RECP")
```

## Examples

### Example 1: REp

``` r
set.seed(123)
library(RECP)
n <- 500; d <- 5; cploc <- 250
x <- rbind(matrix(rt(cploc*d,df = 1), cploc, d), matrix(rt((n-cploc)*d,df = 1)+0.25,n-cploc,d))
print(rank_energy_perm(x))
```

### Example 2: REs

``` r
set.seed(123)
library(RECP)
n <- 500; d <- 5; cploc <- 250
x <- rbind(matrix(rt(cploc*d,df = 1), cploc, d), matrix(rt((n-cploc)*d,df = 1)+0.25,n-cploc,d))
print(rank_energy_eig(x))
```

## Reference

Liu, Y., Li, J., Li, Z., Wang, J., & Wang, M. (2025). Testing the
distribution change in multivariate data using rank energy statistics.
*Journal of Applied Statistics*, 1–19.
<https://doi.org/10.1080/02664763.2025.2578653>
