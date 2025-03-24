# Introduction
This package is the codes for multivariate single change point test by multivariate ranks-based energy statistics

# Main Function
The core function of the package is rank_energy_stat. It performs rank energy statistics.

# Example
x <- datageneration(n=500,d=50,cploc=250,delta=0.05,existcp = TRUE, distri = "Case1")

re_res <- rank_energy_stat(x,499)

print(re_res);
