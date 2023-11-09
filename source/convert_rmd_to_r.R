library(knitr)
# include "purl = FALSE" to exclude code chunks
purl("source/workshop_1_mcmc_en_brms.Rmd",
     output = "source/bayesian_statistics_1_script.R",
     documentation = 1 )
