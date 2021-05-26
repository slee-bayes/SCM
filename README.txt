This repo provides the source codes the base case simulation study in Lee et al. (2021) "A Sequential Choice Model for Multiple Discrete Demand". 

(1) Run_SCM_202105.R
: This is the main R code that (i) simulates a panel dataset of multiple discrete demand by calling "Sim_SCM_202105.R" and (ii) estimates the sequential choice model by calling "Est_SCM_202105.R".

(2) Sim_SCM_202105.R
: This R code simulates a panel dataset of multiple discrete demand by using the sequential search process in Lee et al. (2021).

(3) Sim_OPT_202105.R
: This R code simulates a panel dataset of multiple discrete demand by using the exhaustive grid search, which gaurantees the optimality of the simulated data. This is used for evaluating the optimality of the sequential choice model.

(4) Est_SCM_202105.R
: This R code estimates the sequential choice model by a Bayesian MCMC procedure.
