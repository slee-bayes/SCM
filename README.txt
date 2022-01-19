This repo provides the source codes for the base case simulation study in Lee et al. (2022) "A Sequential Choice Model for Multiple Discrete Demand". 

(1) Run_SCM.R
: This is the main R code that (i) simulates a panel dataset of multiple discrete demand by calling "Sim_SCM.R" and (ii) estimates the sequential choice model by calling "Est_SCM.R".

(2) Sim_SCM.R
: This R code simulates a panel dataset of multiple discrete demand by using the sequential search process in Lee et al. (2022).

(3) Sim_OPT.R
: This R code simulates a panel dataset of multiple discrete demand by using an exhaustive grid search, which gaurantees the optimality of the simulated data. This is used for evaluating the optimality of the sequential choice model.

(4) Est_SCM.R
: This R code estimates the sequential choice model by a Bayesian MCMC procedure.
