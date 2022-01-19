##### A Sequential Choice Model for Multiple Discrete Demand (Lee et al. (2022))
##### Main Code
##### last update 2022.01
rm(list=ls(all=TRUE))

library(bayesm)
library(mvtnorm)

getwd() 



##### Simulation Study #########################################################

### Data Generation

# scale
n_hh=500	# number of households
n_obs=30	# number of observations per household
n_br=2      # number of inside goods

# unit prices
pr_min=rep(2,n_br); pr_max=rep(3,n_br); pr_setting=list(pr_min=pr_min,pr_max=pr_max)

# fixed parameters
psi0=1; offset=c(1,1); e_sd=1
fixparam=list(psi0=psi0,offset=offset,e_sd=e_sd)

# hyper-parameters to be estimated
tpsi=c(0.05,0.07); tpsi_1=log(tpsi); # re-parameterization for 0<psi
tgamma=1; tgamma_1=log(tgamma); # re-parameterization for 0<gamma
tE=40; tE_1=log(tE); # re-parameterization for 0<E

tparam=c(tpsi_1,tgamma_1,tE_1); n_param=length(tparam); n_param
ind_psi=c(1:n_br);  ind_gamma=max(ind_psi)+c(1:length(tgamma)); ind_E=max(ind_gamma)+1; c(ind_psi,ind_gamma,ind_E)

mu_theta=tparam; mu_theta # mean of the heterogeneity distribution
sd_theta=rep(1,n_param); sd_theta # std.dev. of the heterogeneity distribution

### Data Simulation
seed_value=207 # seed value for random generation

# simulation by the sequential choice model
source("Sim_SCM_202105.R") 
set.seed(seed_value)
SIMDATA_SCM=sim_scm(n_hh,n_obs,mu_theta,sd_theta,fixparam,pr_setting)
round(mu_theta-colMeans(SIMDATA_SCM$IndParam),3) 

# simulation by exhaustive search (guarantees the optimal solution)
source("Sim_OPT_202105.r") 
set.seed(seed_value)
SIMDATA_OPT=sim_opt(n_hh,n_obs,mu_theta,sd_theta,fixparam,pr_setting)
round(mu_theta-colMeans(SIMDATA_OPT$IndParam),3)

# compare the two simulated data for evaluating the optimality of the sequential choice model
Compare_byhh=rep(0,n_hh)
for(hh in 1:n_hh){
    Compare_byhh[hh]=mean((rowSums(SIMDATA_SCM$hierdata[[hh]]$X==SIMDATA_OPT$hierdata[[hh]]$X)==n_br))
} # end of for(hh in 1:n_hh)
round(Compare_byhh,3); mean(Compare_byhh)

SIMDATA=SIMDATA_SCM
tIndParam=SIMDATA$IndParam





##### Estimation

### Data
Z=SIMDATA$X_hh; dim(Z) # covariates for mean of distribution of heterogeneity
n_z=ncol(Z); n_z
true=list(mu_theta=mu_theta,sd_theta=sd_theta,tIndParam=tIndParam)
Data=list(hierdata=SIMDATA$hierdata,Z=Z,fixparam=fixparam,true=true)

n_hh=length(Data$hierdata); n_hh
n_obs=nrow(Data$hierdata[[1]]$X); n_obs
n_param=length(tparam); n_param



### Prior
ind_psi; ind_gamma; ind_E
ind_est=c(ind_psi,ind_gamma,ind_E); ind_est; n_param_est=length(ind_est); n_param_est
ind_fix=setdiff(c(1:n_param),ind_est); ind_fix; n_param_fix=length(ind_fix); n_param_fix

n_var=n_param_est; n_var
nu=n_var+3
V=nu*diag(n_var)
Deltabar=matrix(0,nrow=n_z,ncol=n_var,byrow=TRUE)
if(n_z==1) {A=matrix(.01,ncol=1,nrow=1)}
if(n_z>1) {A=diag(rep(.01,n_z))}

Prior=list(nu=nu,V=V,Deltabar=Deltabar,A=A)



### Mcmc
n_edraw=5 # sampling size for the Monte Carlo integration

# initial values
max_spending=rep(0,n_hh)
for(hh in 1:n_hh){max_spending[hh]=max(rowSums(Data$hierdata[[hh]]$X*Data$hierdata[[hh]]$PR))}

init_psi_1=rep(0,n_br); init_theta=c(init_psi_1)
init_gamma_1=rep(0,length(ind_gamma)); init_theta=c(init_theta,init_gamma_1)
init_E_1=log(max(max_spending)+10); init_theta=c(init_theta,init_E_1)
length(init_theta)==n_param # must be true

Thetas0=matrix(init_theta,ncol=n_param,nrow=n_hh,byrow=TRUE); 
Vtheta0=1*diag(n_var)
Delta0=matrix(init_theta[ind_est],ncol=n_var,nrow=n_z,byrow=TRUE)
if(any(ind_fix)){Thetas0[,ind_fix]=tIndParam[,ind_fix]}

# mcmc draw settings
ind_gr1=c(ind_psi,ind_gamma,ind_E); all(ind_gr1%in%ind_est) # Must be TRUE
ind_gr2=NULL; all(ind_gr2%in%ind_est) # Must be TRUE
ind_gr3=NULL; all(ind_gr3%in%ind_est) # Must be TRUE
ind_est_group=list(ind_gr1=ind_gr1,ind_gr2=ind_gr2,ind_gr3=ind_gr3); ind_est_group

# stepsize
ss_gr1=ss_gr2=ss_gr3=0
if(is.null(ind_gr1)==FALSE) {ss_gr1=0.2}; if(is.null(ind_gr2)==FALSE) {ss_gr2=0.2}; if(is.null(ind_gr3)==FALSE) {ss_gr3=0.2}
ss_all=c(ss_gr1,ss_gr2,ss_gr3)

est_setting=list(n_edraw=n_edraw,ss_all=ss_all,ind_est_group=ind_est_group,ind_fix=ind_fix)

R=5000; keep=1; n_report=50; 
Mcmc=list(R=R,keep=keep,n_report=n_report,Thetas0=Thetas0,Vtheta0=Vtheta0,Delta0=Delta0,est_setting=est_setting)

source("Est_SCM_202105.R")
result=est_scm(Data=Data,Prior=Prior,Mcmc=Mcmc)





##### Results ##################################################################################

R2=length(result$Sll); R2  
ind_all=seq(from=1,to=R2,by=1); length(ind_all)
keep_thin=5; keep_thin_1=keep_thin/Mcmc$keep; keep_thin_1
ind_all_thin=seq(from=1,to=R2,by=keep_thin_1); ind_all_thin[1:10]; length(ind_all_thin)
R1=R2/4+1; R1  # burn-in 
ind_valid=seq(from=R1,to=R2,by=1)
ind_valid_thin=seq(from=R1,to=R2,by=keep_thin_1); ind_valid_thin[1:10]; length(ind_valid_thin)

# (1) Likelihood plot
matplot(result$Sll[ind_all_thin],type='l')
matplot(result$Sll[ind_valid_thin],type='l')
logMargDenNR(result$Sll[ind_valid_thin])

# (2) Hyper-parameters
Delta_draw=as.matrix(result$Deltadraw[,n_z*c(1:n_param_est)-(n_z-1)])
matplot(Delta_draw[ind_all_thin,],type='l',main='theta_bar',col=c(1:n_param_est))
abline(h=c(tparam[ind_est]),col=c(1:n_param_est))
mat1=apply(as.matrix(Delta_draw[ind_valid_thin,]),2,quantile,probs=c(0.025,0.5,0.975))
mat1=rbind(tparam[ind_est],colMeans(as.matrix(tIndParam[,ind_est])),mat1); 
rownames(mat1)[1]="true_theta_bar"; rownames(mat1)[2]="average_theta_bar"
print(mat1)

