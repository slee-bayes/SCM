##### A Sequential Choice Model for Multiple Discrete Demand (Lee et al. (2022))
##### Estimation Code
##### last update 2022.01
  
est_scm=function(Data,Prior,Mcmc) {
 
### Required Functions ---------------------------------------------------------
pandterm = function(message) {stop(message, call. = FALSE)}

# function for computing the loglikelihood given stacked data (full vectorization across obs across hhs)
LL_disc_stack=function(X_stack,PR_stack,ID_stack,Thetas_stack,fixparam,n_edraw){
    psi0=fixparam$psi0; e_sd=fixparam$e_sd
    offset=fixparam$offset; os_in=offset[1]; os_out=offset[2]
    
    X=X_stack; PR=PR_stack # n_obs*n_hh x n_br
    n_br=ncol(X); n_obs=nrow(X); 
    
    Psi=matrix(exp(Thetas_stack[,ind_psi]),nrow=n_obs,ncol=length(ind_psi))
    Gamma=matrix(exp(Thetas_stack[,ind_gamma]),nrow=n_obs,ncol=length(ind_gamma)); if(length(ind_gamma)==1) {Gamma=matrix(Gamma,nrow=n_obs,ncol=n_br,byrow=FALSE)}
    E=matrix(exp(Thetas_stack[,ind_E]),nrow=n_obs,ncol=length(ind_E)); EE=matrix(E,nrow=n_obs,ncol=n_br,byrow=FALSE)

    Totexp=rowSums(PR*X); Z_stack=E-Totexp
    ZZ=matrix(Z_stack,nrow=n_obs,ncol=n_br,byrow=FALSE)
    
    ll_vec=rep(-Inf,n_obs)
    IND_pos=(X>0)
    IND_zero=(X==0)
    N_pos=rowSums(IND_pos)
    ind_np=which(N_pos==0 & Z_stack>=0); n_obs_np=length(ind_np)
    ind_p1up=which(N_pos>0 & Z_stack>=0); n_obs_p1up=length(ind_p1up)
    
    LB=matrix(-Inf,nrow=n_obs,ncol=n_br) # lower bound
    UB=matrix(Inf,nrow=n_obs,ncol=n_br) # upper bound
    
    k=1
    ind_UB=(ZZ-k*PR>=0)
    num_1=(ZZ+os_out-k*PR)/(ZZ+os_out); num_1[num_1<0]=0; num=-psi0*log(num_1) # numerator
    denom=log((Gamma*(X+k*1)+os_in)/(Gamma*X+os_in))/Gamma # denominator
    uuu=num/denom 
    uuu[uuu<0]=0; UB[ind_UB]=(log(uuu)-log(Psi))[ind_UB]
    
    k=-1
    ind_LB=(X>0)
    num_1=(ZZ+os_out-k*PR)/(ZZ+os_out); num_1[num_1<0]=0; num=-psi0*log(num_1) # numerator
    denom_1=(Gamma*(X+k*1)+os_in)/(Gamma*X+os_in); denom_1[denom_1<0]=0; denom=log(denom_1)/Gamma # denominator
    lll=num/denom 
    lll[lll<0]=0; LB[ind_LB]=(log(lll)-log(Psi))[ind_LB]
    
    ind_check=which(rowSums(UB>LB)<n_br | Z_stack<0)
    if(length(ind_check)>0){UB[ind_check,]=LB[ind_check,]=0}
    
    ll_exit=rowSums(log(pnorm(UB,0,e_sd)-pnorm(LB,0,e_sd))) # log-likelihood for satisfying the exit condition
    ll_exit[Z_stack<0]=-Inf
    ll_vec=ll_exit
    
    if(n_obs_p1up>0 & n_edraw>0){
        X_p1up=matrix(X[ind_p1up,],nrow=n_obs_p1up); PR_p1up=matrix(PR[ind_p1up,],nrow=n_obs_p1up) # n_obs_p1up*n_hh x n_br
        LB_p1up=matrix(LB[ind_p1up,],nrow=n_obs_p1up); UB_p1up=matrix(UB[ind_p1up,],nrow=n_obs_p1up)
        Psi_p1up=matrix(Psi[ind_p1up,],nrow=n_obs_p1up); Gamma_p1up=matrix(Gamma[ind_p1up,],nrow=n_obs_p1up); EE_p1up=matrix(EE[ind_p1up,],nrow=n_obs_p1up);
        
        n_obs_p1up_stack_all=(n_obs_p1up*n_edraw)
        
        X_p1up_stack=matrix(t(X_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE) # n_obs_p1up*n_hh*n_edraw x n_br
        PR_p1up_stack=matrix(t(PR_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE)
        LB_p1up_stack=matrix(t(LB_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE)
        UB_p1up_stack=matrix(t(UB_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE)
        
        Psi_p1up_stack=matrix(t(Psi_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE)
        Gamma_p1up_stack=matrix(t(Gamma_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE)
        EE_p1up_stack=matrix(t(EE_p1up),nrow=n_obs_p1up_stack_all,ncol=n_br,byrow=TRUE)

        Error_stack_all=matrix(rtrun(mu=rep(0,(n_obs_p1up_stack_all*n_br)),sigma=rep(e_sd,(n_obs_p1up_stack_all*n_br)),a=LB_p1up_stack,b=UB_p1up_stack),nrow=n_obs_p1up_stack_all,ncol=n_br) # n_obs_p1up*n_hh*n_edraw x n_br
        X_stack_all=matrix(0,nrow=n_obs_p1up_stack_all,ncol=n_br)
        
        # 1. starting from all zeros
        # 2. check the bang for the buck for each item and buy or remove one unit of the alternative whose bang for the buck is the best
        # 3. repeat 2. until (i)there is no more alternative with a positive bang for the buck or (ii)constrained by the budget
        n_repeat=1; ind_loop=c(1:n_obs_p1up_stack_all); 
        X_stack=X_stack_all; Error_stack=Error_stack_all
        
        repeat{
            n_obs_p1up_stack=length(ind_loop)
            ZZ_p1up_stack=EE_p1up_stack-matrix(rowSums(PR_p1up_stack*X_stack),nrow=n_obs_p1up_stack,ncol=n_br,byrow=FALSE)
            
            k=1
            bfb_buy=matrix(-Inf,ncol=n_br,nrow=n_obs_p1up_stack)
            ind=((ZZ_p1up_stack-k*PR_p1up_stack)>=0)
            if(any(ind==TRUE)){
                bfb_buy[ind]= (Psi_p1up_stack*exp(Error_stack)/Gamma_p1up_stack)[ind]*log((Gamma_p1up_stack*(X_stack+k)+os_in)[ind]/(Gamma_p1up_stack*X_stack+os_in)[ind]) + psi0*log(((ZZ_p1up_stack-PR_p1up_stack*k+os_out)/(ZZ_p1up_stack+os_out))[ind])
            }
            
            k=-1
            bfb_remove=matrix(-Inf,ncol=n_br,nrow=n_obs_p1up_stack)
            ind=(X_stack>0)
            if(length(ind)>0){
                bfb_remove[ind]= (Psi_p1up_stack*exp(Error_stack)/Gamma_p1up_stack)[ind]*log((Gamma_p1up_stack*(X_stack+k)+os_in)[ind]/(Gamma_p1up_stack*X_stack+os_in)[ind]) + psi0*log(((ZZ_p1up_stack-PR_p1up_stack*k+os_out)/(ZZ_p1up_stack+os_out))[ind])
            }
            
            bfb=cbind(bfb_buy,bfb_remove,rep(0,n_obs_p1up_stack))
            bfb_max=apply(bfb,1,max); 
            ind_max=apply(bfb,1,which.max); # bfb[cbind(seq_along(ind_max), ind_max)]
            
            ind_quit=which(ind_max==(n_br*2+1))
            if(length(ind_quit)>0){
                X_stack_all[ind_loop[ind_quit],]=X_stack[ind_quit,]
                ind_loop=ind_loop[-ind_quit]; if(length(ind_loop)==0) break;
                
                X_stack=X_stack[-ind_quit,]
                PR_p1up_stack=PR_p1up_stack[-ind_quit,]
                Psi_p1up_stack=Psi_p1up_stack[-ind_quit,]
                Gamma_p1up_stack=Gamma_p1up_stack[-ind_quit,]
                EE_p1up_stack=EE_p1up_stack[-ind_quit,]
                Error_stack=Error_stack[-ind_quit,]
                ind_max=ind_max[-ind_quit]
            } # end of if(length(ind_quit)>0)
            
            if(length(ind_loop)==1){X_stack=matrix(X_stack,nrow=length(ind_loop),ncol=n_br); PR_p1up_stack=matrix(PR_p1up_stack,nrow=length(ind_loop),ncol=n_br)}
            ind_move=(ind_max-1)%%n_br+1
            mat_move=matrix(c(rep(1,n_br),rep(-1,n_br),0),nrow=length(ind_move),ncol=(n_br*2+1),byrow=TRUE) # move matrix
            X_stack[cbind(seq_along(ind_move), ind_move)]=X_stack[cbind(seq_along(ind_move), ind_move)]+mat_move[cbind(seq_along(ind_max), ind_max)]
            
            n_repeat=n_repeat+1
        } # end of repeat
        
        ind_match_stack=ifelse(rowSums(X_stack_all==X_p1up_stack)==n_br,1,0) # length=n_obs_p1up*n_edraw
        ind_match=matrix(ind_match_stack,nrow=n_obs_p1up,ncol=n_edraw,byrow=FALSE)
        ll_vec[ind_p1up]=ll_exit[ind_p1up]+log(rowMeans(ind_match));  # ll_i; ll_vec[ind_p1[i]]
    } # end of if(n_obs_p1up>0 & n_edraw>0)
    ll_vec
} 

    

### Main Algorithm -------------------------------------------------------------

## DATA
hierdata=Data$hierdata; 
Z=Data$Z
true_mu_theta=Data$true$mu_theta
true_sd_theta=Data$true$sd_theta
tparam_vec=true_mu_theta; n_param=length(tparam_vec)

fixparam=Data$fixparam; 
psi0=fixparam$psi0; e_sd=fixparam$e_sd
offset=fixparam$offset; os_in=offset[1]; os_out=offset[2]

n_hh=length(hierdata); 
n_br=ncol(hierdata[[1]]$X); n_obs=nrow(hierdata[[1]]$X)
n_z=ncol(Z)

## PRIOR
Deltabar=Prior$Deltabar
A=Prior$A
nu=Prior$nu
V=Prior$V
n_param_est=ncol(V); n_var=n_param_est

## MCMC
# iteration setting
R=Mcmc$R; keep=Mcmc$keep; n_report=Mcmc$n_report

# draw setting
est_setting=Mcmc$est_setting
n_edraw=est_setting$n_edraw
ind_gr1=est_setting$ind_est_group$ind_gr1; ind_gr2=est_setting$ind_est_group$ind_gr2; ind_gr3=est_setting$ind_est_group$ind_gr3; ind_est=sort(unique(c(ind_gr1,ind_gr2,ind_gr3))); ind_est; 
ss_all=(est_setting$ss_all); ss_all

# initial parameter values
Thetas=Mcmc$Thetas0
Delta=Mcmc$Delta0
Vtheta=Mcmc$Vtheta0

#  allocate spaces for saving the draws
Vthetadraw=array(double((R/keep)*n_var*n_var),dim=c(R/keep,n_var*n_var))
Deltadraw=array(double((R/keep)*n_z*n_var),dim=c(R/keep,n_z*n_var))
thetasdraw=array(double((R/keep)*n_hh*n_param),dim=c(n_hh,n_param,R/keep))
Sll=array(0,dim=c(R/keep))
Reject = array(0,dim=c((R/keep),n_hh,length(ss_all)))  # rejection ratio of MH algorithm

# data stacking for faster computing of the model likelihood
X_stack=PR_stack=NULL
n_obs_hh_vec=rep(0,n_hh)
for(hh in 1:n_hh){
    X_hh=Data$hierdata[[hh]]$X; X_stack=rbind(X_stack,X_hh) # n_obs*n_hh x n_br
    PR_hh=Data$hierdata[[hh]]$PR; PR_stack=rbind(PR_stack,PR_hh)
    n_obs_hh_vec[hh]=nrow(X_hh)
}
ID_stack=rep(c(1:n_hh),n_obs_hh_vec) 
Thetas_stack=Thetas[rep(seq(nrow(Thetas)), n_obs_hh_vec),] 
n_obs_stack=sum(n_obs_hh_vec)

# computing the model likelihood at initial values
ll_old_stack=LL_disc_stack(X_stack,PR_stack,ID_stack,Thetas_stack,fixparam,n_edraw)
SLL_old=as.vector(tapply(ll_old_stack, ID_stack, sum)); SLL_old
if(all(SLL_old>-Inf)==FALSE) {
    hh_repeat=which(SLL_old==-Inf); ind_1=which(ID_stack%in%hh_repeat)
    X_stack_1=X_stack[ind_1,]; PR_stack_1=PR_stack[ind_1,]; ID_stack_1=ID_stack[ind_1]; Thetas_stack_1=Thetas_stack[ind_1,]
    ll_old_stack_1=LL_disc_stack(X_stack_1,PR_stack_1,ID_stack_1,Thetas_stack_1,fixparam,n_edraw*10)
    SLL_old_1=as.vector(tapply(ll_old_stack_1, ID_stack_1, sum)); SLL_old[hh_repeat]=SLL_old_1
    if(all(SLL_old>-Inf)==FALSE) {
        hh_repeat=which(SLL_old==-Inf); ind_1=which(ID_stack%in%hh_repeat)
        X_stack_1=X_stack[ind_1,]; PR_stack_1=PR_stack[ind_1,]; ID_stack_1=ID_stack[ind_1]; Thetas_stack_1=Thetas_stack[ind_1,]
        ll_old_stack_1=LL_disc_stack(X_stack_1,PR_stack_1,ID_stack_1,Thetas_stack_1,fixparam,n_edraw*100)
        SLL_old_1=as.vector(tapply(ll_old_stack_1, ID_stack_1, sum)); SLL_old[hh_repeat]=SLL_old_1
        if(all(SLL_old>-Inf)==FALSE) {
            hh_repeat=which(SLL_old==-Inf); ind_1=which(ID_stack%in%hh_repeat)
            X_stack_1=X_stack[ind_1,]; PR_stack_1=PR_stack[ind_1,]; ID_stack_1=ID_stack[ind_1]; Thetas_stack_1=Thetas_stack[ind_1,]
            ll_old_stack_1=LL_disc_stack(X_stack_1,PR_stack_1,ID_stack_1,Thetas_stack_1,fixparam,n_edraw*1000)
            SLL_old_1=as.vector(tapply(ll_old_stack_1, ID_stack_1, sum)); SLL_old[hh_repeat]=SLL_old_1
            if(all(SLL_old>-Inf)==FALSE) {pandterm("Requires Changes in Initial Parameter Values")}
        }
    }
}

# initialize rejection ratio
rej=sumrej=matrix(0,nrow=n_hh,ncol=length(ss_all)) # sum across iterations	

itime=proc.time()[3]
cat("[(n_hh,n_obs,n_br,n_edraw)=(",c(n_hh,n_obs,n_br,n_edraw),")]",fill=TRUE)
cat("[tpsi*=",round(tparam_vec[ind_psi],2),"][tgamma*=",round(tparam_vec[ind_gamma],2),"][tE*=",round(tparam_vec[ind_E],2),"]",fill=TRUE)
cat("MCMC Iteration | Time to End (min)",fill=TRUE)
flush.console()

for (rep in 1:R) {
    
    Thetabar=Z%*%matrix(Delta,ncol=n_var) # mean vector for individual level params  
    rootp=chol(Vtheta)
    rootpi=backsolve(rootp,diag(n_var))
    
    ### (i) draw gr1 params given all the other by MH
    kk=1
    if(ss_all[kk]>0){
        ind=ind_gr1

        Thetas_new=Thetas
        Thetas_new[,ind]=Thetas[,ind] + ss_all[kk]*matrix(rnorm(n_hh*length(ind),mean=0,sd=1),nrow=n_hh,ncol=length(ind))
        Thetas_new_stack=Thetas_new[rep(seq(nrow(Thetas_new)), n_obs_hh_vec),] 
        
        ll_new_stack=LL_disc_stack(X_stack,PR_stack,ID_stack,Thetas_new_stack,fixparam,n_edraw)
        SLL_new=as.vector(tapply(ll_new_stack, ID_stack, sum))
        
        Lpost_old=Lpost_new=rep(0,n_hh)
        Lpost_old=SLL_old+dmvnorm(Thetas[,ind_est]-Thetabar, mean = rep(0,n_param_est), sigma = Vtheta, log = TRUE)
        Lpost_new=SLL_new+dmvnorm(Thetas_new[,ind_est] - Thetabar, mean = rep(0,n_param_est) , sigma = Vtheta, log = TRUE)
        
        Ldiff_vec=Lpost_new-Lpost_old
        Alpha_vec=exp(Ldiff_vec)
        Alpha_vec[Alpha_vec>=1]=1
        
        Unif_vec=runif(n_hh)
        Unif_vec[Alpha_vec>=1]=0
        
        Ind_move=(Unif_vec <= Alpha_vec)
        stay_vec=ifelse(Ind_move,0,1)
        if(any(Ind_move)){
            Thetas[Ind_move,]=Thetas_new[Ind_move,]
            SLL_old[Ind_move]=SLL_new[Ind_move]
        }
        
        rej[,kk]=stay_vec
        sumrej[,kk]=sumrej[,kk]+stay_vec
    }
    
    ### (ii) draw gr2 params given all the other by MH
    kk=2
    if(ss_all[kk]>0){
        ind=ind_gr2

        Thetas_new=Thetas
        Thetas_new[,ind]=Thetas[,ind] + ss_all[kk]*matrix(rnorm(n_hh*length(ind),mean=0,sd=1),nrow=n_hh,ncol=length(ind))
        Thetas_new_stack=Thetas_new[rep(seq(nrow(Thetas_new)), n_obs_hh_vec),] 

        ll_new_stack=LL_disc_stack(X_stack,PR_stack,ID_stack,Thetas_new_stack,fixparam,n_edraw)
        SLL_new=as.vector(tapply(ll_new_stack, ID_stack, sum))
        
        Lpost_old=Lpost_new=rep(0,n_hh)
        Lpost_old=SLL_old+dmvnorm(Thetas[,ind_est]-Thetabar, mean = rep(0,n_param_est), sigma = Vtheta, log = TRUE)
        Lpost_new=SLL_new+dmvnorm(Thetas_new[,ind_est] - Thetabar, mean = rep(0,n_param_est) , sigma = Vtheta, log = TRUE)
        
        Ldiff_vec=Lpost_new-Lpost_old
        Alpha_vec=exp(Ldiff_vec)
        Alpha_vec[Alpha_vec>=1]=1
        
        Unif_vec=runif(n_hh)
        Unif_vec[Alpha_vec>=1]=0
        
        Ind_move=(Unif_vec <= Alpha_vec)
        stay_vec=ifelse(Ind_move,0,1)
        if(any(Ind_move)){
            Thetas[Ind_move,]=Thetas_new[Ind_move,]
            SLL_old[Ind_move]=SLL_new[Ind_move]
        }
        
        rej[,kk]=stay_vec
        sumrej[,kk]=sumrej[,kk]+stay_vec
    }
    
    ### (iii) draw gr3 params given all the other by MH
    kk=3
    if(ss_all[kk]>0){
        ind=ind_gr3

        Thetas_new=Thetas
        Thetas_new[,ind]=Thetas[,ind] + ss_all[kk]*matrix(rnorm(n_hh*length(ind),mean=0,sd=1),nrow=n_hh,ncol=length(ind))
        Thetas_new_stack=Thetas_new[rep(seq(nrow(Thetas_new)), n_obs_hh_vec),] 

        ll_new_stack=LL_disc_stack(X_stack,PR_stack,ID_stack,Thetas_new_stack,fixparam,n_edraw)
        SLL_new=as.vector(tapply(ll_new_stack, ID_stack, sum))
        
        Lpost_old=Lpost_new=rep(0,n_hh)
        Lpost_old=SLL_old+dmvnorm(Thetas[,ind_est]-Thetabar, mean = rep(0,n_param_est), sigma = Vtheta, log = TRUE)
        Lpost_new=SLL_new+dmvnorm(Thetas_new[,ind_est] - Thetabar, mean = rep(0,n_param_est) , sigma = Vtheta, log = TRUE)
        
        Ldiff_vec=Lpost_new-Lpost_old
        Alpha_vec=exp(Ldiff_vec)
        Alpha_vec[Alpha_vec>=1]=1
        
        Unif_vec=runif(n_hh)
        Unif_vec[Alpha_vec>=1]=0
        
        Ind_move=(Unif_vec <= Alpha_vec)
        stay_vec=ifelse(Ind_move,0,1)
        if(any(Ind_move)){
            Thetas[Ind_move,]=Thetas_new[Ind_move,]
            SLL_old[Ind_move]=SLL_new[Ind_move]
        }
        
        rej[,kk]=stay_vec
        sumrej[,kk]=sumrej[,kk]+stay_vec
    }

    ###  draw Vtheta, Delta | {theta_i}
    regout=rmultireg(as.matrix(Thetas[,ind_est]),Z,Deltabar,A,nu,V)
    Vtheta=regout$Sigma
    Delta=regout$B
    
    if(rep%%keep == 0){
        mkeep=rep/keep
        Vthetadraw[mkeep,]=Vtheta
        Deltadraw[mkeep,]=Delta
        thetasdraw[,,mkeep]=Thetas
        Sll[mkeep]=sum(SLL_old)
        Reject[mkeep,,]=rej # n_hh by length(ss_all)
    }
    if(rep%%n_report ==0){
        avgrej=colMeans(sumrej/n_report)
        sumrej=matrix(0,nrow=n_hh,ncol=length(ss_all))
        
        ctime=proc.time()[3]
        timetoend=((ctime-itime)/rep)*(R-rep)
        
        theta_report=colMeans(Thetas)
        cat(" ",rep," (", round(timetoend/60,1),")","SLL=",round(sum(SLL_old),1),"avgrej=",round(c(avgrej[which(ss_all>0)]),2),"ss=",round(c(ss_all[which(ss_all>0)]),3)
            ,"psi*=",round(theta_report[ind_psi],2),"g*=",round(theta_report[ind_gamma],2),"E*=",round(theta_report[ind_E],2),fill=TRUE)
    }
} # end of for (rep in 1:R)
ctime = proc.time()[3]
cat(' Total Time Elapsed: ', round((ctime-itime)/60,2),'\n')  

total_minutes=(ctime-itime)/60
cat('  Simulation Study','\n')
cat("[(n_hh,n_obs,n_br,n_edraw)=(",c(n_hh,n_obs,n_br,n_edraw),")]",fill=TRUE)
cat("[tpsi*=",round(tparam_vec[ind_psi],2),"][tgamma*=",round(tparam_vec[ind_gamma],2),"][tE*=",round(tparam_vec[ind_E],2),"]",fill=TRUE)

OUT=list(Vthetadraw=Vthetadraw,Deltadraw=Deltadraw,thetasdraw=thetasdraw,Sll=Sll,Reject=Reject,total_minutes=total_minutes)
OUT
}


