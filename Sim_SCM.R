##### A Sequential Choice Model for Multiple Discrete Demand (Lee et al. (2022))
##### Data Simulation Code
##### Data Generation by the Sequential Choice Model
    # 1. starting from all zeros
    # 2. check the bang for the buck for each item and "buy or remove" one unit of the alternative whose bang for the buck is the best
    # 3. repeat 2. until (i)there is no more alternative with a positive bang for the buck or (ii)constrained by the budget
##### update: 2022.01

sim_scm=function(n_hh,n_obs,mu_theta,sd_theta,fixparam,pr_setting){

	psi0=fixparam$psi0; e_sd=fixparam$e_sd; offset=fixparam$offset; os_in=offset[1]; os_out=offset[2]
	n_param=length(mu_theta); n_param
	n_br=length(pr_setting$pr_min); n_br
	
	# regressor for heterogeneity distribution
	X_hh=cbind(rep(1,n_hh))
	
	# individual-level parameters
	IndParam=matrix(rnorm(n_hh*n_param,mu_theta,sd_theta),nrow=n_hh,ncol=n_param,byrow=TRUE)
	
	# price
	pr_min=pr_setting$pr_min ; pr_max=pr_setting$pr_max
	PR_array=array(NA,dim=c(n_hh,n_obs,n_br))
	for(k in 1:n_br){PR_array[,,k]=matrix(runif(n_hh*n_obs,min=pr_min[k],max=pr_max[k]),nrow=n_obs,ncol=n_hh)}
	
	hierdata=NULL		
	
	itime=proc.time()[3]
	cat("# of hh | Time to End (min)",fill=TRUE)
	for(hh in 1:n_hh){
	    theta_hh=IndParam[hh,]
	    tpsi_1=theta_hh[ind_psi]; tpsi=exp(tpsi_1)
	    tgamma_1=theta_hh[ind_gamma]; tgamma=exp(tgamma_1); if(length(tgamma)==1) {tgamma=rep(tgamma,n_br)}
	    tE_1=theta_hh[ind_E]; tE=exp(tE_1)

	    PR=PR_array[hh,,]
	    Error=matrix(rnorm(n_obs*n_br,mean=0,sd=e_sd),nrow=n_obs,ncol=n_br)
	    X=matrix(0,nrow=n_obs,ncol=n_br)
	    Z=matrix(0,nrow=n_obs,ncol=1)
	    Umax=matrix(0,nrow=n_obs,ncol=1)
	    
        for(i in 1:n_obs){
			pr_i=PR[i,]; e_i=Error[i,]

            # 1. starting from all zeros
			x_i=rep(0,n_br)
			# 2. check the bang for the buck for each item and buy or remove one unit of the alternative whose bang for the buck is the best
			# 3. repeat 2. until (i)there is no more alternative with a positive bang for the buck or (ii)constrained by the budget
			n_repeat=1
			repeat {
			    k=1
			    bfb_buy=rep(-Inf,n_br)
			    ind=which((tE-sum(pr_i*x_i)-k*pr_i+os_out)>0 & tE-sum(pr_i*x_i)-k*pr_i>=0)
			    if(length(ind)>0){
			        bfb_buy[ind]= tpsi[ind]*exp(e_i[ind])/tgamma[ind]*log((tgamma[ind]*(x_i[ind]+k)+os_in)/(tgamma[ind]*x_i[ind]+os_in)) + psi0*log((tE-sum(pr_i*x_i)-pr_i[ind]*k+os_out)/(tE-sum(pr_i*x_i)+os_out))
			    }
			    
			    k=-1
			    bfb_remove=rep(-Inf,n_br)
			    ind=which(x_i>0)
			    if(length(ind)>0){
			        bfb_remove[ind]= tpsi[ind]*exp(e_i[ind])/tgamma[ind]*log((tgamma[ind]*(x_i[ind]+k)+os_in)/(tgamma[ind]*x_i[ind]+os_in)) + psi0*log((tE-sum(pr_i*x_i)-pr_i[ind]*k+os_out)/(tE-sum(pr_i*x_i)+os_out))
			    }
			    
			    bfb=c(bfb_buy,bfb_remove)
			    if(max(bfb)<=0) break;
			    
			    ind_max=which.max(bfb)
			    if(ind_max<=n_br){x_i[ind_max]=x_i[ind_max]+1}    
			    if(ind_max>n_br){x_i[ind_max-n_br]=x_i[ind_max-n_br]-1}	
			    
			    n_repeat=n_repeat+1
			}
			z_i=tE-sum(pr_i*x_i)
			X[i,]=x_i
			Z[i]=z_i
			Umax[i]=sum(tpsi*exp(e_i)/tgamma*log((tgamma*x_i+os_in))) + psi0*log(z_i+os_out)
        } # end of i loop (1:n_obs)

		hierdata[[hh]]=list(X=X,Z=Z,Umax=Umax,PR=PR)

		if(hh%%100 ==0) {
			ctime=proc.time()[3]
			timetoend=((ctime-itime)/hh)*(n_hh-hh)
			cat(" ",hh," (", round(timetoend/60,1),")",fill=TRUE)
		}
	} # end of hh loop (1:n_hh)

	ctime = proc.time()[3]
	cat(' Total Time Elapsed: ', round((ctime-itime)/60,2),'\n')

	list(hierdata=hierdata,IndParam=IndParam,X_hh=X_hh) 
}
