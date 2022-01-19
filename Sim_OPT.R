##### A Sequential Choice Model for Multiple Discrete Demand (Lee et al. (2022))
##### Data Generation by an Exhaustive Search of the Entire Feasible Set
##### update: 2022.01

sim_opt=function(n_hh,n_obs,mu_theta,sd_theta,fixparam,pr_setting){

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
            max_int_x=tE%/%pr_i
            
            x1_vec=c(0:max_int_x[1]); x2_vec=c(0:max_int_x[2])
            x1_mat=matrix(x1_vec,nrow=length(x2_vec),ncol=length(x1_vec),byrow=TRUE)
            x2_mat=matrix(x2_vec,nrow=length(x2_vec),ncol=length(x1_vec),byrow=FALSE)
            z_mat=tE-pr_i[1]*x1_mat-pr_i[2]*x2_mat; 
            
            u_in_mat=tpsi[1]*exp(e_i[1])/tgamma[1]*log(tgamma[1]*x1_mat+os_in)+tpsi[2]*exp(e_i[2])/tgamma[2]*log(tgamma[2]*x2_mat+os_in)
            u_out_mat=matrix(-Inf,nrow=nrow(z_mat),ncol=ncol(z_mat)); u_out_mat[z_mat>=0]=psi0*log(z_mat[z_mat>=0]+os_out)
            u_mat=u_in_mat+u_out_mat

            ind=which.max(u_mat)
            x_i=c(x1_mat[ind],x2_mat[ind])
            z_i=z_mat[ind]
            umax_i=u_mat[ind]
            
            X[i,]=x_i
            Z[i]=z_i
            Umax[i]=umax_i
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
