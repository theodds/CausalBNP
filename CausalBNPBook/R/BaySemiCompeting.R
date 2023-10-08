BaySemiComp <- function(progression_time,
                        survival_time,
                        status_progress,
                        status_survival,
                        trt,
                        X,
                        num_burn,
                        num_save,
                        num_thin,
                        rho_grid = c(0.2, 0.5, 0.8),
                        time_grid = NULL
                        ) {



  Niter <- num_burn + num_save * num_thin
  burn.in <- num_burn
  lag <- num_thin
  full_data <- cbind(progression_time, survival_time, status_progress, status_survival)
  Z <- cbind(trt, X)

  stopifnot(all(trt %in% c(0,1)))
  stopifnot(all(status_progress %in% c(0,1)))
  stopifnot(all(status_survival %in% c(0,1)))

  out <- main_mcmc(Niter, burn.in, lag, full_data, Z)

  out$time_grid <- time_grid
  if(is.null(time_grid)) out$time_grid <- seq(from = 0, to = max(survival_time) * 1.1, length = 50)
  out$marginal_survival <- Marginal_survival(full_data, Z, out$mcmc0, out$mcmc1, num_save, out$time_grid)
  out$causal_estimand_grid <- vector(mode = 'list', length = length(rho_grid))
  for(i in 1:length(rho_grid)) {
    out$causal_estimand_grid[[i]] <- Estimate_hu(rho_grid[i], full_data, Z, out$mcmc0, out$mcmc1, num_save)
  }

  return(out)
}

## full_data: The first column is observed progression time, the second column is the observed
## survival time, the third column is delta (the censoring indicator for progress),
##  and the last column is the xi (the censoring indicator for survival)
## delta=xi=1, we observe both progression time and death time
## delta=xi=0, we only have the censoring time C
## delta=1, xi=0, we observe progression time, but censored before observing death
## delta=0, xi=1, we observe death only, but neither progression nor censoring

main_mcmc <- function(Niter, burn.in, lag, full_data, Z,
                      seed = sample.int(.Machine$integer.max, 1))
{

  Npat <- nrow(full_data) #the number of patients
  Ncov <- ncol(Z)
  #global variable
  H = 10  #the upper limit for stick breaking of DP.
  lambda0 = 4
  lambda1 = 1 #prior for M
  lambda2 = 1

  id1 = which(Z[,1]==1)
  id0 = which(Z[,1]==0)
  N1 = length(id1)
  N0 = length(id0)
  cov = cbind(rep(1, Npat), Z[,-1])
  cov1 = Z[id1,]
  cov0 = cbind(rep(1, N0), Z[id0,-1])


  temp0 = matrix(0, N0, N0)
  for (i in 1:Ncov)
  {
    temp0 = temp0 + outer(cov0[,i],cov0[,i],'-')^2
  }
  covariance0 = exp(-temp0) + 0.01*diag(N0)
  inv_covariance0 = solve(covariance0)
  temp1 = matrix(0, N1, N1)
  for (i in 1:Ncov)
  {
    temp1 = temp1 + outer(cov1[,i],cov1[,i],'-')^2
  }
  covariance1 = exp(-temp1) + 0.01*diag(N1)
  inv_covariance1 = solve(covariance1)

  var1_betah = solve(t(cov1)%*%inv_covariance1%*%cov1+diag(Ncov))
  var0_betah = solve(t(cov0)%*%inv_covariance0%*%cov0+diag(Ncov))

  partial_cov0 = t(cov0)%*%inv_covariance0
  partial_cov1 = t(cov1)%*%inv_covariance1

  #prior for \Sigma
  Phi1 = cov(full_data[id1,1:2])
  Phi0 = cov(full_data[id0,1:2])


  mcmc1 <- NULL
  mcmc1$M <- rep(NA, Niter)
  mcmc1$Sigma <- array(NA, c(2,2,Niter))
  mcmc1$muh <- array(NA, c(H, 2, N1, Niter))
  mcmc1$betah1 <- array(NA, c(H, Ncov, Niter))
  mcmc1$betah2 <- array(NA, c(H, Ncov, Niter))
  mcmc1$r <- array(NA, c(N1, Niter))
  mcmc1$wh <- array(NA, c(H, Niter))
  mcmc1$imputey <- array(NA, c(N1,2, Niter))

  set.seed(seed)
  ##Initialize
  initial <- init(full_data[id1,1:2], N1, H)
  mcmc1$M[1] = initial$M
  mcmc1$Sigma[,,1] = initial$Sigma
  for (i in 1:N1) mcmc1$muh[,,i,1] = initial$mh
  mcmc1$wh[,1] = initial$wh
  mcmc1$r[,1] = initial$r
  lfit1 = survreg(Surv(full_data[id1,1],full_data[id1,3])~cov1,dist="gaussian")
  lfit2 = survreg(Surv(full_data[id1,2],full_data[id1,4])~cov1,dist="gaussian")
  mcmc1$betah1[,,1] = matrix(rep(lfit1$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc1$betah2[,,1] = matrix(rep(lfit2$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc1$imputey[,,1] = full_data[id1,1:2]
  prior_betah1 = mcmc1$betah1[1,,1]
  prior_betah2 = mcmc1$betah2[1,,1]

  print("MCMC Iteration under treatment:")
  for (iter in 2:Niter)
  {
    print(iter)
    tmp = update_r(mcmc1$Sigma[,,iter-1], mcmc1$muh[,,,iter-1],mcmc1$wh[,iter-1], full_data[id1,], mcmc1$r[,iter-1], N1, H)
    mcmc1$r[,iter] = tmp$r
    mcmc1$imputey[,,iter] = tmp$impute_y
    tmp2 = update_wh_and_M(mcmc1$r[,iter], mcmc1$M[iter-1],lambda1, lambda2, H)
    mcmc1$M[iter] = tmp2$M
    mcmc1$wh[,iter] = tmp2$wh
    mcmc1$muh[,,,iter] = update_muh(mcmc1$wh[,iter-1],mcmc1$r[,iter], mcmc1$Sigma[,,iter-1], mcmc1$muh[,,,iter-1], mcmc1$imputey[,,iter], mcmc1$betah1[,,iter-1],mcmc1$betah2[,,iter-1], N1, cov1,covariance1, inv_covariance1,partial_cov1, H)
    mcmc1$betah1[,,iter] = update_betah1(mcmc1$muh[,,,iter], mcmc1$r[,iter],prior_betah1,var1_betah,partial_cov1, Ncov, H)
    mcmc1$betah2[,,iter] = update_betah2(mcmc1$muh[,,,iter], mcmc1$r[,iter],prior_betah2,var1_betah,partial_cov1,Ncov, H)
    mcmc1$Sigma[,,iter] = update_Sigma(mcmc1$r[,iter], mcmc1$imputey[,,iter], mcmc1$muh[,,,iter], N1, Phi1,lambda0, H)
  }

  mcmc0 <- NULL
  mcmc0$M <- rep(NA, Niter)
  mcmc0$Sigma <- array(NA, c(2,2,Niter))
  mcmc0$muh <- array(NA, c(H, 2, N0, Niter))
  mcmc0$betah1 <- array(NA, c(H, Ncov, Niter))
  mcmc0$betah2 <- array(NA, c(H, Ncov, Niter))
  mcmc0$r <- array(NA, c(N0, Niter))
  mcmc0$wh <- array(NA, c(H, Niter))
  mcmc0$imputey <- array(NA, c(N0,2, Niter))

  set.seed(seed)
  ##Initialize
  initial <- init(full_data[id0,1:2], N0, H)
  mcmc0$M[1] = initial$M
  mcmc0$Sigma[,,1] = initial$Sigma
  for (i in 1:N0) mcmc0$muh[,,i,1] = initial$mh
  mcmc0$wh[,1] = initial$wh
  mcmc0$r[,1] = initial$r
  lfit1 = survreg(Surv(full_data[id0,1],full_data[id0,3])~cov0,dist="gaussian")
  lfit2 = survreg(Surv(full_data[id0,2],full_data[id0,4])~cov0,dist="gaussian")
  mcmc0$betah1[,,1] = matrix(rep(lfit1$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc0$betah2[,,1] = matrix(rep(lfit2$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc0$imputey[,,1] = full_data[id0,1:2]
  prior_betah1 = mcmc0$betah1[1,,1]
  prior_betah2 = mcmc0$betah2[1,,1]

  print("MCMC Iteration under control:")
  for (iter in 2:Niter)
  {
    print(iter)
    tmp = update_r(mcmc0$Sigma[,,iter-1], mcmc0$muh[,,,iter-1],mcmc0$wh[,iter-1], full_data[id0,], mcmc0$r[,iter-1], N0, H)
    mcmc0$r[,iter] = tmp$r
    mcmc0$imputey[,,iter] = tmp$impute_y
    tmp2 = update_wh_and_M(mcmc0$r[,iter], mcmc0$M[iter-1],lambda1, lambda2, H)
    mcmc0$M[iter] = tmp2$M
    mcmc0$wh[,iter] = tmp2$wh
    mcmc0$muh[,,,iter] = update_muh(mcmc0$wh[,iter-1],mcmc0$r[,iter], mcmc0$Sigma[,,iter-1], mcmc0$muh[,,,iter-1], mcmc0$imputey[,,iter], mcmc0$betah1[,,iter-1],mcmc0$betah2[,,iter-1], N0, cov0,covariance0, inv_covariance0,partial_cov0, H)
    mcmc0$betah1[,,iter] = update_betah1(mcmc0$muh[,,,iter], mcmc0$r[,iter],prior_betah1,var0_betah,partial_cov0,Ncov,H)
    mcmc0$betah2[,,iter] = update_betah2(mcmc0$muh[,,,iter], mcmc0$r[,iter],prior_betah2,var0_betah,partial_cov0,Ncov,H)
    mcmc0$Sigma[,,iter] = update_Sigma(mcmc0$r[,iter], mcmc0$imputey[,,iter], mcmc0$muh[,,,iter], N0, Phi0,lambda0, H)
  }

  id <- seq(burn.in+1, Niter, lag)
  mcmc1_save <- NULL
  mcmc1_save$M = mcmc1$M[id]
  mcmc1_save$Sigma = mcmc1$Sigma[,,id]
  mcmc1_save$muh = mcmc1$muh[,,,id]
  mcmc1_save$betah1 = mcmc1$betah1[,,id]
  mcmc1_save$betah2 = mcmc1$betah2[,,id]
  mcmc1_save$r = mcmc1$r[,id]
  mcmc1_save$wh = mcmc1$wh[,id]
  mcmc1_save$imputey = mcmc1$imputey[,,id]

  mcmc0_save <- NULL
  mcmc0_save$M = mcmc0$M[id]
  mcmc0_save$Sigma = mcmc0$Sigma[,,id]
  mcmc0_save$muh = mcmc0$muh[,,,id]
  mcmc0_save$betah1 = mcmc0$betah1[,,id]
  mcmc0_save$betah2 = mcmc0$betah2[,,id]
  mcmc0_save$r = mcmc0$r[,id]
  mcmc0_save$wh = mcmc0$wh[,id]
  mcmc0_save$imputey = mcmc0$imputey[,,id]

  mcmc0 = mcmc0_save
  mcmc1 = mcmc1_save
  #save(mcmc0, mcmc1, file="saved_mcmc.RData")
  return(list(mcmc0=mcmc0, mcmc1=mcmc1))
}

init <- function(response, Num_Patient, H)
{
  hc = hclust(dist(response)^2, "cen")
	r = cutree(hc, k=H)
  wh1 <- table(r)/Num_Patient
  idx <- order(wh1,decreasing=T)
  wh <- wh1[idx]
  tmp0 = split(response,r)
  mh = matrix(NA, H, 2)
  for (h in 1:H)
  {
    mh[h,] = colSums(matrix(tmp0[[h]],ncol=2))/sum(r==h)
  }
  M = 1
  Sigma = matrix(c(0.5, 0.1, 0.1, 0.5),2,2)
  return(list(mh=mh,wh=wh, M=M,Sigma=Sigma, r=r))
}

#wh=mcmc0$wh[,iter-1]; r=mcmc0$r[,iter]; muh=mcmc0$muh[,,,iter-1];data=full_data[id0,];Sigma=mcmc0$Sigma[,,,iter-1];
#cov=cov0; covariance=covariance0; inv_covariance=inv_covariance0; r=mcmc0$r[,iter-1];Num_Patient=N0;
update_r <- function(Sigma, muh, wh, data, r, Num_Patient, H)
{
  impute_y = data[,1:2]
	for (i in 1:Num_Patient)
	{
    if (data[i,3]==1 & data[i,4]==1)
		{
      tmp0 = rep(0, H)
			ph = diag(exp(-0.5*t(data[i,1:2]-t(muh[,,i]))%*%solve(Sigma)%*%(data[i,1:2]-t(muh[,,i]))))*wh
        #dmvnorm(data[i,1:2], muh[2,], Sigma)*wh
		  r[i] = sample(1:H, 1, prob=ph)
		}
		if (data[i,3]==0 & data[i,4]==0)
		{
      tmp0 = rep(0, H)
      for (h in 1:H) tmp0[h] = pmvnorm(lower=data[i,1:2], upper=c(Inf,Inf),mean=muh[h,,i], sigma=Sigma)
			ph = tmp0*wh
      ph[ph<=0] = 0
		  r[i] = sample(1:H, 1, prob=ph)
			impute_y[i,] = rtmvnorm(1, mean=muh[r[i],,i], sigma=Sigma, lower=data[i,1:2],algorithm="gibbs")
			if (impute_y[i,1]==Inf | impute_y[i,2]==Inf) impute_y[i,]=c(data[i,1],data[i,2])
		}
		if (data[i,3]==1 & data[i,4]==0)
		{
		  tmp0 = rep(0, H)
		  for (h in 1:H) tmp0[h] = dtmvnorm.marginal(data[i,1],n=1, mean=muh[h,,i],sigma=Sigma, lower=c(-Inf,data[i,2]),upper=c(Inf,Inf))
		  tmp0[tmp0=="NaN"] = 0
		  tmp0[tmp0=="Inf"] = 0
      ph = tmp0*wh
		  r[i] = sample(1:H, 1, prob=ph)
		  impute_y[i,] = rtmvnorm(1, mean=muh[r[i],,i], sigma=Sigma, lower=c(data[i,1],data[i,2]),upper=c(data[i,1]+1e-6,Inf),algorithm="gibbs")
		  if (impute_y[i,1]==Inf | impute_y[i,2]==Inf) impute_y[i,]=c(data[i,1],data[i,2])
		}
		if (data[i,3]==0 & data[i,4]==1)
		{
		  tmp0 = rep(0, H)
		  for (h in 1:H) tmp0[h] = dtmvnorm.marginal(data[i,2],n=2, mean=muh[h,,i],sigma=Sigma, lower=c(data[i,1],-Inf),upper=c(Inf,Inf))
		  tmp0[tmp0=="NaN"] = 0
		  tmp0[tmp0=="Inf"] = 0
      ph = tmp0*wh
		  r[i] = sample(1:H, 1, prob=ph)
		  impute_y[i,] = rtmvnorm(1, mean=muh[r[i],,i], sigma=Sigma, lower=c(data[i,1],data[i,2]),upper=c(Inf,data[i,2]+1e-6),algorithm="gibbs")
		  if (impute_y[i,1]==Inf | impute_y[i,2]==Inf) impute_y[i,]=c(data[i,1],data[i,2])
		}
	}
	return(list(r=r,impute_y=impute_y))
}

#r<-mcmc$r[,iter]; M=mcmc$M[iter-1]
update_wh_and_M <- function(r, M, lambda1, lambda2, H)
{
   ## returns: wh
    vh <- rep(0,H)  # initialize
    wh <- rep(0,H)
    V <-  1         # record prod_{g<h} (1-vh_h)
    for(h in 1:(H-1)){
      Ah <- which(r==h)
      Bh <- which(r>h)
      vh[h] <-  rbeta(1, 1+length(Ah), M+length(Bh))
      wh[h] <- vh[h]*V
      V <- V*(1-vh[h])
    }
    vh[H] <- 1.0
    wh[H] <- V
    #M <- rgamma(1, lambda1+H-1, lambda2-sum(log(1-vh[1:(H-1)])) )
    M <- rgamma(1, lambda1+H-1, lambda2-max(-100, sum(log(1-vh[1:(H-1)]))) )
    #log(.Machine$double.xmin) = -708.3964
    return(list(wh=wh, M=M))
}

#wh=mcmc0$wh[,iter-1]; r=mcmc0$r[,iter]; muh=mcmc0$muh[,,,iter-1];response=mcmc0$imputey[,,iter];Sigma=mcmc0$Sigma[,,iter-1]; betah1=mcmc0$betah1[,,iter-1];betah2=mcmc0$betah2[,,iter-1]
#cov=cov0; covariance=covariance0; inv_covariance=inv_covariance0; r=mcmc0$r[,iter-1];Num_Patient=N0;partial_cov=partial_cov0
update_muh <- function(wh,r,Sigma, muh, response, betah1,betah2, Num_Patient, cov, covariance, inv_covariance,partial_cov, H)
{
	   for(h in 1:H){
	     if(any(r==h)){      # some data assigned to h-th pointmass
	       Sh <- which(r==h)
	       nh <- length(Sh)
         W = matrix(0, nh, Num_Patient)
         for (i in 1:nh)
         {
           W[i, Sh[i]] = 1
         }
         response1 = response[Sh,1] - Sigma[1,2]/Sigma[2,2]*(response[Sh,2]-muh[h, 2, Sh])
         sigma21 = Sigma[1,1] - Sigma[1,2]*Sigma[2,1]/Sigma[2,2]
	       var1 = chol2inv(chol(inv_covariance+1/sigma21*t(W)%*%W))
	       mu1 = var1%*%(t(W)%*%response1/sigma21 + t(partial_cov)%*%betah1[h,])
	       muh[h,1,] <- rmvn(1, t(mu1), var1)
	       response2 = response[Sh,2] - Sigma[2,1]/Sigma[1,1]*(response[Sh,1]-muh[h, 1, Sh])
	       sigma22 = Sigma[2,2] - Sigma[1,2]*Sigma[2,1]/Sigma[1,1]
	       var2 = chol2inv(chol(inv_covariance+1/sigma22*t(W)%*%W))
	       mu2 = var2%*%(t(W)%*%response2/sigma22 + t(partial_cov)%*%betah2[h,])
	       muh[h,2,] <- rmvn(1, t(mu2), var2)
	     } else {            # no data assinged to h-th pointmass# sample from base measure
	       mu1 <- cov%*%betah1[h,]
	       mu2 <- cov%*%betah2[h,]
	       muh[h,1,] <- rmvn(1, t(mu1), covariance)
	       muh[h,2,] <- rmvn(1, t(mu2), covariance)
	     }
	   }
    return(muh)
}


#muh <- mcmc1$muh[,,,iter];r=mcmc1$r[,iter]
update_betah1 <- function(muh, r,prior_betah1, var,partial_cov, Ncov, H)
{
	betah1 <- matrix(0,H, Ncov)     # initialize
     for(h in 1:H){
      if(any(r==h)){      # some data assigned to h-th pointmass
        mu = var%*%(partial_cov%*%muh[h,1,] + diag(Ncov)%*%prior_betah1 )
        betah1[h,] <- rmvn(1, t(mu), var)
      } else {            # no data assinged to h-th pointmass     # sample from base measure
        betah1[h,] <- rmvn(1, prior_betah1, diag(Ncov))}
    }
    return(betah1)
}





#muh <- mcmc$muh[,,,iter];r=mcmc$r[,iter]
update_betah2 <- function(muh, r,prior_betah2, var, partial_cov, Ncov, H)
{
  betah2 <- matrix(0,H, Ncov)     # initialize
  for(h in 1:H){
    if(any(r==h)){      # some data assigned to h-th pointmass
      mu = var%*%(partial_cov%*%muh[h,2,] + diag(Ncov)%*%prior_betah2 )
      betah2[h,] <- rmvn(1, t(mu), var)
    } else {            # no data assinged to h-th pointmass     # sample from base measure
      betah2[h,] <- rmvn(1, prior_betah2, diag(Ncov))}
  }
  return(betah2)
}


#mcmc=mcmc0; muh=mcmc$muh[,,,iter];r=mcmc$r[,iter]; response=mcmc$imputey[,,iter];Num_Patient=N0; Phi=Phi0
update_Sigma <- function(r, response, muh, Num_Patient, Phi,lambda0, H)
{
	tmp = matrix(0, Num_Patient, 2)
	for (h in 1:H)
	{
		if(any(r==h))
		{
			Sh = which(r==h)
			tmp[Sh,] = t(muh[h,,Sh])
		}
	}
	Sigma <- riwish(lambda0+Num_Patient, Phi + t(response-tmp)%*%(response-tmp))
	return(Sigma)
}

#wh<-mcmc$wh[,iter]; muh=mcmc$muh[,,iter]; sigma2=mcmc$sigma2[iter];r=mcmc$r[,iter];covD=cov;
fmean <- function(wh, muh, theta, cov, covD, inv_covariance, Npat, sigma2)
{
  fx = rep(0, Npat)
  for (i in 1:Npat)
  {
    coeff = t(covD)-cov[i,]
    c1 = exp(-colSums(coeff^2*theta^2))
    fx[i] = exp(sum(wh*c1%*%inv_covariance%*%t(muh)+1/2*wh^2*sigma2))
  }
  return(fx)
}




#xgrid= tim;wh=mcmc0$wh[,iter]; muh=mcmc0$muh[,2,,iter]; betah=mcmc0$betah2[,,iter];sigma2=mcmc0$Sigma[2,2,iter];dpat=covD
#xgrid= tim;wh=mcmc1$wh[,iter]; muh=mcmc1$muh[,1,,iter]; betah=mcmc1$betah1[,,iter];sigma2=mcmc1$Sigma[1,,iter]
#wh=mcmc1$wh[,iter]; muh=mcmc1$muh[,2,,iter]; betah=mcmc1$betah2[,,iter];sigma2=mcmc1$Sigma[2,2,iter]
fmar_survival <- function(xgrid, wh, muh, betah, sigma2, dpat,cov,inv_covariance,H)
{
  fx <- rep(0, length(xgrid))
  for (h in 1:H)
  {
    coeff = t(cov)-dpat
    c1 = exp(-colSums(coeff^2))
    var1 = 1.01 - c1%*%inv_covariance%*%c1
    fx <- fx + wh[h]*pnorm(xgrid, m=betah[h,]%*%dpat+c1%*%inv_covariance%*%(muh[h,]-cov%*%betah[h,]), sd=sqrt(sigma2 + var1))
  }
  return(fx)
}







fbar_iter <- function(mcmc0, mcmc1, u, rho, Npat,partial_cov0,partial_cov1, cov0, cov1, cov,var1, var0, H, nsave)
{
  Mrep = 100
  hu = array(0, c(nsave,Npat,2))
  hu_ave = matrix(0, nsave, 2)
  for (iter in 1:nsave)
  {
    wh0=mcmc0$wh[,iter]; muh0=mcmc0$muh[,,,iter];
    betah10=mcmc0$betah1[,,iter];betah20=mcmc0$betah2[,,iter];
    Sigma0=mcmc0$Sigma[,,iter];
    wh1=mcmc1$wh[,iter]; muh1=mcmc1$muh[,,,iter];
    betah11=mcmc1$betah1[,,iter];betah21=mcmc1$betah2[,,iter];
    Sigma1=mcmc1$Sigma[,,iter];

    mu10_all = betah10%*%t(cov) + t(partial_cov0%*%(t(muh0[,1,])-cov0%*%t(betah10)) )
    mu20_all = betah20%*%t(cov) + t(partial_cov0%*%(t(muh0[,2,])-cov0%*%t(betah20)) )

    mu11_all = betah11%*%t(cov) + t(partial_cov1%*%(t(muh1[,1,])-cov1%*%t(betah11)) )
    mu21_all = betah21%*%t(cov) + t(partial_cov1%*%(t(muh1[,2,])-cov1%*%t(betah21)) )

    for (k in 1:Npat)
    {
      tmp1 = sample(1:H, Mrep, replace=T, prob=wh1)
      mean_tmp1 = t(rbind(mu11_all[,k], mu21_all[,k]))
      mean1 = mean_tmp1[tmp1,]
      tem_sample1 = matrix(0, Mrep, 2)
      tem_sample1 = mgcv::rmvn(Mrep, mean1, Sigma1 + var1[k]*diag(2))
      temp2 = 0
      for (h in 1:H)
      {
        temp2 = temp2 + wh1[h]*pnorm(tem_sample1[,2], mu21_all[h,k], sqrt(Sigma1[2,2] + var1[k]))
      }
      tmp2 = rho*qnorm(temp2) + rnorm(Mrep, 0, sqrt(1-rho^2))
      tmp_value = sum(wh0*pnorm(u, mu20_all[,k], sqrt(Sigma0[2,2] + var0[k])) )
      if (tmp_value>0.999) tmp_value=0.999
      tmp3 = qnorm( tmp_value )
      num1 = sum(tem_sample1[,1] < u & tem_sample1[,2] >=u & tmp2 >= tmp3)

      tmp0 = sample(1:H, Mrep, replace=T, prob=wh0)
      mean_tmp0 = t(rbind(mu10_all[,k], mu20_all[,k]))
      mean0 = mean_tmp0[tmp0,]
      tem_sample0 = matrix(0, Mrep, 2)
      tem_sample0 = mgcv::rmvn(Mrep, mean0, Sigma0 + var0[k]*diag(2))
      temp4 = 0
      for (h in 1:H)
      {
        temp4 = temp4 + wh0[h]*pnorm(tem_sample0[,2], mu20_all[h,k], sqrt(Sigma0[2,2] + var0[k]))
      }
      tmp4 = rho*qnorm(temp4)
      tmp4 = tmp4 + rnorm(Mrep, 0, sqrt(1-rho^2))
      tmp_value = sum(wh1*pnorm(u, mu21_all[,k], sqrt(Sigma1[2,2] + var1[k])) )
      if (tmp_value>0.999) tmp_value=0.999
      tmp5 = qnorm( tmp_value )
      num0 = sum(tem_sample0[,1] < u & tem_sample0[,2] >=u & tmp4 >= tmp5)
      hu[iter,k,1] = num0
      hu[iter,k,2] = num1
    }
  }
    hu_ave[,1] = apply(hu[,,1], 1, mean)
    hu_ave[,2] = apply(hu[,,2], 1, mean)
  return(hu_ave)
}

Marginal_survival <- function(full_data, Z, mcmc0, mcmc1, nsave, tim)
{
  Npat <- nrow(full_data) #the number of patients
  Ncov <- ncol(Z)
  H <- 10


  id1 = which(Z[,1]==1)
  id0 = which(Z[,1]==0)
  N1 = length(id1)
  N0 = length(id0)
  cov = cbind(rep(1, Npat), Z[,-1])
  cov1 = Z[id1,]
  cov0 = cbind(rep(1, N0), Z[id0,-1])

  temp0 = matrix(0, N0, N0)
  for (i in 1:Ncov)
  {
    temp0 = temp0 + outer(cov0[,i],cov0[,i],'-')^2
  }
  covariance0 = exp(-temp0) + 0.01*diag(N0)
  inv_covariance0 = solve(covariance0)
  temp1 = matrix(0, N1, N1)
  for (i in 1:Ncov)
  {
    temp1 = temp1 + outer(cov1[,i],cov1[,i],'-')^2
  }
  covariance1 = exp(-temp1) + 0.01*diag(N1)
  inv_covariance1 = solve(covariance1)

  fmean0 = array(NA, c(Npat,nsave, length(tim)))
  fmean1 = array(NA, c(Npat,nsave, length(tim)))
  for (k in 1:Npat)
  {
    #print(k)
    covD = cov[k,]
    # fgrid0 = NULL
    # fgrid1 = NULL
    for (i in 1:nsave)
    {
      fmean0[k,i,] <- 1-fmar_survival(tim, mcmc0$wh[,i],mcmc0$muh[,2,,i], mcmc0$betah2[,,i], mcmc0$Sigma[2,2,i], covD,cov0,inv_covariance0,H)
      fmean1[k,i,] <- 1-fmar_survival(tim, mcmc1$wh[,i],mcmc1$muh[,2,,i], mcmc1$betah2[,,i], mcmc1$Sigma[2,2,i], covD,cov1,inv_covariance1,H)
    }
  }
  tmp0 = matrix(0, nsave, length(tim))
  tmp1 = matrix(0, nsave, length(tim))
  for (i in 1:nsave)
  {
    tmp0[i,]=apply(fmean0[,i,], 2, mean)
    tmp1[i,]=apply(fmean1[,i,], 2, mean)
  }
  fmean0_ave = apply(tmp0, 2, mean)
  fmean1_ave = apply(tmp1, 2, mean)
  fquantile0 = apply(tmp0, 2, function(x) quantile(x, c(0.025, 0.975)))
  fquantile1 = apply(tmp1, 2, function(x) quantile(x, c(0.025, 0.975)))
  return(list(fmean0_ave=fmean0_ave, fmean1_ave=fmean1_ave, fquantile0=fquantile0, fquantile1=fquantile1))
}

Estimate_hu <- function(rho, full_data, Z, mcmc0, mcmc1, nsave)
{
  Npat <- nrow(full_data) #the number of patients
  Ncov <- ncol(Z)
  id1 = which(Z[,1]==1)
  id0 = which(Z[,1]==0)
  N1 = length(id1)
  N0 = length(id0)
  cov = cbind(rep(1, Npat), Z[,-1])
  cov1 = Z[id1,]
  cov0 = cbind(rep(1, N0), Z[id0,-1])
  c1 = matrix(0, length(id1), Npat)
  c0 = matrix(0, length(id0), Npat)
  H=10
  temp0 = matrix(0, N0, N0)
  for (i in 1:Ncov)
  {
    temp0 = temp0 + outer(cov0[,i],cov0[,i],'-')^2
  }
  covariance0 = exp(-temp0) + 0.01*diag(N0)
  inv_covariance0 = solve(covariance0)
  temp1 = matrix(0, N1, N1)
  for (i in 1:Ncov)
  {
    temp1 = temp1 + outer(cov1[,i],cov1[,i],'-')^2
  }
  covariance1 = exp(-temp1) + 0.01*diag(N1)
  inv_covariance1 = solve(covariance1)
  for (k in 1:Npat)
  {
    covD = cov[k,]
    coeff1 = t(cov1)-covD
    c1[,k] = exp(-colSums(coeff1^2))
    coeff0 = t(cov0)-covD
    c0[,k] = exp(-colSums(coeff0^2))
  }
  var1 = 1.01-diag(t(c1)%*%inv_covariance1%*%c1)
  var0 = 1.01-diag(t(c0)%*%inv_covariance0%*%c0)
  partial_cov1 = t(c1)%*%inv_covariance1
  partial_cov0 = t(c0)%*%inv_covariance0

  u_range = seq(0,6,0.5)
  hu_result <- array(0, c(length(u_range), nsave, 2))
  for (u_index in 1:length(u_range))
    hu_result[u_index,,] = fbar_iter(mcmc0, mcmc1, u_range[u_index], rho, Npat,partial_cov0,partial_cov1, cov0, cov1, cov, var1, var0, H, nsave)
  return(hu_result)
}
