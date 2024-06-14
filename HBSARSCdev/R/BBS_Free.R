library(splines)
# ******************************************************************************
# BBS_Free_v1.R
#   BBS without shape constraints   
# ******************************************************************************
"BBS_Free" <- function(ydata,xdata,vdata=NULL,mcmc=list(),nknots=20,
                       nint=100,xmin=NULL,xmax=NULL) {
  
  # MCMC Iterations
  mcvals=list(nblow=20000,smcmc=10000)
  mcvals[names(mcmc)]=mcmc
  nblow  = mcvals$nblow # number of initial MCMC parameters
  smcmc  = mcvals$smcmc # number of MCMC parameters to save
  nint   = nint         # Number of intervals in grid
  nbasis = nbasis       # Number of cosine functions in each direction
  # ****************************************************************************
  # Parameters for grids of x1 and x2
  if (is.null(xmin)) xmin = min(xdata)  # Minimum of xgrid
  if (is.null(xmax)) xmax = max(xdata)  # Maximum of xgrid
  
  nmcmc = nblow+smcmc
  
  # ****************************************************************************
  # Get Data: vdata (with constant), xdata (x1,x2),ydata
  # ****************************************************************************
  nobs     = length(ydata)
  ynames   = colnames(ydata)
  if (is.null(ynames)) {
    ynames  = "Y"
    ydata   = matrix(ydata)
    colnames(ydata) = ynames
  }
  
  if (is.null(vdata)) {
    vdata           = matrix(1,nobs,1)
    vnames          = "CNST"
    colnames(vdata) = "CNST"
    nparv           = 1
  } else {
    nparv           = ncol(vdata)
    vnames          = colnames(vdata)
    if (is.null(vnames)) {
      vnames = "CNST"
      if (nparv>1) {
        vnames = c(vnames,paste("V",1:(nparv-1),sep=""))
      }
      colnames(vdata) = vnames
    }
  }
  vtv = crossprod(vdata)
  
  # ****************************************************************************
  xrange = xmax-xmin
  # ****************************************************************************
  # Compute xgrid for plotting
  xdelta = xrange/nint
  xgrid  = seq(xmin,xmax,xdelta)
  # ****************************************************************************
  # Get bsplines
  # ****************************************************************************
  # BSpline basis functions
  # Get knots & B-splines based  on quantiles of all data
  phixdata  = bs(xdata,df=nknots,Boundary.knots=c(xmin,xmax))
  phixgrid  = predict(phixdata,xgrid)  # Compute basis functions at xgrid 
  # Find knot location
  nknots    = ncol(phixgrid)
  ntheta    = nknots
  nbasis    = nknots
  xknots    = matrix(0,nknots,1)  # Find max of basis functions: knots?
  for (k in 1:nknots) {
    idk       = which.max(phixgrid[,k])
    xknots[k] = xgrid[idk]
  }
  kall = 1:ntheta
  
  # ****************************************************************************
  # Prior Distributions
  # ****************************************************************************
  # Beta ~ N(b0,V0)
  beta_m0    = matrix(0,nparv,1)
  v0i        = 1/1000
  beta_v0i   = diag(v0i,nparv,nparv)
  beta_v0im0 = beta_v0i %*% beta_m0
  # sigma^2 ~ IG(r0/2,s0/2)
  sigma2_mu = 1   # Prior Mean
  sigma2_sd = 100  # Prior SD
  sigma2_r0 = 2*((sigma2_mu/sigma2_sd)^2 + 2)
  sigma2_s0 = 2*sigma2_mu*(sigma2_r0/2 - 1)
  sigma2_rn = sigma2_r0 + nobs
  # tau^2 ~ IG(r0/2,s0/2)
  tau2_mu = 1   # Prior Mean
  tau2_sd = 100  # Prior SD
  tau2_r0 = 2*((tau2_mu/tau2_sd)^2 + 2)
  tau2_s0 = 2*tau2_mu*(tau2_r0/2 - 1)
  tau2_rn = tau2_r0 + ntheta
  
  # ****************************************************************************
  # Initialize parameters
  # ****************************************************************************
  theta   = matrix(0,ntheta,1)
  betapar = matrix(0,nparv,1)
  sigma   = 1
  sigma2  = sigma^2
  tau     = 1
  tau2    = tau^2
  fxdata  = phixdata %*% theta
  fxgrid  = phixgrid %*% theta
  ptp     = crossprod(phixdata)
  rownames(betapar) = vnames
  
  # ****************************************************************************
  # Save MCMC iterations
  # ****************************************************************************
  thetag   = matrix(0,smcmc,ntheta)
  betag    = matrix(0,smcmc,nparv)
  sigmag   = matrix(0,smcmc,1)
  taug     = matrix(0,smcmc,1)
  fxdatag  = matrix(0,smcmc,nobs)
  fxdatam  = matrix(0,nobs,1)
  fxdatas  = fxdatam
  fxgridg  = matrix(0,smcmc,nint+1)
  fxgridm  = matrix(0,nint+1,1)
  fxgrids  = fxgridm
  yhatm    = matrix(0,nobs,1)
  yhats    = matrix(0,nobs,1)
  colnames(betag)  = vnames
  
  # ****************************************************************************
  # ****************************************************************************
  # Do MCMC ----
  # ****************************************************************************
  cpara = c(nobs,nparv,nblow,nmcmc,smcmc,nint,nbasis)
  
  start_time = Sys.time()
  mcmc_out = get_mcmc_BBS_Free(ydata=ydata,
                               vdata=vdata,
                               vtv=vtv,
                               ptp=ptp,
                               
                               fxdata=fxdata,
                               
                               phixgrid=phixgrid,
                               phixdata=phixdata,
                               
                               kall=kall,
                               
                               tau2_s0=tau2_s0,
                               tau2_rn=tau2_rn,
                               sigma2_s0=sigma2_s0,
                               sigma2_rn=sigma2_rn,
                               beta_v0i=beta_v0i,
                               beta_v0im0=beta_v0im0,
                               
                               cpara=cpara,
                               
                               beta=betapar,
                               theta=theta,
                               tau2=tau2,
                               tau=tau,
                               sigma2=sigma2,
                               sigma=sigma)
  end_time   = Sys.time()
  run_time   = end_time - start_time
  
  fxgridm = mcmc_out$fxgridm
  fxgrids = mcmc_out$fxgrids
  fxdatam = mcmc_out$fxdatam
  fxdatas = mcmc_out$fxdatas
  
  betag   = mcmc_out$mcmcg$betag
  thetag  = mcmc_out$mcmcg$thetag
  sigmag  = mcmc_out$mcmcg$sigmag
  taug    = mcmc_out$mcmcg$taug
  gammag  = mcmc_out$mcmcg$gammag
  fxgridg = mcmc_out$mcmcg$fxgridg
  
  # ****************************************************************************
  # Compute MCMC estimates
  # ****************************************************************************
  # Get HPD Intervals
  qs        = c(0.025,0.975)
  fxgridq   = t(apply(fxgridg,2,quantile,probs=qs))
  
  yhatm       = yhatm/smcmc
  yhats       = sqrt(abs(yhats-smcmc*yhatm^2)/smcmc)
  fxdatam     = fxdatam/smcmc
  fxdatas     = sqrt(abs(fxdatas-smcmc*fxdatam^2)/smcmc)
  fxgridm     = fxgridm/smcmc
  fxgrids     = sqrt(abs(fxgrids-smcmc*fxgridm^2)/smcmc)
  
  outnames    = c("PostMean","PostSD","P(>0)")
  betam       = apply(betag,2,mean)
  betas       = apply(betag,2,sd)
  beta_pg0    = apply(betag>0,2,mean)
  outbeta     = cbind(betam,betas,beta_pg0)
  
  thetam      = apply(thetag,2,mean)
  thetas      = apply(thetag,2,sd)
  theta_pg0   = apply(thetag>0,2,mean)
  outtheta    = cbind(thetam,thetas,theta_pg0)
  
  sigmam      = mean(sigmag)
  sigmas      = sd(sigmag)
  
  taum        = mean(taug)
  taus        = sd(taug)
  
  outsigma    = matrix(c(sigmam,sigmas))
  rownames(outsigma) = c("PostMean","PostSD")
  outtau      = matrix(c(taum,taus))
  rownames(outtau)  = c("PostMean","PostSD")
  
  
  
  colnames(outbeta)  = outnames
  colnames(outtheta) = outnames
  colnames(outsigma) = "Sigma"
  rownames(outbeta)  = vnames
  
  mcmc.out = list()
  mcmc.out$fxgridg = fxgridg
  mcmc.out$betag   = betag
  mcmc.out$thetag  = thetag
  mcmc.out$sigmag  = sigmag
  mcmc.out$taug    = taug
  
  post.est = list()
  post.est$fxgridq     = fxgridq
  post.est$yhatm       = yhatm
  post.est$yhats       = yhats
  post.est$wbetam      = vdata%*%betam
  post.est$fxdatam     = fxdatam
  post.est$fxdatas     = fxdatas
  post.est$fxgridm     = fxgridm
  post.est$fxgrids     = fxgrids
  post.est$outtheta    = outtheta
  post.est$outtau      = outtau
  post.est$outbeta     = outbeta
  post.est$outsigma    = outsigma
  
  out=list()
  out$mcmc.out=mcmc.out
  out$post.est=post.est
  
  out$y           = ydata
  out$v           = vdata
  out$x           = xdata
  out$yresid      = ydata-post.est$wbetam
  out$nknots      = nknots
  out$xgrid       = xgrid
  out$phixdata    = phixdata
  out$phixgrid    = phixgrid
  out$mcmc.time   = run_time
  
  return(out)
}