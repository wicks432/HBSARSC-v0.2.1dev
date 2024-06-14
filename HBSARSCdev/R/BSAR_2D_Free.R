# ******************************************************************************
# BSAR_2D_Free_v1.R
#   BSAR 2D without shape constraints
#   INPUT
#       save(ydata,xdata,vdata,file="Data2d.Rdata")
#      ydata = vector of dependent variables
#      xdata = n x 2 matrix of x values
#      vdata = fixed effects
#              NOTE vdata includes a column of ones for the intercept
#   y = vdata*beta + f(x1,x2) + epsilon
#   f(x1,x2)  = sum_j sum_k theta_{j,k}phi_j(x1)phi_k(x2)
#   theta_{j,k} = N(0,tau^2*exp(-(j+k)*gamma))
#   epsilon     = N(0,sigma2)
#   Put constant in vdata
#
#   Simulation program is gBSAR_D2_v1.R
#   If a simulation, parameter values overide nbasis, xgrid, etc
# ******************************************************************************
# Avoids arrays by vectorization
# Indexes theta where the fastest running index is leftmost, as in R 
#    1, 1
#    2, 1
#    3, 1
#    1, 2
#    2, 2
#    3, 2
#    1, 3
#    2, 3
#    3, 3
#    This vectorizes theta by columns in matrix(theta)
# ******************************************************************************
#  if nbasis = 3
#  theta_mat = [
#      theta_11, theta_12, theta_13,
#      theta_21, theta_22, theta_23,
#      theta_31, theta_32, theta_33
#  ]
#  theta = matrix(theta_mat) 
#  theta_11
#  theta_21
#  theta_31
#  theta_12
#  theta_22 
#  theta_32
#  theta_13
#  theta_23
#  theta_33
# ******************************************************************************
# The bivariate basis is the outer product of the univariate basis
# This gives dim of theta as nbasis^2, which can be large
# If f is smooth, then the coefficients for phi_j*phi_k are close
# to zero for moderate j & k.
# Program reduces the number of parameters by excluding 
# pairs of j and k that are large frequencies
# See rad and freq_all below
# ******************************************************************************
"BSAR_2D_Free" <- function(ydata,xdata,vdata=NULL,mcmc=list(),nbasis=20,
                           nint=100,xmin=NULL,xmax=NULL,iflagJoint=2) {
  
  # MCMC Iterations
  mcvals=list(nblow=20000,smcmc=10000)
  mcvals[names(mcmc)]=mcmc
  nblow  = mcvals$nblow # number of initial MCMC parameters
  smcmc  = mcvals$smcmc # number of MCMC parameters to save
  nint   = nint         # Number of intervals in grid
  nbasis = nbasis       # Number of cosine functions in each direction
  # ****************************************************************************
  # Parameters for grids of x1 and x2
  if (is.null(xmin)) xmin = apply(xdata, 2, min)  # Minimum of xgrid
  if (is.null(xmax)) xmax = apply(xdata, 2, max)  # Maximum of xgrid
  
  #     about 25% of all (j,k)
  kmax       = nbasis       # j + k <= kmax when computing bivarate basis
  ntheta1D   = nbasis +1    # Includes constant
  iflagJoint = iflagJoint   # Effects combinations of (j,k) in theta_{j,k}
  # 1 = all (j,k): Only if nbasis is smallish
  # 2 = Lower half of (j,k)
  # 3 = all (j,k) minus those within a circle
  
  nmcmc      = nblow+smcmc
  
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
  nparx    = ncol(xdata)
  if (nparx != 2) cat("\nERROR: xdata should be n x 2: nparx = ",nparx)
  xnames   = colnames(xdata)
  if (is.null(xnames)) {
    xnames = c("X1","X2")
    colnames(xdata) = xnames
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
  xgrid  = matrix(0,nint+1,nparx)
  xdelta = matrix(0,nparx,1)
  # Loop over dimensions of x to compute xgrid & xdata
  for (j in 1:nparx) {
    xdelta[j] = xrange[j]/nint
    xgrid[,j] = seq(xmin[j],xmax[j],xdelta[j])
  }  # end compute xgrid and univerate phi
  
  # ****************************************************************************
  # Get indicies for theta_{j,k,...}, vectorized xgrid
  kall       = 1:nbasis
  kall0      = matrix(c(0,kall))
  freq_all   = kall0
  for (j in 2:nparx) {
    f0       = matrix(1,ntheta1D,1) %x% freq_all 
    f1       = kall0 %x% matrix(1,nrow(freq_all),1)
    freq_all = cbind(f0,f1)
  }  # End loop over dimensions 2 to nparx
  # Remove constant  
  freq_all   = freq_all[-1,]              # Remove constant
  freq_sum   = apply(freq_all,1,sum)
  # ****************************************************************************
  # Remover higher-order interactions of basis functions
  # All (j,k) minus those inside a circle
  if (iflagJoint==3) {
    rad      = sqrt((nbasis-freq_all[,1])^2+(nbasis-freq_all[,2])^2)
    idin     = which(rad >= kmax)
    freq_sum = freq_sum[idin]
    freq_all = freq_all[idin,]
  }
  # Lower half of (j,k)
  if (iflagJoint==2) {
    idin       = which(freq_sum<=kmax)
    freq_sum   = freq_sum[idin]
    freq_all   = freq_all[idin,]
  }
  ntheta      = nrow(freq_all)
  theta_names = NULL
  for (kk in 1:ntheta) {
    j  = freq_all[kk,1]
    k  = freq_all[kk,2]
    theta_names = c(theta_names,paste("theta(",j,",",k,")"))
  }
  
  
  
  # ****************************************************************************
  # Compute Cosine functions for xdata
  # phi0(x)  = matrix(1/sqrt(xrange),nint+1,1)
  # phij(x)  = sqrt(2/xrange)*cos(k*pi*((x-xmin)/xrange))
  # fxdata   = phixdata %*% theta where theta is a vector
  # ****************************************************************************
  # ****************************************************************************
  # Suppose nbasis = 3.  Each row of phixdata is
  # phi_1(x) * phi_1(y)
  # phi_2(x) * phi_1(y)
  # phi_3(x) * phi_1(y)
  # phi_1(x) * phi_2(y)
  # phi_2(x) * phi_2(y)
  # phi_3(x) * phi_2(y)
  # phi_1(x) * phi_3(y)
  # phi_2(x) * phi_3(y)
  # phi_3(x) * phi_3(y)
  # for each observed pair (x,y)
  # ****************************************************************************
  #  theta = matrix(theta_mat) 
  #  theta_11
  #  theta_21
  #  theta_31
  #  theta_12
  #  theta_22 
  #  theta_32
  #  theta_13
  #  theta_23
  #  theta_33
  # ****************************************************************************
  phix0    = 1/sqrt(xrange[1])
  phiy0    = 1/sqrt(xrange[2])
  phixdata = matrix(0,nobs,ntheta)
  for (i in 1:nobs) {
    xi        = xdata[i,1]
    yi        = xdata[i,2]
    # Loop over indices where i+j <= kmax
    for (kk in 1:ntheta) {
      k1      = freq_all[kk,1]   # Frequency for xi
      k2      = freq_all[kk,2]   # Frequency for yi
      if(k1==0){
        phixi = phix0
      } else {
        phixi = sqrt(2/xrange[1])*cos(k1*pi*(xi-xmin[1])/xrange[1])
      }
      if (k2==0) {
        phiyi = phiy0
      } else {
        phiyi = sqrt(2/xrange[2])*cos(k2*pi*(yi-xmin[2])/xrange[2])
      }
      phixdata[i,kk] = phixi*phiyi
    }  # end loop over indices
  } # End lover over data
  
  # ****************************************************************************
  # Compute Cosine functions for xdata
  # phi0(x) = matrix(1/sqrt(xrange),nint+1,1)
  # phij(x)  = sqrt(2/xrange)*cos(k*pi*((x-xmin)/xrange))
  # fxdata   = phixdata %*% theta where theta is a vector
  # ****************************************************************************
  phixgrid = matrix(0,(nint+1)^2,ntheta)
  rcpp_get_phixgrid_2D(phixgrid, xgrid, xrange, xmin, freq_all, phix0, phiy0, ntheta, nint)
  
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
  # gamma ~ Exponential(w0)
  w0   = 1
  wk   = sum(freq_sum)/2
  
  # ****************************************************************************
  # Initialize parameters
  # ****************************************************************************
  theta   = matrix(0,ntheta,1)
  betapar = matrix(0,nparv,1)
  sigma   = 1
  sigma2  = sigma^2
  tau     = 1
  tau2    = tau^2
  gampar  = 1
  fxdata  = phixdata %*% theta
  fxgrid  = phixgrid %*% theta
  fxgridm_mat = matrix(fxgrid,nint+1,nint+1)
  ptp     = crossprod(phixdata)
  rownames(theta)   = theta_names
  rownames(betapar) = vnames
  egamkall          = exp(-freq_sum*gampar)
  # Worry about egamkall getting too small
  iz0 = which(egamkall < 1E-15)   # Guard against numerical errors
  if(length(iz0)>0) egamkall[iz0] = 1E-15
  
  # ****************************************************************************
  # Save MCMC iterations
  # ****************************************************************************
  thetag   = matrix(0,smcmc,ntheta)
  betag    = matrix(0,smcmc,nparv)
  sigmag   = matrix(0,smcmc,1)
  taug     = matrix(0,smcmc,1)
  gammag   = matrix(0,smcmc,1)
  fxdatam  = matrix(0,nobs,1)
  fxdatas  = fxdatam
  fxgridm  = matrix(0,(nint+1)^2,1)
  fxgrids  = fxgridm
  yhatm    = matrix(0,nobs,1)
  yhats    = matrix(0,nobs,1)
  colnames(thetag) = theta_names
  colnames(betag)  = vnames
  
  # ****************************************************************************
  # ****************************************************************************
  # Do MCMC ----
  # ****************************************************************************
  cpara = c(nobs,nparv,nblow,nmcmc,smcmc,nint,ntheta)
  
  start_time = Sys.time()
  mcmc_out = get_mcmc_2DBSAR_Free(ydata=ydata,
                                  vdata=vdata,
                                  vtv=vtv,
                                  ptp=ptp,
                                  
                                  fxdata=fxdata,
                                  
                                  phixgrid=phixgrid,
                                  phixdata=phixdata,
                                  
                                  wk=wk,
                                  egamkall=egamkall,
                                  freq_sum=freq_sum,
                                  
                                  tau2_s0=tau2_s0,
                                  tau2_rn=tau2_rn,
                                  sigma2_s0=sigma2_s0,
                                  sigma2_rn=sigma2_rn,
                                  beta_v0i=beta_v0i,
                                  beta_v0im0=beta_v0im0,
                                  
                                  cpara=cpara,
                                  
                                  beta=betapar,
                                  theta=theta,
                                  gampar=gampar,
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
  
  betag  = mcmc_out$mcmcg$betag
  thetag = mcmc_out$mcmcg$thetag
  sigmag = mcmc_out$mcmcg$sigmag
  taug   = mcmc_out$mcmcg$taug
  gammag = mcmc_out$mcmcg$gammag
  
  # ****************************************************************************
  # Compute MCMC estimates
  # ****************************************************************************
  yhatm       = yhatm/smcmc
  yhats       = sqrt(abs(yhats-smcmc*yhatm^2)/smcmc)
  fxdatam     = fxdatam/smcmc
  fxdatas     = sqrt(abs(fxdatas-smcmc*fxdatam^2)/smcmc)
  fxgridm     = fxgridm/smcmc
  fxgrids     = sqrt(abs(fxgrids-smcmc*fxgridm^2)/smcmc)
  fxgridm_mat = matrix(fxgridm,nint+1,nint+1)
  fxgrids_mat = matrix(fxgrids,nint+1,nint+1)
  
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
  
  gammam      = mean(gammag)
  gammas      = sd(gammag)
  outsigma    = matrix(c(sigmam,sigmas))
  rownames(outsigma) = c("PostMean","PostSD")
  outsmooth   = c(taum,gammam,taus,gammas)
  outsmooth   = matrix(outsmooth,2,2)
  rownames(outsmooth)  = c("Tau","Gamma")
  colnames(outsmooth)  = c("PostMean","PostSD")
  
  
  
  colnames(outbeta)  = outnames
  colnames(outtheta) = outnames
  colnames(outsigma) = "Sigma"
  rownames(outbeta)  = vnames
  rownames(outtheta) = theta_names
  
  mcmc.out = list()
  mcmc.out$betag  = betag
  mcmc.out$thetag = thetag
  mcmc.out$sigmag = sigmag
  mcmc.out$taug   = taug
  mcmc.out$gammag = gammag
  
  post.est = list()
  post.est$yhatm       = yhatm
  post.est$yhats       = yhats
  post.est$wbetam      = vdata%*%betam
  post.est$fxdatam     = fxdatam
  post.est$fxdatas     = fxdatas
  post.est$fxgridm     = fxgridm
  post.est$fxgrids     = fxgrids
  post.est$fxgridm_mat = fxgridm_mat
  post.est$fxgrids_mat = fxgrids_mat
  post.est$outtheta    = outtheta
  post.est$outsmooth   = outsmooth
  post.est$outbeta     = outbeta
  post.est$outsigma    = outsigma
  
  x1mat = xgrid[,1] %*% matrix(1,1,nint+1)
  x2mat = matrix(1,nint+1,1) %*% t(xgrid[,2])
  
  out=list()
  out$mcmc.out=mcmc.out
  out$post.est=post.est
  
  out$y           = ydata
  out$v           = vdata
  out$x           = xdata 
  out$x1mat       = x1mat
  out$x2mat       = x2mat
  out$nbasis      = nbasis
  out$xgrid       = xgrid
  out$mcmc.time   = run_time
  
  return(out)
}