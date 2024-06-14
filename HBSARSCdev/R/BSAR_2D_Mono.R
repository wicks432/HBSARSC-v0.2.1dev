# ******************************************************************************
# BSAR_2D_Mono_v1.R
#   Monotone in 2D
#   INPUT
#      save(ydata,xdata,vdata,file="Data_Mono2D.Rdata")
#      ydata  = n x 1 dependent variable
#      xdata  = n x 2 independent variables
#      vdata  = Covariates
#               NOTE: vdata should have a column of ones
#   y = vdata*beta + f(x1,x2) + epsilon
#   z(x1,x2)  = theta_{0,0} + sum_j sum_k theta_{j,k}phi_j(x1)phi_k(x2)
#   Functions ANOVA 
#       f(x1,x2) = f1(x1) + f2(x2) + f12(x1,x2)
#   f(x1,x2)  = d*int_x1 int_x2 g(z(s1,s2)) ds1 ds2
#   d         = plus or minus one
#   g(.)      = exp or square
#   theta_{j,k} = N(0,tau^2*exp(-(j+k)*gamma))
#   epsilon     = N(0,sigma2)
#   Put constant in vdata
#
# ******************************************************************************
#   Note: compute f(data) by approximation from the f(grid)
# ******************************************************************************
#
#   Simulation program is gBSAR_D2_Mono_v1.R
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
"BSAR_2D_Mono" <- function(ydata,xdata,vdata=NULL,mcmc=list(),nbasis=20,
                           nint=100,xmin=NULL,xmax=NULL,iflagJoint=2,
                           iflagZ=0,iflagCenter=1,iflagpn=1,iflagMainpn=c(1,1)) {
  
  # MCMC Iterations
  mcvals=list(nblow=20000,smcmc=10000,nskip=1)
  mcvals[names(mcmc)]=mcmc
  nblow  = mcvals$nblow # number of initial MCMC parameters
  smcmc  = mcvals$smcmc # number of MCMC parameters to save
  nskip  = mcvals$nskip
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
  iflagZ      = iflagZ      # 0 = Model using Z^2 or 1 = exp(Z)
  iflagCenter = iflagCenter # make int f = 0
  
  iflagpn     = iflagpn     # 1 or -1 for interaction effects
  iflagMainpn = iflagMainpn # 1 or -1 for main effects
  
  nmcmc       = nblow+nskip*smcmc
  # shape constraints: use adaptive Metropolis
  maxmodmet = 10          # Maximum times to modify Metropolis
  nblow0    = 1000         # Initial number of initial MCMC parmeters
  nmcmcall  = maxmodmet*nblow0 + nmcmc
  
  # ****************************************************************************
  # Get Data: vdata with constant, xdata (x1,x2),ydata
  # ****************************************************************************
  nobs     = length(ydata)
  ynames   = colnames(ydata)
  if(is.null(ynames)){
    ynames  = "Y"
    ydata   = matrix(ydata)
    colnames(ydata) = ynames
  }
  nparx    = ncol(xdata)
  if(nparx != 2) cat("\nERROR: xdata should be n x 2: nparx = ",nparx)
  xnames   = colnames(xdata)
  if(is.null(xnames)){
    xnames = c("X1","X2")
    colnames(xdata) = xnames
  }
  
  if(is.null(vdata)){
    vdata           = matrix(1,nobs,1)
    vnames          = "CNST"
    colnames(vdata) = "CNST"
    nparv           = 1
  }else{
    nparv           = ncol(vdata)
    vnames          = colnames(vdata)
    if(is.null(vnames)){
      vnames = "CNST"
      if(nparv>1){
        vnames = c(vnames,paste("V",1:(nparv-1),sep=""))
      }
      colnames(vdata)  = vnames
    }
  }
  vtv             = crossprod(vdata)
  
  # ****************************************************************************
  # ****************************************************************************
  # Real data: Compute grids, phix, etc
  xmin   = floor(apply(xdata,2,min))
  xmax   = ceiling(apply(xdata,2,max))
  xrange = xmax-xmin
  # ****************************************************************************
  # Compute xgrid for plotting
  xgrid     = matrix(0,nint+1,nparx)
  xdelta    = matrix(0,nparx,1)
  # Loop over dimensions of x to compute xgrid & xdata
  for(j in 1:nparx){
    xdelta[j] = xrange[j]/nint
    xgrid[,j] = seq(xmin[j],xmax[j],xdelta[j])
  }  # end compute xgrid and univerate phi
  
  # ****************************************************************************
  # Get indicies for theta_{j,k,...}, vectorized xgrid
  kall         = 1:nbasis
  kall0        = matrix(c(0,kall))
  freq_all     = kall0
  xgrid_vec    = matrix(xgrid[,1])
  for(j in 2:nparx){
    f0           = matrix(1,ntheta1D,1) %x% freq_all 
    f1           = kall0 %x% matrix(1,nrow(freq_all),1)
    freq_all     = cbind(f0,f1)
    x0           = matrix(1,nint+1,1) %x% xgrid_vec
    x1           = xgrid[,j] %x% matrix(1,nrow(xgrid_vec),1) 
    xgrid_vec    = cbind(x0,x1)
  }  # End loop over dimensions 2 to nparx
  # Remove constant  
  #freq_all = freq_all[-1,]              # Remove constant
  freq_sum = apply(freq_all,1,sum)
  # ****************************************************************************
  # Remover higher-order interactions of basis functions
  # All (j,k) minus those inside a circle
  if(iflagJoint==3){
    rad        = sqrt((nbasis-freq_all[,1])^2+(nbasis-freq_all[,2])^2)
    idin       = which(rad >= kmax )
    freq_sum   = freq_sum[idin]
    freq_all   = freq_all[idin,]
  }
  # Lower half of (j,k)
  if(iflagJoint==2){
    idin       = which(freq_sum<=kmax)
    freq_sum   = freq_sum[idin]
    freq_all   = freq_all[idin,]
  }
  ntheta     = nrow(freq_all)
  theta_names = NULL
  for(kk in 1:ntheta){
    j  = freq_all[kk,1]
    k  = freq_all[kk,2]
    theta_names  = c(theta_names,paste("theta(",j,",",k,")"))
  }
  
  
  # ****************************************************************************
  # Compute Cosine functions for xgrid
  # phi0(x)  = matrix(1/sqrt(xrange),nint+1,1)
  # phij(x)  = sqrt(2/xrange)*cos(k*pi*((x-xmin)/xrange))
  # ****************************************************************************
  # ****************************************************************************
  # Suppose nbasis = 3.  Each row of phixgrid is
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
  zxg     = (xgrid[,1] - xmin[1])/xrange[1]
  phix1   = sqrt(2/xrange[1])*cos(pi*zxg %*% t(kall))
  zxg     = (xgrid[,2] - xmin[2])/xrange[2]
  phix2   = sqrt(2/xrange[2])*cos(pi*zxg %*% t(kall))
  # Add constant
  phix1   = cbind(matrix(1/sqrt(xrange[1]),nint+1,1),phix1)
  phix2   = cbind(matrix(1/sqrt(xrange[2]),nint+1,1),phix2)
  phixgrid   = matrix(0,(nint+1)^2,ntheta)
  icount = 0
  for(k in 1:(nint+1)){
    phik   = matrix(phix2[k,])
    phik   = matrix(phik[freq_all[,2]+1])  # Row is freq + 1
    for(j in 1:(nint+1)){
      icount = icount + 1
      phij   = matrix(phix1[j,])
      phij   = matrix(phij[freq_all[,1]+1])  # Row is freq + 1
      phijk   = matrix(phik * phij)
      phixgrid[icount,] = t(phijk)
    }  # End loop over x1
  }  # End loop over x2
  
  
  # ****************************************************************************
  # Basis functions for "main" effects
  kall      = matrix(1:nbasis)
  phiMain   = array(0, dim=c(nint+1,ntheta1D,nparx))
  for(j in 1:nparx){
    zx        = matrix((xgrid[,j]-xmin[j])/xrange[j])   
    phix      = sqrt(2/xrange[j])*cos(pi*zx %*% t(kall))
    phix      = matrix(phix,nint+1,nbasis)
    # Add constants
    phix     = cbind(matrix(sqrt(1/xrange[j]),nint+1,1),phix)
    phiMain[,,j] = phix
  }
  
  theta1_names = paste("Theta1(",c(0,kall),")",sep="")
  theta2_names = paste("Theta2(",c(0,kall),")",sep="")
  
  # ****************************************************************************
  # Locate indices that define data (x1,x2)
  #   Instead of computing f(x_i) at observations,
  #   find where x_i falls in xgrid and use fxgrid
  #   Find location of xdata in xgrid
  #   xgrid[xinx[i,1],1] < xi1 <= xgrid[xinx[i,1]+1,1]
  #   xgrid[xinx[i,2],2] < xi2 <= xgrid[xinx[i,2]+1,2]
  # ****************************************************************************
  # ****************************************************************************
  # Interaction function f12(x1,x2)
  xinx = matrix(0,nobs,nparx)   # indices 
  # Loop over observations
  for (i in 1:nobs) {  
    # Loop over dimensions
    for (k in 1:nparx) {  
      xik   = xdata[i,k]
      if (xik == xmax[k]) {
        xinx[i,k] = nint+1
      }else{
        xinx[i,k] = max(which(xik >= xgrid[,k]))
      }
    }# End loop over dimensions
  } # End loop over observations
  
  xinxMain  = array(0,c(nobs,2,nparx))
  xoverMain = array(0,c(nobs,nparx))   # Relative excess over grid
  for (k in 1:nparx) {
    # Loop over observations
    for (i in 1:nobs) { 
      xi  = xdata[i,k]
      idx = which(xi == xgrid[,k])  # Is it exactly on xgrid?
      if (length(idx) > 0 ) {
        xinxMain[i,1,k] = idx
        xinxMain[i,2,k] = idx
        xoverMain[i,k]  = 0 
      } else {
        idx = max(which(xi > xgrid[,k]))
        xinxMain[i,1,k] = idx
        xinxMain[i,2,k] = idx+1
        xoverMain[i,k]  = (xi - xgrid[idx,k])/(xdelta[k])
      }
    } # End loop over observations for xinx and xover
  }
  
  # Get excess for that xdata is over boundary of xgrid based on xinx
  # Used in numerical integration
  xover = matrix(0,nobs,nparx)   # excess of data over grid
  for(k in 1:nparx){
    xover[,k] = xdata[,k] - xgrid[xinx[,k],k]   # excess of data over grid
    xover[,k] = xover[,k]/xdelta[k]             # proportion excess relative to grid
  }
  # Approximation is weighted average data points
  
  f12apx_weights = cbind(
    (1-xover[,1])*(1-xover[,2]),      # Weights for f00  
    xover[,1]*(1-xover[,2]),          # weights for f10
    (1-xover[,1])*xover[,2],          # weights for f01
    xover[,1]*xover[,2]               # weights for f11
  )
  
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
  wkMain = sum(kall)/2
  
  # ****************************************************************************
  # Initialize parameters
  # ****************************************************************************
  thetaMain     = matrix(0,ntheta1D,2)  # Theta for main effects function
  thetaMain[1,] = 1
  theta         = matrix(0,ntheta,1)    # Theta for interaction function
  theta[1]      = 1
  xipar         = 0                   # theta_0 = exp(xipar) in Chi-Square model 
  xiparMain     = matrix(0,2,1)
  
  betapar  = matrix(0,nparv,1)
  sigma    = 1
  sigma2   = sigma^2
  tau      = 1
  tau2     = tau^2
  gampar   = 1
  tauMain  = matrix(1,nparx,1)
  tau2Main = tauMain^2
  gamMain  = matrix(1,nparx,1)
  
  
  # Compute f at xgrid and data
  # Interaction function
  zxgrid  = phixgrid %*% theta
  if (iflagZ == 1) {
    ggrid   = exp(zxgrid)
  } else {
    ggrid   = zxgrid^2
  }
  ggrid_mat    = matrix(ggrid,nint+1,nint+1)
  fx12grid_mat = iflagpn*rcpp_CumInt2D(ggrid_mat,xdelta)
  # Is f center?
  if (iflagCenter == 1) {
    c            = rcpp_Int2D(fx12grid_mat,xdelta) # Total integral of f over xgrid
    c            = c/(xrange[1]*xrange[2])
    fx12grid_mat = fx12grid_mat - c
  }
  # Main effects
  zMain        = matrix(0,nint+1,nparx)
  gMain        = zMain
  fxMain       = matrix(0,nint+1,nparx)
  fxMain_mat   = array(0, dim=c(nint+1,nint+1,nparx))
  for (j in 1:nparx) {
    zMain[,j]  = phiMain[,,j] %*% thetaMain[,j]
    if(iflagZ == 1){
      gMain[,j] = exp(zMain[,j])
    }else{
      gMain[,j] = zMain[,j]^2
    }
    fxMain[,j]  = iflagMainpn[j]*matrix(rcpp_CumTrap(gMain[,j],xdelta[j]))
    if (iflagCenter == 1) {
      cj         = rcpp_trapint(fxMain[,j],xdelta[j])
      fxMain[,j] = fxMain[,j] - cj/xrange[j]
    }  # End Center
    if (j == 1) {
      fxMain_mat[,,j] = gMain[,j] %*% matrix(1,1,nint+1)
    } else {
      fxMain_mat[,,j] = matrix(1,nint+1,1) %*% t(fxMain[,j])
    }
  }  # End over main effects 
  fxgrid_mat   = fx12grid_mat + fxMain_mat[,,1] + fxMain_mat[,,2]
  fx12data     = rcpp_fdata_apx(fx12grid_mat,f12apx_weights,xinx-1,xover,nint)   # f12 at data
  fxMaindata   = matrix(0,nobs,2)
  fxdata       = fx12data
  for (j in 1:nparx) {
    fxMaindata[,j] = rcpp_GetSCfxobs(xinxMain[,,j]-1,xoverMain[,j],fxMain[,j])   #  f1 at data
    fxdata         = fxdata + fxMaindata[,j]
  }
  
  rownames(theta)   = theta_names
  rownames(betapar) = vnames
  egamkall          = exp(-freq_sum*gampar)
  egamkallMain      = cbind(exp(-kall0*gamMain[1]),exp(-kall0*gamMain[2]))
  
  
  # ****************************************************************************
  # Save MCMC iterations
  # ****************************************************************************
  thetag           = matrix(0,smcmc,ntheta)
  thetaMaing       = array(0,dim=c(smcmc,ntheta1D,2))
  betag            = matrix(0,smcmc,nparv)
  sigmag           = matrix(0,smcmc,1)
  taug             = matrix(0,smcmc,1)
  tauMaing         = matrix(0,smcmc,2)
  gammag           = matrix(0,smcmc,1)
  gamMaing         = matrix(0,smcmc,2)
  fxdatam          = matrix(0,nobs,1)
  fxdatas          = fxdatam
  fxgridm_mat      = matrix(0,(nint+1),nint+1)
  fxgrids_mat      = fxgridm_mat
  yhatm            = matrix(0,nobs,1)
  yhats            = matrix(0,nobs,1)
  colnames(thetag) = theta_names
  colnames(betag)  = vnames
  
  # ****************************************************************************
  # ****************************************************************************
  # Use Metropolis-Hastings with shape constraints
  # Random walk Metropolis-Hastings with t-errors for fmodel > 1
  # pmet = P(accept candidate theta)
  # If pmet is > 0.7, then step size is too small: increase metm 
  # If pmet is < 0.3, then step size is too large: decrease metm 
  # Keep mets > metm  (std dev > mean)
  # Why?  If mets > metm, then alpha < 3
  # Also IG has a dead zone near zeros if alpha > 2.25
  # So, you want to keep 2 < alpha < 2.25.
  # Fix alpha.  Focus on metm to get mets and beta
  # ****************************************************************************
  # Adaptive Metropolis Hastings for theta
  theta_met     = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=ntheta)
  theta1_met    = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=ntheta1D)
  theta2_met    = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=ntheta1D)
  thetaMain_met = cbind(theta1_met,theta2_met)
  
  # ****************************************************************************
  # ****************************************************************************
  # Do MCMC ----
  # ****************************************************************************
  cpara = c(nobs,nparv,nblow,nblow0,maxmodmet,nskip,nmcmcall,smcmc,nint,ntheta,
            nbasis,nparx,ntheta1D,iflagpn,iflagCenter,iflagZ)
  met_para = list(thetaMain_met=thetaMain_met,theta_met=theta_met)
  start_time = Sys.time()
  mcmc_out = get_mcmc_2DBSAR_Mono(ydata=ydata,
                                  vdata=vdata,
                                  vtv=vtv,
                                  iflagMainpn=iflagMainpn,
                                  
                                  xinxMain=xinxMain-1,
                                  xinx=xinx-1,
                                  xoverMain=xoverMain,
                                  xover=xover,
                                  xdelta=xdelta,
                                  xrange=xrange,
                                  
                                  fxdata=fxdata,
                                  fxMaindata=fxMaindata,
                                  fx12data=fx12data,
                                  fx12grid_mat=fx12grid_mat,
                                  fxMain_mat=fxMain_mat,
                                  f12apx_weights=f12apx_weights,
                                  
                                  phixgrid=phixgrid,
                                  phiMain=phiMain,
                                  
                                  wk=wk,
                                  wkMain=wkMain,
                                  egamkallMain=egamkallMain,
                                  egamkall=egamkall,
                                  freq_sum=freq_sum,
                                  kall0=kall0,
                                  
                                  tau2_s0=tau2_s0,
                                  tau2_rn=tau2_rn,
                                  sigma2_s0=sigma2_s0,
                                  sigma2_rn=sigma2_rn,
                                  beta_v0i=beta_v0i,
                                  beta_v0im0=beta_v0im0,
                                  
                                  cpara=cpara,
                                  met_para=met_para,
                                  
                                  beta=betapar,
                                  thetaMain=thetaMain,
                                  theta=theta,
                                  xiparMain=xiparMain,
                                  xipar=xipar,
                                  gamMain=gamMain,
                                  gampar=gampar,
                                  tau2Main=tau2Main,
                                  tau2=tau2,
                                  sigma2=sigma2,
                                  sigma=sigma)
  
  end_time = Sys.time()
  run_time = end_time - start_time
  
  
  fxgridm_mat = mcmc_out$fxgridm_mat
  fxgrids_mat = mcmc_out$fxgrids_mat
  yhatm       = mcmc_out$yhatm
  yhats       = mcmc_out$yhats
  fxdatam     = mcmc_out$fxdatam
  fxdatas     = mcmc_out$fxdatas
  
  betag      = mcmc_out$mcmcg$betag
  thetag     = mcmc_out$mcmcg$thetag
  thetaMaing = mcmc_out$mcmcg$thetaMaing
  sigmag     = mcmc_out$mcmcg$sigmag
  taug       = mcmc_out$mcmcg$taug
  tauMaing   = mcmc_out$mcmcg$tauMaing
  gammag     = mcmc_out$mcmcg$gammag
  gamMaing   = mcmc_out$mcmcg$gamMaing
  
  theta_met = get_mat_to_list_met(mcmc_out$metg$theta_met)
  thetaMain_met_tmp = vector("list", length = 2)
  thetaMain_met_tmp[[1]] = get_mat_to_list_met(mcmc_out$metg$thetaMain_met[,1])
  thetaMain_met_tmp[[2]] = get_mat_to_list_met(mcmc_out$metg$thetaMain_met[,2])
  thetaMain_met = thetaMain_met_tmp
  
  # ****************************************************************************
  # Compute MCMC estimates
  # ****************************************************************************
  yhatm       = yhatm/smcmc
  yhats       = sqrt(abs(yhats-smcmc*yhatm^2)/smcmc)
  fxdatam     = fxdatam/smcmc
  fxdatas     = sqrt(abs(fxdatas-smcmc*fxdatam^2)/smcmc)
  fxgridm_mat = fxgridm_mat/smcmc
  fxgrids_mat = sqrt(abs(fxgrids_mat-smcmc*fxgridm_mat^2)/smcmc)
  theta_met$pmet = theta_met$pmet/smcmc
  thetaMain_met[[1]]$pmet = thetaMain_met[[1]]$pmet/smcmc
  thetaMain_met[[2]]$pmet = thetaMain_met[[2]]$pmet/smcmc
  
  met.out = list()
  met.out$theta_met     = theta_met
  met.out$thetaMain_met = thetaMain_met
  
  outnames    = c("PostMean","PostSD","P(>0)")
  betam       = apply(betag,2,mean)
  betas       = apply(betag,2,sd)
  beta_pg0    = apply(betag>0,2,mean)
  outbeta     = cbind(betam,betas,beta_pg0)
  
  thetam      = apply(thetag,2,mean)
  thetas      = apply(thetag,2,sd)
  theta_pg0   = apply(thetag>0,2,mean)
  outtheta    = cbind(thetam,thetas,theta_pg0)
  
  
  theta1_names = paste("Theta1(",kall0,")",sep="")
  theta1m      = matrix(apply(thetaMaing[,,1],1,mean))
  theta1s      = matrix(apply(thetaMaing[,,1],1,sd))
  theta1_pg0   = matrix(apply(thetaMaing[,,1]>0,1,mean))
  outtheta1    = cbind(theta1m,theta1s,theta1_pg0)
  colnames(outtheta1) = c("PostMean","PostSD","P(>0)")
  rownames(outtheta1) = theta1_names
  
  theta2_names = paste("Theta2(",kall0,")",sep="")
  theta2m      = matrix(apply(thetaMaing[,,2],1,mean))
  theta2s      = matrix(apply(thetaMaing[,,2],1,sd))
  theta2_pg0   = matrix(apply(thetaMaing[,,2]>0,1,mean))
  outtheta2    = cbind(theta2m,theta2s,theta2_pg0)
  colnames(outtheta2) = c("PostMean","PostSD","P(>0)")
  rownames(outtheta2) = theta2_names
  
  sigmam      = mean(sigmag)
  sigmas      = sd(sigmag)
  
  taum        = mean(taug)
  taus        = sd(taug)
  gammam      = mean(gammag)
  gammas      = sd(gammag)
  tauMainm    = apply(tauMaing,2,mean)
  tauMains    = apply(tauMaing,2,sd)
  gamMainm    = apply(gamMaing,2,mean)
  gamMains    = apply(gamMaing,2,sd)
  a12         = matrix(c(taum,taus),1,2)
  b12         = matrix(c(gammam,gammas),1,2)
  a1          = cbind(tauMainm,tauMains)
  b1          = cbind(gamMainm,gamMains)
  outsmooth   = rbind(a1,a12,b1,b12)
  colnames(outsmooth) = c("PostMean","PostSTD")
  rownames(outsmooth) = c("tau_1","tau_2","tau_12",
                          "gamma_1","gamma_2","gamma_12")
  
  outsigma = matrix(c(sigmam,sigmas))
  rownames(outsigma) = c("PostMean","PostSD")
  
  
  
  
  colnames(outbeta)  = outnames
  colnames(outtheta) = c("PostMean","PostSD","P>0")
  colnames(outsigma) = "Sigma"
  rownames(outbeta)  = vnames
  rownames(outtheta) = theta_names
  
  
  mcmc.out = list()
  mcmc.out$betag      = betag
  mcmc.out$thetag     = thetag
  mcmc.out$thetaMaing = thetaMaing
  mcmc.out$sigmag     = sigmag
  mcmc.out$taug       = taug
  mcmc.out$tauMaing   = tauMaing
  mcmc.out$gammag     = gammag
  mcmc.out$gamMains   = gamMains
  
  post.est = list()
  post.est$yhatm       = yhatm
  post.est$yhats       = yhats
  post.est$wbetam      = vdata%*%betam
  post.est$fxdatam     = fxdatam
  post.est$fxdatas     = fxdatas
  post.est$fxgridm_mat = fxgridm_mat
  post.est$fxgrids_mat = fxgrids_mat
  post.est$outtheta    = outtheta
  post.est$outsmooth   = outsmooth
  post.est$outbeta     = outbeta
  post.est$outsigma    = outsigma
  
  x1mat = xgrid[,1] %*% matrix(1,1,nint+1)
  x2mat = matrix(1,nint+1,1) %*% t(xgrid[,2])
  
  out=list()
  out$met.out=met.out
  out$mcmc.out=mcmc.out
  out$post.est=post.est
  
  out$y           = ydata
  out$v           = vdata
  out$x           = xdata 
  out$x1mat       = x1mat
  out$x2mat       = x2mat
  out$nbasis      = nbasis
  out$iflagZ      = iflagZ
  out$iflagCenter = iflagCenter
  out$iflagMainpn = iflagMainpn
  out$xgrid       = xgrid
  out$mcmc.time   = run_time
  
  return(out)
}