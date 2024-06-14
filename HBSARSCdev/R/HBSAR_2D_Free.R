# ******************************************************************************
# HBSAR_2D_v1
#   
#   HB Spectral Analysis Regression 
#   2D function f(x1,x2)
#   OUPUT DATA
#      Make a .Rdata file
#           save(ydata,xdata,vdata,wdata,zdata,id_group, file="HBSARM.Rdata")
#      ydata    = vector of dependent variable
#      xdata    = 2 columns for x1 and x2 
#      vdata    = OPTIONAL. covariaties for fixed effects
#      wdata    = matrix for random effects
#                 NOTE: wdata includes vector of ones for the random intercepts
#      zdata    = matrix for upper level model for group-level variables
#                 NOTE: zdata includes vector of ones for intercepts
#      id_group = vector that gives group number for each observations
#   Lower-level model for group j
#         Yj         = v*alpha + w*beta_j  + f_j(x1,x2) + epsilon
#
#         Bivariate Spectral representation in each group
#         Z_j(x1,x2)  = sum_k1 sum_k2 theta_{j,k1,k2} phi_k1(x1) phi_k2(x2)
#         Free f:
#              f_j(x1,x2) = Z_j(x1,x2)
#         Monotone f:
#              f_j(x1,x2) = int_a1^x1 int_a2^x2 g(Z_(x1,x2)) dx1 dx2
#              where g(z) = z^2 (Chi-Squared Model) or = exp(z)  (Log-Normal Model)
#         epsilon_ij ~ N(0,sigma^2) or N(0,sigma_j^2)
#         alpha is fixed effect, beta_j is random effect
#
#     HB model for spectral coefficients
#       Lower Level Model for j = 1, .., J
#          theta_{j,k1,k2} = theta_{0,k1,ki2} + N(0,tau_{k1,k2}^2*exp(-(k1+k2)*gamma_j))
#          Priors:  tau_{k1,k2} ~ IG
#                   gamma_j ~ G(xi_gamma,zeta_gamma)
#                   In program xi_gamma is gamma_alpha and zeta_gamma is gamma_beta 
#       Upper Level Model
#          theta_{0,k1,k2} = N(0,N(0,eta^2(1+(k1+k2)/zeta_gamma)^{-eta_gamma}))
#          where eta_gamma and zeta_gamma are prior parameters for gamma_j
#
#      HB Model for beta_j
#         beta_j     = Phi'z_j + delta_j  were delta_j ~ N(0,Lambda)
#      HB Error variance 
#         sigma_j^2 ~ IG(sigma2_alpha/2,sigma2_beta/2)
#
# ******************************************************************************
# ******************************************************************************
"HBSAR_2D_Free" <- function(ydata,xdata,vdata=NULL,wdata=NULL,zdata=NULL,id_group=NULL,
                            mcmc=list(), nbasis=20, nint=100, xmin=NULL, xmax=NULL,
                            iflagJoint   = 2, # Effects combinations of (j,k) in theta_{j,k}
                            # 1 = all (j,k): Only if nbasis is smaller
                            # 2 = Lower half of (j,k)
                            # 3 = all (j,k) minus those within a circle
                            iflagCenter  = 1, # Used with shape constraints: 1 if int f = 0 and 0 if f(xmin) = 0
                            iflagZ       = 0, # 0 = Model using Z^2 or 1 = exp(Z)
                            iflagHBsigma = 1, # 1 for HB sigma and 0 else
                            nExtreme     = 2,
                            shape=c('Free','Increasing','Decreasing')) {
  
  if (!is.matrix(ydata)) ydata = as.matrix(ydata)
  if (ncol(xdata) != 2) stop()  
  
  # MCMC Iterations
  mcvals=list(nblow=20000,smcmc=10000,nskip=1)
  mcvals[names(mcmc)]=mcmc
  nblow     = mcvals$nblow        # number of initial MCMC parameters
  smcmc     = mcvals$smcmc        # number of MCMC parameters to save
  nskip     = mcvals$nskip        # Save very nskip iterations
  nmcmc     = nblow + nskip*smcmc # total number of MCMC parameters
  nbasis    = nbasis              # Number of knots
  kmax      = nbasis              # Include phi_j*phi_k if j + k < xmax
  
  # ****************************************************************************
  # Parameters for grids of x1 and x2
  if (is.null(xmin)) xmin = apply(xdata, 2, min)  # Minimum of xgrid
  if (is.null(xmax)) xmax = apply(xdata, 2, max)  # Maximum of xgrid
  
  
  # ****************************************************************************
  # Set Parameters
  # ifalgsc gives free or shaped constraint
  # 0 is free f
  # 1 is monotone
  # ****************************************************************************
  if (shape == "Free") {
    iflagsc = 0; iflagpn = 1
  } else if (shape == "Increasing") {
    iflagsc = 1; iflagpn = 1
  } else if (shape == "Decreasing") {
    iflagsc = 1; iflagpn = -1
  } else {
    stop(paste("error input shape=", shape, 
               "\n choose shape=c(\"Free \", \"Increasing\")",
    ))
  }
  
  
  # ****************************************************************************
  # Bivariate basis is tensor product of univarite basis
  # Smoothing prior declines like exp(-(k1+k2)*gamma)
  # so coefficients theta_{j,k1,k2} shrink to zero very rapidly
  # iflagJoint provides 2 methods to eliminate large k1+k2
  iflagJoint  = iflagJoint   # Effects combinations of (j,k) in theta_{j,k}
  # 1 = all (j,k): Only if nbasis is smaller
  # 2 = Lower half of (j,k)
  # 3 = all (j,k) minus those within a circle
  
  # ****************************************************************************
  # Adaptive MCMC in pre burn-in period
  maxmodmet = 0
  nblow0    = 0
  # ****************************************************************************
  nmcmcall  = maxmodmet*nblow0 + nmcmc
  
  # ****************************************************************************
  # Get data
  # ****************************************************************************
  yname           = colnames(ydata)
  if(is.null(yname)){
    yname           = "Y"
    ydata           = matrix(ydata)
    colnames(ydata) = yname
  } 
  ntot            = nrow(ydata)         # Total number of observations
  nparx           = ncol(xdata)         # Should be 2
  xnames          = colnames(xdata)
  if(is.null(xnames)){
    xnames = c("X1","X2")
    colnames(xdata) = xnames
  } 
  # Find groups
  ugroups    = sort(unique(id_group))  # Unique group numbers
  gnames     = paste("G",ugroups,sep="_")
  ngroups    = length(ugroups)
  # id_group may not be sequential from 1 ... ngroups
  # Need an index to get correct betam
  id_beta  = matrix(0,ntot,1)
  for(j in 1:ngroups){
    id = which(id_group == ugroups[j])
    id_beta[id] = matrix(j,length(id),1)
  }
  #
  # vdata, wdata, and zdata
  # vdata is optional fixed effects
  if(is.null(vdata)){
    vnames = "NoFixedEffects"
  }else{
    # fixed effects
    nparv = ncol(vdata)
    vnames = colnames(vdata)
  }
  # zdata at minimum has constant 1
  if(is.null(wdata)){
    wdata            = matrix(1,ntot,1)
    wnames           = "CNST"
    colnames(wdata)  = wnames
    nparw            = 1
  }else{
    if (is.null(colnames(wdata))) colnames(wdata) <- paste0("w",1:ncol(wdata))
    oname = colnames(wdata)
    wdata = cbind(1, wdata)
    colnames(wdata) <- c("CNST", oname)
    nparw = ncol(wdata)
  }
  # zdata at minimum has constant 1
  if(is.null(zdata)){
    zdata            = matrix(1,ngroups,1)
    znames           = "CNST"
    colnames(zdata)  = znames
    nparz            = 1
  }else{
    if (is.null(colnames(zdata))) colnames(zdata) <- paste0("Z",1:ncol(zdata))
    oname = colnames(zdata)
    zdata = cbind(1, zdata)
    nparz = ncol(zdata)
  }
  
  ztz   = crossprod(zdata)
  
  # ****************************************************************************
  # Xgrid
  xrange     = xmax-xmin
  xgrid      = matrix(0,nint+1,nparx)
  xdelta     = xrange/nint
  for (j in 1:nparx) {
    xgrid[,j] = seq(xmin[j],xmax[j],xdelta[j])
  }
  
  # ****************************************************************************
  # Get indices for theta_{j,k,...}, vectorized xgrid
  kall         = 1:nbasis
  kall0        = matrix(c(0,kall))
  freq_all     = kall0
  ntheta1D     = nbasis + 1 
  for (j in 2:nparx) {
    f0           = matrix(1,ntheta1D,1) %x% freq_all 
    f1           = kall0 %x% matrix(1,nrow(freq_all),1)
    freq_all     = cbind(f0,f1)
  }  # End loop over dimensions 2 to nparx
  # Remove constant free f
  freq_all   = freq_all[-1,]     
  freq_sum   = apply(freq_all,1,sum)
  # ****************************************************************************
  # Remover higher-order interactions of basis functions
  # All (j,k) minus those inside a circle
  if (iflagJoint==3) {
    rad        = sqrt((nbasis-freq_all[,1])^2+(nbasis-freq_all[,2])^2)
    idin       = which(rad >= kmax )
    freq_sum   = freq_sum[idin]
    freq_all   = freq_all[idin,]
  }
  # Lower half of (j,k)
  if (iflagJoint==2) {
    idin       = which(freq_sum<=kmax)
    freq_sum   = freq_sum[idin]
    freq_all   = freq_all[idin,]
  }
  
  
  
  ntheta = nrow(freq_all)
  
  # ****************************************************************************
  theta_names = NULL
  for (kk in 1:ntheta) {
    j  = freq_all[kk,1]
    k  = freq_all[kk,2]
    theta_names  = c(theta_names,paste("theta(",j,",",k,")"))
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
  phixdata = matrix(0,ntot,ntheta)
  for (i in 1:ntot) {
    xi       = xdata[i,1]
    yi       = xdata[i,2]
    # Loop over indices where i+j <= kmax
    for(kk in 1:ntheta){
      k1      = freq_all[kk,1]   # Frequency for xi
      k2      = freq_all[kk,2]   # Frequency for yi
      if(k1==0){
        phixi  = phix0
      }else{
        phixi  = sqrt(2/xrange[1])*cos(k1*pi*(xi-xmin[1])/xrange[1])
      }
      if(k2==0){
        phiyi  = phiy0
      }else{
        phiyi  = sqrt(2/xrange[2])*cos(k2*pi*(yi-xmin[2])/xrange[2])
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
  phixgrid   = matrix(0,(nint+1)^2,ntheta)
  icount = 0
  for(j in 1:(nint+1)){
    yj = xgrid[j,2]
    for(i in 1:(nint+1)){
      icount = icount + 1
      xi    = xgrid[i,1]
      # Loop over indices where j+k <= kmax
      for(kk in 1:ntheta){
        k1      = freq_all[kk,1]   # Frequency for xi
        k2      = freq_all[kk,2]   # Frequency for yj
        if(k1==0){
          phixi  = phix0
        }else{
          phixi  = sqrt(2/xrange[1])*cos(k1*pi*(xi-xmin[1])/xrange[1])
        }
        if(k2==0){
          phiyj  = phiy0
        }else{
          phiyj  = sqrt(2/xrange[2])*cos(k2*pi*(yj-xmin[2])/xrange[2])
        }
        phixgrid[icount,kk] = phixi*phiyj
      }  # end loop over indices
    }  # End loop over xi
  } # End lover over yj
  
  
  # ****************************************************************************
  outpt  = MakeHBPointer(id_group)
  nobs1  = outpt$nobs
  iptHB1 = outpt$iptd
  iptHB0 = vector("list", length=length(iptHB1))
  for (j in 1:ngroups) {
    iptHB0[[j]] = iptHB1[[j]] - 1
  }
  
  # ****************************************************************************
  # Compute vtv in different groups
  # Note: alpha is fixed effect that is common to all groups
  #       but error variance depends on groups.  
  vtv     = NULL
  if(nparv > 0){
    if(iflagHBsigma==0){
      vtv  = crossprod(vdata)  # Homogeneous variances
    }else{
      # HB variances
      vtv   = array(0,dim=c(nparv,nparv,ngroups))
      for(j in 1:ngroups){
        vj       = vdata[iptHB1[[j]],]
        vtv[,,j] = crossprod(vj)
      } # End loop over groups
    } # Homogeneous vs HB variances
  } # End compute vtv 
  
  # ****************************************************************************
  # Compute wtw in each group
  # beta_i is random effects
  wtw     = 0
  if(nparw > 0){
    wtw     = array(0,dim=c(nparw,nparw,ngroups))
    for(j in 1:ngroups){
      if(nparw > 1){
        wj       = wdata[iptHB1[[j]],]
      }else{
        wj       = wdata[iptHB1[[j]]]
      }
      wtw[,,j] = crossprod(wj)
    } # End loop over groups
  }  # End wtw
  
  # ****************************************************************************
  # Compute phix'phix for each group
  phi2     = array(0,dim=c(ntheta,ntheta,ngroups))
  for(j in 1:ngroups){
    phij   = phixdata[iptHB1[[j]],]
    phi2[,,j] = crossprod(phij)
  } # End loop over groups
  
  
  # ****************************************************************************
  # ****************************************************************************
  # Make lists of parameters for adaptive Metropolis
  sigma2_met = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=1) # Adaptive Metropolis for HB Sigma2 hyperparameters
  gamma_met  = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=1) # Adaptive Metropolis for HB gamma  hyperparameters
  
  
  # ****************************************************************************
  # Initialize priors
  # Fixed effect alpha ~ N(m0,V0)
  if (nparv>0) {
    alpha_m0    = matrix(0,nparv,1)
    alpha_v0    = 1000*diag(nparv)
    alpha_v0i   = diag(nparv)/1000
    alpha_v0im0 = alpha_v0i %*% alpha_m0
  } else {  # no fixed effects: use placeholders
    alpha_m0    = matrix(0, 1,1)
    alpha_v0    = matrix(1000, 1,1)
    alpha_v0i   = matrix(1/1000, 1,1)
    alpha_v0im0 = matrix(0, 1,1)
  }
  
  # Error variance: homogeneous or HB
  if (iflagHBsigma == 0) {
    # Sigma2 ~ IG(r0/2,s0/2)
    sigma2_r0   = 4
    sigma2_s0   = 2
    sigma       = matrix(1,ngroups,1)
    sigma2      = sigma^2
    sigma2_rn   = sigma2_r0 + ntot
    
    sigma2_alpha = sigma2_beta = sigma2_alpha_m0 = sigma2_alpha_v02 = sigma2_beta_r0 = sigma2_beta_s0 = 0
    
  } else {
    # HB model for variance: sigma_j^2 ~ IG(alpha/2,beta/2)
    # Hyperpriors:   alpha ~ N(m0,v02)I(alpha>2) and beta~G(r0/2,s0/2)
    sigma2_alpha     = 8
    sigma2_beta      = 2
    sigma2_mu        = sigma2_beta/(sigma2_alpha - 2)
    sigma2_alpha_m0  = 8
    sigma2_alpha_v02 = 25
    sigma2_beta_r0   = 8
    sigma2_beta_s0   = 2
    sigma2           = matrix(1,ngroups,1)
    sigma            = sqrt(sigma2)
    
    sigma2_rn = 0
    sigma2_r0 = 0
    sigma2_s0 = 0
  }
  
  # HB Regression beta = zdata*Phi + delta
  phidim    = nparw*nparz
  phi_b0    = matrix(0,phidim,1)    # prior mean 
  phi_v0i   = diag(phidim)/10000    # prior precision
  phi_v0ib0 = matrix(0,phidim,1)
  
  # Var(delta) = Lambda ~ IW
  lambda_f0  = 5
  lambda_fn  = lambda_f0 + ngroups
  lambda_g0  = 10*diag(nparw)
  lambda_g0i = diag(nparw)/10
  
  # theta_jk ~ N(theta_0k, tau_k^2*exp(-k*gamma))
  # tau_k^2  ~ IG(r0/2,s0/2)
  tau2_mu        = 1   # Note: Keep tau2_mu = 1
  tau2_sigma     = 10
  tau2_r0        = 2*((tau2_mu/tau2_sigma)^2 + 2)
  tau2_s0        = 2*tau2_mu*(tau2_r0/2 - 1)
  
  tau2_rn    = tau2_r0 + ngroups
  
  # ****************************************************************************
  w0   = 1
  wk   = sum(freq_sum)/2
  
  
  
  # ****************************************************************************
  # Smoothing parameter  gamma
  gammall         = matrix(1,ngroups,1)       # lower level smoothing
  
  # Initialize parameters for gamma_j ~ G(alpha,beta)I(gamma < gmax)
  # Reparametrized to achieve better mixing
  # mu = alpha/beta and sigma2 = alpha/beta2
  # Priors: beta ~ Gamma(r0,s0) and alpha ~ N(alpha_m0,alpha_v02)I(0<alpha)
  gmax             = 3       # Maximum value for gamma
  gamma_mu         = 1
  gamma_sigma      = 1
  gamma_sigma2     = gamma_sigma^2
  gamma_beta       = gamma_mu/gamma_sigma2
  gamma_alpha      = gamma_mu*gamma_beta
  gamma_prob       = pgamma(gmax,shape=gamma_alpha,rate=gamma_beta)
  # hyper-parameters for gamma_mu and gamma_sigma
  # gamma_mu ~ Truncated N(m0,v02)
  gamma_mu_m0      = 1
  gamma_mu_v02     = 4
  mu_bot           = .1
  mu_top           = gmax
  
  # gamma_sigma ~ Truncated N(m0,v02)
  gamma_sigma_m0   = 1
  gamma_sigma_v02  = 16
  sigma_bot        = .5
  sigma_top        = 10
  
  # ****************************************************************************
  tau2all         = matrix(1,ntheta,1)   # upper level variances for theta_jk
  tauall          = sqrt(tau2all)
  
  # ****************************************************************************
  # Upper-level smoothing parameters
  # theta_0k ~ N(0,eta02*(1+k/gamma_beta)^{-gamma_alpha})
  eta0            = 20           # Variance parameter
  eta02           = eta0^2
  # Estimate eta02? eta02 ~ IG(r0/2,s0/2)
  eta02_m0        = 1
  eta02_r0        = 5
  eta02_s0        = eta02_m0*(eta02_r0 - 2)
  eta02_rn        = eta02_r0 + ntheta
  
  
  
  # gamma0   ~ Gamma(alpha,beta)
  gamma0          = 1                                       
  gam0vec         = (1+freq_sum/gamma_beta)^(-gamma_alpha)
  th0v            = eta02*gam0vec
  
  
  # ****************************************************************************
  # Initialize parameters
  if (nparv>0) {
    alpha = matrix(0,nparv,1)
  } else {
    alpha = matrix(0,1,1)
  }
  betall          = matrix(0,ngroups,nparw)
  phi             = matrix(0,nparz,nparw)
  lambda          = diag(1,nparw,nparw)
  lambdai         = diag(1,nparw,nparw)
  theta0          = matrix(0,ntheta,1)
  thetall         = matrix(0,ntheta,ngroups)    # lower level spectral coefficients
  
  if (nparv>0) {
    valpha = vdata%*%alpha
  } else {
    valpha = matrix(0,ntot,1)
  }
  
  wbeta   = matrix(0,ntot,1)
  for (j in 1:ngroups) {   # Compute w*beta
    wj    = as.matrix(wdata[iptHB1[[j]],])
    bj    = as.matrix((betall[j,]))
    wbeta[iptHB1[[j]],1] = wj %*% bj
  }
  
  
  
  # ****************************************************************************
  # Matrics for saving MCMC iterations
  if(nparv>0){
    alphag  = matrix(0,smcmc,nparv)
  }else{
    alphag  = 0
  }
  # Error variance: homogeneous or HB
  if(iflagHBsigma == 0){
    sigmag    = matrix(0,smcmc,1)
  }else{
    sigmag    = matrix(0,smcmc,ngroups)
    sigma2_alphag = matrix(0,smcmc,1)
    sigma2_betag  = matrix(0,smcmc,1)
  }
  tauallg       = matrix(0,smcmc,ntheta)      # StdDev for spectral coefficients
  eta0g         = matrix(0,smcmc,1)           # STDEV for upper level spectral coefficients
  gamma_mugg    = matrix(0,smcmc,1)           # Smoothing parameter for upper level model
  gammallg      = matrix(0,smcmc,ngroups)      # Smoothing parameters for lower level model
  thetam        = matrix(0,ntheta,ngroups)  
  thetas        = thetam
  thetallg      = array(0,dim=c(ntheta,ngroups,smcmc))
  theta0g       = matrix(0,smcmc,ntheta)
  fxgridm       = array(0,dim=c(nint+1,nint+1,ngroups))
  fxgrids       = array(0,dim=c(nint+1,nint+1,ngroups))
  f0xgridm      = matrix(0,nint+1,nint+1)
  f0xgrids      = matrix(0,nint+1,nint+1)
  betallg       = array(0, c(ngroups, nparw, smcmc))
  phig          = matrix(0,smcmc,phidim)
  nparw2        = nparw*nparw
  lambdag       = matrix(0,smcmc,nparw2)
  betam         = betall
  betas         = betall
  
  gamma_alphag  = matrix(0,smcmc,1)
  gamma_betag   = matrix(0,smcmc,1)
  
  
  fxgridall       = array(0,dim=c(nint+1,nint+1,ngroups))
  fxobs           = matrix(0,ntot)
  
  f0xgrid         = matrix(0,nint+1,nint+1)
  fxobsm          = matrix(0,ntot,1)  # Compute average over MCMC of fj(x) at observations
  
  # ****************************************************************************
  # prepare cpp----
  ## convert
  data_all = list(ydata=ydata, valpha=valpha, wbeta=wbeta, vdata=vdata, wdata=wdata, zdata=zdata)
  cpara = c(ntot, ngroups, nparv, nparw, nparw2, nparz, phidim, nbasis, nblow, nblow0, maxmodmet,
            nskip, nmcmcall, smcmc, nint, ntheta, iflagHBsigma)
  
  met_para = list(gamma_met=gamma_met, sigma2_met=sigma2_met)
  
  gamma_para = list(gmax=gmax, gamma_mu=gamma_mu, gamma_sigma=gamma_sigma, gamma_sigma2=gamma_sigma2,
                    gamma_beta=gamma_beta, gamma_alpha=gamma_alpha, gamma_prob=gamma_prob,
                    gamma_mu_m0=gamma_mu_m0, gamma_mu_v02=gamma_mu_v02, mu_bot=mu_bot, mu_top=mu_top,
                    gamma_sigma_m0=gamma_sigma_m0, gamma_sigma_v02=gamma_sigma_v02, sigma_bot=sigma_bot,
                    sigma_top=sigma_top, gam0vec=gam0vec, wk=wk)
  
  tau_para = list(tau2_s0=tau2_s0, tau2_rn=tau2_rn)
  eta_para = list(eta02_s0=eta02_s0, eta02_rn=eta02_rn, eta02=eta02, eta0=eta0, th0v=th0v)
  
  sigma_para = c(sigma2_s0, sigma2_rn, sigma2_alpha, sigma2_beta, sigma2_alpha_m0,
                 sigma2_alpha_v02, sigma2_beta_r0, sigma2_beta_s0)
  
  v_para = list(alpha_v0i=alpha_v0i, alpha_v0im0=alpha_v0im0)
  z_para = list(phi_v0i=phi_v0i, lambda_fn=lambda_fn, lambda_g0i=lambda_g0i)
  
  basis_para = list(kall=kall,freq_sum=freq_sum)
  phi_para = list(phixgrid=phixgrid, phixdata=phixdata)
  
  fx_para = list(fxgridall=fxgridall, f0xgrid=f0xgrid, fxdata=fxobs)
  
  # ****************************************************************************
  # MCMC Loop ----
  start_time = Sys.time()
  mcmc_out = get_mcmc_2DHBSAR(data_all=data_all,
                              nobs=nobs1,
                              iptHB1=iptHB0,
                              ztz=ztz,
                              vtv=vtv,
                              wtw=wtw,
                              phi2=phi2,
                              
                              cpara=cpara,
                              gamma_para=gamma_para,
                              tau_para=tau_para,
                              eta_para=eta_para,
                              met_para=met_para,
                              sigma_para=sigma_para,
                              v_para=v_para,
                              z_para=z_para,
                              phi_para=phi_para,
                              fx_para=fx_para,
                              basis_para=basis_para,
                              
                              alpha=alpha,
                              betall=betall,
                              phi=phi,
                              lambda=lambda,
                              lambdai=lambdai,
                              thetall=thetall,
                              theta0=theta0,
                              gammall=gammall,
                              tau2all=tau2all,
                              tauall=tauall,
                              sigma2=sigma2,
                              sigma=sigma,
                              
                              betallg=betallg,
                              thetallg=thetallg)
  end_time   = Sys.time()
  run_time   = end_time - start_time
  
  alphag        = mcmc_out$mcmcg$alphag
  sigmag        = mcmc_out$mcmcg$sigmag
  sigma2_alphag = mcmc_out$mcmcg$sigma2_alphag
  sigma2_betag  = mcmc_out$mcmcg$sigma2_betag
  tauallg       = mcmc_out$mcmcg$tauallg
  eta0g         = mcmc_out$mcmcg$eta0g
  gamma_mugg    = mcmc_out$mcmcg$gamma_mugg
  gammallg      = mcmc_out$mcmcg$gammallg
  theta0g       = mcmc_out$mcmcg$theta0g
  phig          = mcmc_out$mcmcg$phig
  lambdag       = mcmc_out$mcmcg$lambdag
  
  
  fxgridm    = mcmc_out$fxgridm
  fxgrids    = mcmc_out$fxgrids
  f0xgridm   = mcmc_out$f0xgridm
  f0xgrids   = mcmc_out$f0xgrids
  fxdatam    = mcmc_out$fxdatam
  fxdatas    = mcmc_out$fxdatas
  
  betam         = mcmc_out$betam
  betas         = mcmc_out$betas
  thetam        = mcmc_out$thetam
  thetas        = mcmc_out$thetas
  
  met.out = list()
  met.out$gamma_met  = get_mat_to_list_met(mcmc_out$metg$gamma_met)
  met.out$sigma2_met = get_mat_to_list_met(mcmc_out$metg$sigma2_met)
  
  # Compute summary statistics
  met.out$sigma2_met$pmet = met.out$sigma2_met$pmet/smcmc
  met.out$gamma_met$pmet  = met.out$gamma_met$pmet/smcmc
  
  fxgridm     = fxgridm/smcmc
  fxgrids     = sqrt(abs(fxgrids - smcmc*fxgridm^2)/smcmc)
  fxgridmm    = apply(fxgridm,1:2,mean)   # Mean of f over populations
  f0xgridm    = f0xgridm/smcmc
  f0xgrids    = sqrt(abs(f0xgrids-smcmc*f0xgridm^2)/smcmc)
  fxobsm      = fxobsm/smcmc
  theta0m     = apply(theta0g,2,mean)
  theta0s     = apply(theta0g,2,sd)
  thetam      = thetam/smcmc
  thetas      = sqrt(abs(thetas-smcmc*thetam^2)/smcmc)
  colnames(thetam) = gnames
  colnames(thetas) = gnames
  rownames(thetam) = theta_names
  rownames(thetas) = theta_names
  eta0m       = apply(eta0g,2,mean)
  eta0s       = apply(eta0g,2,sd)
  
  
  if(nparv > 0){
    alpham    = matrix(apply(alphag,2,mean))
    alphas    = matrix(apply(alphag,2,sd))
    alpha_pg0 = matrix(apply(alphag>0,2,mean))
    valpham   = vdata %*% alpham
  }else{
    alpham    = 0
    alphas    = 1
    valpham   = matrix(0,ntot,1)
  }
  if(iflagHBsigma == 0){
    sigmam   = mean(sigmag)
    sigmas   = sd(sigmag)
  }else{
    sigmam      = apply(sigmag,2,mean)
    sigmas      = apply(sigmag,2,sd)
    sigma2_alpham = mean(sigma2_alphag)
    sigma2_alphas = sd(sigma2_alphag)
    sigma2_betam  = mean(sigma2_betag)
    sigma2_betas  = sd(sigma2_betag)
    sigmm         = mean(sigmam)
    sigms         = sd(sigmam)
    sigsm         = mean(sigmas)
    sigss         = sd(sigmas)
    sigmmin       = min(sigmam)
    sigmmax       = max(sigmam)
  }
  
  tauallm     = apply(tauallg,2,mean)
  taualls     = apply(tauallg,2,sd)
  
  gammallm      = apply(gammallg,2,mean) 
  gammalls      = apply(gammallg,2,sd) 
  gamma_mugm    = mean(gamma_mugg)
  gamma_mugs    = sd(gamma_mugg)
  gamma_alpham  = mean(gamma_alphag)
  gamma_alphas  = sd(gamma_alphag)
  gamma_betam   = mean(gamma_betag)
  gamma_betas   = sd(gamma_betag)
  
  phim       = matrix(apply(phig,2,mean),nparz,nparw)
  phis       = matrix(apply(phig,2,sd),  nparz,nparw)
  phi_pg0    = matrix(apply(phig>0,2,mean),nparz,nparw)
  lambdam    = matrix(apply(lambdag,2,mean),nparw,nparw)
  lambdas    = matrix(apply(lambdag,2,sd),  nparw,nparw)
  betam      = betam/smcmc
  betas      = sqrt(abs(betas - smcmc*betam^2)/smcmc)
  
  
  wbetam     = apply(wdata*betam[id_beta,],1,sum)
  yresid     = ydata - valpham - wbetam
  
  # ****************************************************************************
  # Plot individual-level f_0
  x1mat = xgrid[,1] %*% matrix(1,1,nint+1)
  x2mat = matrix(1,nint+1,1) %*% t(xgrid[,2])
  # ****************************************************************************
  
  mcmc.out = list()
  if (nparv>0) mcmc.out$alphag = alphag
  mcmc.out$sigmag = sigmag
  if (iflagHBsigma == 1) {
    mcmc.out$sigma2_alphag = sigma2_alphag
    mcmc.out$sigma2_betag  = sigma2_betag
  }
  mcmc.out$tauallg      = tauallg
  mcmc.out$gammallg     = gammallg
  mcmc.out$gamma_mugg   = gamma_mugg
  mcmc.out$gamma_alphag = gamma_alphag
  mcmc.out$gamma_betag  = gamma_betag
  mcmc.out$eta0g        = eta0g
  mcmc.out$theta0g      = theta0g
  mcmc.out$thetallg     = thetallg
  mcmc.out$lambdag      = lambdag
  mcmc.out$phig         = phig
  
  post.est = list()
  post.est$fxgridm  = fxgridm
  post.est$fxgrids  = fxgrids
  post.est$f0xgridm = f0xgridm
  post.est$f0xgrids = f0xgrids
  post.est$fxdatam  = fxobsm
  post.est$theta0m  = theta0m
  post.est$theta0s  = theta0s
  post.est$thetam   = thetam
  post.est$thetas   = thetas
  post.est$eta0m    = eta0m
  post.est$eta0s    = eta0s
  
  if(nparv > 0){
    post.est$alpham    = alpham
    post.est$alphas    = alphas
    post.est$alpha_pg0 = alpha_pg0
    post.est$valpham   = valpham
  } else {
    post.est$valpham   = matrix(0,ntot,1)
  }
  post.est$wbetam   = wbetam
  post.est$sigmam   = sigmam
  post.est$sigmas   = sigmas
  if(iflagHBsigma == 1) {
    post.est$sigma2_alpham = sigma2_alpham
    post.est$sigma2_alphas = sigma2_alphas
    post.est$sigma2_betam  = sigma2_betam
    post.est$sigma2_betas  = sigma2_betas
    post.est$sigmm         = sigmm
    post.est$sigms         = sigms
    post.est$sigsm         = sigsm
    post.est$sigss         = sigss
    post.est$sigmmin       = sigmmin
    post.est$sigmmax       = sigmmax
  }
  
  post.est$tauallm = tauallm
  post.est$taualls = taualls
  
  post.est$gammallm     = gammallm
  post.est$gammalls     = gammalls
  post.est$gamma_mugm   = gamma_mugm
  post.est$gamma_mugs   = gamma_mugs
  post.est$gamma_alpham = gamma_alpham
  post.est$gamma_alphas = gamma_alphas
  post.est$gamma_betam  = gamma_betam
  post.est$gamma_betas  = gamma_betas
  
  post.est$phim    = phim
  post.est$phis    = phis
  post.est$phi_pg0 = phi_pg0
  post.est$lambdam = lambdam
  post.est$lambdas = lambdas
  post.est$betam   = betam
  post.est$betas   = betas
  
  out=list()
  out$met.out=met.out
  out$mcmc.out=mcmc.out
  out$post.est=post.est
  
  out$valpham = valpham
  out$wbetam  = wbetam
  out$fxdatam = fxobsm
  
  out$y           = ydata
  out$v           = vdata
  out$x           = xdata 
  out$w           = wdata
  out$z           = zdata
  out$x1mat       = x1mat
  out$x2mat       = x2mat
  out$shape       = shape
  out$nbasis      = nbasis
  out$xgrid       = xgrid
  out$group       = id_group
  out$ngroup      = ngroups
  out$iflagCenter = iflagCenter
  out$yresid      = yresid
  out$mcmc.time   = run_time
  
  return(out)
}