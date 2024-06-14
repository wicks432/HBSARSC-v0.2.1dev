library(splines)
# ******************************************************************************
# HBBS
#   HB B-Spline
#      Common knots at quantiles for all data
#   DATA INPUT
#      Create <mydata.Rdata> file 
#      Please search & set infile = "mydata" without extension
#      <mydata.Rdata> has the following data structures:
#      ydata    = vector of dependent variable
#      xdata    = vector of independent variables for f
#      vdata    = NULL or matrix of variables for fixed effects alpha
#                 NULL if none. 
#                 Generally, does not include '1' for intercept since 
#                 random intercept makes more sense than common intercept,
#                 but user can decide
#      wdata    = matrix of variables for random effects beta_j
#                 Generally, includes column of '1' for intercept
#                 Column of '1' is in ONE of either vdata or wdata
#      zdata    = matrix of variables for upper level
#                 Must include column of '1' for intercept
#      id_groups = Vector that gives group number
#   Lower-level model
#      Observation i from group j
#         Y_ij       = v_ij*alpha + w_ij*beta_j + f_j(x_ij) + epsilon_ij
#         epsilon_ij ~ N(0,sigma^2) or N(0,sigma_j^2)
#         alpha is fixed effect, beta_j is random effect
#   Upper-level 
#      HB Model for beta_j
#         beta_j     = Phi'z_j + delta_j  were delta_j ~ N(0,Lambda)
#      HB Error variance 
#         sigma_j^2 ~ IG(sigma2_alpha/2,sigma2_beta/2)
#
#      f_j gets complex depending on shape restrictions
#      Unconstrained:
#         f_j(x) = sum_{k=1}^K theta_{j,k} phi_k(x)
#      Shape constraints: 
#         Derivatives of f are functions of h(x) and g(z)
#         g(Z) = Z^2    for square-root Gaussian process
#         g(Z) = exp(Z) for log Gaussian process
#         h_j(x) = 1 or softmax of square wave for U & S shapes
#      HB Smoothing Prior for Spectral Coefficients
#         iflagZcnst = 0 if free f, and iflagZcnst = 1 for shape constraints
#         All models:  
#            theta_j,k ~ N(theta_0,k,tau_k^2)       for k >  iflagZcnst
#         Distribution of constant term (with shape constraint) depends on g
#         log(theta_{j,0}) ~ N(log(theta_{0,0}),tau_0^2) for square-root Gaussian
#         theta_{j,0}      ~ N(theta_{0,0},tau_0^2)      for log Guassian
#
#      Upper Level Model for theta_{0,k} for k > iflagZcnst
#          theta_{0,k} ~ N(0,eta^2)
#
#      S-shaped & Multi-modal functions use softmax version of
#         of square wave with zeros at inflection point.  
#         Squish function for one inflection point
#      h_j(x) = [1-exp(psi_j*(x-omega_j))]/[1+exp(psi_j*(x-omega_j))]
#      HB model
#         psi_j or nu_j ~ log normal
#         omega_j = xmin + xrange*(1/(1+exp(-rho_j)))  for 1 inflection point
#         omega =xmin + xrange*(Ordered logit for more than one omega)
#         rho_j ~ normal
# ******************************************************************************
# ******************************************************************************
"HBBS" <- function(ydata,xdata,vdata=NULL,wdata=NULL,zdata=NULL,id_group=NULL,
                   mcmc=list(), nknots=20, nint=100, xmin=NULL, xmax=NULL,
                   iflagCenter  = 1, # Used with shape constraints: 1 if int f = 0 and 0 if f(xmin) = 0
                   iflagZ       = 0, # 0 = Model using Z^2 or 1 = exp(Z)
                   iflagHBsigma = 1, # 1 for HB sigma and 0 else
                   iflagACE     = 0, # 1 for Auto Correlate Error model
                   iflaglm      = 0, # 0 = Normal, 1=Binary Probit, 2=Ordinal DV
                   iflagpsi     = 1, # 1 = estimate slope of squish function for S and multi-modal)
                   nExtreme     = 2,
                   shape=c('Free','Increasing','Decreasing','IncreasingConvex','DecreasingConcave',
                           'IncreasingConcave','DecreasingConvex','IncreasingS','DecreasingS',
                           'IncreasingRotatedS','DecreasingRotatedS','InvertedU','Ushape',
                           'IncMultiModal','DecMultiModal')) {
  
  if (!is.matrix(ydata)) ydata = as.matrix(ydata)
  if (!is.matrix(xdata)) xdata = as.matrix(xdata)
  
  # MCMC Iterations
  mcvals=list(nblow=20000,smcmc=10000,nskip=1)
  mcvals[names(mcmc)]=mcmc
  nblow     = mcvals$nblow        # number of initial MCMC parameters
  smcmc     = mcvals$smcmc        # number of MCMC parameters to save
  nskip     = mcvals$nskip        # Save very nskip iterations
  nmcmc     = nblow + nskip*smcmc # total number of MCMC parameters
  nknots    = nknots              # Number of knots
  # ******************************************************************************
  # Mesh for computing and plotting f
  xmin      = ifelse(is.null(xmin), min(xmin), xmin)     # Minimum of xgrid
  xmax      = ifelse(is.null(xmax), min(xmax), xmax)     # Maximum of xgrid
  xrange    = xmax-xmin
  xmid      = (xmin+xmax)/2
  nint      = nint  # mesh size for computing & plotting f
  xdelta    = (xmax-xmin)/nint
  xgrid     = matrix(seq(xmin,xmax,xdelta))
  # ****************************************************************************
  # Set Parameters
  # ifalgsc gives free or shaped constraint
  # 0 is free f
  # 1 is monotone
  # 2 is increasing and convex or decreasing and concave
  # 3 is increasing and concave or decreasing and convex
  # 4 is U+Concave or inverse U+Convex:  
  #   If iflagpn = 1,  then U and convex
  #   if iflagpn = -1, then Inverted U and concave
  # 5 is S increasing convex-to-concave or decreasing concave-to-convex
  # 6 is S increasing concave-to-convex or decreasing convex-to-concave
  # 7 is multi-modal
  #   iflagpn = 1 =>  U, inverted U, U, ... or min, max, min, ...
  #   iflagpn = -1 => inverted U, U, ...   or max, min, max, min
  # ****************************************************************************
  if (shape == "Free") {
    iflagsc = 0; iflagpn = 1
  } else if (shape == "Increasing") {
    iflagsc = 1; iflagpn = 1
  } else if (shape == "Decreasing") {
    iflagsc = 1; iflagpn = -1
  } else if (shape == "IncreasingConvex") {
    iflagsc = 2; iflagpn = 1
  } else if (shape == "DecreasingConcave") {
    iflagsc = 2; iflagpn = -1
  } else if (shape == "IncreasingConcave") {
    iflagsc = 3; iflagpn = 1
  } else if (shape == "DecreasingConvex") {
    iflagsc = 3; iflagpn = -1
  } else if (shape == "InvertedU") {
    iflagsc = 4; iflagpn = 1
  } else if (shape == "Ushape") {
    iflagsc = 4; iflagpn = -1
  } else if (shape == "IncreasingS") {        # Increasing S convex -to-concave
    iflagsc = 5; iflagpn = 1
  } else if (shape == "DecreasingS") {        # Decreasing S concave-to-convex
    iflagsc = 5; iflagpn = -1
  } else if (shape == "IncreasingRotatedS") { # Increasing S concave-to-convex
    iflagsc = 6; iflagpn = 1
  } else if (shape == "DecreasingRotatedS") { # Decreasing S convex -to-concave
    iflagsc = 6; iflagpn = -1
  } else if (shape == "IncMultiModal") {
    iflagsc = 7; iflagpn = 1
  } else if (shape == "DecMultiModal") {
    iflagsc = 7; iflagpn = -1
  } else {
    stop(paste("error input shape=", shape, 
               "\n choose shape=c(\"Free \", \"Increasing\", \"Decreasing\", \"IncreasingConvex\", 
               \"DecreasingConvex\",\"IncreasingS\", \"DecreasingS\", \"IncreasingRotatedS\", 
               \"DecreasingRotatedS\",  \"InvertedU\", \"Ushape\", \"IncMultiModal\", \"DecMultiModal\")",
    ))
  }
  
  iflagZcnst=0
  if(iflagsc>0) iflagZcnst = 1  # Include constant 
  if((iflagsc==5)|(iflagsc==6)) iflagZ = 0  # Identification issue with squish & Log-normal
  if(iflagsc < 5) iflagpsi   = 0    # Doesn't use psi
  if(iflaglm == 1) iflagHBsigma = 0  # Probit sets sigma = 1 for identification
  # ****************************************************************************
  # Give # of stationary and inflection points
  # U-shaped functions have one stationary point
  nstat = 0; nflex = 0; nomega = 0; id_stat=0; id_flex=0
  # If S or multi-modal, give inflection/stationary points
  id_flex      = 0  # Indices for inflection points in combined set
  id_stat      = 0  # Indices for stationary points in combined set
  omega0_input = 0
  zeta0_input  = 0
  if (iflagsc == 4) {
    nstat  = 1
    nflex  = 0
    nomega = 1
  }
  
  # S-shaped functions have one inflection point & one omega
  if((iflagsc == 5)|(iflagsc == 6)){
    nstat   = 0 
    nflex   = 1
    nomega  = 1
    id_flex = 1
    omega0_input = xmid  # Initialize to mid point
    # Get zeta0_input from omega0_input
    z           = (omega0_input-xmin)/xrange
    zeta0_input = log(z/(1-z))  # logit of omega0_input
  }
  # Multi-modal functions: number of stationary & inflection points are given
  if(iflagsc == 7){
    nstat   = nExtreme      # number of stationary points
    nflex   = nstat-1       # number of inflection points
    nomega  = nstat+nflex   # Total number of omegas
    id_stat = seq(1,by=2,length.out=nstat)  # Identifies stationary points from theta
    id_flex = seq(2,by=2,length.out=nflex)  # Identifies inflextion points from theta
    # Give initial values:  Equally spacing
    odelta        = xrange/(nomega+1)
    omega0_input  = seq(xmin+odelta,xmax-odelta,odelta)
    # omega0_input = matrix(c(3,5,6,7,8))
    # Find corresponding  zeta0_input for ordinal logit
    z           = (omega0_input-xmin)/xrange
    x           = z
    x[-1]       = z[-1] - z[-nomega]
    x0          = 1-sum(x)
    zeta0_input = log(x/x0)
  }
  # ****************************************************************************
  # Don't use exp with S shape because of weak identification
  if((iflagsc==5)|(iflagsc==6)) iflagZ = 0
  # ****************************************************************************
  # ****************************************************************************
  
  ntot    = nrow(ydata)   # Total number of observations
  ugroups = sort(unique(id_group))  # Unique group numbers
  ngroups = length(ugroups)   # Total number of groups
  
  if(is.null(vdata)){
    nparv  = 0
    vnames = "NoVdata"
    vdata  = matrix(0, ntot, 1)
    vtv    = 0
  }else{
    nparv   = ncol(vdata)  # Data matrix for fixed effects
    vnames  = colnames(vdata)
  }
  
  if (is.null(wdata)) {
    wdata = matrix(1, ntot)
    colnames(wdata) <- "CNST"
  } else {
    if (is.null(colnames(wdata))) colnames(wdata) <- paste0("W",1:ncol(wdata))
    oname = colnames(wdata)
    wdata = cbind(1, wdata)
    colnames(wdata) <- c("CNST", oname)
  }
  
  if (is.null(zdata)) {
    zdata = matrix(1, ngroups)
    colnames(zdata) <- "CNST"
  } else {
    if (is.null(colnames(zdata))) colnames(zdata) <- paste0("Z",1:ncol(zdata))
    oname = colnames(zdata)
    zdata = cbind(1, zdata)
    colnames(zdata) <- c("CNST", oname)
  }
  
  # Data Parameters
  nparw      = ncol(wdata)  # Data matrix for random effects
  nparw2     = nparw*nparw
  nparz      = ncol(zdata)  # Data matrix for upper-level model
  yname      = colnames(ydata)
  xname      = colnames(xdata)
  vnames     = colnames(vdata)
  wnames     = colnames(wdata)
  znames     = colnames(zdata)
  gnames     = paste("G",ugroups,sep="_")
  ztz        = crossprod(zdata)
  
  # ******************************************************************************
  # Locate indices that define data (x1,x2)
  #   Used with shape constraints.
  #   Instead of computing f(x_i) at observations,
  #   find where x_i falls in xgrid and use fxgrid
  #   Find location of xdata in xgrid
  #   xgrid[xinx[i,1]] <= xi <= xgrid[xinx[i,2]]
  #
  # ******************************************************************************
  xinx  =  matrix(0,ntot,2)
  xover  = matrix(0,ntot,1)   # Relative excess over grid
  # Loop over observations
  for (i in 1:ntot) {
    xi  = xdata[i]
    idx = which(xi == xgrid)  # Is it exactly on xgrid?
    if (length(idx) > 0 ) {
      xinx[i,1]  = idx
      xinx[i,2]  = idx
      xover[i]   = 0 
    } else {
      idx       = max(which(xi > xgrid))
      xinx[i,1] = idx
      xinx[i,2] = idx+1
      xover[i]  = (xi - xgrid[idx])/(xdelta)
    }
  } # End loop over observations for xinx and xover
  
  # ******************************************************************************
  # Adaptive MCMC in pre burn-in period
  maxmodmet = 0
  nblow0    = 0
  if (iflagsc>0) {
    # shape constraints: use adaptive Metropolis
    maxmodmet = 10           # Maximum times to modify Metropolis
    nblow0    = 1000         # Initial number of initial MCMC parmeters
  }
  # ******************************************************************************
  
  nmcmcall  = maxmodmet*nblow0 + nmcmc
  
  
  yodata = ydata
  # 0/1 data for probit model
  if (iflaglm == 1) { # Probit Model
    yodata = ydata           # yodata are the observed 0/1 variables
    # Initialize the latent y, which is N(mean,1)
    ydata = -1 + 2*yodata   # ydata are the latent variables
    id_probit1 = which(yodata==1)  # index for ones
    id_probit0 = which(yodata==0)  # index for zeros
    ydatam     = matrix(0,ntot,1)
  } else {
    id_probit1=id_probit0=0
  }
  
  # Ordinal data
  if (iflaglm == 2) { 
    # ydata becomes latent variables and yodata is observed data
    yodata = ydata  # yodata is observed data, ydata will be updated
    # cutpoints
    ymin = min(ydata)  # should be 1
    if(ymin != 1) cat("\n DATA ERROR Min(Y) is not 1\n")
    ymax = max(ydata)  
    
    cutt = seq(from = ymin+.5, to = ymax-.5, by = 1)
    cutt = c(ymin-10,cutt,ymax+10)  # Add really big lower and upper bounds for convenience
    ncut = length(cutt)
    id_cut_free = 3:(ncut-2)    # Rows of free cutpoints
    # For each observation, get cutpoints that bound it below and above
    # cutt[id_cut_bot] < ydata < cutt[id_cut_top]
    id_cut_bot  = ydata
    id_cut_top  = ydata + 1
    ybot        = cutt[id_cut_bot]          # Lower bounds for data
    ytop        = cutt[id_cut_top]          # Upper bounds for data
    cuttall     = matrix(cutt,ncut,ngroups) # Cut points for each group
    cuttallm    = matrix(0,ncut,ngroups)    # Estimated mean
    cuttalls    = cuttallm                  # estimated sd
    ydatam      = matrix(0,ntot,1)
  } else {
    cutt = ncut = id_cut_free = id_cut_bot = id_cut_top = ybot = ytop = 0
    cuttall = cuttallm = cuttalls = matrix(0, 1,1)
  }
  
  
  
  # ****************************************************************************
  outpt    = MakeHBPointer(id_group)
  nobs     = outpt$nobs
  iptHB    = outpt$iptd
  id_first = outpt$id_first  # id for first observation in each group
  id_beta  = outpt$id_beta
  iptHB1   = vector("list", length=length(iptHB))
  for (j in 1:ngroups) {
    iptHB1[[j]] = iptHB[[j]] - 1
  }
  
  # Compute vtv in each group
  vtv     = 0
  if(nparv > 0){
    vtv   = array(0,dim=c(nparv,nparv,ngroups))
    for (j in 1:ngroups) {
      if (nparv > 1) {
        vj  = vdata[iptHB[[j]],]
      } else {
        vj  = vdata[iptHB[[j]]]
      }
      vtv[,,j] = crossprod(vj)
      
    }  # End loop over groups
  } else {  # End vtv
    vtv = array(0,dim=c(1,1,ngroups))
  }
  # Compute wtw in each group
  wtw   = 0
  if(nparw > 0) {
    wtw = array(0,dim=c(nparw,nparw,ngroups))
    for (j in 1:ngroups) {
      if (nparw > 1) {
        wj = wdata[iptHB[[j]],]
      } else {
        wj = wdata[iptHB[[j]]]
      }
      wtw[,,j] = crossprod(wj)
    } # End loop over groups
  }  # End wtw
  
  
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
  for(k in 1:nknots){
    idk       = which.max(phixgrid[,k])
    xknots[k] = xgrid[idk]
  }
  kall      = 1:ntheta
  
  # Include constant function
  if(iflagZcnst == 1){
    phixgrid  = cbind(matrix(1,nint+1,1),phixgrid)
    kall      = c(0,kall)
    ntheta    = ntheta + 1
  }
  # Precompute phi'phi for each group
  phi2      = 0
  if(iflagsc==0){
    phi2    = array(0,dim=c(nbasis,nbasis,ngroups))
    for(j in 1:ngroups){
      phixj     = phixdata[iptHB[[j]],]
      phi2[,,j] = crossprod(phixj)
    } # End loop over groups
  }else{
    # Don't use phixdata if there are shape constraints
    phixdata = matrix(0,1,1)
    phi2 = array(0,dim=c(1,1,1))    
  }
  
  theta_names = paste("Theta_",kall,sep="")
  
  
  
  # ******************************************************************************
  # ******************************************************************************
  # Make lists of parameters for adaptive Metropolis
  theta_met = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=ntheta)  # Adaptive Metropolis parameter for generic theta
  theta_met = matrix(rep(theta_met, ngroups), ncol=ngroups) # Make ngroups copies
  # use theta_met[[j]]$m0 to get m0 for jth population. 
  mslope_log_met = GetMetVec(w=0.5, m0=0.01, alpha=4.0, ndim=1)
  mslope_log_met = matrix(rep(mslope_log_met, ngroups), ncol=ngroups)
  
  
  psi_met = GetMetVec(w=0.5, m0=0.1, alpha=4.0, ndim=1)  # Adaptive Metropolis parameter for generic psi
  psi_met = matrix(rep(psi_met, ngroups), ncol=ngroups) # Make ngroups copies
  zeta_met = GetMetVec(w=0.5, m0=0.01, alpha=4.0, ndim=ifelse(nomega==0,1,nomega))  # Adaptive Metropolis parameter for generic zeta
  zeta_met = matrix(rep(zeta_met, ngroups), ncol=ngroups) # Make ngroups copies
  
  sigma2_met = GetMetVec(w=0.5, m0=0.001, alpha=4.0, ndim=1) # Adaptive Metropolis for HB Sigma2 hyperparameters
  
  # ******************************************************************************
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
  
  # theta_jk ~ N(theta_0k, tau_k^2)
  # tau_k^2  ~ IG(r0/2,s0/2)
  tau2_mu    = 1   # Note: Keep tau2_mu = 1
  tau2_sigma = 10
  tau2_r0    = 2*((tau2_mu/tau2_sigma)^2 + 2)
  tau2_s0    = 2*tau2_mu*(tau2_r0/2 - 1)
  
  tau2_rn    = tau2_r0 + ngroups
  
  # ****************************************************************************
  # Initialize squish fuction
  hfunall   = NULL
  hfunm     = NULL
  hfunj_new = NULL
  
  # ******************************************************************************
  tau2all         = matrix(1,ntheta,1)   # upper level variances for theta_jk
  tauall          = sqrt(tau2all)
  
  # ****************************************************************************
  # Upper-level smoothing parameters
  # theta_0k ~ N(0,eta02)
  eta0            = 1           # Variance parameter
  eta02           = eta0^2
  # Estimate eta02? eta02 ~ IG(r0/2,s0/2)
  eta02_m0        = 1
  eta02_r0        = 5
  eta02_s0        = eta02_m0*(eta02_r0 - 2)
  eta02_rn        = eta02_r0 + ntheta
  
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
  xiparall        = 0
  xipar0          = 0
  # If Z^2, model theta_0j = exp(xi_j) and theta_00 = exp(xi_0)
  if (iflagZcnst == 1 & iflagZ==0) {
    xiparall      = matrix(0,ngroups,1)      # xipar_j0 = log(theta_j0)
    xipar0        = 0                        # xipar_00 = log(theta_00)
    thetall[1,]   = matrix(exp(xiparall),1,ngroups) 
    theta0[1]     = exp(xipar0)
  }
  
  if (nparv > 0) {
    valpha = vdata%*%alpha
  } else {
    valpha = matrix(0,ntot,1)
  }
  
  wbeta = matrix(0,ntot,1)
  for (j in 1:ngroups) {   # Compute w*beta
    wj    = as.matrix(wdata[iptHB[[j]],])
    bj    = as.matrix((betall[j,]))
    wbeta[iptHB[[j]],1] = wj %*% bj
  }
  
  # ****************************************************************************
  # Initialize rho for ACE model
  rhoall    = matrix(0,ngroups,1)
  rho_m0    = 0    # Prior mean
  rho_s0    = .5    # Prior sd
  rho_v0    = rho_s0^2
  rho_v0i   = 1/rho_v0
  rho_v0im0 = rho_v0i*rho_m0
  
  
  # Monotone convex/concave:  Add linear term
  mslope0              = 1
  mslopeall            = matrix(1,ngroups,1)
  mslope_new           = NULL
  if ((iflagsc==2)|(iflagsc==3)) {
    mslope0              = 1
    mslopeall            = matrix(1,ngroups,1)
    mslope0_log          = log(mslope0)
    mslopeall_log        =  matrix(mslope0_log,ngroups,1)
    mslope_log_m0        = 0
    mslope_log_v0        = 1000
    mslope_log_v0i       = 1/1000
    mslope_log_v0im0     = mslope_log_v0i*mslope_log_m0
    mslope_log_sigma2    = 1
    mslope_log_sigma     = 1   
    mslope_log_sigma2_r0 = 4
    mslope_log_simga2_rn = 4 + ngroups
    mslope_log_sigma2_s0 = 2
    mslope_logg          = matrix(0,smcmc,ngroups)
    mslope0_logg         = matrix(0,smcmc,1)
    mslopeg              = mslope_logg
    mslope0g             = mslope0_logg
    mslope_log_sigmag    = matrix(0,smcmc,1)
  } else {
    mslope0_log=mslopeall_log=mslope_log_v0i=mslope_log_v0im0=mslope_log_sigma2=mslope_log_simga2_rn=mslope_log_sigma2_s0=0
    
    mslope_logg          = matrix(0,smcmc,ngroups)
    mslope0_logg         = matrix(0,smcmc,1)
    mslopeg              = mslope_logg
    mslope0g             = mslope0_logg
    mslope_log_sigmag    = matrix(0,smcmc,1)
  }
  
  # Data matrices for omega
  idball       = matrix(0,1,ngroups)
  omegall      = matrix(0,1,ngroups)
  
  # U-shaped
  if (iflagsc == 4) {
    nomega        = 1
    zeta0         = 0    
    zeta_sigma2   = 1
    zeta_sigma    = 1
    zeta_m0       = 0
    zeta_v02      = 100
    zeta_v0       = sqrt(zeta_v02)
    zeta_r0       = 4
    zeta_s0       = 2
    zeta_rn       = zeta_r0 + ngroups
    zetall        = matrix(0,1,ngroups)
    sez           = matrix(apply(exp(zetall),2,cumsum),nomega,ngroups)
    omegall       = xmin + xrange*sez/(1+sez[nomega,])
    omegall       = matrix(omegall,1,ngroups)
    omega0        = mean(omegall)
  } else {
    zeta0 = zeta_sigma2 = zeta_sigma = zeta_m0 = zeta_v02 = zeta_v0 = zeta_r0 = zeta_s0 = zeta_rn = zeta_r0 = sez = omega0 = 0
    zetall = matrix(0,1,ngroups)
  }
  # S and Multi-modal
  if ((iflagsc > 4)) {
    # **************************************************************************
    # S-shaped & Multi-modal has has squish function 
    # For one inflection point: logistic function 
    #    h(x) = (1-exp(psi*(x-omega)))/(1+exp(psi*(x-omega))
    # For multiple inflection points:
    #    Softmax version of square wave based on logistic function
    # HB model for omega is ordered logit
    # omegaj = xmin + xrange*(sum_k<=j exp(zetak)/(1+sum_all exp(zetak))
    # zetaj       ~ N(zeta0,zetaj_sigma2)
    # zeta0     ~ N(zeta0_m0,zeta_v02)
    # zetaj_sigma2 ~ IG(zeta_r0/2,zeta_s0/2)
    zeta0         = zeta0_input    
    zeta_sigma2   = matrix(1,nomega,1)
    zeta_sigma    = sqrt(zeta_sigma2)
    zeta_m0       = zeta0_input
    zeta_v02      = matrix(16,nomega,1)
    zeta_v0       = sqrt(zeta_v02)
    zeta_r0       = 4
    zeta_s0       = 2
    zeta_rn       = zeta_r0 + ngroups
    zetall        = zeta0 %*% matrix(1,1,ngroups) 
    sez           = matrix(apply(exp(zetall),2,cumsum),nomega,ngroups)
    omegall       = xmin + xrange*sez/(1+sez[nomega,])
    omega0        = omega0_input
    
    # **************************************************************************
    # HB model for log(psi_j) (if iflagpsi = 1)
    # psij_log       ~ N(psi_log_mu,psi_log_sigma2)
    # psi_log_mu     ~ N(psi_log_m0,psi_log_v02)
    # psi_log_sigma2 ~ IG(psi_log_r0/2,psi_log_s0/2)
    psi_fixed      = 1/xdelta           # Use psi_fixed if iflagpsi = 0
    psi_fixed_log  = log(psi_fixed)
    psi0           = psi_fixed
    psi0_log       = psi_fixed_log
    psiall         = matrix(psi_fixed,ngroups,1)  # All values of psi in lower leve model
    psiall_log     = log(psiall)                 # Do HB & metropolis on log(psi)
    psi_log_mu     = psi_fixed_log
    psi_log_sigma2 = 1
    psi_log_sigma  = sqrt(psi_log_sigma2)
    psi_log_m0     = psi_fixed_log
    psi_log_v02    = 9
    psi_log_v0     = sqrt(psi_log_v02)
    psi_log_r0     = 5
    psi_log_s0     = 1
    psi_log_rn     = psi_log_r0 + ngroups
    
    # Compute softmax square wave at inflection points
    hfunall  = matrix(0,nint+1,ngroups)
    idball   = matrix(0,nflex+2,ngroups)
    for (j in 1:ngroups) {
      omegaj      = matrix(omegall[,j])
      hfunall[,j] = rcpp_GetSquish(omegaj[id_flex-1],psiall[j],xgrid)
      # multi-modal
      if (iflagsc >= 7) {
        # Find different intervals
        breakers = sort(c(xmin,omegaj[id_flex],xmax))
        breakers = floor(breakers/xdelta)*xdelta  # round to xgrid gaps
        idb      = match(breakers,xgrid)          # find location in xgrid
        idball[,j] = idb                          # Indidices for 
      } # end mult-modal
    } # end loop over groups
    
  } else {  # End models with stationary and/or inflection points
    psi_fixed      = 0
    psi_fixed_log  = 0
    psi0           = 0
    psi0_log       = 0
    psiall         = 0
    psiall_log     = 0
    psi_log_mu     = 0
    psi_log_sigma2 = 0
    psi_log_sigma  = 0
    psi_log_m0     = 0
    psi_log_v02    = 0
    psi_log_v0     = 0
    psi_log_r0     = 0
    psi_log_s0     = 0
    psi_log_rn     = 0
    hfunall        = matrix(0,nint+1,ngroups)
    idball         = matrix(0,nflex+2,ngroups)
  }
  
  # ****************************************************************************
  # ****************************************************************************
  # Matrics for saving MCMC iterations
  if (nparv>0) {
    alphag = matrix(0, smcmc, nparv)
  } else {
    alphag = matrix(0, smcmc, 1)
  }
  if (iflagHBsigma==1) {
    sigmag = matrix(0, smcmc, ngroups)
  } else {
    sigmag = matrix(0, smcmc, 1)
  }
  sigma2_alphag = matrix(0, smcmc, 1)
  sigma2_betag  = matrix(0, smcmc, 1)
  
  thetallg     = array(0, c(ntheta, ngroups, smcmc))
  theta0g      = array(0, c(smcmc, ntheta))
  tauallg      = array(0, c(smcmc, ntheta)) # StdDev for spectral coefficients
  eta0g        = array(0, c(smcmc, 1)) # STDEV for upper level spectral coefficients
  
  betallg = array(0, c(ngroups, nparw, smcmc))
  phig    = array(0, c(smcmc, phidim))
  lambdag = array(0, c(smcmc, nparw2))
  
  # ACE
  rhoallg = matrix(0, smcmc, ngroups)
  
  # squish
  zetag   = array(0, c(nomega, ngroups, smcmc))
  omegag  = array(0, c(nomega, ngroups, smcmc))
  zeta0g  = array(0, c(smcmc, nomega))
  omega0g = array(0, c(smcmc, nomega))
  psiallg = array(0, c(smcmc, ngroups))
  psi0g   = array(0, c(smcmc, 1))
  
  # slope
  mslope_logg       = matrix(0, smcmc, ngroups)
  mslope0_logg      = matrix(0, smcmc, 1)
  mslopeg           = matrix(0, smcmc, ngroups)
  mslope0g          = matrix(0, smcmc, 1)
  mslope_log_sigmag = matrix(0, smcmc, 1)
  
  fxgridallg = array(0, c(nint+1, ngroups, smcmc))
  f0xgridg   = array(0, c(smcmc, nint+1))
  
  fxgridall       = matrix(0,nint+1,ngroups)
  fxdata          = matrix(0,ntot)
  fxdatam         = matrix(0,ntot,1)  # Compute average over MCMC of fj(x) at observations
  # Compute f at theta
  for (j in 1:ngroups) {
    
    fpars  = list(iflagsc=iflagsc, 
                  iflagpn=iflagpn, 
                  iflagCenter=iflagCenter, 
                  iflagZ=iflagZ, 
                  theta=thetall[,j], 
                  phix=phixgrid, 
                  delta=xdelta, 
                  range=xrange,
                  xmin=xmin, 
                  xgrid=xgrid, 
                  nstat=nstat,
                  hfun=hfunall[,j],
                  omega=omegall[,j],
                  id_stat=id_stat-1,
                  idb=idball[,j]-1, 
                  mslope=mslopeall[j])
    
    fxj           = rcpp_Getfx(fpars)
    fxgridall[,j] = fxj 
    
    # **************************************************************************
    # Compute fj(x_obs)
    # fxobs = GetUpfxobs(xgrid,xinx,xover,fxgrid)
    a         = xinx[iptHB[[j]],]-1
    b         = matrix(xover[iptHB[[j]]])
    c         = matrix(fxgridall[,j])
    fxdata[iptHB[[j]]] = rcpp_GetSCfxobs(a,b,c)
  }
  # Estimate upper level function with mean of lower level functions
  f0xgrid = matrix(apply(fxgridall,1,mean))
  
  # ****************************************************************************
  # prepare cpp----
  ## convert
  data_all = list(ydata=ydata, yodata=yodata, valpha=valpha, wbeta=wbeta, vdata=vdata, wdata=wdata, zdata=zdata)
  cpara = c(ntot, ngroups, nparv, nparw, nparw2, nparz, phidim, nbasis, nblow, nblow0, maxmodmet,
            nskip, nmcmcall, smcmc, nint, ntheta, nstat, nflex, nomega,
            iflagHBsigma, iflagpsi, iflagCenter, iflagACE, iflagZ, iflaglm, iflagpn, iflagsc, iflagZcnst)
  
  probit_para  = list(id_probit1=id_probit1-1, id_probit0=id_probit0-1)
  ordinal_para = list(id_cut_free=id_cut_free-1, id_cut_bot=id_cut_bot-1, id_cut_top=id_cut_top-1,
                      ybot=ybot, ytop=ytop, cuttall=cuttall, cuttallm=cuttallm, cuttalls=cuttalls)
  
  ace_para = list(rhoall=rhoall, rho_v0i=rho_v0i, rho_v0im0=rho_v0im0)
  
  met_para = list(theta_met=theta_met, psi_met=psi_met, zeta_met=zeta_met,
                  sigma2_met=sigma2_met, mslope_log_met=mslope_log_met)
  
  
  tau_para = list(tau2_s0=tau2_s0, tau2_rn=tau2_rn)
  eta_para = list(eta02_s0=eta02_s0, eta02_rn=eta02_rn, eta02=eta02, eta0=eta0)
  
  squish_para1 = list(zeta_m0=zeta_m0, zeta_v02=zeta_v02, zeta_s0=zeta_s0, zeta_rn=zeta_rn, psi_log_m0=psi_log_m0,
                      psi_log_v02=psi_log_v02, psi_log_s0=psi_log_s0, psi_log_rn=psi_log_rn, psi_fixed=psi_fixed,
                      psi_fixed_log=psi_fixed_log, psi0=psi0, psi_log_mu=psi_log_mu, psi_log_sigma2=psi_log_sigma2,
                      psi_log_sigma=psi_log_sigma)
  
  squish_para2 = list(zeta_sigma2=zeta_sigma2, zeta_sigma=zeta_sigma, zeta0=zeta0, omega0=omega0, zetall=zetall, 
                      psiall=psiall, psiall_log=psiall_log, omegall=omegall, hfunall=hfunall, idball=idball-1, id_stat=id_stat-1,id_flex=id_flex-1)
  
  sigma_para = c(sigma2_s0, sigma2_rn, sigma2_alpha, sigma2_beta, sigma2_alpha_m0,
                 sigma2_alpha_v02, sigma2_beta_r0, sigma2_beta_s0)
  
  v_para = list(alpha_v0i=alpha_v0i, alpha_v0im0=alpha_v0im0)
  z_para = list(phi_v0i=phi_v0i, lambda_fn=lambda_fn, lambda_g0i=lambda_g0i)
  
  slope_para = list(mslopeall=mslopeall, mslope0_log=mslope0_log, mslopeall_log=mslopeall_log,
                    mslope_log_v0i=mslope_log_v0i, mslope_log_v0im0=mslope_log_v0im0, 
                    mslope_log_sigma2=mslope_log_sigma2, mslope_log_simga2_rn=mslope_log_simga2_rn,
                    mslope_log_sigma2_s0=mslope_log_sigma2_s0)
  basis_para = list(kall=kall)
  phi_para = list(phixgrid=phixgrid, phixdata=phixdata)
  
  fx_para = list(fxgridall=fxgridall, f0xgrid=f0xgrid, fxdata=fxdata)
  dpara = c(xdelta,xrange,xmin,xmax,xipar0)
  
  # ****************************************************************************
  # MCMC Loop ----
  start_time = Sys.time()
  mcmc_out = get_mcmc_HBBS(data_all=data_all, 
                           dpara=dpara,
                           xgrid=xgrid,
                           xinx=xinx-1, 
                           xover=xover, 
                           nobs=nobs, 
                           iptHB1=iptHB1,
                           ztz=ztz, 
                           vtv=vtv,
                           wtw=wtw, 
                           phi2=phi2, 
                           
                           cpara=cpara, 
                           ace_para=ace_para,
                           probit_para=probit_para, 
                           ordinal_para=ordinal_para,
                           tau_para=tau_para,
                           eta_para=eta_para, 
                           met_para=met_para, 
                           squish_para1=squish_para1, 
                           squish_para2=squish_para2,
                           sigma_para=sigma_para, 
                           v_para=v_para, 
                           z_para=z_para,
                           phi_para=phi_para, 
                           fx_para=fx_para, 
                           basis_para=basis_para, 
                           slope_para=slope_para, 
                           
                           alpha=alpha, 
                           betall=betall,
                           phi=phi, 
                           lambda=lambda, 
                           lambdai=lambdai, 
                           thetall=thetall, 
                           theta0=theta0,
                           xiparall=xiparall,
                           tau2all=tau2all, 
                           tauall=tauall, 
                           sigma2=sigma2, 
                           sigma=sigma,
                           
                           betallg=betallg,
                           rhoallg=rhoallg,
                           
                           thetallg=thetallg,
                           
                           zetag=zetag,
                           omegag=omegag,
                           zeta0g=zeta0g, 
                           omega0g=omega0g,
                           psiallg=psiallg,
                           psi0g=psi0g,
                           
                           mslope_logg=mslope_logg,
                           mslope0_logg=mslope0_logg,
                           mslopeg=mslopeg,
                           mslope0g=mslope0g,
                           mslope_log_sigmag=mslope_log_sigmag,
                           
                           fxgridallg=fxgridallg,
                           f0xgridg=f0xgridg)
  end_time   = Sys.time()
  run_time   = end_time - start_time
  
  alphag        = mcmc_out$mcmcg$alphag
  sigmag        = mcmc_out$mcmcg$sigmag
  sigma2_alphag = mcmc_out$mcmcg$sigma2_alphag
  sigma2_betag  = mcmc_out$mcmcg$sigma2_betag
  tauallg       = mcmc_out$mcmcg$tauallg
  eta0g         = mcmc_out$mcmcg$eta0g
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
  rhoallm       = mcmc_out$rhoallm
  rhoalls       = mcmc_out$rhoalls
  cuttallm      = mcmc_out$cuttallm
  cuttalls      = mcmc_out$cuttalls
  
  met_list = list()
  theta_met_tmp = vector("list", length = ngroups)
  for (j in 1:ngroups) {
    theta_met_tmp[[j]] = get_mat_to_list_met(mcmc_out$metg$theta_met[,j])
  }
  met_list$theta_met = theta_met_tmp
  
  if ((iflagsc==2)|(iflagsc==3)) {
    mslope_log_met_tmp = vector("list", length = ngroups)
    for (j in 1:ngroups) {
      mslope_log_met_tmp[[j]] = get_mat_to_list_met(mcmc_out$metg$mslope_log_met[,j])
    }
    met_list$mslope_log_met = mslope_log_met_tmp
  }
  
  if (iflagsc>=4) {
    
    if (iflagpsi>0) {
      psi_met_tmp = vector("list", length = ngroups)
      for (j in 1:ngroups) {
        psi_met_tmp[[j]] = get_mat_to_list_met(mcmc_out$metg$psi_met[,j])
      }
      met_list$psi_met = psi_met_tmp
    }
    
    zeta_met_tmp = vector("list", length = ngroups)
    for (j in 1:ngroups) {
      zeta_met_tmp[[j]] = get_mat_to_list_met(mcmc_out$metg$zeta_met[,j])
    }
    met_list$zeta_met = zeta_met_tmp
  }
  
  met_list$sigma2_met = get_mat_to_list_met(mcmc_out$metg$sigma2_met)
  
  # Compute summary statistics
  met_list$sigma2_met$pmet = met_list$sigma2_met$pmet/smcmc
  for (j in 1:ngroups) {
    met_list$theta_met[[j]]$pmet = met_list$theta_met[[j]]$pmet/smcmc
    if (iflagsc>=4) {
      met_list$zeta_met[[j]]$pmet = met_list$zeta_met[[j]]$pmet/smcmc
      if (iflagpsi==1) met_list$psi_met[[j]]$pmet = met_list$psi_met[[j]]$pmet/smcmc
    }
  }
  
  # Get HPD Intervals
  qs      = c(0.025,0.975)
  fxgridq = array(0,dim=c(nint+1,2,ngroups))
  for(j in 1:ngroups){
    fxjg         = fxgridallg[,j,]
    fxjq         = t(apply(fxjg,1,quantile,probs=qs))
    fxgridq[,,j] = fxjq
  }
  
  fxgridm     = fxgridm/smcmc
  fxgrids     = sqrt(abs(fxgrids - smcmc*fxgridm^2)/smcmc)
  f0xgridm    = f0xgridm/smcmc
  f0xgrids    = sqrt(abs(f0xgrids - smcmc*f0xgridm^2)/smcmc)
  fxdatam     = fxdatam/smcmc
  fxdatas     = sqrt(abs(fxdatas - smcmc*fxdatam^2)/smcmc)
  theta0m     = apply(theta0g,2,mean)
  theta0s     = apply(theta0g,2,sd)
  thetam      = thetam/smcmc
  thetas      = sqrt(abs(thetas-smcmc*thetam^2)/smcmc)
  colnames(thetam) = gnames
  colnames(thetas) = gnames
  rownames(thetam) = theta_names
  rownames(thetas) = theta_names
  eta0m      = matrix(apply(eta0g,2,mean))
  eta0s      = matrix(apply(eta0g,2,sd))
  rho_out    = NULL
  rhoallm    = NULL
  if (iflagACE == 1) {
    rhoallm    = matrix(apply(rhoallg,2,mean))
    rhoalls    = matrix(apply(rhoallg,2,sd))
    rhoall_pg0 = matrix(apply(rhoallg>0,2,mean))
    rho_out    = cbind(rhoallm,rhoalls,rhoall_pg0)
    colnames(rho_out)  = c("PostMean","PostSD","P(>0)")
    rownames(rho_out)  = gnames
  }
  
  
  alpha_out = NULL
  if(nparv > 0){
    alpham    = matrix(apply(alphag,2,mean))
    alphas    = matrix(apply(alphag,2,sd))
    alpha_pg0 = matrix(apply(alphag>0,2,mean))
    alpha_out = cbind(alpham,alphas,alpha_pg0)
    colnames(alpha_out) = c("PostMean","PostSD","P(>0)")
    rownames(alpha_out) = vnames
    valpham             = vdata %*% alpham
  }else{
    alpham  = 0
    alphas  = 1
    valpham = matrix(0,ntot,1)
  }
  if(iflagHBsigma == 0){
    sigmam = mean(sigmag)
    sigmas = sd(sigmag)
  }else{
    sigmam        = apply(sigmag,2,mean)
    sigmas        = apply(sigmag,2,sd)
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
  
  
  phim       = matrix(apply(phig,2,mean),nparz,nparw)
  phis       = matrix(apply(phig,2,sd),  nparz,nparw)
  phi_pg0    = matrix(apply(phig>0,2,mean),nparz,nparw)
  lambdam    = matrix(apply(lambdag,2,mean),nparw,nparw)
  lambdas    = matrix(apply(lambdag,2,sd),  nparw,nparw)
  betam      = betam/smcmc
  betas      = sqrt(abs(betas - smcmc*betam^2)/smcmc)
  beta_out   = list(betam,betas)
  
  # Estimated latents for Probit and Ordered Probit
  if(iflaglm>=1){
    ydatam = ydatam/smcmc
  }
  
  # Slope for monotone convex/concave
  mslope_out  = NULL
  mslope0_out = NULL
  if((iflagsc==2)|(iflagsc==3)){
    mslope_logm           = apply(mslope_logg,2,mean)
    mslope_logs           = apply(mslope_logg,2,sd)
    mslope_log_pg0        = apply(mslope_logg>0,2,mean)
    
    mslope0_logm           = mean(mslope0_logg)
    mslope0_logs           = sd(mslope0_logg)
    mslope0_log_pg0        = mean(mslope0_logg>0)
    
    mslopem               = apply(mslopeg,2,mean)
    mslopes               = apply(mslopeg,2,sd)
    mslope_pg0            = apply(mslopeg>0,2,mean)
    mslope_out            = matrix(cbind(mslopem,mslopes,mslope_pg0),ngroups,3)
    colnames(mslope_out)  = c("PostMean","PostSD","P>0")
    rownames(mslope_out)  = gnames 
    
    mslope0m              = mean(mslope0g)
    mslope0s              = sd(mslope0g)
    mslope0_pg0           = mean(mslope0g>0)
    mslope0_out           = matrix(c(mslope0m,mslope0s,mslope0_pg0),3,1)
    colnames(mslope0_out) = "Slope for Monotone"
    rownames(mslope0_out) = c("PostMean","PostSD","P>0")
    
    mslope_log_sigmam     = apply(mslope_log_sigmag,2,mean)
    mslope_log_sigmas     = apply(mslope_log_sigmag,2,sd)
    
  }  # End Estiators for slope of monotone convex/concave
  
  omega_out   = NULL
  omega0_out  = NULL
  
  if(iflagsc >= 4){
    if(nomega == 1){
      onames = "omega"
    }else{
      ostat_names = paste("_Stat",1:nstat,sep="")
      oflex_names = paste("_Flex",1:nflex,sep="")
      onames      = paste("Omega",1:nomega,sep="")
      onames[id_stat] = paste(onames[id_stat],ostat_names,sep="_")
      onames[id_flex] = paste(onames[id_flex],oflex_names,sep="_")
    }
    
    if(nomega == 1){
      omega0m    = mean(omega0g)
      omega0s    = sd(omega0g)
      omega0_pg0 = mean(omega0g>0) 
      omega0_out = matrix(c(omega0m,omega0s,omega0_pg0))
      rownames(omega0_out)   =  c("PostMean","PostSD","P(>0)")
      colnames(omega0_out)   = "Omega0"
      
      omegam     = apply(omegag, c(1,2), mean)
      omegas     = apply(omegag, c(1,2), sd)
      omegam     = matrix(omegam)
      omegas     = matrix(omegas)
      rownames(omegam)    = gnames
      rownames(omegas)    = gnames
      omega_out = cbind(omegam,omegas)
      colnames(omega_out) = c("Omega_PostMean","PostSD")
      rownames(omega_out) = gnames
    }else{
      omega0m    = matrix(apply(omega0g,2,mean))
      omega0s    = matrix(apply(omega0g,2,sd))
      omega0_pg0 = matrix(apply(omega0g>0,2,mean)) 
      omega0_out = cbind(omega0m,omega0s,omega0_pg0)
      colnames(omega0_out)   = c("Omega0_PostMean","PostSD","P(>0)")
      rownames(omega0_out)   = onames
      
      omegam    = apply(omegag, c(1,2), mean)
      omegas    = apply(omegag, c(1,2), sd)
      omegam    = t(matrix(omegam,nomega,ngroups))
      omegas    = t(matrix(omegas,nomega,ngroups))
      colnames(omegam)    = onames
      colnames(omegas)    = onames
      rownames(omegam)    = gnames
      rownames(omegas)    = gnames
      omega_out = list(omegam,omegas)
    }
    
    if((iflagpsi>0)|(iflagsc>4)){
      psim    = apply(psiallg,2,mean)
      psis    = apply(psiallg,2,sd)
      psi0m   = mean(psi0g)
      psi0s   = sd(psi0g)
    }
  }
  
  if(iflaglm == 2){
    cuttallm   = cuttallm/smcmc
    cuttalls   = sqrt(abs(cuttalls - smcmc*cuttallm^2)/smcmc)
    cuttallm   = data.frame(cuttallm)
    cuttalls   = data.frame(cuttalls)
    colnames(cuttallm) = gnames
    colnames(cuttalls) = gnames
    cutnames          = paste("Cut",0:(length(cutt)-1))
    rownames(cuttallm) = cutnames
    rownames(cuttalls) = cutnames
    cutmm              = apply(cuttallm,1,mean)
    cutms              = apply(cuttallm,1,sd)
    cutsm              = apply(cuttalls,1,mean)
    cutss              = apply(cuttalls,1,sd)
    cutmin             = apply(cuttallm,1,min)
    cutmax             = apply(cuttallm,1,max)
    cutv               = apply(cuttalls^2,1,mean) + cutms^2
    cutsd              = sqrt(cutv)
    cutout             = cbind(cutmm,cutsm,cutsm,cutss,cutmin,cutmax)
    colnames(cutout)   = c("AVG_PostM","AVG_PostSD","STD_PostM","STD_PostSD","Min_PostM","Max_PostM")
    rownames(cutout)   = cutnames
    
  }
  
  wbetam     = apply(wdata*betam[id_beta,],1,sum)
  
  yhat     = valpham + wbetam 
  if (iflaglm == 0) {
    yresid = ydata - yhat
  } else {
    yresid = ydatam - yhat
  }
  
  if (iflagACE==1) {
    yresidACE = yresid
    for (j in 1:ngroups) {
      
      fxdj = fxdata[iptHB[[j]]]
      nj   = nobs[j]
      
      # normal  include Y residuals
      yresidj     = yresid[iptHB[[j]]]
      yresidj[-1] = yresidj[-1] - rhoallm[j]*(yresidj[-nj]-fxdj[-nj])
      yresidACE[iptHB[[j]]] = yresidj
      
    } # End ACE adjustment
  }
  
  
  
  if(iflagACE==1){
    rhoutm   = matrix(apply(rho_out[,1:2],2,mean))
    rhouts   = matrix(apply(rho_out[,1:2],2,sd))
    rhoutmin = matrix(apply(rho_out[,1:2],2,min))
    rhoutmax = matrix(apply(rho_out[,1:2],2,max))
  }
  
  
  if(nomega>0){
    # Get summary stats for lower level model
    ooutmm  = apply(omegam,2,mean)
    ooutsm  = apply(omegas,2,mean)
    ooutms  = apply(omegam,2,sd)
    ooutss  = apply(omegas,2,sd)
    oout    = cbind(ooutmm,ooutsm,ooutms,ooutss)
  }
  
  mcmc.out = list()
  if (nparv>0) mcmc.out$alphag = alphag
  mcmc.out$sigmag = sigmag
  if (iflagHBsigma == 1) {
    mcmc.out$sigma2_alphag = sigma2_alphag
    mcmc.out$sigma2_betag  = sigma2_betag
  }
  mcmc.out$fxgridallg   = fxgridallg
  mcmc.out$f0xgridg     = f0xgridg
  mcmc.out$tauallg      = tauallg
  mcmc.out$eta0g        = eta0g
  mcmc.out$theta0g      = theta0g
  mcmc.out$thetallg     = thetallg
  mcmc.out$lambdag      = lambdag
  mcmc.out$phig         = phig
  if (iflagACE==1) mcmc.out$rhoallg = rhoallg
  if ((iflagsc==2)||(iflagsc==3)) {
    mcmc.out$mslope_logg       = mslope_logg
    mcmc.out$mslope0_logg      = mslope0_logg
    mcmc.out$mslopeg           = mslopeg
    mcmc.out$mslope0g          = mslope0g
    mcmc.out$mslope_log_sigmag = mslope_log_sigmag
  }
  if (iflagsc >= 4) {
    mcmc.out$omegag  = omegag
    mcmc.out$omega0g = omega0g
    mcmc.out$zetag   = zetag
    mcmc.out$zeta0g  = zeta0g
    if (iflagpsi > 0) {
      mcmc.out$psiallg = psiallg
      mcmc.out$psi0g   = psi0g
    }
  }
  
  post.est = list()
  post.est$fxgridm  = fxgridm
  post.est$fxgridq  = fxgridq
  post.est$fxgrids  = fxgrids
  post.est$f0xgridm = f0xgridm
  post.est$f0xgrids = f0xgrids
  post.est$fxdatam  = fxdatam
  post.est$fxdatas  = fxdatas
  post.est$theta0m  = theta0m
  post.est$theta0s  = theta0s
  post.est$thetam   = thetam
  post.est$thetas   = thetas
  post.est$eta0m    = eta0m
  post.est$eta0s    = eta0s
  
  if(iflagACE == 1){
    post.est$rho_out = rho_out
  }
  if(nparv > 0){
    post.est$alpha_out = alpha_out
    post.est$valpham   = vdata %*% alpham
  } else {
    post.est$valpham   = matrix(0,ntot,1)
  }
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
  
  
  post.est$phim    = phim
  post.est$phis    = phis
  post.est$phi_pg0 = phi_pg0
  post.est$lambdam = lambdam
  post.est$lambdas = lambdas
  post.est$betam   = betam
  post.est$betas   = betas
  
  # Estimated latents for Probit and Ordered Probit
  if(iflaglm>=1){
    post.est$ydatam = ydatam
  }
  
  # Slope for monotone convex/concave
  if((iflagsc==2)|(iflagsc==3)){
    post.est$mslope_out  = mslope_out
    post.est$mslope0_out = mslope0_out
  }  # End Estiators for slope of monotone convex/concave
  
  if(iflagsc >= 4){
    
    post.est$omega_out  = omega_out
    post.est$omega0_out = omega0_out
    
    if(iflagpsi>0){
      post.est$psim    = psim
      post.est$psis    = psis
      post.est$psi0m   = psi0m
      post.est$psi0s   = psi0s
    }
  }
  
  if(iflaglm == 2){
    post.est$cutout   = cutout
  }
  
  wbetam = apply(wdata*betam[id_beta,],1,sum)
  if(iflagACE==1){
    post.est$rhoutm   = rhoutm
    post.est$rhouts   = rhouts
    post.est$rhoutmin = rhoutmin
    post.est$rhoutmax = rhoutmax
  }
  
  
  if(nomega>0){
    post.est$oout = oout
  }
  
  # *******************************************************
  out=list()
  out$met.out=met_list
  out$mcmc.out=mcmc.out
  out$post.est=post.est
  
  out$valpham = valpham
  out$wbetam  = wbetam
  out$fxdatam = fxdatam
  
  out$y           = ydata
  out$v           = vdata
  out$x           = xdata 
  out$w           = wdata
  out$z           = zdata
  out$shape       = shape
  out$nknots      = nknots
  out$xgrid       = xgrid
  out$group       = id_group
  out$ngroup      = ngroups
  out$iflaglm     = iflaglm
  out$iflagCenter = iflagCenter
  out$yresid      = yresid
  out$mcmc.time   = run_time
  if(iflagACE==1){
    out$yresidACE = yresidACE
  }
  
  return(out)
}
