# ***************************************************
# MakeHBPointer
#    Make a pointer into the matrics for groups
#    Input
#        id   = identifies groups
#    Output
#        nobs     = number of observations in each group
#        iptd     = list of pointers into group
#        id_beta  = indices to expanded beta to match wdata
#        id_first = index of first observation in each group
#                   used with ACE
#    Useage
#        xj   = xdata[iptd[[j]],]  gives you the x data for group j
#        Compute wdata_j*beta_j for each group
#        wdata*beta = matrix(apply(wdata*bpar[id_beta,],1,sum))
#    Call
#        outpt = MakeHBPointer(id)
# ***************************************************
MakeHBPointer <- function(id){
  uid = sort(unique(id))   # unique populations
  ng  = length(uid)        # number of unique groups
  id_first = matrix(0,ng,1)
  
  # id_group may not be sequential from 1 ... ngroups
  # Need an index to get correct betam
  
  nobs = matrix(0,ng,1)          # number of observations in each group
  id_beta  = matrix(0,ng,1)      # match rows of beta to wdata
  for(j in 1:ng){                # loop over groups
    jrow = which(id == uid[j])   # identify rows for group j
    id_first[j] = jrow[1]        # first observation
    nj   = length(jrow)          # Number of observations in group j
    id_beta[jrow] = matrix(j,nj,1)
    
    nobs[j] = length(jrow)       # number of observations in group j
    if(j == 1){
      iptd = list(jrow)          # start the list iptd
    }else{
      iptd[[j]] = jrow           # add new member to list
    }
  }  # end loop over groups
  return(list("nobs"=nobs,"iptd"=iptd,"id_first"=id_first,"id_beta"=id_beta))
} # End MakeHBPointer


# *************************************************************
# Use Metropolis-Hastings with shape constraints
# * Random walk Metropolis-Hastings with t-errors for fmodel > 1
# *  pmet = P(accept candidate theta)
# * If pmet is > 0.7, then step size is too small: increase metm 
# * If pmet is < 0.3, then step size is too large: decrease metm 
# * Keep mets > metm  (std dev > mean)
# * Why?  If mets > metm, then alpha < 3
# * Also IG has a dead zone near zeros if alpha > 2.25
# * So, you want to keep 2 < alpha < 2.25.
# * Fix alpha.  Focus on metm to get mets and beta
# **************************************************************
# GetMetList
#    Constructor function to initialize list of parameters for
#    adaptive metropolis
#   xxxx_met = GetMetList(w,m0,alpha,ndim)
GetMetList <- function(w=.5,m0=.001,alpha=4,ndim=1){
  beta        = m0*(alpha-1)
  m_AM        = matrix(m0,ndim,1)
  met_var     = matrix(m0,ndim,1)
  pmet        = 0
  icount      = 0
  beta_AM     = matrix(beta,ndim,1)
  x = list("w"=w,"m0"=m0,"alpha"=alpha,"beta"=beta, "ndim" = ndim,
           "m_AM"=m_AM,"beta_AM"= beta_AM,
           "met_var"=met_var,"pmet"=0,"icount"=0);
  return(x)
}  # End constructor function for adaptive Metropolis parameters


# ***************************************************
# Generate variance parameter for adaptive Metropolis
#  var_met = AdaptMetVar(xxxx_met)
AdaptMetVar <- function(x){
  met_var  = matrix(x$beta_AM/rgamma(x$ndim,shape=x$alpha))
  return(met_var)
} # End generate variance parameter for adaptive Metropolis


# ***************************************************
# AdaptMetUpdate
#   Update parameters in adaptive Metropolis
AdaptMetUpdate <- function(x){
  x$icount   = x$icount + 1
  x$m_AM     = x$m_AM + (x$met_var - x$m_AM)/x$icount
  mnew       = x$w*x$met_var + (1-x$w)*x$m_AM
  x$beta_AM  = (x$alpha-1)*mnew
  return(x)
}  # End adaptive Metropolis update


# ***************************************************
# UpdateMet
#   Adjust mean of t-distribution for Metropolis based on acceptance rate
#   Used warm-up iterations before regular MCMC iterations
#   xxxx_met = UpdateMet(xxx_met,nb=nblow0,fct=10)
UpdateMet <- function(x,nb=nblow0,fct=10){
  pmet   = x$pmet/nb
  ndim   = x$ndim
  if (pmet > 0.6){  #// Too many acceptances, increase m0  
    m0 = x$m0*fct
    x  = GetMetList(m0=m0,ndim=ndim)
  }else if (pmet < 0.3){  #// Too few acceptances, decrease metm 
    m0 = x$m0/fct;
    if(m0 < 1e-8) m0 = 1e-8
    x  = GetMetList(m0=m0,ndim=ndim)
    #cat(paste(j,'th group', " pmet = ",pmet[j]/nblow0," < 0.3. Reduce metm and redo MCMC loop",'\n', sep=''));
  }else{  # Seems to be a good rejection rate: reset pmet 
    x$pmet   = 0
  }  # End adjust m0 for group j for metropolis based on rejection rates
  return(x)
}


# ***************************************************
# GetSquish
#     Compute squish function
#     Softmax modification of square wave 
#     Used with S-shapes and multi-modal functiona to 
#     flip second derivatives at inflection points
# INPUT
#     omega     = inflection points
#     psi       = slope
# OUTPUT
#     hfun      = squish function of xgrid
# hfun = GetSquish(omega,psi)
# ***************************************************
GetSquish = function(omega,psi,xgrid,nint){
  nflex   = length(omega)
  deven   = as.numeric(nflex%%2 == 0)
  hfun    = matrix(deven,nint+1,1)    # 1 + sum_k  ... if even
  # Loop over inflection points
  for(j in 1:nflex){
    z      = psi*(xgrid-omega[j])
    # Worry about underflow and over flow
    # if z < -50, then hfun = 1
    # if z >  50, then hfun = -1
    idz     = which(z< -50)
    if(length(idz>0)) z[idz] = -50
    idz     = which(z>50)
    if(length(idz>0)) z[idz] = 50
    hj     = (1-exp(z))/(1+exp(z))
    hfun     = hfun + (-1)^(j-1)*hj
  }  # End loop over inflection points
  return(hfun)
}  # End Squish


# ***************************************************
# GetMonofxgrid
#   Compute increasing f = int_a^x exp[Z(s)] ds on xgrid
# INPUT
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************
GetMonofxgrid  = function(fpars){
  z     = fpars$phixgrid %*% fpars$theta
  if(fpars$iflagZ == 1){
    dfx  = exp(z) 
  }else{
    dfx  = z^2
  }
  
  fx         = fpars$iflagpn*CumTrap(dfx,fpars$xdelta)
  if(fpars$iflagCenter==0){  # do not center f
    return(fx)
  }else{  # center f
    # Compute mean
    c         = trapint(fx,fpars$xdelta)
    fx0       = fx - c/fpars$xrange
    return(fx0)
  }
}  # End GetMonofxgrid


# ***************************************************
# GetMonoSamefxgrid
#   First and second derivative have the same sign
#   delta = 1  => increasing and convex
#   delta = -1 => decreasing and concave
#   f(x) = int_a^x int_a^s g[Z(t)] dt ds on xgrid
# INPUT
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************
GetMonoSamefxgrid  = function(fpars){
  z     = fpars$phixgrid %*% fpars$theta
  if(fpars$iflagZ == 1){
    d2fx  = exp(z)
  }else{
    d2fx  = z^2
  }
  
  dfx         = CumTrap(d2fx,fpars$xdelta)
  fx          = CumTrap(dfx,fpars$xdelta) + (fpars$xgrid-fpars$xmin)*fpars$mslope
  fx          = iflagpn*fx
  if(fpars$iflagCenter==0){  # do not center f
    return(fx)
  }else{  # center f
    # Compute mean
    c         = trapint(fx,fpars$xdelta)
    return(fx-c/fpars$xrange)
  }
}  # End GetMonoSamefxgrid


# ***************************************************
# GetMonoDiffxgrid
#   First and second derivative have different signs
#   delta = 1  => increasing and convex
#   delta = -1 => decreasing and concave
#   f(x) = int_a^(b-x) int_a^s g[Z(t)] dt ds on xgrid
# INPUT
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************
GetMonoDiffxgrid  = function(fpars){
  z     = fpars$phixgrid %*% fpars$theta
  if(fpars$iflagZ == 1){
    d2fx  = exp(z)
  }else{
    d2fx  = z^2
  }
  dfx         = rev(CumTrap(rev(d2fx),fpars$xdelta))  # int_s^b f(t)dt
  fx          = CumTrap(dfx,fpars$xdelta) + (fpars$xgrid-fpars$xmin)*fpars$mslope # int_a^x
  fx          = fpars$iflagpn*fx
  if(fpars$iflagCenter==0){  # do not center f
    return(fx)
  }else{  # center f
    # Compute mean
    c         = trapint(fx,fpars$xdelta)
    return(fx-c/fpars$xrange)
  }
}  # End GetMonoDiffxgrid


# ***************************************************
# GetSfxgrid
#   S-shaped
#   iflagpn =  1 and d3 =  1 is increasing convex-to-concave
#   iflagpn = -1 and d3 =  1 is decreasing concave-to-convex
#   iflagpn =  1 and d3 = -1 is increasing concave-to-convex
#   iflagpn = -1 and d3 = -1 is decreasing and convex-to-concave
#   f(x) = int_a^x int_a^s h(t)g[Z(t)] dt ds on xgrid
# INPUT
#    Input 
#       fpars is a list
#       d3 is 1 or -1
#   Output
#       fx             = function computed on xgrid
# ***************************************************
GetSfxgrid  = function(fpars,d3){
  z     = fpars$phixgrid %*% fpars$theta
  if(fpars$iflagZ == 1){
    d2fx  = exp(z)
  }else{
    d2fx  = z^2
  }
  # d3 =  1 -> hfun is "squish down" for "S"
  # d3 = -1 -> hfun is "squish up" for "rotated S"
  d2fx        = d3*fpars$hfun*d2fx
  dfx         = CumTrap(d2fx,fpars$xdelta)
  c2          = min(c(0,dfx))
  fx          = CumTrap(dfx,fpars$xdelta)
  fx          = iflagpn*(fx - c2*(fpars$xgrid-fpars$xmin))
  if(fpars$iflagCenter==0){  # do not center f
    return(fx)
  }else{  # center f
    # Compute mean
    c         = trapint(fx,fpars$xdelta)
    fx        = fx - c/fpars$xrange
    return(fx)
  }
}  # End GetSfxgrid


# ***************************************************
# Get convex, U-shaped or concave, inverted U shaped on segement x
# INPUT
#    x     = grid between neighboring inflection points
#    f     = h(x)*g(Z(x)) computed on x
#    omega = stationary point
# ***************************************************
GetUseg   = function(x,gz,omega,xdelta){
  fx        = CumTrap(gz,xdelta)
  fx        = CumTrap(fx,xdelta)
  idbw      = which(x <= omega)      # indicies for x <= omega
  zg1       = gz[idbw]               # g(Z) from xmin to omega
  zg1       = trapint(zg1,xdelta)    # int_a^omega g(Z)
  fx        = (fx - zg1*(x-x[1]))
  
  return(fx)
}  # End get U on segement


# ***************************************************
# GetMMfxgrid
#   Multi-modal   functions
#   iflagpn =  1 U, inverted U, ... or min, max, min
#   iflagpn = -1 Inverted U, U, ... or max, min, max
# INPUT
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************
GetMMfxgrid  = function(fpars){
  z     = fpars$phixgrid %*% fpars$theta
  omega_s = fpars$omega[fpars$id_stat]
  idb     = fpars$idb
  if(fpars$iflagZ == 1){
    d2fx  = fpars$hfun*exp(z)
  }else{
    d2fx  = fpars$hfun*z^2
  }
  fx       = 0  # Note: f(xmin) = 0.
  for(k in 1:fpars$nstat){
    idk    = idb[k]:idb[k+1]
    xk     = fpars$xgrid[idk];
    d2fxk  = d2fx[idk];
    fxk    = GetUseg(xk,d2fxk,omega_s[k],fpars$xdelta)
    fxk    = fxk[-1] + fx[length(fx)]              
    fx     = matrix(c(fx,fxk))
  } # End get f on segments
  fx      = -fpars$iflagpn*fx
  if(fpars$iflagCenter==0){  # do not center f
    return(fx)
  }else{  # center f
    # Compute mean
    c         = trapint(fx,fpars$xdelta)
    fx        = fx - c/fpars$xrange
    return(fx)
  }
}  # End GetMMfxgrid


# ***************************************************
# GetUfxgrid
# U Convex or Inverted U Concave
# iflagpn = -1 for U Convex
# iflagpn =  1 for Inverted U Concave
# INPUT
#   hfun    = squish function
#   omega   = stationary point
#   fpara   = theta parameters for f
#   phix    = basis functions on xgrid
#   xdelta  = interval size
#   xrange  = xmax - xmin is range of x data
#   iflagpn = 1 Inverted U concave and -1 U concave
# Output
#   fx      = f evaluated at x
#   fpars   = list(theta,hfun,omega)
# ***************************************************
GetUfxgrid  = function(fpars){
  z          = fpars$phixgrid %*% fpars$theta
  if(fpars$iflagZ == 1){
    dfx  = exp(z)
  }else{
    dfx  = z^2
  }
  idbw      = which(fpars$xgrid <= as.numeric(fpars$omega))  # indiciec for x <= omega
  zg1       = dfx[idbw]      # g(Z) from xmin to omega
  zintaw    = trapint(zg1,fpars$xdelta)  # int_a^omega g(Z)
  fx        = CumTrap(dfx,fpars$xdelta)
  fx        = CumTrap(fx,fpars$xdelta)
  fx        = fpars$iflagpn*(fx - zintaw*(fpars$xgrid-fpars$xmin))
  if(fpars$iflagCenter==0){  # do not center f
    return(fx)
  }else{  # center f
    # Compute mean
    c         = trapint(fx,fpars$xdelta)
    return(fx-c/fpars$xrange)
  }
}  # End GetUfxgrid


# ***************************************************
# Getfx
#    Compute function on x
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************
Getfx = function(fpars){
  fx     = 0
  # free function
  if(fpars$iflagsc == 0)  return(fpars$phixgrid %*% fpars$theta)
  # ***************************************************
  # monotone function
  if(fpars$iflagsc == 1){
    fx = GetMonofxgrid(fpars)
    return(fx)
  } 
  
  # Increasing convex or decreasing and concave
  if(fpars$iflagsc == 2){
    fx = GetMonoSamefxgrid(fpars)
    return(fx)
  }  
  # Increasing concave or decreasing and convex
  if(fpars$iflagsc == 3){
    fx = GetMonoDiffxgrid(fpars)
    return(fx)
  }  
  # ***************************************************
  # U-shaped function
  if(fpars$iflagsc == 4){
    fx = GetUfxgrid(fpars)
    return(fx)
  }  
  # S-shaped (increasing,convex-to-concave) or (decreasing,concave-to-convex)
  if(fpars$iflagsc == 5){
    fx  = GetSfxgrid(fpars,1)
    return(fx)
  }
  # S-shaped (increasing,concave-to-convex) or (decreasing,convex-to-concave)
  if(fpars$iflagsc == 6){
    fx  = GetSfxgrid(fpars,-1)
    return(fx)
  }
  
  # Multi-modal functions
  if(fpars$iflagsc == 7){
    fx  = GetMMfxgrid(fpars)
  }
  
}  # End Getfx


# ***************************************************
# GetSCfxobs
#   Compute f at observations
#   INPUT
#      xgrid   = grid for x
#      xinx    = gives index for xi in xgrid
#                xgrid[xinx[i]] < xdata[i] <= xgrid[xinx[i]+1]
#      xover   = xdata - xgrid[xinx] is excess over boundary
#      fxgrid  = f computed on xgrid
#                Note: f is increasing on xgrid
#   OUTPUT
#      fxobs   = f computed at observation
# fxobs = GetSCfxobs(xgrid,xinx,xover,fxgrid)
# ***************************************************
GetSCfxobs = function(xgrid,xinx,xover,fxgrid){
  n      = length(xinx)
  fxobs  = matrix(0,n,1)
  f0     = fxgrid[xinx[,1]]
  f1     = fxgrid[xinx[,2]]
  fxobs  = f0 + (f1-f0)*xover
  return(fxobs)
}


# ***************************************************
# ***************************************************
# More Functions
# ***************************************************
# TrapInt does trapidzoid integration
# Int f(x)*dx = xdelta*(f1/2 + f2 + ... + fn-1 + fn/2)
# ***************************************************
trapint <- function(f,xdelta){
  nr = length(f)
  a  = sum(f)
  a  = a - (f[1] + f[nr])/2
  return(a*xdelta)
}


# ***************************************************
# ***************************************************
# Int2D
#  Integral over 2 x 2 grid
#   int_a^x1 int_b^x2 f(s,t) dx1*dx2 
#   for (x1,x2) in grid of values
#   Simple Trapiziod.  
#   Compute heights in each cell = average corner falues of f
#   Cumsum up and down grid
#
#   INPUT
#      fmat   = f computed over grid
#      xdelta = vector of delta
#
#   OUTPUT
#      fint   = int_a^x1 int_b^x2 f(s,t) dx1*dx2 at grid
#   CALL
#      fint   = Int2D(fmat,xdelta)
# ***************************************************
Int2D = function(fmat,xdelta) {
  ny   = ncol(fmat)
  nx   = nrow(fmat)
  # Compute average of the 4 fmat values at each cell
  fint = fmat[,-ny] + fmat[,-1]
  fint = fint[-nx,] + fint[-1,]
  fint = sum(fint)/4   
  fint = fint*(xdelta[1]*xdelta[2])  # Multiply by length and width of grids.
  return(fint)
}  # End Int2D
# ***************************************************
# ***************************************************


# ***************************************************
# ***************************************************
# CumInt2D
#   Cumulative integral over 2 x 2 grid
#   int_a^x1 int_b^x2 f(s,t) dx1*dx2 
#   for (x1,x2) in grid of values
#   Simple Trapiziod.  
#   Compute heights in each cell = average corner falues of f
#   Cumsum up and down grid
#
#   INPUT
#      fmat   = f computed over grid
#      xdelta = vector of delta
#
#   OUTPUT
#      fint   = int_a^x1 int_b^x2 f(s,t) dx1*dx2 at grid
#   CALL
#      fint   = CumInt2D(fmat,xdelta)
# ***************************************************
CumInt2D = function(fmat,xdelta){
  ny   = ncol(fmat)
  nx   = nrow(fmat)
  # Compute average of the 4 fmat values at each cell
  fint = fmat[,-ny] + fmat[,-1]
  fint = fint[-nx,] + fint[-1,]
  fint = fint/4   
  # Add 0 to first row and first column
  fint = cbind(matrix(0,nx-1),fint)
  fint = rbind(matrix(0,1,ny),fint)
  # Cumsum volumns of blocks in rows and columns
  fint = t(apply(fint,1,cumsum))  # CAUTION  apply & cumsum dumps results in columns
  fint = apply(fint,2,cumsum)
  fint = fint*(xdelta[1]*xdelta[2])  # Multiply by length and width of grids.
  return(fint)
} # End CumInt2D


# ***************************************************
# CUMTrap computes cumulative integral of function
# int_a^x f(s)ds over a grid of points a, a+delta, a+2*delta, ..., b
# Based on Trapizod rule.
# INPUT
#   f     = function on grid
#   delta = size of grid interval
# OUTPUT
#   fcum  = cumulative function corresonding to f
# ***************************************************
CumTrap = function(f,delta){
  nr    = NROW(f)
  fcum	= cumsum(f[1:nr-1]) + cumsum(f[2:nr])
  fcum	= fcum*delta/2
  fcum	= matrix(c(0,fcum),nr,1)
  return(fcum)
}


# ***************************************************
#   RNDTNA generates form a truncated above normal distribution
#   Truncated on X < xtop
#    INPUT
#		mu		= vector of means
#		sigma		= vector of stds
#		xbot		= lower limit
#   OUTPUT
#		x		= truncated normal
# ***************************************************
rndtna = function(mu,sigma,xtop){
  u = runif(n=length(mu))
  fa = 0
  fb = pnorm(xtop,mean=mu,sd=sigma)  # normal cdf
  p  = fa + u*(fb-fa)
  # Worry about p = 0 or p = 1
  z2 = p <= 0.00001
  if(sum(z2)>0) p  = z2*(0.00001 - p) + p
  z2 = p >= 0.99999
  if(sum(z2)>0) p  = z2*(0.99999 - p) + p
  x  = qnorm(p,mean=mu,sd=sigma)  # inverse normal cdf
  # make sure that x < xtop
  z  = x > xtop
  if(sum(z)>0) x = z*(xtop-x) + x
  return(x)
}


# ***************************************************
#   RNDTNB generates form a truncated normal distribution
#   Truncated on [xbot, infinity)
#    INPUT
#		mu		= vector of means
#		sigma		= vector of stds
#		xbot		= lower limit
#   OUTPUT
#		x		= truncated normal
# ***************************************************
rndtnb = function(mu,sigma,xbot){
  u = runif(n=length(mu))
  fa = pnorm(xbot,mean=mu,sd=sigma)  # normal cdf
  fb = 1
  p  = fa + u*(fb-fa)
  # Worry about p = 0 or p = 1
  z2 = p <= 0.00001
  if(sum(z2)>0) p  = z2*(0.00001 - p) + p
  z2 = p >= 0.99999
  if(sum(z2)>0) p  = z2*(0.99999 - p) + p
  x  = qnorm(p,mean=mu,sd=sigma)  # inverse normal cdf
  # make sure that x > xbot
  z  = x < xbot
  if(sum(z)> 0) x = z*(xbot-x) + x
  return(x)
}


# ***************************************************
#   RNDTNAB generates form a truncated normal distribution
#   Truncated on xbot < X < xtop
#    INPUT
#		mu		= vector of means
#		sigma		= vector of stds
#		xbot		= lower limit
#   OUTPUT
#		x		= truncated normal
# ***************************************************
rndtnab = function(mu,sigma,xbot,xtop){
  u = runif(n=length(mu))
  fa = pnorm(xbot,mean=mu,sd=sigma)  # normal cdf
  fb = pnorm(xtop,mean=mu,sd=sigma)  # normal cdf
  # Need to gaurd against bad falues for fa and fb
  x0 = (xbot+xtop)/2                 # Good guess that stasifies xbot < x < xtop
  z  = fa < fb                       # Idenitfy correct ordering for z
  delta = z*(fb-fa)                  # delta = fb-fa if fb>fa
  # delta = 0     if fb<=fa
  p  = fa + u*delta                 
  # Worry about p = 0 or p = 1
  z2 = p <= 0.00001
  if(sum(z2)>0) p  = z2*(0.00001 - p) + p
  z2 = p >= 0.99999
  if(sum(z2)>0) p  = z2*(0.99999 - p) + p
  x  = qnorm(p,mean=mu,sd=sigma)     # inverse normal cdf
  if(sum(z)> 0) x  = z*(x-x0) + x0   # x = x if fa < fb
  # x = x0 if fa >= fb
  # force x between a and b
  z  = x < xbot
  if(sum(z) > 0) x = z*(xbot-x) + x
  z  = x > xtop
  if(sum(z)>0)   x = z*(xtop-x) + x
  
  
  return(x)
}


"get_mat_to_list_met" <- function(met_mat) {
  ndim = met_mat[1]
  met_list = list(
    ndim    = ndim,
    pmet    = met_mat[2],
    icount  = met_mat[3],
    alpha   = met_mat[4],
    beta    = met_mat[5],
    w       = met_mat[6],
    m0      = met_mat[7],
    beta_AM = met_mat[8          : (8 +   ndim-1)],
    m_AM    = met_mat[(8+ndim)   : (8 + 2*ndim-1)],
    met_var = met_mat[(8+2*ndim) : (8 + 3*ndim-1)])
  
  return(met_list)
}

# ******************************************************************************
# ******************************************************************************
# GetSmooth
# Generate smoothing parameters tau2 and gamma
# INPUT
#    tau2      = Current value variance parameter
#    gampar    = Current value of rate of decline
#    theta     = Fourier coefficients
#    egamkall  = exp(-gampar*k)
#    wk        = sum(kall)/2 - w0 where kall = 1, 2, ..., nbasis
#    sigma2    = error variance
# OUTPUT LIST
#    tau2
#    gampar
#    egamkall
# CALL
# newsmooth = getsmooth(tau2,gampar,theta,egamkall,wk)
getsmooth = function(tau2,gampar,theta,egamkall,wk,tau2_s0,tau2_rn,freq_sum){
  # Full conditional of tau2 is IG(rn/2,sn/2)
  #  rgamma(n, shape, rate = 1, scale = 1/rate)
  #  The Gamma distribution with parameters shape = a and scale = s has density
  #      f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
  nb          = length(egamkall)
  theta2      = theta^2
  tau2_sn     = tau2_s0 + sum(theta2/(egamkall))   # Posterior rate parameters
  tau2        = tau2_sn/(2*rgamma(1,tau2_rn/2))
  tau2i       = 1/tau2
  tau         = sqrt(tau2)
  
  # ****************************************************************************
  # Generate gamma with slice sampling
  ck          = matrix(theta2/(2*tau2),nb,1)
  # Worry about theta_k = 0.  Big hassel, but it matters because we take logs
  idtheta0  = which(ck == 0)
  ntheta0   = length(idtheta0)
  if(ntheta0 == nb){
    gampar  = .1 # all rows = 0. punt!
  }else{
    if(ntheta0>0) ck[idtheta0] = matrix(1,ntheta0,1)  # set component to 1
    u       = runif(nb)                           # Uniform rvs
    bn 		  = gampar + (log(ck - log(u)*egamkall) - log(ck))/freq_sum
    #if(ntheta0>0) bn		  = bn[-idtheta0]			                      # Remove entries with ck = 0 @
    bmin    = min(bn)
    u  		  = runif(1)
    gampar  = bmin + log(u + (1-u)*exp(-wk*bmin))/wk
  }
  egamkall  = matrix(exp(-gampar*freq_sum))
  iz0 = which(egamkall < 1E-15)   # Guard against numerical errors
  if(length(iz0)>0) egamkall[iz0] = 1E-15
  return(list("tau2"=tau2,"gampar"=gampar,"egamkall"=egamkall))
}  # End getsmooth
# ******************************************************************************
# ******************************************************************************
# GetSmoothSC
# Generate smoothing parameters tau2 and gamma
# Assumpes shape constraints
# INPUT
#    tau2      = Current value variance parameter
#    gampar    = Current value of rate of decline
#    theta     = Fourier coefficients
#    freqs     = Frequencies of cosines
#    egamkall  = exp(-gampar*k)
#    wk        = sum(kall)/2 - w0 where kall = 1, 2, ..., nbasis
#    sigma2    = error variance
# OUTPUT LIST
#    tau2
#    gampar
#    egamkall
# CALL
# newsmooth = getsmoothSC(tau2,gampar,theta,freqs,egamkall,wk)
getsmoothSC = function(tau2,gampar,theta,freqs,egamkall,wk,tau2_s0,tau2_rn,freq_sum,iflagZ){
  # Full conditional of tau2 is IG(rn/2,sn/2)
  #  rgamma(n, shape, rate = 1, scale = 1/rate)
  #  The Gamma distribution with parameters shape = a and scale = s has density
  #      f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
  nb          = length(egamkall)
  #  Need to worry about log(theta_0) if using Z^2
  theta2   = theta^2
  if(iflagZ == 0){
    xipar     = log(theta[1])
    theta2[1] = xipar^2
  }
  tau2_sn     = tau2_s0 + sum(theta2/(egamkall))   # Posterior rate parameters
  tau2        = tau2_sn/(2*rgamma(1,tau2_rn/2))
  tau2i       = 1/tau2
  tau         = sqrt(tau2)
  
  # ****************************************************************************
  # Generate gamma with slice sampling
  # Drop theta[1], which is the constant
  theta2    = theta2[-1]
  egamkall  = egamkall[-1]
  nb        = nb-1
  ck          = matrix(theta2/(2*tau2),nb,1)
  # Worry about theta_k = 0.  Big hassle, but it matters because we take logs
  idtheta0  = which(ck == 0)
  ntheta0   = length(idtheta0)
  if(ntheta0 == nb){
    gampar  = .1 # all rows = 0. punt!
  }else{
    if(ntheta0>0) ck[idtheta0] = matrix(1,ntheta0,1)  # set component to 1
    u       = runif(nb)                           # Uniform rvs
    bn 		  = gampar + (log(ck - log(u)*egamkall) - log(ck))/freqs[-1]
    #if(ntheta0>0) bn		  = bn[-idtheta0]			                      # Remove entries with ck = 0 @
    bmin    = min(bn)
    u  		  = runif(1)
    gampar  = bmin + log(u + (1-u)*exp(-wk*bmin))/wk
  }
  egamkall  = matrix(exp(-gampar*freqs))
  iz0 = which(egamkall < 1E-15)   # Guard against numerical errors
  if(length(iz0)>0) egamkall[iz0] = 1E-15
  return(list("tau2"=tau2,"gampar"=gampar,"egamkall"=egamkall))
}  # End getsmoothSC
# ******************************************************************************
# ******************************************************************************
# fdata_apx 
#   Approximate f at data by using f computed on xgrid
#   Weighted average of corner points for grid cell that data are in
# INPUT
#   fmat        = f computed on xgrid
#   fapx_weight = weights for the four corners based
#   xinx        = identifies indices of grid where data belongs
# RETURN
#   fdata       = approxmation of f at data
# CALL
#   fdata  = fdata_apx(fmat,f12apx_weights,xinx)
# ******************************************************************************
fdata_apx = function(fmat,f12apx_weights,xinx,nint){
  nobs = nrow(xinx)       # number of observations
  # Approximate fdata_anal by using grid
  fdata  = matrix(0,nobs,1)
  for(i in 1:nobs){
    i1   = xinx[i,1]      # Grid index for x1
    i2   = xinx[i,2]      # Grid index for x2
    fi   = 0              # Weighted average of neighbors
    # Not on right & not on top boundary
    if((i1 < nint+1)&(i2 < nint+1)){
      fi   = fmat[i1,i2]*f12apx_weights[i,1] + 
        fmat[i1+1,i2]*f12apx_weights[i,2] +
        fmat[i1,i2+1]*f12apx_weights[i,3] + 
        fmat[i1+1,i2+1]*f12apx_weights[i,4]
    }
    # Observation on right boundary and not top
    if((i1 == nint+1)&(i2 < nint+1)){
      fi  = fmat[i1,i2]*(1-xover[i,2]) + fmat[i1,i2+1]*xover[i,2]
    }
    # Observation on top boundary and not right
    if((i1 < nint+1)&(i2 == nint+1)){
      fi  = fmat[i1,i2]*(1-xover[i,1]) + fmat[i1+1,i2]*xover[i,1]
    }
    # Observation on right & top boundary: fmat is exact
    if((i1==nint+1) & i2==nint+1)  fi = fmat[i1,i2]
    fdata[i] = fi
  }
  return(fdata)
}  # End fdata_apx to approximate fdata with f at grid

# ******************************************************************************
# ******************************************************************************
# GetMainfxobs
#   Compute f at observations for main effects
#   INPUT
#      xgrid   = grid for x
#      xinx    = gives index for xi in xgrid
#                xgrid[xinx[i]] < xdata[i] <= xgrid[xinx[i]+1]
#      xover   = (xdata - xgrid[xinx])/xdelta is excess over boundary
#      idn0    = index of xdata that is not on grid
#      fxgrid  = f computed on xgrid
#                Note: f is increasing on xgrid
#   OUTPUT
#      fxobs   = f computed at observation
# fxobs = GetMainfxobs(xgrid,xdelta,xinx,xover,fxgrid)
# ******************************************************************************
GetMainfxobs = function(xgrid,xdelta,xinx,xover,fxgrid){
  n      = length(xinx)
  fxobs  = matrix(0,n,1)
  idn0   = which(xover>0)
  if(sum(idn0)<n){  # Got data on grid points
    fxobs[-idn0] = fxgrid[xinx[-idn0]]
  }else{  # overshoot
    y0     = fxgrid[xinx[idn0]]
    y1     = fxgrid[(xinx[idn0]+1)]
    fxobs[idn0] = y0 + (y1-y0)*xover[idn0]  # Linearly interpolate
  }
  #fxob   = iflagpn*fxobs
  return(fxobs)
}  # End GetMainfxobs