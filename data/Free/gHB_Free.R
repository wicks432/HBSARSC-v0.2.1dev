#-------------------------------------------------------------------------------
#   Generate HB Data
#   Search "Getf" for function specification
#   Lower-level model
#      Observation i from group j
#         Y_ij       = v_ij*alpha + w_ij*beta_j + f_j(x_ij) + epslion_ij
#         epsilon_ij ~ N(0,sigma^2) or N(0,sigma_j^2)
#         alpha is fixed effect, beta_j is random effect
#   Upper-level 
#      HB Model for beta_j
#         beta_j     = Phi'z_j + delta_j  were delta_j ~ N(0,Lambda)
#      HB Error variance 
#         sigma_j^2 ~ IG(sigma2_alpha/2,sigma2_beta/2)
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
outfile      ="HB"

iflagHBsigma = 1    # 1 for HB model.
iflagACE     = 0    # 1 for ACE errors
iflagCenter  = 1    # 1 for int f = 0 and 0 for f(0) = 0

ntot         = 2000  # Total number of observations
ngroups      = 100   # Number of groups


nparv       = 1      # Number of v variables: nparv = 0 if no fixed effects
nparw       = 2      # number of W variables, including the intercept
nparz       = 1      # number of z varibables, including constant

xmin        = 0             # Minimum value for x
xmax        = 10            # Maximum value for x
nint        = 100           # Number of grid points


# Generate error variance sigma: homogeneous or heterogeneous model
sigmat  = 1
if(iflagHBsigma == 1){
  # HB error variance
  # sigma_j^2 is IG(alpha/2,beta/2)
  sigma2_mu     = sigmat^2           # Mean of sigma2
  sigma2_v2     = sigma2_mu/4        # Variance of sigma2
  sigma2_alphat = 2*(sigma2_mu^2/sigma2_v2 + 1)
  sigma2_betat  = 2*sigma2_mu*(sigma2_alphat/2 -1 )
  sigma2t       = sigma2_betat/(2*rgamma(ngroups,shape=sigma2_alphat/2,rate=1))
  sigmat        = sqrt(sigma2t)     # error std dev
}

# tau_k^2 ~ IG(tau2_r/2, tau2_s/2)
tau2_mu    = 9     # Mean
tau2_sd    = 2     # Standard deviation
tau2_r     = 2*(2+(tau2_mu/tau2_sd)^2)
tau2_s     = (tau2_r-4)*tau2_mu

#---------------------------------------------------
# MakeHBPointer
#    Make a pointer into the matrics for groups
#    Input
#        id   = identifies groups
#    Output
#        nobs = number of observations in each group
#        iptd = list of pointers into group
#        id_beta = indices to expaned beta to match wdata
#    Useage
#        xj   = xdata[iptd[[j]],]  gives you the x data for group j
#        Compute wdata_j*beta_j for each group
#        wdata*beta = matrix(apply(wdata*bpar[id_beta,],1,sum))
#    Call
#        outpt = MakeHBPointer(id)
#---------------------------------------------------
MakeHBPointer <- function(id){
  uid = sort(unique(id))   # unique populations
  ng  = length(uid)        # number of unique groups
  
  # id_group may not be sequential from 1 ... ngroups
  # Need an index to get correct betam
  
  nobs = matrix(0,ng,1)          # number of observations in each group
  id_beta  = matrix(0,ng,1)      # match rows of beta to wdata
  for(j in 1:ng){                # loop over groups
    jrow = which(id == uid[j])   # identify rows for group j
    nj   = length(jrow)          # Number of observations in group j
    id_beta[jrow] = matrix(j,nj,1)
    
    nobs[j] = length(jrow)       # number of observations in group j
    if(j == 1){
      iptd = list(jrow)          # start the list iptd
    }else{
      iptd[[j]] = jrow           # add new member to list
    }
  }  # end loop over groups
  return(list("nobs"=nobs,"iptd"=iptd,"id_beta"=id_beta))
}  # End MakeHBPointer
#---------------------------------------------------

#-----------------------------------------------------
# TrapInt does trapidzoid integration
# Int f(x)*dx = xdelta*(f1/2 + f2 + ... + fn-1 + fn/2)
#-----------------------------------------------------
trapint <- function(f,xdelta){
  nr = length(f)
  a  = sum(f)
  a  = a - (f[1] + f[nr])/2
  return(a*xdelta)
}


#----------------------------------------------------------
# CumTrap computes cumulative integral of function
# int_a^x f(s)ds over a grid of points a, a+delta, a+2*delta, ..., b
# INPUT
#   f     = function on grid
#   delta = size of grid interval
# OUTPUT
#   fcum  = cumulative function corresonpding to f
#----------------------------------------------------------
CumTrap = function(f,delta){
  nr    = nrow(f)
  fcum	= cumsum(f[1:nr-1]) + cumsum(f[2:nr])
  fcum	= fcum*delta/2
  fcum	= matrix(c(0,fcum),nr,1)
  return(fcum)
}  # End Cumtrap 
#----------------------------------------------------------

# Make groups for HB
nobs     = floor(ntot/ngroups)
ntot     = nobs*ngroups
id_group = matrix(1:ngroups) %x% matrix(1,nobs,1)

uid         = sort(unique(id_group))
if(ngroups != length(uid)) print("Something wrong with HB groups")
# Lists to extract group data
idout       = MakeHBPointer(id_group)   # Pointers into groups
nobs        = idout$nobs                # Vector: number of observations in each group
iptHB       = idout$iptd                # xj  = xdata[[iptHB[j]]] = observations in group j
id_beta     = idout$id_beta             # wdata*beta = matrix(apply(wdata*betat[id_beta,],1,sum))

# Got fixed effects
if(nparv>0){
  # Generate v data
  vdata            = matrix(rnorm(ntot*nparv),ntot,nparv)
  vnames           = paste("V",1:nparv,sep="")
  colnames(vdata)  = vnames
  alphat           = matrix(floor(2+.5*rnorm(nparv)*1000)/1000)
  alphat[1]        = -1
  rownames(alphat) = vnames
  valphat          = vdata %*% alphat
}else{
  valphat          = NULL
  alphat           = NULL
  vdata            = NULL
  vnames           = "Error!"
}


# Generate Z data
zdata           = matrix(1,ngroups,nparz)
znames          = "CNST"
if(nparz>1){
  zdata[,2:nparz] = rnorm(ngroups*(nparz-1))
  zk              = paste("Z",1:(nparz-1),sep="")
  znames          = c(znames,zk)
}


colnames(zdata) = znames

# Generate Wdata 
wdata    = matrix(1,ntot,nparw)
if(nparw>1){
  for(j in 1:(nparw-1)){
    wdata[,(j+1)] = rnorm(ntot)
    wnames    = paste("W",1:(nparw-1),sep="")
    wnames    = c("CNST",wnames)
    colnames(wdata) = wnames
  }
}else{
  wnames           = "CNST"
  colnames(wdata) = wnames
}


#-----------------------------------------------------------
# Upper level model for random effects beta_j
# beta_j = phi'z_j + delta_j and delta_j ~ N(0,lambda)
phit = 2 + .5*rnorm(nparw*nparz)
phit = floor(phit*100)/100
phit = matrix(phit,nparz,nparw)
phit[1,1]  = 10
phit[1,2]  = 2

colnames(phit)  = wnames
row.names(phit) = znames
# Generate error variance
lamb0      = .5    # base variance
if(nparw==1){
  lambdat   = lamb0
  lambdat12 = sqrt(lambdat)
}else{
  v        = (lamb0 + (lamb0/2)*runif(nparw))
  lambdat  = diag(v)
  lambdat12 = chol(lambdat)  # Upper diagonal: t(lambdat12)  %*% lambdat12 = lambdat
  colnames(lambdat) = wnames
}


# Generate beta
berrors = matrix(rnorm(ngroups*nparw),ngroups,nparw) %*% lambdat12
betat   = zdata%*%phit + berrors
colnames(betat) = wnames

# ACE correlations
rhoallt = matrix(0,ngroups,1)
if(iflagACE == 1)  rhoallt =  matrix(rbeta(ngroups,30,20))  # Correlation parameters


# Generate X data on [xmin, xmax]
xdata     = matrix(runif(ntot,xmin,xmax))

#-----------------------------------------------------
# Grid on [xmim,xmax] for plotting functions
# Also, grid for numerical integration
xdelta = (xmax-xmin)/nint
xgrid  = seq(xmin,xmax,xdelta)
xgrid  = matrix(xgrid,nint+1,1)
xrange = xmax-xmin
xmid   = (xmin+xmax)/2
zgrid  = (xgrid-xmid)/xrange  # Standardized

#------------------------------------------------------
# Locate indices that define data (x1,x2)
#   Instead of computing f(x_i) at observations,
#   find where x_i falls in xgrid and use fxgrid
#   Find location of xdata in xgrid
#   xgrid[xinx[i,1]] <= xi <= xgrid[xinx[i,2]]
#------------------------------------------------------
xinx  =  matrix(0,ntot,2)
xover  = matrix(0,ntot,1)   # Relative excess over grid
# Loop over observations
for(i in 1:ntot){  
  xi   = xdata[i]
  idx   = which(xi == xgrid)  # Is it exactly on xgrid?
  if(length(idx) > 0 ){
    xinx[i,1]  = idx
    xinx[i,2]  = idx
    xover[i]   = 0 
  }else{
    idx    = max(which(xi > xgrid))
    xinx[i,1]  = idx
    xinx[i,2]  = idx+1
    xover[i]   = (xi - xgrid[idx])/(xdelta)
  }
} # End loop over observations for xinx and xover

#--------------------------------------------
# Compute true f. Comment or uncomment lines
#--------------------------------------------
Getf = function(x,fpars){
  
  #----------------------------------------------------
  # Multimodal: linear + normal + linear*normal
  #----------------------------------------------------
  f = fpars$b1 *(x) +
    fpars$b2*dnorm(x,fpars$b3,fpars$b4) +
    fpars$b5*(x-fpars$b8)*dnorm(x,fpars$b6,fpars$b7)
  
    
  return(f)
}
b10 = 1.2
b20 = 10
b30 = 3
b40 = 1
#-----------------------
b50 = 30
b60 = 6
b70 = 1
b80 = b60-b70/4
fpars    = list(b1=b10,b2=b20,b3=b30,b4=b40,b5=b50,b6=b60,b7=b70,b8=b80)
f0xgridt = Getf(xgrid,fpars)
# Make \int f = 0
if(iflagCenter==1){
  c = trapint(f0xgridt,xdelta)
  f0xgridt = f0xgridt - c/xrange
}


#----------------------------------------------
# Generate group-level functions
#----------------------------------------------
fxgridt     = matrix(0,nint+1,ngroups)
fxgridb0t   = fxgridt
fxdatat     = matrix(0,ntot,1)
for(j in 1:ngroups){
  xj = xdata[iptHB[[j]]]       # xdata for group j
  # Perturb coefficients b0 with random errros
  b1 = rnorm(1,b10,.1)
  b2 = rnorm(1,b20,.5)
  b3 = rnorm(1,b30,.25)
  b4 = runif(1,.8*b40,1.5*b40)
  #-----------------------
  b5 = rnorm(1,b50,2)
  b6 = rnorm(1,b60,.25)
  b7 = runif(1,.8*b70,1.5*b70)
  b8 = b6-b7/4
  fpars = list(b1=b1,b2=b2,b3=b3,b4=b4,b5=b5,b6=b6,b7=b7,b8=b8)
  fxj   = Getf(xgrid,fpars)    # Function on grid
  fxjd  = Getf(xj,fpars)       # Function on xdata for group j
  # Make \int f = 0
  if(iflagCenter==1){
    c    = trapint(fxj,xdelta)
    fxj  = fxj  - c/xrange
    fxjd = fxjd - c/xrange
  }
  fxgridt[,j]         = fxj
  fxdatat[iptHB[[j]]] = fxjd
}  # End get fj for each group

#---------------------------------------------------------
# Reset f0 to mean of fj
f0xgridt = apply(fxgridt,1,mean)   

#----------------------------------------------
# Compute Wi*beta_i
wbetat = matrix(apply(wdata*betat[id_beta,],1,sum)) 
# Get errors
rerr   = matrix(0,ntot,1)      # Random error for each group
for(j in 1:ngroups){ # Loop over populations
  nj  = nobs[j]
  if(iflagHBsigma == 0){
    sigj = sigmat
  }else{
    sigj = sigmat[j]
  }
  # Normal errors
  e0j = sigj*rnorm(nj)
  errj = e0j
  if(iflagACE == 1){
    rhoj     = rhoallt[j]
    errj[1]  = e0j[1]/sqrt(1-rhoj^2)    # Stationary starting value
    for(k in 2:nj) errj[k] = e0j[k] + rhoj*errj[k-1] # AR errors
  }
  rerr[iptHB[[j]]] = errj
} # End generate errors

#-------------------------------------------------
# Generate Y
ym     = valphat + wbetat + fxdatat  # mean of Y
ydata  = ym +  rerr 
yresid = ydata - valphat - wbetat   # Parametric residuals for plotting

# Plot curves and residuals
ymin  = floor(min(c(f0xgridt,fxgridt,yresid)))
ymax  = ceiling(max(c(f0xgridt,fxgridt,yresid)))
matplot(xgrid,cbind(f0xgridt,fxgridt),type="l",ylim=c(ymin,ymax),
        xlab="X",ylab="f")
points(xdata,yresid,pch=16)


boxplot(betat,main="Random Coefficients")
abline(h=0)
if(iflagHBsigma == 1) hist(sigmat,main="HB SD sigma")
if(iflagACE == 1)  hist(rhoallt,main="Auto-correlation rho")


colnames(ydata) = "Y"
colnames(xdata) = "X"
colnames(wdata) = wnames
colnames(zdata) = znames




# Save model information used in estimation
save(
  iflagHBsigma, # 1 for HB sigma2 and 0 for homo sigma2
  iflagCenter,  # 1 for center f
  iflagACE,     # 1 for ACE errors
  nparv,        # number of fixed coefficients
  nparw,        # number of random coefficients
  nparz,        # number of upper level covariates
  xmin,         # minimum of range for X
  xmax,         # maximum of range for X
  xrange,       # range of the data
  xdelta,       # dx
  nint,         # number of grid points
  xgrid,        # grid for plotting f(x)
  file=paste(outfile,"_model_pars.Rdata",sep="")
)

# Save model parameters
save(
  fxgridt,
  f0xgridt,
  alphat,
  betat,
  sigmat,
  rhoallt,
  phit,
  lambdat,
  file=paste(outfile,"_true_pars.Rdata",sep="")
)
# save data
save(ydata,xdata,vdata,wdata,zdata,id_group,
     file=paste(outfile,"_data.Rdata",sep=""))
# NOTE: wdata and zdata includes a column of ones

