rm(list=ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library("HBSARSCdev")

"intsim" <- function(f,delta) { # integration using Simpson's rule
  r=length(f)
  if (r==2*floor(r/2)) {
    stop('ERROR: Even number of rows for simpson\'s integration')
  } else if(r==3) {	# Simple Simpson's integration
    t=c(1,4,1)
    fint=sum(t*f)*delta/3
  } else {	# Composite Simpson's integration
    t=c(1,rep(c(4,2),times=(r-3)/2),4,1)
    fint=sum(t*f)*delta/3
  }
  fint
}

"calc_RISE" <- function(ft,fhat,xgrid) {
  sse  = (ft - fhat)^2
  RISE = sqrt(intsim(sse, xgrid[2]-xgrid[1]))
  return(RISE)
}

"calc_RISE_j" <- function(ft,fhat,xgrid) {
  ngroup = ncol(ft)
  RISE   = numeric(ngroup)
  for (j in 1:ngroup) {
    sse     = (ft[,j] - fhat[,j])^2
    RISE[j] = sqrt(intsim(sse, xgrid[2]-xgrid[1]))
  }
  
  return(RISE)
}

load("data/Free/HB_data.Rdata")

xmin=0 
xmax=10
ngroups = length(unique(id_group))
iflagCenter  = 1 # Used with shape constraints: 1 if int f = 0 and 0 if f(xmin) = 0
iflagZ       = 0 # 0 = Model using Z^2 or 1 = exp(Z)
iflagHBsigma = 1 # 1 for HB sigma and 0 else
iflagACE     = 1 # 1 for Auto Correlate Error model
iflaglm      = 0 # 0 = Normal, 1=Binary Probit, 2=Ordinal DV
iflagpsi     = 1 # 1 = estimate slope of squish function for S and multi-modal)

nbasis=nknots=20
mcmc=list(nblow=50000,smcmc=10000)
shape = "Free"
set.seed(1)
fout_HBSAR = HBSAR(ydata=ydata, xdata=xdata, vdata=vdata, wdata=wdata[,-1,drop=F], 
                   zdata=NULL, id_group=id_group, nbasis=nbasis,
                   xmin=xmin, xmax=xmax, mcmc=mcmc,
                   iflagZ=iflagZ, iflagHBsigma=iflagHBsigma, iflagACE=iflagACE, 
                   iflaglm=iflaglm, shape=shape)

set.seed(1)
fout_HBBS  = HBBS(ydata=ydata, xdata=xdata, vdata=vdata, wdata=wdata[,-1,drop=F], 
                  zdata=NULL, id_group=id_group, nknots=nknots,
                  xmin=xmin, xmax=xmax, mcmc=mcmc,
                  iflagZ=iflagZ, iflagHBsigma=iflagHBsigma, iflagACE=iflagACE, 
                  iflaglm=iflaglm, shape=shape)


fout_BSAR = fout_BBS = vector("list", length=ngroups)
nbasis2=nknots2=10
for (j in 1:ngroups) {
  bId = id_group==j
  yj  = ydata[bId]
  xj  = xdata[bId]
  vj  = cbind(wdata, vdata)[bId,]
  
  set.seed(1); 
  fout_BSAR[[j]] = BSAR_Free(ydata=yj, xdata=xj, vdata=vj, nbasis=nbasis2, mcmc=mcmc, nint=100, xmin=xmin, xmax=xmax)
  set.seed(1); 
  fout_BBS[[j]] = BBS_Free(ydata=yj, xdata=xj, vdata=vj, nknots=nknots2, mcmc=mcmc, nint=100, xmin=xmin, xmax=xmax)
}


load("data/Free/HB_true_pars.Rdata")
load("data/Free/HB_model_pars.Rdata")

par(mfrow=c(1,3))
f0t = phit[1] + f0xgridt
fjt = matrix(betat[,1],length(xgrid),ngroups,byrow=T)+fxgridt
plot(xgrid, f0t, type="l", col="red", lwd=2, xlab="x", ylab="f+b0", main="(a) True",
     ylim = range(c(f0t,fjt)))
for (j in 1:ngroups) {
  lines(xgrid, fjt[,j])
}
lines(xgrid, f0t, col="red", lwd=2)

f0_HBBS = fout_HBBS$post.est$phim[1]+fout_HBBS$post.est$f0xgridm
fj_HBBS = matrix(fout_HBBS$post.est$betam[,1],length(xgrid),ngroups,byrow=T)+fout_HBBS$post.est$fxgridm
plot(xgrid, f0_HBBS, type="l", col="red", lwd=2, xlab="x", ylab="f+b0", main="(b) HBBS",
     ylim = range(c(f0_HBBS,fj_HBBS)))
for (j in 1:ngroups) {
  lines(xgrid, fj_HBBS[,j])
}
lines(xgrid, f0_HBBS, col="red", lwd=2)

f0_HBSAR = fout_HBSAR$post.est$phim[1]+fout_HBSAR$post.est$f0xgridm
fj_HBSAR = matrix(fout_HBSAR$post.est$betam[,1],length(xgrid),ngroups,byrow=T)+fout_HBSAR$post.est$fxgridm
plot(xgrid, f0_HBSAR, type="l", col="red", lwd=2, xlab="x", ylab="f+b0", main="(c) HBSAR",
     ylim = range(c(f0_HBSAR,fj_HBSAR)))
for (j in 1:ngroups) {
  lines(xgrid, fj_HBSAR[,j])
}
lines(xgrid, f0_HBSAR, col="red", lwd=2)

