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

load("data/MC/HB_data.Rdata")

xmin=0 
xmax=10
ngroups = length(unique(id_group))
iflagCenter  = 0 # Used with shape constraints: 1 if int f = 0 and 0 if f(xmin) = 0
iflagZ       = 0 # 0 = Model using Z^2 or 1 = exp(Z)
iflagHBsigma = 1 # 1 for HB sigma and 0 else
iflagACE     = 0 # 1 for Auto Correlate Error model
iflaglm      = 0 # 0 = Normal, 1=Binary Probit, 2=Ordinal DV
iflagpsi     = 1 # 1 = estimate slope of squish function for S and multi-modal)

nbasis=nknots=20
mcmc=list(nblow=40000,smcmc=10000)
shape = "IncreasingConcave"


set.seed(1);
fout_HBBS  = HBBS(ydata=ydata, xdata=xdata, vdata=vdata, wdata=wdata[,-1,drop=F], 
                          zdata=NULL, id_group=id_group, nknots=nknots,
                          xmin=xmin, xmax=xmax, mcmc=mcmc, iflagCenter=iflagCenter,
                          iflagZ=iflagZ, iflagHBsigma=iflagHBsigma, iflagACE=iflagACE, 
                          iflaglm=iflaglm, shape=shape)

set.seed(1); 
fout_HBSAR = HBSAR(ydata=ydata, xdata=xdata, vdata=vdata, wdata=wdata[,-1,drop=F], 
                           zdata=NULL, id_group=id_group, nbasis=nbasis,
                           xmin=xmin, xmax=xmax, mcmc=mcmc, iflagCenter=iflagCenter,
                           iflagZ=iflagZ, iflagHBsigma=iflagHBsigma, iflagACE=iflagACE, 
                           iflaglm=iflaglm, shape=shape)


#save.image("make_figure7.Rdata")

load("data/MC/HB_true_pars.Rdata")
load("data/MC/HB_model_pars.Rdata")

"fig_fall" <- function(fout, main) {
  xgrid   = fout$xgrid
  f0_HBBS = fout$post.est$f0xgridm
  fj_HBBS = fout$post.est$fxgridm
  plot(xgrid, f0_HBBS, type="l", col="red", lwd=2, xlab="x", ylab="f", main=main,
       ylim = range(c(f0_HBBS,fj_HBBS)))
  for (j in 1:ngroups) {
    lines(xgrid, fj_HBBS[,j])
  }
  lines(xgrid, f0_HBBS, col="red", lwd=2)
}

#pdf("Figure7.pdf", width = 15, height = 10)
par(mfrow=c(1,3))
f0t = f0xgridt
fjt = fxgridt
plot(xgrid, f0t, type="l", col="red", lwd=2, xlab="x", ylab="f", main="(a) True",
     ylim = range(c(f0t,fjt)))
for (j in 1:ngroups) {
  lines(xgrid, fjt[,j])
}
lines(xgrid, f0t, col="red", lwd=2)

fig_fall(fout_HBBS,"(b) HB B-Spline")
fig_fall(fout_HBSAR,"(c) HBSAR Cosine")
#dev.off()



