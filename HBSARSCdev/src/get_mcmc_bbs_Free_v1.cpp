// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include "misc.h"

#define use_progress

#ifdef use_progress
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif


using namespace arma;

// [[Rcpp::export]]
Rcpp::List get_mcmc_BBS_Free(arma::colvec& ydata,
                             arma::mat& vdata,
                             arma::mat& vtv,
                             arma::mat& ptp,
                             
                             arma::colvec& fxdata,
                             
                             arma::mat& phixgrid,
                             arma::mat& phixdata,
                             
                             arma::colvec& kall,
                             
                             double tau2_s0,
                             double tau2_rn,
                             double sigma2_s0,
                             double sigma2_rn,
                             arma::mat& beta_v0i,
                             arma::colvec& beta_v0im0,
                             
                             arma::vec& cpara,
                             
                             arma::colvec& beta,
                             arma::mat& theta,
                             double tau2,
                             double tau,
                             double sigma2,
                             double sigma) {
  
  
  unsigned int ntot      = cpara(0);
  unsigned int nparv     = cpara(1);
  unsigned int nblow     = cpara(2);
  unsigned int nmcmcall  = cpara(3);
  unsigned int smcmc     = cpara(4);
  unsigned int nint      = cpara(5);
  unsigned int ntheta    = cpara(6);
  unsigned int ngrid     = nint+1;
  
  colvec fxgrid  = zeros<colvec>(ngrid);
  colvec fxgridm = zeros<colvec>(ngrid);
  colvec fxgrids = zeros<colvec>(ngrid);
  colvec fxdatam = zeros<colvec>(ntot);
  colvec fxdatas = zeros<colvec>(ntot);
  
  // Matrics for saving MCMC iterations
  mat thetag  = zeros<mat>(smcmc,ntheta);
  mat betag   = zeros<mat>(smcmc, nparv);
  vec sigmag  = zeros<mat>(smcmc);    // Error variance
  vec taug    = zeros<vec>(smcmc);    // StdDev for spectral coefficients
  mat fxgridg = zeros<mat>(smcmc, ngrid);
  
  colvec yhat   = zeros<colvec>(ntot);
  colvec yhatm  = zeros<colvec>(ntot);
  colvec yhats  = zeros<colvec>(ntot);
  
  // local
  colvec yresid = zeros<colvec>(ntot);
  colvec vbeta  = zeros<colvec>(ntot);
  colvec theta2 = zeros<colvec>(ntheta);
  colvec bn;
  mat rMat, vni, vni12, vbi, vbi12, vbni, vbni12;
  
  // sigma 
  double sigma2_sn;
  
  // theta
  double thvar, thvari;
  
  // tau
  double tau2_sn;
  
#ifdef use_progress
  Progress p(nmcmcall, true);
#endif
  
  unsigned int isave = 0;
  for (unsigned int imcmc=0; imcmc<nmcmcall; imcmc++) {
    
#ifdef use_progress
    if ( Progress::check_abort() ) { // check for user abort
      break;
    }
    p.increment();
#else
    R_CheckUserInterrupt(); // void return type (not bool!!)
#endif
    // Generate Beta
    yresid = ydata - fxdata;
    bn     = (vdata.t() * yresid)/sigma2 + beta_v0im0;
    vbi    = vtv/sigma2 + beta_v0i;
    vbi12  = chol(vbi);
    rMat   = randn(nparv, 1);
    beta   = solve(vbi, bn + vbi12.t()*rMat);
    vbeta  = vdata * beta;
    
    // Generate Theta
    yresid      = ydata - vbeta;
    thvar       = tau2;
    thvari      = 1/thvar;
    bn          = (phixdata.t() * yresid)/sigma2;
    vni         = ptp/sigma2;
    vni.diag() += thvari;
    vni12       = chol(vni);
    rMat        = randn(ntheta, 1);
    theta       = solve(vni, bn + vni12.t()*rMat);
    theta2      = square(theta);
    fxdata      = phixdata * theta;
    
    // Generate sigma2
    yresid    = ydata - vbeta - fxdata;
    sigma2_sn = sigma2_s0 + sum(square(yresid));
    sigma2    = sigma2_sn/(2 * randg(distr_param(sigma2_rn / 2, 1.0)));
    sigma     = sqrt(sigma2);
    
    // Generate Smoothing Parameter tau2
    tau2_sn = tau2_s0 + accu(theta2); // Posterior rate parameters
    tau2    = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
    tau     = sqrt(tau2);
    
    // Save MCMC iterations
    if (imcmc >= nblow) {
      yhat               = vbeta + fxdata;
      
      yhatm             += yhat;
      yhats             += square(yhat);
      fxdatam           += fxdata;
      fxdatas           += square(fxdata);
      fxgrid             = phixgrid*theta;
      fxgridm           += fxgrid;
      fxgrids           += square(fxgrid);
      fxgridg.row(isave) = fxgrid.t();
      
      betag.row(isave)  = beta.t();
      thetag.row(isave) = theta.t();
      sigmag(isave)     = sigma;
      taug(isave)       = tau;
      
      isave++;
    } // End save MCMC
  } // end of mcmc loop
  
  cout << "MCMC is done!" << endl;
  
  Rcpp::List mcmcg = Rcpp::List::create(Rcpp::Named("betag")   = betag,
                                        Rcpp::Named("thetag")  = thetag,
                                        Rcpp::Named("sigmag")  = sigmag,
                                        Rcpp::Named("taug")    = taug,
                                        Rcpp::Named("fxgridg") = fxgridg);
  
  // maximum element of list is 20.
  return Rcpp::List::create(Rcpp::Named("fxgridm")  = fxgridm,        //  1
                            Rcpp::Named("fxgrids")  = fxgrids,        //  2
                            Rcpp::Named("fxdatam")  = fxdatam,        //  3
                            Rcpp::Named("fxdatas")  = fxdatas,        //  4
                            Rcpp::Named("mcmcg")    = mcmcg);         //  5
  
}
