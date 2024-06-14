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
Rcpp::List get_mcmc_2DBSAR_Mono(arma::colvec& ydata,
                                arma::mat& vdata,
                                arma::mat& vtv,
                                arma::vec& iflagMainpn,
                                
                                arma::ucube& xinxMain,
                                arma::umat& xinx,
                                arma::mat& xoverMain,
                                arma::vec& xover,
                                arma::vec& xdelta,
                                arma::vec& xrange,
                                
                                arma::colvec& fxdata,
                                arma::mat& fxMaindata,
                                arma::colvec& fx12data,
                                arma::mat& fx12grid_mat,
                                arma::cube& fxMain_mat,
                                arma::mat& f12apx_weights,
                                
                                arma::mat& phixgrid,
                                arma::cube& phiMain,
                                
                                double wk,
                                double wkMain,
                                arma::mat& egamkallMain,
                                arma::colvec& egamkall,
                                arma::colvec& freq_sum,
                                arma::colvec& kall0,
                                
                                double tau2_s0,
                                double tau2_rn,
                                double sigma2_s0,
                                double sigma2_rn,
                                arma::mat& beta_v0i,
                                arma::colvec& beta_v0im0,
                                
                                arma::vec& cpara,
                                Rcpp::List& met_para,
                                
                                arma::colvec& beta,
                                arma::mat& thetaMain,
                                arma::colvec& theta,
                                arma::colvec& xiparMain,
                                double xipar,
                                arma::colvec& gamMain,
                                double gampar,
                                arma::colvec& tau2Main,
                                double tau2,
                                double sigma2,
                                double sigma) {
  
  
  
  unsigned int ntot      = cpara(0);
  unsigned int nparv     = cpara(1);
  unsigned int nblow     = cpara(2);
  unsigned int nblow0    = cpara(3);
  unsigned int maxmodmet = cpara(4);
  unsigned int nskip     = cpara(5);
  unsigned int nmcmcall  = cpara(6);
  unsigned int smcmc     = cpara(7);
  unsigned int nint      = cpara(8);
  unsigned int ntheta    = cpara(9);
  unsigned int nbasis    = cpara(10);
  unsigned int nparx     = cpara(11);
  unsigned int ntheta1D  = cpara(12);
  
  unsigned int ngrid     = (nint+1)*(nint+1);
  
  int iflagpn      = cpara(13); 
  bool bflagCenter = cpara(14);
  bool bflagZ      = cpara(15);
  
  colvec theta_met  = met_para["theta_met"];
  mat thetaMain_met = met_para["thetaMain_met"];
  
  double tau = sqrt(tau2);
  arma::colvec tauMain = sqrt(tau2Main);
  
  mat x1mat = ones<mat>(1,nint+1);
  mat x2mat = ones<mat>(nint+1,1);
  
  mat fxgrid_mat  = zeros<mat>(nint+1,nint+1);
  mat fxgridm_mat = zeros<mat>(nint+1,nint+1);
  mat fxgrids_mat = zeros<mat>(nint+1,nint+1);
  colvec fxdatam = zeros<colvec>(ntot);
  colvec fxdatas = zeros<colvec>(ntot);
  
  // Matrics for saving MCMC iterations
  mat thetag      = zeros<mat>(smcmc,ntheta);
  cube thetaMaing = zeros<cube>(ntheta1D,nparx,smcmc);
  mat betag       = zeros<mat>(smcmc, nparv);
  vec sigmag      = zeros<vec>(smcmc);           // Error variance
  vec taug        = zeros<vec>(smcmc);           // StdDev for spectral coefficients
  mat tauMaing    = zeros<mat>(smcmc,nparx);   // StdDev for spectral coefficients
  vec gammag      = zeros<vec>(smcmc);           // Smoothing parameters for lower level model
  mat gamMaing    = zeros<mat>(smcmc,nparx); // Smoothing parameters for lower level model
  
  
  colvec yhatm = zeros<colvec>(ntot);
  colvec yhats = zeros<colvec>(ntot);
  
  // local
  colvec ym           = zeros<colvec>(ntot);
  colvec yhat         = zeros<colvec>(ntot);
  colvec ymean_new    = zeros<colvec>(ntot);
  colvec ymean_old    = zeros<colvec>(ntot);
  colvec yresid_new   = zeros<colvec>(ntot);
  colvec yresid_old   = zeros<colvec>(ntot);
  colvec fx12data_new = zeros<colvec>(ntot);
  colvec yresid       = zeros<colvec>(ntot);
  colvec vbeta        = zeros<colvec>(ntot);
  colvec theta2 = zeros<colvec>(ntheta);
  colvec bn;
  unsigned int xid;
  mat rMat, vni, vni12, vbi, vbi12, vbni, vbni12;
  Rcpp::List newsmooth;
  vec vec_tmp;
  
  // sigma 
  double sigma2_sn;
  
  // theta
  colvec thv, ths;
  mat ggrid_mat_new    = zeros<mat>(nint+1,nint+1);
  mat fx12grid_mat_new = zeros<mat>(nint+1,nint+1);
  colvec zgrid_new   = zeros<colvec>(ngrid);
  colvec ggrid_new   = zeros<colvec>(ngrid);
  colvec fxMainj_new = zeros<colvec>(nint+1);
  colvec fxMaindata_new = zeros<colvec>(ntot);
  colvec met_var_new, met_var, met_std, theta_old, thetak_old, thetak2_old, t2_old,
  theta_new, thetak_new, thetak2_new, t2_new, t2Main_old, t2Main_new;
  uvec z;
  double testp, ck, xipar_old, xipar2_old, xipar_new, xipar2_new, sse_new, sse_old,
  t0resid2_old;
  unsigned int nz2;
  
  
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
    // *************************************************************************
    // Generate Beta
    // *************************************************************************
    yresid = ydata - fxdata;
    bn     = (vdata.t() * yresid)/sigma2 + beta_v0im0;
    vbi    = vtv/sigma2 + beta_v0i;
    vbi12  = chol(vbi);
    rMat   = randn(nparv, 1);
    beta   = solve(vbi, bn + vbi12.t()*rMat);
    vbeta  = vdata * beta;
    
    // *************************************************************************
    // Generate Main Effects Functions
    // *************************************************************************
    for (unsigned int j=0; j<nparx; j++) {
      //xid = j==0 ? 1 : 0;
      if (j==0) {
        ym = vbeta + fxMaindata.col(1) + fx12data;  // Mean without main effect
      } else {
        ym = vbeta + fxMaindata.col(0) + fx12data;  // Mean without main effect
      }
      //variance of theta
      thv = tau2Main(j)*egamkallMain.col(j);
      ths = sqrt(thv);
      // get variances for t-dist random walk
      met_var_new = AdaptMetVar(thetaMain_met.col(j));
      ck          = 5.66; // Constant from Harito
      met_var     = ck*met_var_new;
      met_std     = sqrt(met_var);
      // Sometimes met_std goes to 0.  Shake it up
      z           = find(met_std < 1e-10);
      nz2         = z.n_elem;
      if (nz2 > 0) {
        colvec sub_met_std = zeros<colvec>(nz2); 
        for (unsigned int jj=0; jj<nz2; jj++) {
          if (randu() < 0.1) {
            sub_met_std(jj) = 1e-10;
          }
        } 
        met_std(z) = sub_met_std;
        met_var(z) = pow(met_std(z), 2);
      }
      theta_old   = thetaMain.col(j);
      // Random walk from for normals
      if (!bflagZ) {
        // With E(Z^2), constraint theta_j0 > 0 with theta_j0 = exp(xi_j)
        xipar_old     = xiparMain(j);                        // theta_j0 = exp(xi_j)
        thetak_old    = theta_old.tail(nbasis); //Skip theta_j0
        xipar2_old    = pow(xipar_old,2);  // squared residual  k=0
        thetak2_old   = square(thetak_old); // squared residuals k>0
        t2Main_old    = zeros<colvec>(ntheta1D);
        t2Main_old(0) = t0resid2_old;
        t2Main_old.tail(nbasis) = thetak2_old;
        
        // Generate theta_new for Z^2
        xipar_new     = xipar_old + met_std(0)*ths(0)*randn();
        thetak_new    = thetak_old + met_std.tail(nbasis)%ths.tail(nbasis)%randn(nbasis,1);
        theta_new     = zeros<colvec>(ntheta1D);
        theta_new(0)  = exp(xipar_new);
        theta_new.tail(nbasis) = thetak_new;
        xipar2_new    = pow(xipar_new,2);
        thetak2_new   = square(thetak_new);
        t2Main_new    = zeros<colvec>(ntheta1D);
        t2Main_new(0) = xipar2_new;
        t2Main_new.tail(nbasis) = thetak2_new;
        
      } else {
        // Generate theta_new for exp(z)
        thetak_new = theta_old + met_std%ths%randn(ntheta1D,1);
        t2Main_old = square(theta_old);
        t2Main_new = square(theta_new);
      } // End z^2 or exp(z) for generating theta
      
      // Compute f at new theta
      zgrid_new   = phiMain.slice(j) * theta_new;
      if (!bflagZ) {
        ggrid_new = square(zgrid_new);
      } else {
        ggrid_new = exp(zgrid_new);
      }
      fxMainj_new = iflagMainpn(j)*rcpp_CumTrap(ggrid_new,xdelta(j));
      // Is f center?
      if (bflagCenter) {
        fxMainj_new  -= rcpp_trapint(fxMainj_new,xdelta(j))/xrange(j);
      }
      vec_tmp = xoverMain.col(j);
      fxMaindata_new = rcpp_GetSCfxobs(xinxMain.slice(j),vec_tmp,fxMainj_new);
      // Get metropolis test probabilities
      ymean_new      = ym + fxMaindata_new;
      ymean_old      = ym + fxMaindata.col(j);
      yresid_new     = ydata - ymean_new;
      yresid_old     = ydata - ymean_old;
      sse_new        = accu(square(yresid_new));
      sse_old        = accu(square(yresid_old));
      testp          = -sse_new/(2*sigma2) + sse_old/(2*sigma2);
      testp          = testp - accu(t2Main_new/thv)/2 + accu(t2Main_old/thv)/2; // Priors for theta
      
      // Accept new draw of theta
      if (log(randu()) < testp) { // Accept candidate
        if (!bflagZ) xiparMain(j) = xipar_new;
        thetaMain.col(j)      = theta_new;
        fxMaindata.col(j)     = fxMaindata_new;
        if (j == 0) {
          fxMain_mat.slice(j) = fxMainj_new * x1mat;
        } else {
          fxMain_mat.slice(j) = x2mat * fxMainj_new.t();
        }
        
        thetaMain_met.col(j).tail(thetaMain_met(0,j)) = met_var_new;
        thetaMain_met(1,j)++; // pmet
      } // End MH
      fxgrid_mat = fxMain_mat.slice(0) + fxMain_mat.slice(1) + fx12grid_mat;
      fxdata     = fxMaindata.col(0) + fxMaindata.col(1) + fx12data;
      // Update adaptive metropolis mean and beta
      thetaMain_met.col(j) = AdaptMetUpdate2(thetaMain_met.col(j));
    } // End Loop to generate Main Effects
    
    
    // *************************************************************************
    // Generate Theta for interactions f12(x1,x2)
    // *************************************************************************
    ym  = vbeta + fxMaindata.col(0) + fxMaindata.col(1); // mean without f12(x1,x2)
    // HB variance of theta
    thv = tau2*egamkall;
    ths = sqrt(thv);
    // get variances for t-dist random walk
    met_var_new = AdaptMetVar(theta_met);
    ck          = 5.66; // Constant from Harito
    met_var     = ck*met_var_new;
    met_std     = sqrt(met_var);
    // Sometimes met_std goes to 0.  Shake it up
    z           = find(met_std < 1e-10);
    nz2         = z.n_elem;
    if (nz2 > 0) {
      colvec sub_met_std = zeros<colvec>(nz2); 
      for (unsigned int jj=0; jj<nz2; jj++) {
        if (randu() < 0.1) {
          sub_met_std(jj) = 1e-10;
        }
      } 
      met_std(z) = sub_met_std;
      met_var(z) = pow(met_std(z), 2);
    }
    theta_old   = theta;
    // Random walk from for normals
    if (!bflagZ) {
      // With E(Z^2), constraint theta_j0 > 0 with theta_j0 = exp(xi_j)
      xipar_old    = xipar;                        // theta_j0 = exp(xi_j)
      thetak_old   = theta_old.tail(ntheta-1); //Skip theta_j0
      xipar2_old   = pow(xipar_old,2);  // squared residual  k=0
      thetak2_old  = square(thetak_old); // squared residuals k>0
      t2_old       = zeros<colvec>(ntheta);
      t2_old(0)    = t0resid2_old;
      t2_old.tail(ntheta-1) = thetak2_old;
      
      // Generate theta_new for Z^2
      xipar_new    = xipar_old + met_std(0)*ths(0)*randn(); 
      thetak_new   = thetak_old + met_std.tail(ntheta-1)%ths.tail(ntheta-1)%randn(ntheta-1,1);
      theta_new    = zeros<colvec>(ntheta);
      theta_new(0) = exp(xipar_new);
      theta_new.tail(thetak_new.n_elem) = thetak_new;
      xipar2_new   = pow(xipar_new,2);
      thetak2_new  = square(thetak_new);
      t2_new       = zeros<colvec>(ntheta);
      t2_new(0)    = xipar2_new;
      t2_new.tail(ntheta-1) = thetak2_new;
    } else {
      // Generate theta_new for exp(z)
      thetak_new = theta_old + met_std%ths%randn(ntheta,1);
      t2_old     = square(theta_old);
      t2_new     = square(theta_new);
    } // End z^2 or exp(z) for generating theta
    
    // Compute f at new theta
    zgrid_new   = phixgrid * theta_new;
    if (!bflagZ) {
      ggrid_new = square(zgrid_new);
    } else {
      ggrid_new = exp(zgrid_new);
      
    }
    ggrid_mat_new    = reshape(ggrid_new, nint+1,nint+1);
    fx12grid_mat_new = iflagpn*rcpp_CumInt2D(ggrid_mat_new,xdelta);
    // Is f center?
    if (bflagCenter) {
      fx12grid_mat_new  -= rcpp_Int2D(fx12grid_mat_new,xdelta)/(xrange(0)*xrange(1));
    }
    fx12data_new = rcpp_fdata_apx(fx12grid_mat_new,f12apx_weights,xinx,xover,nint);
    // Get metropolis test probabilities
    ymean_new      = ym + fx12data_new;
    ymean_old      = ym + fx12data;
    yresid_new     = ydata - ymean_new;
    yresid_old     = ydata - ymean_old;
    sse_new        = accu(square(yresid_new));
    sse_old        = accu(square(yresid_old));
    testp          = -sse_new/(2*sigma2) + sse_old/(2*sigma2);
    testp          = testp - accu(t2_new/thv)/2 + accu(t2_old/thv)/2; // Priors for theta
    
    // Accept new draw of theta
    if (log(randu()) < testp) { // Accept candidate
      if (!bflagZ) xipar = xipar_new;
      theta        = theta_new;
      fx12data     = fx12data_new;
      fx12grid_mat = fx12grid_mat_new;
      theta_met.tail(theta_met(0)) = met_var_new;
      theta_met(1)++; // pmet
      fxgrid_mat   = fxMain_mat.slice(0) + fxMain_mat.slice(1) + fx12grid_mat;
      fxdata       = fxMaindata.col(0) + fxMaindata.col(1) + fx12data;
    } // End MH
    // Update adaptive metropolis mean and beta
    theta_met = AdaptMetUpdate2(theta_met);
    
    // Generate sigma2
    yresid    = ydata - vbeta - fxdata;
    sigma2_sn = sigma2_s0 + sum(square(yresid));
    sigma2    = sigma2_sn/(2 * randg(distr_param(sigma2_rn / 2, 1.0)));
    sigma     = sqrt(sigma2);
    
    // *************************************************************************
    // Generate Smoothing Parameter tau2 and gamma for Main functions
    // *************************************************************************
    for (unsigned int j=0; j<nparx; j++) {
      newsmooth = rcpp_getsmoothSC(tau2Main(j),gamMain(j), thetaMain.col(j), square(thetaMain.col(j)),
                                   kall0, egamkallMain.col(j),wkMain,tau2_s0,tau2_rn,bflagZ);
      tau2Main(j)         = newsmooth["tau2"];
      gamMain(j)          = newsmooth("gampar");
      colvec tmp1         = newsmooth("egamkall");
      egamkallMain.col(j) = tmp1;
    }
    tauMain = sqrt(tau2Main);
    // *************************************************************************
    // Generate Smoothing Parameter tau2 and gamma for interaction function
    // *************************************************************************
    newsmooth = rcpp_getsmoothSC(tau2,gampar,theta,theta2,freq_sum,egamkall,wk,
                                 tau2_s0,tau2_rn,bflagZ);
    tau2        = newsmooth["tau2"];
    tau         = sqrt(tau2);
    gampar      = newsmooth("gampar");
    colvec tmp2 = newsmooth("egamkall");
    egamkall    = tmp2;
    
    // Save MCMC iterations
    if ((imcmc >= nblow+nblow0*maxmodmet) && (imcmc % nskip == 0)) {
      yhat         = vbeta + fxdata;
      
      yhatm       += yhat;
      yhats       += square(yhat);
      fxdatam     += fxdata;
      fxdatas     += square(fxdata);
      fxgridm_mat += fxgrid_mat;
      fxgrids_mat += square(fxgrid_mat);
      
      betag.row(isave)  = beta.t();
      thetag.row(isave) = theta.t();
      thetaMaing.slice(isave) = thetaMain;
      
      sigmag(isave)       = sigma;
      taug(isave)         = tau;
      gammag(isave)       = gampar;
      tauMaing.row(isave) = tauMain.t();
      gamMaing.row(isave) = gamMain.t();
      
      isave++;
    } // End save MCMC
  } // end of mcmc loop
  
  cout << "MCMC is done!" << endl;
  
  Rcpp::List metg = Rcpp::List::create(Rcpp::Named("theta_met")     = theta_met,
                                       Rcpp::Named("thetaMain_met") = thetaMain_met);
  
  
  Rcpp::List mcmcg = Rcpp::List::create(Rcpp::Named("betag")      = betag,
                                        Rcpp::Named("thetag")     = thetag,
                                        Rcpp::Named("thetaMaing") = thetaMaing,
                                        Rcpp::Named("sigmag")     = sigmag,
                                        Rcpp::Named("taug")       = taug,
                                        Rcpp::Named("tauMaing")   = tauMaing,
                                        Rcpp::Named("gammag")     = gammag,
                                        Rcpp::Named("gamMaing")   = gamMaing);
  
  // maximum element of list is 20.
  return Rcpp::List::create(Rcpp::Named("fxgridm_mat") = fxgridm_mat,  //  1
                            Rcpp::Named("fxgrids_mat") = fxgrids_mat,  //  2
                            Rcpp::Named("yhatm")       = yhatm,        //  3
                            Rcpp::Named("yhats")       = yhats,        //  4
                            Rcpp::Named("fxdatam")     = fxdatam,      //  5
                            Rcpp::Named("fxdatas")     = fxdatas,      //  6
                            Rcpp::Named("mcmcg")       = mcmcg,        //  7
                            Rcpp::Named("metg")        = metg);        //  8
  
}
