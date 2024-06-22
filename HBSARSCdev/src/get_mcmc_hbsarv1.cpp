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
Rcpp::List get_mcmc_HBSAR(Rcpp::List& data_all,
                          arma::colvec& dpara,
                          arma::colvec& xgrid,
                          arma::umat& xinx,
                          arma::vec& xover,
                          arma::vec& nobs,
                          Rcpp::List& iptHB1,
                          arma::mat& ztz,
                          arma::cube& vtv,
                          arma::cube& wtw,
                          arma::cube& phi2,
                          
                          arma::vec& cpara,
                          Rcpp::List& ace_para,
                          Rcpp::List& probit_para,
                          Rcpp::List& ordinal_para,
                          Rcpp::List& gamma_para,
                          Rcpp::List& tau_para,
                          Rcpp::List& eta_para,
                          Rcpp::List& met_para,
                          Rcpp::List& squish_para1,
                          Rcpp::List& squish_para2,
                          arma::vec& sigma_para,
                          Rcpp::List& v_para,
                          Rcpp::List& z_para,
                          Rcpp::List& phi_para,
                          Rcpp::List& fx_para,
                          Rcpp::List& basis_para,
                          Rcpp::List& slope_para,
                          
                          arma::colvec& alpha,
                          arma::mat& betall,
                          arma::mat& phi,
                          arma::mat& lambda,
                          arma::mat& lambdai,
                          arma::mat& thetall,
                          arma::colvec& theta0,
                          arma::colvec& xiparall,
                          arma::colvec& gammall,
                          arma::colvec& tau2all,
                          arma::colvec& tauall,
                          arma::colvec& sigma2,
                          arma::colvec& sigma,
                          
                          arma::cube& betallg,
                          arma::mat& rhoallg,
                          
                          arma::cube& thetallg,
                          
                          arma::cube& zetag,
                          arma::cube& omegag,
                          arma::mat& zeta0g, 
                          arma::mat& omega0g,
                          arma::mat& psiallg,
                          arma::colvec& psi0g,
                          
                          arma::mat& mslope_logg,
                          arma::colvec& mslope0_logg,      
                          arma::mat& mslopeg,     
                          arma::colvec& mslope0g,
                          arma::colvec& mslope_log_sigmag,
                          
                          arma::cube& fxgridallg,
                          arma::mat& f0xgridg) {
  

  
  unsigned int ntot      = cpara(0);
  unsigned int ngroup    = cpara(1);
  unsigned int nparv     = cpara(2);
  unsigned int nparw     = cpara(3);
  unsigned int nparw2    = cpara(4);
  unsigned int nparz     = cpara(5);
  unsigned int phidim    = cpara(6);
  unsigned int nbasis    = cpara(7);
  unsigned int nblow     = cpara(8);
  unsigned int nblow0    = cpara(9);
  unsigned int maxmodmet = cpara(10);
  unsigned int nskip     = cpara(11);
  unsigned int nmcmcall  = cpara(12);
  unsigned int smcmc     = cpara(13);
  unsigned int nint      = cpara(14);
  unsigned int ntheta    = cpara(15);
  unsigned int nstat     = cpara(16);
  unsigned int nflex     = cpara(17); 
  unsigned int nomega    = cpara(18);
  
  bool bflagHBsigma = cpara(19);
  bool bflagpsi     = cpara(20);  
  bool bflagCenter  = cpara(21);
  bool bflagACE     = cpara(22);
  bool bflagZ       = cpara(23);
  
  int iflaglm          = cpara(24);
  int iflagpn          = cpara(25);   
  int iflagsc          = cpara(26);
  int iflagZcnst       = cpara(27);
  bool save_large_mcmc = cpara(28);
  
  double xdelta = dpara(0);
  double xrange = dpara(1);
  double xmin   = dpara(2);  
  double xmax   = dpara(3);
  double xipar0 = dpara(4);
  
  int ingroup = ngroup;
  
  colvec ydata  = data_all["ydata"];
  colvec yodata = data_all["yodata"];
  colvec valpha = data_all["valpha"];
  colvec wbeta  = data_all["wbeta"];
  mat vdata     = data_all["vdata"];
  mat wdata     = data_all["wdata"];
  mat zdata     = data_all["zdata"];
  
  
  colvec rhoall    = ace_para["rhoall"];
  double rho_v0i   = ace_para["rho_v0i"];
  double rho_v0im0 = ace_para["rho_v0im0"];
  
  uvec id_probit1 = probit_para["id_probit1"];
  uvec id_probit0 = probit_para["id_probit0"];
  
  uvec id_cut_free = ordinal_para["id_cut_free"];
  uvec id_cut_bot  = ordinal_para["id_cut_bot"];
  uvec id_cut_top  = ordinal_para["id_cut_top"];
  colvec ybot      = ordinal_para["ybot"]; 
  colvec ytop      = ordinal_para["ytop"];
  mat cuttall      = ordinal_para["cuttall"];
  mat cuttallm     = ordinal_para["cuttallm"];
  mat cuttalls     = ordinal_para["cuttalls"];
  
  mat theta_met      = met_para["theta_met"];
  colvec gamma_met   = met_para["gamma_met"];
  mat psi_met        = met_para["psi_met"];
  mat zeta_met       = met_para["zeta_met"];
  colvec sigma2_met  = met_para["sigma2_met"];
  mat mslope_log_met = met_para["mslope_log_met"];
  
  double gmax            = gamma_para["gmax"];
  double gamma_mu        = gamma_para["gamma_mu"];
  double gamma_sigma     = gamma_para["gamma_sigma"];
  double gamma_beta      = gamma_para["gamma_beta"];
  double gamma_alpha     = gamma_para["gamma_alpha"];
  double gamma_prob      = gamma_para["gamma_prob"];
  double gamma_mu_m0     = gamma_para["gamma_mu_m0"];
  double gamma_mu_v02    = gamma_para["gamma_mu_v02"];
  double mu_bot          = gamma_para["mu_bot"];
  double mu_top          = gamma_para["mu_top"];
  double gamma_sigma_m0  = gamma_para["gamma_sigma_m0"];
  double gamma_sigma_v02 = gamma_para["gamma_sigma_v02"];
  double sigma_bot       = gamma_para["sigma_bot"];
  double sigma_top       = gamma_para["sigma_top"];
  colvec gam0vec         = gamma_para["gam0vec"];
  double wk              = gamma_para["wk"];
  
  double tau2_s0 = tau_para["tau2_s0"];
  double tau2_rn = tau_para["tau2_rn"];
  
  double eta02_s0 = eta_para["eta02_s0"];
  double eta02_rn = eta_para["eta02_rn"];
  double eta02    = eta_para["eta02"];
  double eta0     = eta_para["eta0"];
  colvec th0v     = eta_para["th0v"];
  
  colvec mslopeall            = slope_para["mslopeall"];
  double mslope0_log          = slope_para["mslope0_log"];
  colvec mslopeall_log        = slope_para["mslopeall_log"];
  double mslope_log_v0i       = slope_para["mslope_log_v0i"];
  double mslope_log_v0im0     = slope_para["mslope_log_v0im0"];
  double mslope_log_sigma2    = slope_para["mslope_log_sigma2"];
  double mslope_log_simga2_rn = slope_para["mslope_log_simga2_rn"];
  double mslope_log_sigma2_s0 = slope_para["mslope_log_sigma2_s0"];
  
  colvec zeta_m0        = squish_para1["zeta_m0"];
  colvec zeta_v02       = squish_para1["zeta_v02"];
  double zeta_s0        = squish_para1["zeta_s0"];
  double zeta_rn        = squish_para1["zeta_rn"];
  double psi_log_m0     = squish_para1["psi_log_m0"];
  double psi_log_v02    = squish_para1["psi_log_v02"];
  double psi_log_s0     = squish_para1["psi_log_s0"];
  double psi_log_rn     = squish_para1["psi_log_rn"];
  double psi_fixed      = squish_para1["psi_fixed"];
  double psi_fixed_log  = squish_para1["psi_fixed_log"];
  double psi0           = squish_para1["psi0"];
  double psi_log_mu     = squish_para1["psi_log_mu"];
  double psi_log_sigma2 = squish_para1["psi_log_sigma2"];
  double psi_log_sigma  = squish_para1["psi_log_sigma"];
  
  colvec zeta_sigma2    = squish_para2["zeta_sigma2"];
  colvec zeta_sigma     = squish_para2["zeta_sigma"];
  colvec zeta0          = squish_para2["zeta0"];
  colvec omega0         = squish_para2["omega0"];
  mat    zetall         = squish_para2["zetall"];
  colvec psiall         = squish_para2["psiall"];
  colvec psiall_log     = squish_para2["psiall_log"];
  mat    omegall        = squish_para2["omegall"];
  mat    hfunall        = squish_para2["hfunall"];
  umat   idball         = squish_para2["idball"];
  uvec   id_stat        = squish_para2["id_stat"];
  uvec   id_flex        = squish_para2["id_flex"];
  
  
  double sigma2_s0        = sigma_para(0);
  double sigma2_rn        = sigma_para(1);
  
  double sigma2_alpha     = sigma_para(2);
  double sigma2_beta      = sigma_para(3);
  double sigma2_alpha_m0  = sigma_para(4);
  double sigma2_alpha_v02 = sigma_para(5);
  double sigma2_beta_r0   = sigma_para(6);
  double sigma2_beta_s0   = sigma_para(7);
  
  mat alpha_v0i      = v_para["alpha_v0i"];
  colvec alpha_v0im0 = v_para["alpha_v0im0"];
  
  mat phi_v0i      = z_para["phi_v0i"];
  double lambda_fn = z_para["lambda_fn"];
  mat lambda_g0i   = z_para["lambda_g0i"];
  
  colvec kall = basis_para["kall"];
  
  mat phixgrid = phi_para["phixgrid"];
  mat phixdata = phi_para["phixdata"];
  
  mat fxgridall  = fx_para["fxgridall"];
  colvec f0xgrid = fx_para["f0xgrid"];
  colvec fxdata  = fx_para["fxdata"];
  colvec fxdatam = fxdata;
  
  mat thetam        = zeros<mat>(ntheta,ngroup);  
  mat thetas        = zeros<mat>(ntheta,ngroup);
  mat betam         = zeros<mat>(size(betall));
  mat betas         = zeros<mat>(size(betall));
  
  mat fxgridm       = zeros<mat>(nint+1,ngroup);
  mat fxgrids       = zeros<mat>(nint+1,ngroup);
  colvec fxgrid_old = zeros<colvec>(nint+1);
  colvec f0xgridm   = zeros<colvec>(nint+1);
  colvec f0xgrids   = zeros<colvec>(nint+1);
  colvec f0Bxgridm  = zeros<colvec>(nint+1);
  colvec f0Bxgrids  = zeros<colvec>(nint+1);
  colvec fxdatas    = zeros<colvec>(ntot);
  
  // Matrics for saving MCMC iterations
  unsigned int nparv_save = nparv > 0 ? nparv : 1;
  mat alphag = zeros<mat>(smcmc, nparv_save);
  
  unsigned int nsigma_save = bflagHBsigma ? ngroup : 1;
  mat sigmag        = zeros<mat>(smcmc,nsigma_save);    // Error variance
  vec sigma2_alphag = zeros<vec>(smcmc);
  vec sigma2_betag  = zeros<vec>(smcmc);
  mat tauallg       = zeros<mat>(smcmc,ntheta);    // StdDev for spectral coefficients
  vec eta0g         = zeros<vec>(smcmc);           // STDEV for upper level spectral coefficients
  vec gamma_mugg    = zeros<vec>(smcmc);           // Smoothing parameter for upper level model
  mat gammallg      = zeros<mat>(smcmc,ngroup);    // Smoothing parameters for lower level model
  vec gamma_alphag  = zeros<vec>(smcmc);
  vec gamma_betag   = zeros<vec>(smcmc);
  mat theta0g       = zeros<mat>(smcmc,ntheta);
  mat phig          = zeros<mat>(smcmc,phidim);
  mat lambdag       = zeros<mat>(smcmc,nparw2);
  
  
  // slope parameter
  colvec mresid, met_var_mslope_new;
  double mslope_log_old, met_var_mslope, met_std_mslope, mslope_log_new, mslope_new, 
  sse_slopes, mslope0, mslope_log_sigma2_sn, mslope_log_sigma, slope_bn;
  
  // local
  colvec yresid = zeros<colvec>(ntot);
  colvec rstar  = zeros<colvec>(ntot);
  colvec ymean  = zeros<colvec>(ntot);
  colvec yresidj, rj, rjlead, rjlag, yjstar, v, rVec, bn, bj;
  mat rMat, ranMat, vbi, vbi12, vbni, vbni12;
  double sse, ck, vni, vn;
  unsigned int nj;
  
  // ACE
  double ace_vbn, ace_sbn, ace_bn, fa, fb, u, prob, rhoj, srj, rmin, rmax;
  
  // alpha
  mat vtvj, vj;
  
  // beta
  mat wj, zbreg, wtwj;
  colvec zbregj, v0ib0, wjbj;
  
  
  // sigma 
  //if non-homogeneous, sigma2(j) has same value among all j
  double s2, sigma2_sn;
  
  // sigma_alpha, sigma_beta 
  double sigma2_beta_sn, sigma2_beta_rn, testpr, var_met_sigma, std_met_sigma, sigma2_alpha_new, lsum, pnew, pold;
  
  // phi
  mat phi_vbni   = zeros<mat> (nparz*nparw,nparz*nparw);
  mat phi_vbni12 = zeros<mat> (nparz*nparw,nparz*nparw);
  mat ztb        = zeros<mat> (nparz,nparz);
  colvec phi_vec = zeros<colvec> (nparz*nparw);
  colvec phi_bn  = zeros<colvec> (nparz*nparw);
  colvec a       = zeros<colvec> (nparz);
  mat resid_mat  = zeros<mat> (ngroup,nparw);
  mat lambda_gni = zeros<mat> (nparw, nparw);
  mat lambda_gn  = zeros<mat> (nparw, nparw);
  
  // theta
  mat vbn, vbn12, phixj, phi2j, gamjvec, thetaj, thv, thvi, bmat;
  double sse_new, sse_old, gammaj, gpar, nz, nz2, t0resid2_new, t0resid2_old, testp;
  uvec zc, zvi, zv, z, zi;
  colvec tkresid2_new = zeros<colvec>(ntheta-1);
  colvec tresid2_new  = zeros<colvec>(ntheta);
  colvec tkresid2_old = zeros<colvec>(ntheta-1);
  colvec tresid2_old  = zeros<colvec>(ntheta);
  Rcpp::List fpars;
  colvec resid_new, resid_old, fxgrid_new, thvk, fxj, theta0k, fxdata_new, fxdata_star_new,
  fxdata_old, fxdata_star_old, theta_old, thetak_old, gamvec, ths, met_var_new, met_var,
  met_std, hfunj_new;
  double xiparj_old=0;
  double xiparj_new=0;
  colvec thetak_new = zeros<colvec>(ntheta-1);
  colvec theta_new  = zeros<colvec>(ntheta);
  
  // squish
  double psij_log_old, std_met_psi,
  psi_log_v2, psi_log_v, psi_log_mn, zeta_sn, psi_log_sn,
  x, x0, psij_old, psij_log_new, psij_new;
  
  vec breakers;
  colvec resid, zeta_vn2, zeta_vn, zeta_mn, zetaj_old, 
  var_met_psi_new, zvec, xvec, met_var_psi_new;
  
  rowvec resid_zeta        = zeros<rowvec> (ngroup);
  colvec resid_zeta_new    = zeros<colvec> (nomega);
  colvec resid_zeta_old    = zeros<colvec> (nomega);
  colvec zetaj_new         = zeros<colvec> (nomega);
  colvec var_met_zeta      = zeros<colvec> (nomega);
  colvec std_met_zeta      = zeros<colvec> (nomega);
  colvec met_var_zeta_new  = zeros<colvec> (nomega);
  colvec omegaj_new        = zeros<colvec> (nomega);
  colvec omegaj_old        = zeros<colvec> (nomega);
  colvec met_var_omega_new = zeros<colvec> (nomega);
  colvec sez               = zeros<colvec> (nomega);
  
  unsigned int kbot, tbot;
  colvec vpar, thvb;
  
  // theta0
  colvec t0var;
  double xi_vn, xi_bn;
  
  unsigned int tk;
  
  colvec resid2, gamk;
  double tau2_sn;
  
  // tau2
  colvec resid_tau  = zeros<colvec>(ngroup);
  colvec resid2_tau = zeros<colvec>(ngroup);
  
  // eta02
  colvec resid2_eta;
  double eta02_sn;
  
  // gamma
  colvec u1, gam0vec_new, th0v_new, theta02;
  double u2, bmin, bmax, gamma_sum, lgamma_sum, var_met_g0, std_met_g0, gamma_mu_new,
  gamma_sigma_new, gamma_sigma2_new, gamma_beta_new, gamma_alpha_new, gamma_prob_new,
  testpa;
  
  // ordinal
  colvec yk, yok, mk, cutt, ytopk, ybotk, ord_fa, ord_fb, ord_p, colvec_tmp;
  double sk, jk, cutjk, cbot, ctop;
  uvec id_cut_botk, id_cut_topk, idz, zjbot, zjtop;
  vec vec_tmp;
  
  colvec ydatam = zeros<colvec>(ntot);
  
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
    // cout << imcmc << endl;
    // Generate fixed effect alpha
    if (nparv > 0) {
      yresid = ydata - wbeta - fxdata;   // residuals
      vbi    = alpha_v0i;
      bn     = alpha_v0im0;
      for (unsigned int j = 0; j < ngroup; j++) {
        uvec gidx = iptHB1[j];
        rj   = yresid.elem(gidx);
        vtvj = vtv.slice(j);
        vj   = vdata.rows(gidx);
        // Correction for ACE 
        if (bflagACE) {
          rhoj = rhoall(j);
          // First observation has stationary variance 
          srj        = 1-rhoj;
          rj(0)     *= srj;
          vj.row(0) *= srj;
          // Observations 2 to n are adjusted for rho
          rj.tail(rj.n_elem-1) -= rhoj*rj.head(rj.n_elem-1);
          vj.tail_rows(vj.n_rows-1) -= rhoj*vj.head_rows(vj.n_rows-1);
          vtvj = vj.t() * vj;
        }
        
        vbi += vtvj/sigma2(j);
        bn  += (vj.t() * rj)/sigma2(j);
      }
      
      vbi12  = chol(vbi);
      rMat   = randn(nparv, 1);
      alpha  = solve(vbi, bn + vbi12.t()*rMat);
      valpha = vdata * alpha;
    } // End generate alpha
    
    // Do HB parametric model
    yresid = ydata - fxdata - valpha;
    //--------------------------------------------
    // Generate Beta_j
    zbreg = zdata * phi;
    sse = 0;
    for (unsigned int j = 0; j < ngroup; j++) {
      uvec gidx = iptHB1[j];
      zbregj = zbreg.row(j).t();
      v0ib0  = lambdai * zbregj;
      wj     = wdata.rows(gidx);
      rj     = yresid.elem(gidx);
      wtwj   = wtw.slice(j);
      if (bflagACE) {
        rhoj = rhoall(j);
        // First observation has stationary variance 
        srj = 1-rhoj;
        rj(0) *= srj;
        wj.row(0) *= srj;
        // Observations 2 to n are adjusted for rho
        rj.tail(rj.n_elem-1) -= rhoj*rj.head(rj.n_elem-1);
        wj.tail_rows(wj.n_rows-1) -= rhoj*wj.head_rows(wj.n_rows-1);
        wtwj = wj.t() * wj;
      }
      
      vbni   = wtwj / sigma2(j) + lambdai;
      bn     = (wj.t() * rj) / sigma2(j) + v0ib0;
      vbni12 = chol(vbni);   // Upper triangluar: t(vbni12) %*% vbni12 = vbni
      // Generate draw by using solve: avoides inverse of vni
      ranMat        = randn(nparw, 1);
      bj            = solve(vbni, bn + vbni12.t()*ranMat);
      betall.row(j) = bj.t();
      wj            = wdata.rows(gidx);
      wjbj          = wj * bj;
      wbeta.rows(gidx) = wjbj;
      // Compute residual
      rj            = yresid.elem(gidx) - wjbj;
      // Correct for ACE
      if (bflagACE) {
        rj(0)  *= srj;
        rj.tail(rj.n_elem-1) -= rhoj*rj.head(rj.n_elem-1);
      }
      
      if (bflagHBsigma) {
        if (iflaglm == 1) {
          // probit model has fixed error variance
          sigma2(j) = 1;
          sigma(j)  = 1;
        } else{
          // HB error variance
          sigma2_rn = sigma2_alpha + nobs(j);
          sigma2_sn = sigma2_beta  + accu(square(rj));
          s2        = sigma2_sn / (2 * randg(distr_param(sigma2_rn / 2, 1.0)));
          sigma2(j) = s2;
          sigma(j)  = sqrt(s2);
        } // End probit or other
      } // End generate sigma2 in HB model
    } // End Loop over groups to generate beta_j and sigma2_j
    // HB model for sigma if not probit
    if (bflagHBsigma && iflaglm != 1) {
      //-----------------------------------------------------------------
      // Generate sigma2_alpha and sigma2_beta
      // HB model for variance: sigma_j^2 ~ IG(alpha/2,beta/2)
      // Hyperpriors:   alpha ~ N(m0,v02)I(alpha>2) and beta~G(r0/2,s0/2)
      // Generate sigma2_beta
      sigma2_beta_sn = sigma2_beta_s0 + accu(1/sigma2);
      sigma2_beta_rn = sigma2_beta_r0 + ingroup*sigma2_alpha;
      sigma2_beta    = randg(distr_param(sigma2_beta_rn/2,2/sigma2_beta_sn));
      // Generate sigma2_alpha
      met_var_new      = AdaptMetVar(sigma2_met);
      ck               = 5.66;   // Constant from Harito
      var_met_sigma    = ck*met_var_new(0);
      std_met_sigma    = sqrt(var_met_sigma);
      sigma2_alpha_new = rcpp_rndtnb(sigma2_alpha,std_met_sigma,2); // Generate candidate from truncated normal
      
      lsum        = accu(log(sigma2));
      //  IG Likelihood
      testpr      = ingroup*log(sigma2_beta/2)*(sigma2_alpha_new - sigma2_alpha)/2;
      testpr      = testpr - ingroup*(lgamma(sigma2_alpha_new/2)-lgamma(sigma2_alpha/2));
      testpr      = testpr - lsum*(sigma2_alpha_new-sigma2_alpha)/2;
      // Prior
      testpr      = testpr - pow(sigma2_alpha_new - sigma2_alpha_m0, 2)/(2*sigma2_alpha_v02);
      testpr      = testpr + pow(sigma2_alpha     - sigma2_alpha_m0, 2)/(2*sigma2_alpha_v02);
      // Generating function
      pnew        = 1-normcdf(2.0,sigma2_alpha,std_met_sigma);
      pold        = 1-normcdf(2.0,sigma2_alpha_new,std_met_sigma);
      testpr      = testpr + log(pnew) - log(pold);
      if (log(randu())<testpr) {
        sigma2_alpha = sigma2_alpha_new;
        //sigma2_mu             = sigma2_beta/(sigma2_alpha -2);
        sigma2_met(1)++; // pmet
        sigma2_met.tail(sigma2_met(0)) = met_var_new;
      }
      AdaptMetUpdate(sigma2_met);
    }
    
    // Update mean
    ymean  = valpha + wbeta + fxdata;
    yresid = ydata - ymean;
    
    // Generate sigma2 for homogenous error variance
    if (!bflagHBsigma) {
      //Probit or else
      if (iflaglm == 1) {
        for (unsigned int j=0; j<ngroup; j++) {
          sigma2(j) = 1;
          sigma(j)  = 1;
        }
      } else {
        rstar = yresid;
        // Adjust for ACE
        if (bflagACE) {
          // Loop over groups
          for (unsigned int j=0; j<ngroup; j++) {
            uvec gidx = iptHB1[j];
            rhoj = rhoall(j);
            srj  = 1-rhoj;
            rj   = yresid.elem(gidx);
            nj   = nobs(j);
            // Correct for ACE
            if (bflagACE) {
              rj(0) *= srj;
              rj.tail(rj.n_elem-1) -= rhoj*rj.head(rj.n_elem-1);
            }
            rstar.elem(gidx) = rj;
          }  // End loop over groups
        } // End adjust for ACE
        // Generate homogeneous sigma
        sigma2_sn = sigma2_s0 + sum(square(rstar));
        sigma2    = sigma2_sn/(2 * randg(distr_param(sigma2_rn / 2, 1.0))) * ones<colvec>(ngroup);
        sigma     = sqrt(sigma2);
      } // End probit or else
    }  // End generate homogeneous sigma
    
    
    //------------------------------------------
    // Generate Phi (Beta = Zdata*Phi + delta)
    phi_vbni   = kron(lambdai, ztz) + phi_v0i;
    phi_vbni12 = chol(phi_vbni);
    ztb        = zdata.t() * betall;
    for (unsigned int j = 0; j < nparw; j++) {
      if (nparw > 1) {
        a = ztb * lambdai.col(j);
      }
      else {
        a = ztb * lambdai(0, 0);
      }
      
      phi_bn.subvec(nparz*j, nparz*(j+1)-1) = a;
      
    }
    rMat = randn(phidim, 1);
    phi_vec  = solve(phi_vbni, phi_bn + phi_vbni12.t()*rMat);
    phi      = reshape(phi_vec, nparz,nparw);
    //------------------------------------------
    // Generate Lambda from Beta = Zdata*breg + delta  
    resid_mat = betall - zdata * phi;
    if (nparw > 1) {
      lambda_gni = lambda_g0i + resid_mat.t()*resid_mat;
      lambda_gn  = inv(lambda_gni);
      lambdai    = wishrnd(lambda_gn, lambda_fn); // random draw from a Wishart distribution
      lambda     = inv(lambdai);
    } else {
      sse       = accu(square(resid_mat));
      lambda_gn = lambda_g0i + sse;
      lambda    = lambda_gn / (2 * randg(distr_param(lambda_fn / 2, 1.0)));
      lambdai   = 1/lambda;
    }
    
    // Generate rho for ACE model
    ymean  = valpha + wbeta + fxdata;
    yresid = ydata - ymean;
    if (bflagACE) {
      // Loop over groups
      for (unsigned int j=0; j<ngroup; j++) {
        uvec gidx = iptHB1[j];
        rj       = yresid.elem(gidx);
        s2       = sigma2(j);
        rjlead   = rj.tail(rj.n_elem-1); // r2, r3, ..., rn
        rjlag    = rj.head(rj.n_elem-1); // r1, r2, ..., rn-1
        
        // Posterior mean and variance 
        ace_vbn  = 1/(rho_v0i + accu(square(rjlead))/s2);
        vec_tmp  = rjlag.t()*rjlead;
        ace_bn   = ace_vbn*(rho_v0im0 + vec_tmp(0)/s2);
        
        // Slice sampling
        v      = sqrt(1-rhoj*rhoj)*randu();
        v      = sqrt(1-square(v));
        vec_tmp.zeros(v.n_elem+1);
        vec_tmp.head(v.n_elem) = v;
        vec_tmp.tail(1) = .99;
        rmin    = max(-vec_tmp);
        rmax    = min(vec_tmp);
        ace_sbn = sqrt(ace_vbn);
        fa   = normcdf(-1.0, ace_bn, ace_sbn);
        fb   = normcdf( 1.0, ace_bn, ace_sbn);
        u    = randu();
        prob = fa + (fb-fa)*u;
        if (prob>0.9999) {
          // Hits top 
          rhoj = .9;
        } else {
          if (prob<0.0001) {
            // Hits bottom
            rhoj = -0.9;
          } else {
            // Generate between bounds
            rhoj = R::qnorm(prob, ace_bn, ace_sbn, true, false);
          }
        }
        rhoall(j) = rhoj;
      } // End loop over groups
    } // End generate ACE
    
    //------------------------------------------
    // Do HB nonparametric model
    yresid = ydata - valpha - wbeta;   // Residuals without fxdata
    // Free model: do not need numerical integration
    if (iflagsc==0) {
      // Loop over groups to generate theta_j
      for (unsigned int j = 0; j < ngroup; j++) { // Loop over groups to generate theta_j
        uvec gidx = iptHB1[j];
        rj    = yresid.elem(gidx);
        phixj = phixdata.rows(gidx);
        phi2j = phi2.slice(j); // Get phi_j'phi_j
        
        if (bflagACE) {
          rhoj  = rhoall(j);
          srj   = 1-rhoj;
          // Adjust first observation for stationary variance
          rj(0)  *= srj;
          phixj.row(0) *= srj;
          // Adjust observations 2 to n for rho
          rj.tail(rj.n_elem-1) -= rhoj*rj.head(rj.n_elem-1);
          phixj.tail_rows(phixj.n_rows-1) -= rhoj*phixj.head_rows(phixj.n_rows-1);
          phi2j = phixj.t() * phixj;
        } // End adjusting for ACE
        thetaj  = thetall.col(j);
        gammaj  = gammall(j);
        gamjvec = exp(-kall*gammaj);
        thv     = tau2all%gamjvec;
        // Worry about rounding and numerical stability for small variances
        uvec id = find(thv < 1e-10);
        if (id.n_elem > 0) {
          thv.elem( id ) = 1e-10*ones<colvec>(id.n_elem);
        }
        thvi        = 1/thv;
        mat bmat    = phi2j/sigma2(j);
        mat vbni    = bmat;
        vbni.diag() += thvi;
        rMat = randn(ntheta, 1);
        // Issue with poorly conditioned vbni
        if (rcond(vbni) > 1e-20) {
          vbni12 = chol(vbni); // Upper triangluar: t(vbni12) %*% vbni12 = vbni
          bn     = (phixj.t() * rj)/sigma2(j) + theta0/thv;
          // Generate draw by using solve: avoides inverse of vni
          
          thetaj = solve(vbni, bn + vbni12.t()*rMat);
        } else {
          colvec eigval = zeros<colvec>(nbasis);
          colvec eigvec = zeros<colvec>(nbasis,nbasis);
          eig_sym( eigval, eigvec, vbni );
          uvec id0 = find(eigval < 1e-10);
          eigval.elem(id0) = 1e-10*ones<vec>(id0.n_elem);
          vbni   = eigvec * eigval.diag() * eigvec.t();
          vbni12 = chol(vbni);   // Upper triangluar: t(vbni12) %*% vbni12 = vbni
          bn     = vbn * ((phixj.t() * rj)/sigma2(j) + theta0/thv);
          // Generate draw by using solve: avoides inverse of vni
          thetaj = solve(vbni, bn + vbni12.t()*rMat);
        }
        
        thetall.col(j) = thetaj;
        
        // Compute new fj
        fxgridall.col(j)  = phixgrid * thetaj;  // fj computed on xgrid
        phixj             = phixdata.rows(gidx);
        fxj               = phixj * thetaj;     // fj computed at observations
        fxdata.rows(gidx) = fxj;
      } // End loop over groups to generate theta_j
    } else {
      // Got shape constraint.  Do metropolis algorithm to generate theta_{j,k}
      // Only S and multi-modal used omega (stationary/infection) and squish function
      // Initialize parameters that are not used in all models
      // Reset pmet if beyound adaptation period
      if (imcmc == nblow + maxmodmet*nblow0) {
        sigma2_met(1) = 0;
        gamma_met(1)  = 0;
        for (unsigned int j = 0; j < ngroup; j++) {
          theta_met(1,j)      = 0; 
          zeta_met(1,j)       = 0;
          psi_met(1,j)        = 0;
          mslope_log_met(1,j) = 0;
        }  // end loop over groups to set pmet = 0
      } // end set pmet to zero at start of SMCMC
      
      //-------------------------------------------------------------------
      // Worry about g(Z) = Z^2 or g(Z) =exp(Z)
      // theta_{j,0} has different models for Z^2 and exp(Z)
      // If Z^2, then theta_{j,0} is log-normal
      // If exp(Z), then theta_{j,0} is normal.
      // df/dx = Z^2 so theta_j,0 = exp(xipar_j)
      // Generate theta with shape constraints
      //theta00 = theta0(0);                   // Upper-level Fourier coefficient for constant
      theta0k = theta0.tail(theta0.n_elem-1);  // Upper-level Fourier coefficient for consines
      // Loop over groups to generate theta_j
      // kbot = iflagZcnst + 1; // bottom index
      for (unsigned int j = 0; j < ngroup; j++) {
        uvec gidx = iptHB1[j];
        yresidj   = yresid(gidx);    // Parametric residuals for group j
        
        fxdata_old = fxdata(gidx);   // fJ(x) at current values of theta
        theta_old  = thetall.col(j); // Old value for theta_j
        
        gammaj = gammall(j);
        gamvec = exp(-gammaj*kall);
        
        thv    = tau2all % gamvec;
        // Worry about zero variances: which creates testp=NaN.
        // If Var(theta_jk) is too small, set theta_jk = theta_0k
        zv     = find(ths < 1e-10);
        zvi    = find(ths >= 1e-10);
        nz     = zv.n_elem;
        if (nz > 0) {
          thv(zv) = 1e-20*ones<colvec>(nz);
        }
        ths    = sqrt(thv);
        //-----------------------------------------------------------------
        // get variances for t-dist random walk
        met_var_new = AdaptMetVar(theta_met.col(j));
        
        ck        = 5.66;   // Constant from Harito
        met_var   = ck*met_var_new;
        met_std   = sqrt(met_var);
        // Somtimes met_std goes to 0.  Shake it up
        z         = find(met_std < 1e-10);
        nz2       = z.n_elem;
        if(nz2 > 0){
          colvec sub_met_std = zeros<colvec>(nz2); 
          for (unsigned int jj=0; jj<nz2; jj++) {
            if (randu() < 0.1) {
              sub_met_std(jj) = 1e-10;
            }
          } 
          met_std(z) = sub_met_std;
          met_var(z) = pow(met_std(z), 2);
        }
        
        // Random walk from for normals
        if (!bflagZ) {
          // With E(Z^2), constraint theta_j0 > 0 with theta_j0 = exp(xi_j)
          xiparj_old     = xiparall(j);                        // theta_j0 = exp(xi_j)
          thetak_old     = theta_old.tail(theta_old.n_elem-1); //Skip theta_j0
          t0resid2_old   = pow(xiparj_old - xipar0,2);  // squared residual  k=0
          tkresid2_old   = square(thetak_old - theta0k); // squared residuals k>0
          tresid2_old    = zeros<colvec>(ntheta);
          tresid2_old(0) = t0resid2_old;
          tresid2_old.tail(tkresid2_old.n_elem) = tkresid2_old;
          
          // Generate theta_new for Z^2
          xiparj_new   = xiparj_old + met_std(0)*ths(0)*randn(); 
          thetak_new   = thetak_old + met_std.tail(met_std.n_elem-1)%ths.tail(ths.n_elem-1)%randn(nbasis,1);
          theta_new(0) = exp(xiparj_new);
          theta_new.tail(thetak_new.n_elem) = thetak_new;
          
          // If variance is too small, replace theta_jk with theta_0k
          if (nz>0) theta_new(zv) = theta0(zv);
          t0resid2_new   = pow(xiparj_new - xipar0,2);
          tkresid2_new   = square(thetak_new - theta0k);
          tresid2_new    = zeros<colvec>(ntheta);
          tresid2_new(0) = t0resid2_new;
          tresid2_new.tail(tkresid2_new.n_elem) = tkresid2_new;
          
          
        } else {
          // Generate theta_new for exp(z)
          theta_new = theta_old + met_std%ths%randn(ntheta);
          // If variance is too small, replace theta_jk with theta_0k
          if(nz>0) theta_new(zv) = theta0(zv);
          tresid2_old = square(theta_old - theta0);
          tresid2_new = square(theta_new - theta0);
        }
        // HB Prior
        // Remove entries where thv = 0 & thetak_new = theta_0k
        if (nz>0) {
          tresid2_new = tresid2_new(zvi);
          tresid2_old = tresid2_old(zvi);
          thv         = thv(zvi);
        }
        
        testp = 0;
        testp = testp - accu(tresid2_new/thv)/2 + accu(tresid2_old/thv)/2;
        
        // Get fJ(x) for theta_new
        // Get function for group j
        fpars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                   Rcpp::Named("iflagpn")     = iflagpn,
                                   Rcpp::Named("iflagCenter") = bflagCenter,
                                   Rcpp::Named("iflagZ")      = bflagZ,
                                   Rcpp::Named("theta")       = theta_new,
                                   Rcpp::Named("phix")        = phixgrid,
                                   Rcpp::Named("delta")       = xdelta,
                                   Rcpp::Named("range")       = xrange,
                                   Rcpp::Named("xmin")        = xmin,
                                   Rcpp::Named("xgrid")       = xgrid,
                                   Rcpp::Named("nstat")       = nstat,
                                   Rcpp::Named("hfun")        = hfunall.col(j),
                                   Rcpp::Named("omega")       = omegall.col(j),
                                   Rcpp::Named("idb")         = idball.col(j),
                                   Rcpp::Named("id_stat")     = id_stat,
                                   Rcpp::Named("mslope")      = mslopeall(j));
        
        
        fxgrid_new = rcpp_Getfx(fpars);
        //---------------------------------------------------
        // Compute fj(x_obs)
        // fxdata = GetUpfxdata(xgrid,xinx,xover,fxgrid)
        umat a     = xinx.rows(gidx);
        vec  b     = xover(gidx);
        fxdata_new = rcpp_GetSCfxobs(a, b, fxgrid_new);
        //---------------------------------------------------
        // Metropolis Test
        // Compute Likelihood
        // Adjust for ACE
        yresidj          = yresid.elem(gidx);
        yjstar           = yresidj;
        fxdata_star_new  = fxdata_new;
        fxdata_star_old  = fxdata_old;
        if (bflagACE) {
          rhoj = rhoall(j);
          srj  = 1 - rhoj;
          nj   = nobs(j);
          // Adjust first observation for stationary variance
          yjstar(0)          *= srj;
          fxdata_star_new(0) *= srj;
          fxdata_star_old(0) *= srj;
          // AR model for observations 2 to n
          yjstar.tail(nj-1) = yjstar.tail(nj-1) - rhoj*yjstar.head(nj-1);
          fxdata_star_new.tail(nj-1) = fxdata_star_new.tail(nj-1) - rhoj*fxdata_star_new.head(nj-1);
          fxdata_star_old.tail(nj-1) = fxdata_star_old.tail(nj-1) - rhoj*fxdata_star_old.head(nj-1);
        } // End ACE adjustment
        // ******************************************************************
        // Add log-likelihood to testp
        // ******************************************************************
        resid_new = yjstar - fxdata_star_new;  // Note: fxdata_new_j only for group j
        sse_new   = accu(square(resid_new));
        resid_old = yjstar - fxdata_star_old;
        sse_old   = accu(square(resid_old));
        testp = testp - (sse_new - sse_old)/(2*sigma2(j)); // Likelihood
        
        if (log(randu()) < testp) { // Accept candidate
          if (!bflagZ) xiparall(j) = xiparj_new;
          thetall.col(j)                        = theta_new;
          fxdata.rows(gidx)                     = fxdata_new;
          fxgridall.col(j)                      = fxgrid_new;
          theta_met.col(j).tail(theta_met(0,j)) = met_var_new;
          theta_met(1,j)++; // pmet
          // Use sse_old in updating Squish function
        }
        // ********************************************************************
        // Monotone convex/concave:  Generate linear term ----
        // ********************************************************************
        if ((iflagsc==2)||(iflagsc==3)) {
          testp              = 0;
          fxdata_old         = fxdata.elem(gidx);
          fxgrid_old         = fxgridall.col(j);
          mslope_log_old     = mslopeall_log(j);
          met_var_mslope_new = AdaptMetVar(mslope_log_met.col(j));
          ck                 = 5.66; // Constant from Harito
          met_var_mslope     = ck*met_var_mslope_new(0);
          met_std_mslope     = sqrt(met_var_mslope);
          mslope_log_new     = mslope_log_old + met_std_mslope*randn();
          mslope_new         = exp(mslope_log_new);
          sse_slopes         = - pow(mslope_log_new-mslope0_log,2)/(2*mslope_log_sigma2) + pow(mslope_log_old-mslope0_log,2)/(2*mslope_log_sigma2);
          testp              = testp + sse_slopes;
          
          // Get fJ(x) for theta_new
          // Get function for group j
          fpars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                     Rcpp::Named("iflagpn")     = iflagpn,
                                     Rcpp::Named("iflagCenter") = bflagCenter,
                                     Rcpp::Named("iflagZ")      = bflagZ,
                                     Rcpp::Named("theta")       = thetall.col(j),
                                     Rcpp::Named("phix")        = phixgrid,
                                     Rcpp::Named("delta")       = xdelta,
                                     Rcpp::Named("range")       = xrange,
                                     Rcpp::Named("xmin")        = xmin,
                                     Rcpp::Named("xgrid")       = xgrid,
                                     Rcpp::Named("nstat")       = nstat,
                                     Rcpp::Named("hfun")        = hfunall.col(j),
                                     Rcpp::Named("omega")       = omegall.col(j),
                                     Rcpp::Named("idb")         = idball.col(j),
                                     Rcpp::Named("id_stat")     = id_stat,
                                     Rcpp::Named("mslope")      = mslope_new);
          
          
          fxgrid_new = rcpp_Getfx(fpars);
          
          //---------------------------------------------------
          // Compute fj(x_obs)
          // fxdata = GetUpfxdata(xgrid,xinx,xover,fxgrid)
          umat a     = xinx.rows(gidx);
          vec  b     = xover(gidx);
          fxdata_new = rcpp_GetSCfxobs(a, b, fxgrid_new);
          
          //---------------------------------------------------
          // Metropolis Test
          // Compute Likelihood
          // Adjust for ACE
          yresidj          = yresid.elem(gidx);
          yjstar           = yresidj;
          fxdata_star_new  = fxdata_new;
          fxdata_star_old  = fxdata_old;
          if (bflagACE) {
            rhoj = rhoall(j);
            srj  = 1- rhoj;
            nj   = nobs(j);
            // Adjust first observation for stationary variance
            yjstar(0)          *= srj;
            fxdata_star_new(0) *= srj;
            fxdata_star_old(0) *= srj;
            // AR model for observations 2 to n
            yjstar.tail(nj-1) = yjstar.tail(nj-1) - rhoj*yjstar.head(nj-1);
            fxdata_star_new.tail(nj-1) = fxdata_star_new.tail(nj-1) - rhoj*fxdata_star_new.head(nj-1);
            fxdata_star_old.tail(nj-1) = fxdata_star_old.tail(nj-1) - rhoj*fxdata_star_old.head(nj-1);
          } // End ACE adjustment
          // ******************************************************************
          // Add log-likelihood to testp
          // ******************************************************************
          resid_new = yjstar - fxdata_star_new;  // Note: fxdata_new_j only for group j
          sse_new   = accu(square(resid_new));
          resid_old = yjstar - fxdata_star_old;
          sse_old   = accu(square(resid_old));
          testp = testp - (sse_new - sse_old)/(2*sigma2(j)); // Likelihood
          if (log(randu()) < testp) { // Accept candidate
            fxdata.rows(gidx) = fxdata_new;
            fxgridall.col(j)  = fxgrid_new;
            mslopeall_log(j)  = mslope_log_new;
            mslopeall(j)      = mslope_new;
            mslope_log_met.col(j).tail(mslope_log_met(0,j)) = met_var_mslope_new(0);
            mslope_log_met(1,j)++; // pmet
          } // End accept candidate theta_new
          // End Generate New f for slopes for Monotone convex/concave
        } // End generate slope for Monotone convex/concave
        // ********************************************************************
        // U-shaped functions: generate stationary point. No Squish function
        // ********************************************************************
        if (iflagsc == 4) {
          testp=0;
          fxdata_old  = fxdata.elem(gidx);
          fxgrid_old  = fxgridall.col(j);
          // ******************************************************************
          // Generate zeta and omega is logit  
          zetaj_old        = zetall.col(j);
          omegaj_old       = omegall.col(j);
          met_var_zeta_new = AdaptMetVar(zeta_met.col(j));
          ck               = 5.66;   // Constant from Harito
          var_met_zeta(0)  = ck*met_var_zeta_new(0);
          std_met_zeta(0)  = sqrt(var_met_zeta(0));
          // 5% of time generate omega from uniform on xmin to xmax
          // 95% of time generate zeta from random walk metropolis
          if (randu()<.05) {
            // Randomly reinitialize zeta & omega to avoid getting stuck
            omegaj_new(0) = xmin + (xmax-xmin)*randu();
            // Find corresponding  zeta0_input for ordinal logit
            x            = (omegaj_new(0)-xmin)/xrange;
            x0           = 1-sum(x);
            zetaj_new(0) = log(x/x0);
          } else {
            zetaj_new(0)  = zetaj_old(0) + std_met_zeta(0)*randn();
            sez(0)        = exp(zetaj_new(0));
            omegaj_new(0) = xmin + xrange*sez(0)/(1+sez(0));
          }  // End generation of zeta and omega
          // Get prior distributions fro zeta ~ N(zeta0,zeta_sigma2)
          resid_zeta_new(0) = (zetaj_new(0) - zeta0(0))/zeta_sigma(0);
          resid_zeta_old(0) = (zetaj_old(0) - zeta0(0))/zeta_sigma(0);
          sse_new           = accu(pow(resid_zeta_new(0), 2));
          sse_old           = accu(pow(resid_zeta_old(0), 2));
          testp             = testp - sse_new/2 + sse_old/2;  // Priors
          // Get function for group j: just S or multi-modal
          fpars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                     Rcpp::Named("iflagpn")     = iflagpn,
                                     Rcpp::Named("iflagCenter") = bflagCenter,
                                     Rcpp::Named("iflagZ")      = bflagZ,
                                     Rcpp::Named("theta")       = thetall.col(j),
                                     Rcpp::Named("phix")        = phixgrid,
                                     Rcpp::Named("delta")       = xdelta,
                                     Rcpp::Named("range")       = xrange,
                                     Rcpp::Named("xmin")        = xmin,
                                     Rcpp::Named("xgrid")       = xgrid,
                                     Rcpp::Named("nstat")       = nstat,
                                     Rcpp::Named("hfun")        = hfunall.col(j),
                                     Rcpp::Named("omega")       = omegaj_new,
                                     Rcpp::Named("idb")         = idball.col(j),
                                     Rcpp::Named("id_stat")     = id_stat,
                                     Rcpp::Named("mslope")      = mslopeall(j));
          
          fxgrid_new = rcpp_Getfx(fpars);
        
          //---------------------------------------------------
          // Compute fj(x_obs)
          // fxdata = GetUpfxdata(xgrid,xinx,xover,fxgrid)
          umat a     = xinx.rows(gidx);
          vec  b     = xover(gidx);
          fxdata_new = rcpp_GetSCfxobs(a, b, fxgrid_new);
          //---------------------------------------------------
          // Metropolis Test
          // Compute Likelihood
          // Adjust for ACE
          yresidj          = yresid.elem(gidx);
          yjstar           = yresidj;
          fxdata_star_new  = fxdata_new;
          fxdata_star_old  = fxdata_old;
          if (bflagACE) {
            rhoj = rhoall(j);
            srj  = 1- rhoj;
            nj   = nobs(j);
            // Adjust first observation for stationary variance
            yjstar(0)          *= srj;
            fxdata_star_new(0) *= srj;
            fxdata_star_old(0) *= srj;
            // AR model for observations 2 to n
            yjstar.tail(nj-1) = yjstar.tail(nj-1) - rhoj*yjstar.head(nj-1);
            fxdata_star_new.tail(nj-1) = fxdata_star_new.tail(nj-1) - rhoj*fxdata_star_new.head(nj-1);
            fxdata_star_old.tail(nj-1) = fxdata_star_old.tail(nj-1) - rhoj*fxdata_star_old.head(nj-1);
          } // End ACE adjustment
          // ******************************************************************
          // Add log-likelihood to testp
          // ******************************************************************
          resid_new = yjstar - fxdata_star_new;  // Note: fxdata_new_j only for group j
          sse_new   = accu(square(resid_new));
          resid_old = yjstar - fxdata_star_old;
          sse_old   = accu(square(resid_old));
          testp = testp - (sse_new - sse_old)/(2*sigma2(j)); // Likelihood
          if (log(randu()) < testp) { // Accept candidate
            fxgridall.col(j)  = fxgrid_new;
            fxdata.rows(gidx) = fxdata_new;
            omegall.col(j)    = omegaj_new(0);
            zetall.col(j)     = zetaj_new(0);
            // Fixed or random slope
            zeta_met.col(j).tail(zeta_met(0,j)) = met_var_zeta_new(0);
            zeta_met(1,j)++;
          } // End accept candidate  for
        } // End generating New U-shaped function
        
        // ********************************************************************
        // Generate New S or Multi-modal Parameters and Squish functions ----
        // ********************************************************************
        if (iflagsc > 4) {
          testp       = 0;
          fxdata_old  = fxdata.rows(gidx);
          fxgrid_old  = fxgridall.col(j);
          // HB model for log(psi)
          // Generate new psi_log
          if (bflagpsi) {
            psij_old      = psiall(j);
            psij_log_old  = psiall_log(j);
            // Random walk metropolis on zeta and psi_log
            met_var_psi_new = AdaptMetVar(psi_met.col(j));
            ck              = 5.66;   // Constant from Harito
            var_met_psi_new = ck*met_var_psi_new;
            std_met_psi     = sqrt(var_met_psi_new(0));
            psij_log_new    = psij_log_old + std_met_psi*randn();
            psij_new        = exp(psij_log_new);
            testp = testp - (pow(psij_log_new-psi_log_mu,2) - pow(psij_log_old-psi_log_mu,2))/(2*psi_log_sigma2);
          } else{
            //# Constant value for psi = psi_fixed
            psij_log_new   = psi_fixed_log;
            psij_new       = psi_fixed;
          }
          // End generate new psi_log
          //# ******************************************************************
          // Generate zeta and omega is ordered logit  
          zetaj_old        = zetall.col(j);
          omegaj_old       = omegall.col(j);
          met_var_zeta_new = AdaptMetVar(zeta_met.col(j));
          ck               = 5.66;   // Constant from Harito
          var_met_zeta     = ck*met_var_zeta_new;
          std_met_zeta     = sqrt(var_met_zeta);
          // 5% of time generate omega from uniform on xmin to xmax
          // 95% of time generate zeta from random walk metropolis
          if (randu()<.05) {
            // Randomly reinitialize zeta & omega to avoid getting stuck
            omegaj_new = sort(xmin + (xmax-xmin)*randu(nomega,1));
            // Find corresponding  zeta0_input for ordinal logit
            zvec      = (omegaj_new-xmin)/xrange;
            xvec      = zvec;
            if (nomega > 1) {
              xvec.tail(nomega-1) = zvec.tail(nomega-1) - zvec.head(nomega-1);
              x0 = 1-sum(xvec);
            }
            zetaj_new = log(xvec/x0);
          } else {
            if (nomega == 1) {
              zetaj_new(0)  = zetaj_old(0) + std_met_zeta(0)*randn();
              sez(0)        = sum(exp(zetaj_new(0)));
              omegaj_new(0) = xmin + xrange*sez(0)/(1+sez(0));
            } else {
              zetaj_new   = zetaj_old + std_met_zeta%randn(nomega);
              sez         = cumsum(exp(zetaj_new));
              omegaj_new  = xmin + xrange*sez/(1+sez(nomega-1));
            }
          }  // End generation of zeta and omega
          // Get prior distributions fro zeta ~ N(zeta0,zeta_sigma2)
          resid_new = (zetaj_new - zeta0)/zeta_sigma;
          resid_old = (zetaj_old - zeta0)/zeta_sigma;
          sse_new   = accu(square(resid_new));
          sse_old   = accu(square(resid_old));
          testp     = testp - sse_new/2 + sse_old/2;  // Priors
          // Find intervals based on inflection points
          vec_tmp    = zeros<vec>(id_flex.n_elem+2);
          vec_tmp(0) = xmin;
          vec_tmp.subvec(1,vec_tmp.n_elem-2) = omegaj_new(id_flex);
          vec_tmp(vec_tmp.n_elem-1) = xmax;
          
          breakers = sort(vec_tmp);
          breakers = floor(breakers/xdelta)*xdelta; // round to xgrid gaps
          uvec idb_new = zeros<uvec>(breakers.n_elem);
          // for (unsigned int jj = 0; jj < breakers.n_elem; jj++) {
          //   for (unsigned int jjj= 0; jjj < nint+1; jjj++) {
          //     if (xgrid(jjj)==breakers(jj)) idb_new(jj) = jjj; // find location in xgrid
          //   }
          // }
          for (unsigned int jj = 0; jj < breakers.n_elem; jj++) {
            colvec tmp  = abs(xgrid - breakers(jj));
            idb_new(jj) = tmp.index_min();
          }
          
          // Compute softmax square wave at inflection points
          hfunj_new = rcpp_GetSquish(omegaj_new(id_flex),psij_new,xgrid);
          // Get fJ(x) for theta_new
          // Get function for group j: just S or multi-modal
          fpars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                     Rcpp::Named("iflagpn")     = iflagpn,
                                     Rcpp::Named("iflagCenter") = bflagCenter,
                                     Rcpp::Named("iflagZ")      = bflagZ,
                                     Rcpp::Named("theta")       = thetall.col(j),
                                     Rcpp::Named("phix")        = phixgrid,
                                     Rcpp::Named("delta")       = xdelta,
                                     Rcpp::Named("range")       = xrange,
                                     Rcpp::Named("xmin")        = xmin,
                                     Rcpp::Named("xgrid")       = xgrid,
                                     Rcpp::Named("nstat")       = nstat,
                                     Rcpp::Named("hfun")        = hfunj_new,
                                     Rcpp::Named("omega")       = omegaj_new,
                                     Rcpp::Named("idb")         = idb_new,
                                     Rcpp::Named("id_stat")     = id_stat,
                                     Rcpp::Named("mslope")      = mslopeall(j));
          fxgrid_new = rcpp_Getfx(fpars);
          //---------------------------------------------------
          // Compute fj(x_obs)
          // fxdata = GetUpfxdata(xgrid,xinx,xover,fxgrid)
          umat a     = xinx.rows(gidx);
          vec  b     = xover(gidx);
          fxdata_new = rcpp_GetSCfxobs(a, b, fxgrid_new);
          
          //---------------------------------------------------
          // Metropolis Test
          // Compute Likelihood
          // Adjust for ACE
          yresidj          = yresid.elem(gidx);
          yjstar           = yresidj;
          fxdata_star_new  = fxdata_new;
          fxdata_star_old  = fxdata_old;
          if (bflagACE) {
            rhoj = rhoall(j);
            srj  = 1 - rhoj;
            nj   = nobs(j);
            // Adjust first observation for stationary variance
            yjstar(0)          *= srj;
            fxdata_star_new(0) *= srj;
            fxdata_star_old(0) *= srj;
            // AR model for observations 2 to n
            yjstar.tail(nj-1) = yjstar.tail(nj-1) - rhoj*yjstar.head(nj-1);
            fxdata_star_new.tail(nj-1) = fxdata_star_new.tail(nj-1) - rhoj*fxdata_star_new.head(nj-1);
            fxdata_star_old.tail(nj-1) = fxdata_star_old.tail(nj-1) - rhoj*fxdata_star_old.head(nj-1);
          } // End ACE adjustment
          // ******************************************************************
          // Add log-likelihood to testp
          // ******************************************************************
          resid_new = yjstar - fxdata_star_new;  // Note: fxdata_new_j only for group j
          sse_new   = accu(square(resid_new));
          resid_old = yjstar - fxdata_star_old;
          sse_old   = accu(square(resid_old));
          testp = testp - (sse_new - sse_old)/(2*sigma2(j)); // Likelihood
          if (log(randu()) < testp) { // Accept candidate
            hfunall.col(j)    = hfunj_new;
            fxgridall.col(j)  = fxgrid_new;
            fxdata.rows(gidx) = fxdata_new;
            omegall.col(j)    = omegaj_new;
            zetall.col(j)     = zetaj_new;
            idball.col(j)     = idb_new;
            // Fixed or random slope
            if (bflagpsi) {
              psiall(j)     = psij_new;
              psiall_log(j) = psij_log_new;
              psi_met.col(j).tail(psi_met(0,j)) = met_var_psi_new;
              psi_met(1,j)++;
            }  // End Fixed or random slope of squish function
            zeta_met.col(j).tail(zeta_met(0,j)) = met_var_zeta_new;
            zeta_met(1,j)++;
          } // End accept candidate  for 
        } // End generating new S or Multi-modal
        // ********************************************************************
        // Update adaptive metropolis mean and beta
        theta_met.col(j) = AdaptMetUpdate2(theta_met.col(j));
        if (iflagsc == 2 || iflagsc == 3) {
          mslope_log_met.col(j) = AdaptMetUpdate2(mslope_log_met.col(j));
        }
        
        if (iflagsc > 4) {
          if (bflagpsi) psi_met.col(j) = AdaptMetUpdate2(psi_met.col(j));
          zeta_met.col(j) = AdaptMetUpdate2(zeta_met.col(j));
        } // Generate Squish function
      } //  End loop over groups
    } // End shape constraints
    // ---------------------------------------------------------------
    //  Generate upper-level spectral coefficients
    // -----------------------------------------------
    //  Generate upper-level model theta_0k
    //  theta_jk ~ N(theta_0k,tau_k^2*eta_j^2*exp(-k*gamma_j)) for k>0
    //  theta_0k ~ N(0,eta0^2*(1+k/gamma_beta)^(-gamma_alpha))
    // Generate theta_0k.  Loop over k
    for (unsigned int k = 0; k < ntheta; k++) {
      // Get variance of theta_jk
      // Note: ktheta = c(rep(0,iflagZcnst),1:nbasis)
      // theta_jk ~ N(theta_0k,tau_k^2) if 1 <= iflagZcnst
      // theta_jk ~ N(theta_0k,tau_k^2*eta_j^2*exp(-k*gamma_j) if k > iflagZcnst
      thvk = tau2all(k)*exp(-kall(k)*gammall);  // Var(theta_jk) for j = 1, ..., J
      
      if (iflagZcnst == 1 && !bflagZ) {
        // Got intercept, which is positive with theta_00 = exp(xi)
        xi_vn       = 1/(sum(1/thvk) + 1/th0v(k)); // variance
        if (k == 0) {
          xi_bn     = xi_vn*accu(xiparall)/tau2all(k); // mean
          xipar0    = xi_bn + sqrt(xi_vn)*randn();
          theta0(0) = exp(xipar0);
        } else {
          xi_bn     = xi_vn*accu(thetall.row(k).t()/thvk); // mean
          theta0(k) = xi_bn + sqrt(xi_vn)*randn();
        } // End got intercept
      } else {
        // No intercepts
        xi_vn       = 1/(sum(1/thvk) + 1/th0v(k)); // variance
        xi_bn       = xi_vn*accu(thetall.row(k).t()/thvk); // mean
        theta0(k)   = xi_bn + sqrt(xi_vn)*randn();
      }
    } // End loop to generate theta_0k
    
    // Generate upper-level model for slope parameter Convex/concave monotone
    if((iflagsc==2)|(iflagsc==3)){
      vni                  = ingroup/mslope_log_sigma2 + mslope_log_v0i;
      vn                   = 1/vni;
      slope_bn             = vn*(sum(mslopeall_log)/mslope_log_sigma2 + mslope_log_v0im0);
      mslope0_log          = slope_bn + sqrt(vn)*randn();
      mslope0              = exp(mslope0_log);
      // Variance
      mresid               = mslopeall_log - mslope0_log;
      mslope_log_sigma2_sn = mslope_log_sigma2_s0 + accu(square(mresid));
      mslope_log_sigma2    = mslope_log_sigma2_sn/(2*randg(distr_param(mslope_log_simga2_rn/2, 1.0)));
      mslope_log_sigma     = sqrt(mslope_log_sigma2);
    }
    
    //--------------------------------------------------------
    // Generate upper-level model for squish function parameters
    if (iflagsc>=4) {
      // Generate zeta0
      zeta_vn2 = 1/(ingroup/zeta_sigma2 + 1/zeta_v02);
      zeta_vn  = sqrt(zeta_vn2);
      zeta_mn  = zeta_vn2%(sum(zetall,1)/zeta_sigma2 + zeta_m0/zeta_v02);
      zeta0    = zeta_mn + zeta_vn%randn(nomega);
      sez      = cumsum(exp(zeta0));
      omega0   = xmin + xrange*sez/(1+sez(nomega-1));
      // ********************************************************************
      // Generate zeta_sigma2
      for (unsigned int k = 0; k < nomega; k++) {
        resid_zeta     = zetall.row(k) - zeta0(k);
        zeta_sn        = zeta_s0 + accu(square(resid_zeta));
        zeta_sigma2(k) = zeta_sn/(2*randg(distr_param(zeta_rn/2, 1.0)));
      } // End generate zeta_sigma2
      // ********************************************************************
      // Generate psi for upper-level model
      if (iflagsc>4 && bflagpsi) {
        // Generate psi_log_mu
        psi_log_v2  = 1/(ingroup/psi_log_sigma2 + 1/psi_log_v02);
        psi_log_v   = sqrt(psi_log_v2);
        psi_log_mn  = psi_log_v2*(accu(psiall_log)/psi_log_sigma2 + psi_log_m0/psi_log_v02);
        psi_log_mu  = psi_log_mn + psi_log_v*randn();
        // Generate psi_log_sigma2
        resid          = psiall_log - psi_log_mu;
        sse            = accu(square(resid));
        psi_log_sn     = psi_log_s0 + sse;
        psi_log_sigma2 = psi_log_sn/(2*randg(distr_param(psi_log_rn/2, 1.0)));
        psi_log_sigma  = sqrt(psi_log_sigma2);
        
        psi0 = exp(psi_log_mu);
      }  // End Generate psi0
    }  // End generate upper-level squish parameters
    //--------------------------------------------------------------------------
    // Compute upper-level f0
    if (iflagsc==0) {
      f0xgrid  = phixgrid * theta0;
    } else {  // Shape constraints: Need to do bias correction
      f0xgrid  = mean(fxgridall,1);
    } // End compute f0
    
    //-----------------------------------------------------------------
    // End nonparametric model
    //-----------------------------------------------------
    ymean  = valpha + wbeta + fxdata;  // mean of Y
    yresid = ydata - ymean;
    
    // Generate tau2all for HB Smoothing Prior
    tk = 0;   // Index to handle theta_j0 = exp(xi_j) with shape constraints
    // If shape constraint, theta_j0 = exp(xipar_j) & xipar_j ~ N(xipar_0,tau_0^2)
    if (iflagZcnst == 1) {
      resid_tau  = xiparall - xipar0;
      resid2_tau = square(resid_tau);
      tau2_sn    = tau2_s0 + accu(resid2_tau);
      tau2all(0) = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
      tk         = 1;   //  Increment index to skip tau2_0
    }
    for (unsigned int k=tk; k<ntheta; k++) {
      resid_tau  = thetall.row(k).t() - theta0(k);
      resid2_tau = square(resid_tau);
      gamk       = exp(-kall(k)*gammall);
      tau2_sn    = tau2_s0 + accu(resid2_tau/gamk);
      tau2all(k) = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
    }
    // Get std dev
    tauall = sqrt(tau2all);
    //--------------------------------------------------------------------
    // Generate eta0^2 where theta_0k ~ N(0,eta0^2*(1+kall/gamma_beta)^(-gamma_alpha))
    if (iflagZcnst == 0) {
      resid2_eta = square(theta0)/gam0vec;
    } else {
      resid2_eta = square(theta0.tail(theta0.n_elem-1))/gam0vec.tail(gam0vec.n_elem-1);
    }
    
    sse        = accu(resid2_eta);
    eta02_sn   = eta02_s0 + sse;
    eta02      = eta02_sn/(2*randg(distr_param(eta02_rn/2, 1.0)));
    eta0       = sqrt(eta02);
    th0v       = eta02*gam0vec;
    
    //--------------------------------------------------------------------
    // Generate gamma for each population
    tbot = iflagZcnst;
    for (unsigned int j=0; j<ngroup; j++) {  // Loop over groups to generate gammaj with slice sampling
      gammaj    = gammall(j);
      gamvec    = exp(-kall.subvec(tbot,ntheta-1)*gammaj);
      resid     = thetall(span(tbot,ntheta-1),j) - theta0.subvec(tbot,ntheta-1);
      resid2    = square(resid);
      colvec ck = resid2/(2*tau2all.subvec(tbot,ntheta-1));
      
      // Worry about ck = 0
      z         = find(ck == 0);   // Find index where ck == 0
      zi        = find(ck != 0);   // Find index where ck == 0
      nz        = z.n_elem;
      if (nz == nbasis) {  // all theta's are zeros!
        // Reset gamma_j 
        gammaj     = 1;
        gammall(j) = gammaj;
      } else {
        
        if(nz>0) ck(z) = ones<colvec>(nz);   // Set zeros to 1 for the time being
        u1   = randu(nbasis, 1);
        bn   = gammaj + (log(ck - log(u1)%gamvec) - log(ck))/kall.subvec(tbot,ntheta-1);
        if (nz>0) bn = bn(zi);   // drop the z's that are 0
        
        // Slice sampling for gamma^(alpha-1)
        u2   = randu();
        bmin = gammaj*(pow(u2, 1/(gamma_alpha - 1)));
        colvec bg = zeros<colvec>(bn.n_elem+1);
        bg.head(bn.n_elem) = bn;
        bg.tail(1) = gmax;
        bmax = min(bg);  // Sometimes gamma wanders off to a large number.  Limit it to 5
        if(bmin>bmax){
          gammaj  = 1.0;
        }else{
          u2   = randu();
          gpar = wk - gamma_beta;
          gammaj = bmax + log(u2 + (1.0-u2)*exp((bmin-bmax)*gpar))/gpar;
        }  
        gammall(j) = gammaj;
        
      } // End if ck all zeros
    } // End generate gammaj for each poulation
    // Summary stats used in generating hyperparameters
    gamma_sum  = accu(gammall);
    lgamma_sum = accu(log(gammall));
    
    //------------------------------------------------------------------------
    // Generate alpha and beta for gamma_j ~ G(alpha,beta)
    // Generate gamma_mu = gamma_alpha/gamma_beta and gamma_sigma2 = gamma_mu/gamma_beta
    // from truncated normal distribution
    kbot             = iflagZcnst;
    met_var_new      = AdaptMetVar(gamma_met); // Variance for Metropolis
    var_met_g0       = 5.66*met_var_new(0);
    std_met_g0       = sqrt(var_met_g0);      // STD DEV for Metropolis
    
    gamma_mu_new     = rcpp_rndtnab(gamma_mu,std_met_g0,mu_bot,mu_top);
    gamma_sigma_new  = rcpp_rndtnab(gamma_sigma,std_met_g0,sigma_bot,sigma_top);
    gamma_sigma2_new = pow(gamma_sigma_new,2);
    gamma_beta_new   = gamma_mu_new/gamma_sigma2_new;
    gamma_alpha_new  = gamma_beta_new*gamma_mu_new;
    gamma_prob_new   = R::pgamma(gmax,gamma_alpha_new,1/gamma_beta_new, 1, 0);
    
    // Likelihood for gamma_j ~ G(alpha,beta)I(gamma_j < gmax)
    testpa = 0;
    testpa = testpa - ingroup*(log(gamma_prob_new) - log(gamma_prob));  // Normalizing constant
    testpa = testpa + ingroup*gamma_alpha_new*log(gamma_beta_new);
    testpa = testpa - ingroup*gamma_alpha*log(gamma_beta);
    testpa = testpa - ingroup*lgamma(gamma_alpha_new);
    testpa = testpa + ingroup*lgamma(gamma_alpha);
    testpa = testpa + (gamma_alpha_new - gamma_alpha)*lgamma_sum;
    testpa = testpa - (gamma_beta_new  - gamma_beta)*gamma_sum;
    // Likelihood theta_0j ~ N(0,eta0^2*(1+k/gamma_beta)^(-gamma_alpha))
    gam0vec_new    = pow(1+kall/gamma_beta_new, -gamma_alpha_new);
    th0v_new       = eta02*gam0vec_new;
    theta02        = square(theta0.subvec(kbot, ntheta-1));
    testpa         = testpa - accu(theta02/(th0v_new.subvec(kbot, ntheta-1)))/2 + accu(theta02/(th0v.subvec(kbot, ntheta-1)))/2;
    testpa         = testpa - accu(log(th0v_new.subvec(kbot, ntheta-1)))/2 + accu(log(th0v.subvec(kbot, ntheta-1)))/2;
    
    // Prior for gamma_mu is truncated normal
    testpa = testpa - pow(gamma_mu_new - gamma_mu_m0, 2)/(2*gamma_mu_v02);
    testpa = testpa + pow(gamma_mu     - gamma_mu_m0, 2)/(2*gamma_mu_v02);
    // Prior for sigma is truncated normal
    testpa = testpa - pow(gamma_sigma_new - gamma_sigma_m0, 2)/(2*gamma_sigma_v02);
    testpa = testpa + pow(gamma_sigma     - gamma_sigma_m0, 2)/(2*gamma_sigma_v02);
    
    // Truncated random walk for mu
    pnew = normcdf(mu_top,gamma_mu,std_met_g0)     - normcdf(mu_bot,gamma_mu,std_met_g0);      // Normalizing constant for gamma_mu_new
    pold = normcdf(mu_top,gamma_mu_new,std_met_g0) - normcdf(mu_bot,gamma_mu_new,std_met_g0);  // Normalizing constant for gamma_mu_old
    testpa = testpa + log(pnew) - log(pold);
    // Truncated random walk for sigma
    pnew = normcdf(sigma_top,gamma_sigma,std_met_g0)     - normcdf(sigma_bot,gamma_sigma,std_met_g0);            // Normalizing constant for gamma_sigma_new
    pold = normcdf(sigma_top,gamma_sigma_new,std_met_g0) - normcdf(sigma_bot,gamma_sigma_new,std_met_g0);    // Normalizing constant for gamma_sigma_old
    testpa = testpa + log(pnew) - log(pold);  
    if(log(randu())<testpa){
      gamma_alpha  = gamma_alpha_new;
      gamma_beta   = gamma_beta_new;
      gamma_mu     = gamma_mu_new;
      gamma_sigma  = gamma_sigma_new;
      //gamma_sigma2 = gamma_sigma2_new;
      gamma_prob   = gamma_prob_new;
      gam0vec      = gam0vec_new;
      th0v         = th0v_new;
      //gamma0       = gamma_mu;
      gamma_met(1)++; // pmet
      gamma_met.tail(gamma_met(0)) = met_var_new;
    }  // End generate gamma_alpha and gamma_beta
    // End HB Smoothing Prior
    AdaptMetUpdate(gamma_met);
    
    //-------------------------------------------------------
    // If probit model, generate latent ydata where ydata0 are observed 0/1
    if (iflaglm == 1) {
      ydata(id_probit1) = rcpp_rndtnb2(ymean(id_probit1),1.0,0.0);  // Y > 0
      ydata(id_probit0) = rcpp_rndtna2(ymean(id_probit0),1.0,0.0);  // Y < 0
    }
    
    // If ordinal data, generate latents ydata & cutpoints.  
    // y0data are observed data
    if (iflaglm == 2) {
      // Loop over groups
      for (unsigned int k=0; k<ngroup; k++) {
        uvec gidx   = iptHB1[k];
        yk          = ydata(gidx);    // latents for group k
        yok         = yodata(gidx);   // Observed data
        mk          = ymean(gidx);    // means for group k
        sk          = sigma(k);
        cutt        = cuttall.col(k);   // cutpoints 
        id_cut_botk = id_cut_bot(gidx); // indicies in cutt for lower bound
        id_cut_topk = id_cut_top(gidx); // indicies in cutt for upper bound
        ybotk       = cutt(id_cut_botk);
        ytopk       = cutt(id_cut_topk);
        ord_fa      = zeros<colvec>(yk.n_elem);
        for (unsigned int i=0; i<yk.n_elem; i++) {
          ord_fa(i) = normcdf(ybotk(i),mk(i),sk);
        }
        ord_fb      = zeros<colvec>(yk.n_elem);
        for (unsigned int i=0; i<yk.n_elem; i++) {
          ord_fb(i) = normcdf(ytopk(i),mk(i),sk);
        }
        ord_p       = ord_fa + (ord_fb-ord_fa)%randu<vec>(yk.n_elem);
        // Check for extreme values of p
        idz         = find(ord_p > 0.99999*ones<colvec>(ord_p.n_elem));
        if (idz.n_elem>0) ord_p(idz) = 0.99999*ones<colvec>(idz.n_elem);
        idz         = find(ord_p < 0.00001*ones<colvec>(ord_p.n_elem));
        if (idz.n_elem>0) ord_p(idz) = 0.00001*ones<colvec>(idz.n_elem);
        // Generate new latent values
        for (unsigned int i=0; i<yk.n_elem; i++) {
          yk(i)     = R::qnorm(ord_p(i), mk(i), sk, true, false);  
        }
        // Check bounds on yk
        idz         = find(yk > ytopk);
        if (idz.n_elem>0) yk(idz)  = ytopk(idz) - 0.00001*ones<colvec>(idz.n_elem);
        idz         = find(yk < ybotk);
        if (idz.n_elem>0) yk(idz) = ybotk(idz) + 0.00001*ones<colvec>(idz.n_elem);
        // Store it
        ydata(gidx) = yk;
        // Generate free cutpoints
        for (unsigned int j=0; j<id_cut_free.n_elem; j++) {
          jk    = id_cut_free(j);   // ID of free cutpoints
          cutjk = cutt(jk);         // Current value of free cutpoint
          zjbot = find(yk < cutjk); // IDs for latents less than cutjk
          zjtop = find(yk > cutjk); // IDs for latents greater than cutjk
          if (zjbot.n_elem> 0) {
            vec_tmp.zeros(zjbot.n_elem+1);
            vec_tmp.head(zjbot.n_elem) = yk(zjbot);
            vec_tmp.tail(1) = cutt(jk-1);
            cbot = max(vec_tmp);
          } else {
            cbot = cutt(jk-1);
          }
          if (zjtop.n_elem>0) {
            vec_tmp.zeros(zjtop.n_elem+1);
            vec_tmp.head(zjtop.n_elem) = yk(zjtop);
            vec_tmp.tail(1) = cutt(jk+1);
            ctop = min(vec_tmp);
          } else {
            ctop = cutt(jk+1);
          }
          cutjk    = cbot + (ctop-cbot)*randu();  // Generate cutpoint
          cutt(jk) = cutjk;                       // store it
        } // End generate free cutpoints
        // Store results
        cuttall.col(k) = cutt;
        ydata(gidx)    = yk;
        ybot(gidx)     = ybotk;
        ytop(gidx)     = ytopk;
      }  // End loop over groups
    }  // End generate latents and cutpoints for ordinal data
    
    // In pre-burnin period, adjust mean for adaptive Metropolis
    if (iflagsc>0) {
      // shape constrained: modify metropolis
      if ((imcmc < nblow0*maxmodmet) && (imcmc == floor((double)imcmc/(double)nblow0)*nblow0) && imcmc > 0) {
        for (unsigned int j=0; j<ngroup; j++) {
          theta_met.col(j) = UpdateMet2(theta_met.col(j),nblow0,10.0);
          if (iflagsc == 2 || iflagsc == 3) {
            mslope_log_met.col(j) = UpdateMet2(mslope_log_met.col(j),nblow0,10.0);
          }
          // Got Squish parameters
          if (iflagsc>=4) {
            if (bflagpsi) psi_met.col(j) = UpdateMet2(psi_met.col(j),nblow0,10.0);
            zeta_met.col(j) = UpdateMet2(zeta_met.col(j),nblow0,10.0);
          } // End got Squish Parameters
        } // End loop over groups
        
        if (bflagHBsigma) UpdateMet(sigma2_met,nblow0,10.0);
        UpdateMet(gamma_met,nblow0,10.0);
      }// End test in adjustment period
      // -------------------------------------------------------
    } // Got shape constraints
    //-------------------------------------------------------
    //-------------------------------------------------------
    
    // Save MCMC iterations
    if ((imcmc >= nblow+nblow0*maxmodmet) && (imcmc % nskip == 0)) {
      
      if (nparv > 0) alphag.row(isave) = alpha.t();
      if (!bflagHBsigma) {
        sigmag.row(isave)    = sigma(0);
      } else {
        sigmag.row(isave)    = sigma.t();
        sigma2_alphag(isave) = sigma2_alpha;
        sigma2_betag(isave)  = sigma2_beta;
      }
      
      tauallg.row(isave)       = tauall.t();
      eta0g(isave)             = eta0;
      gamma_mugg(isave)        = gamma_mu;
      gamma_alphag(isave)      = gamma_alpha;
      gamma_betag(isave)       = gamma_beta;
      gammallg.row(isave)      = gammall.t();
      thetam                  += thetall;
      thetas                  += square(thetall);
      if (save_large_mcmc) {
        thetallg.slice(isave)  = thetall;
      }
      theta0g.row(isave)       = theta0.t();
      fxdatam                 += fxdata;
      fxdatas                 += square(fxdata);
      fxgridm                 += fxgridall;
      fxgrids                 += square(fxgridall);
      if (save_large_mcmc) {
        fxgridallg.slice(isave) = fxgridall;
      }
      f0xgridg.row(isave)      = f0xgrid.t();
      f0xgridm                += f0xgrid;
      f0xgrids                += square(f0xgrid);
      phig.row(isave)          = phi_vec.t();
      lambdag.row(isave)       = vectorise(lambda).t();
      betam                   += betall;
      betas                   += square(betall);
      betallg.slice(isave)     = betall;
      
      if (bflagACE) {
        rhoallg.row(isave) = rhoall.t();
      }
      
      // Save latents for probit or ordinal probit
      if (iflaglm >= 1) {
        ydatam += ydata;
      }
      
      if (iflaglm == 2) {
        cuttallm += cuttall;
        cuttalls += square(cuttall);
      }
      
      if ((iflagsc==2)||(iflagsc==3)) {
        mslope_logg.row(isave)   = mslopeall_log.t();
        mslope0_logg(isave)      = mslope0_log;
        mslopeg.row(isave)       = mslopeall.t();
        mslope0g(isave)          = mslope0;
        mslope_log_sigmag(isave) = mslope_log_sigma;
      } // Save slope parameters for monotone convex/concave
      
      // U-shaped has squish function
      if (iflagsc >= 4) {
        zetag.slice(isave)  = zetall;
        omegag.slice(isave) = omegall;
        zeta0g.row(isave)   = zeta0.t();
        omega0g.row(isave)  = omega0.t();
        if (iflagsc>4 && bflagpsi) {
          psiallg.row(isave) = psiall.t();
          psi0g(isave)       = psi0;
        }  // Estimated psi for squish function
      }  // Got squish function
      
      isave++;
    } // End save MCMC
  } // end of mcmc loop
  
  cout << "MCMC is done!" << endl;
  
  
  // Compute summary statistics
  sigma2_met(1) /= (nskip*smcmc);
  gamma_met(1)  /= (nskip*smcmc);
  for (unsigned int j=0; j<ngroup; j++) {
    theta_met(1,j) /= (nskip*smcmc);
    if(iflagsc>=4){
      zeta_met(1,j) /= (nskip*smcmc); 
      if (bflagpsi) {
        psi_met(1,j) /= (nskip*smcmc);
      }
    }
  }
  
  Rcpp::List metg = Rcpp::List::create(Rcpp::Named("sigma2_met")     = sigma2_met,
                                       Rcpp::Named("gamma_met")      = gamma_met,
                                       Rcpp::Named("theta_met")      = theta_met,
                                       Rcpp::Named("zeta_met")       = zeta_met,
                                       Rcpp::Named("psi_met")        = psi_met,
                                       Rcpp::Named("mslope_log_met") = mslope_log_met);
  
  Rcpp::List mcmcg = Rcpp::List::create(Rcpp::Named("alphag")        = alphag,
                                        Rcpp::Named("sigmag")        = sigmag,
                                        Rcpp::Named("sigma2_alphag") = sigma2_alphag,
                                        Rcpp::Named("sigma2_betag")  = sigma2_betag,
                                        Rcpp::Named("tauallg")       = tauallg,
                                        Rcpp::Named("eta0g")         = eta0g,
                                        Rcpp::Named("gamma_mugg")    = gamma_mugg,
                                        Rcpp::Named("gammallg")      = gammallg,
                                        Rcpp::Named("theta0g")       = theta0g,
                                        Rcpp::Named("phig")          = phig,
                                        Rcpp::Named("lambdag")       = lambdag);
  
  // maximum element of list is 20.
  return Rcpp::List::create(Rcpp::Named("metg")     = metg,           //  1
                            Rcpp::Named("fxgridm")  = fxgridm,        //  2
                            Rcpp::Named("fxgrids")  = fxgrids,        //  3
                            Rcpp::Named("f0xgridm") = f0xgridm,       //  4
                            Rcpp::Named("f0xgrids") = f0xgrids,       //  5
                            Rcpp::Named("fxdatam")  = fxdatam,        //  6
                            Rcpp::Named("fxdatas")  = fxdatas,        //  7
                            Rcpp::Named("betam")    = betam,          //  8
                            Rcpp::Named("betas")    = betas,          //  9
                            Rcpp::Named("thetam")   = thetam,         // 10
                            Rcpp::Named("thetas")   = thetas,         // 11
                            Rcpp::Named("cuttallm") = cuttallm,       // 12
                            Rcpp::Named("cuttalls") = cuttalls,       // 13
                            Rcpp::Named("mcmcg")    = mcmcg);         // 14
  
}
