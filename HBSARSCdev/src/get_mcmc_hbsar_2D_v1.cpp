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
Rcpp::List get_mcmc_2DHBSAR(Rcpp::List& data_all,
                            arma::vec& nobs,
                            Rcpp::List& iptHB1,
                            arma::mat& ztz,
                            arma::cube& vtv,
                            arma::cube& wtw,
                            arma::cube& phi2,
                            
                            arma::vec& cpara,
                            Rcpp::List& gamma_para,
                            Rcpp::List& tau_para,
                            Rcpp::List& eta_para,
                            Rcpp::List& met_para,
                            arma::vec& sigma_para,
                            Rcpp::List& v_para,
                            Rcpp::List& z_para,
                            Rcpp::List& phi_para,
                            Rcpp::List& fx_para,
                            Rcpp::List& basis_para,
                            
                            arma::colvec& alpha,
                            arma::mat& betall,
                            arma::mat& phi,
                            arma::mat& lambda,
                            arma::mat& lambdai,
                            arma::mat& thetall,
                            arma::colvec& theta0,
                            arma::colvec& gammall,
                            arma::colvec& tau2all,
                            arma::colvec& tauall,
                            arma::colvec& sigma2,
                            arma::colvec& sigma,
                            
                            arma::cube& betallg,
                            arma::cube& thetallg) {
  
  
  
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
  
  bool bflagHBsigma = cpara(16);
  
  int ingroup = ngroup;
  
  colvec ydata  = data_all["ydata"];
  colvec valpha = data_all["valpha"];
  colvec wbeta  = data_all["wbeta"];
  mat vdata     = data_all["vdata"];
  mat wdata     = data_all["wdata"];
  mat zdata     = data_all["zdata"];
  
  colvec gamma_met   = met_para["gamma_met"];
  colvec sigma2_met  = met_para["sigma2_met"];
  
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
  colvec freq_sum = basis_para["freq_sum"];
  
  mat phixgrid = phi_para["phixgrid"];
  mat phixdata = phi_para["phixdata"];
  
  cube fxgridall = fx_para["fxgridall"];
  mat f0xgrid    = fx_para["f0xgrid"];
  colvec fxdata  = fx_para["fxdata"];
  colvec fxdatam = fxdata;
  
  mat thetam     = zeros<mat>(ntheta,ngroup);  
  mat thetas     = zeros<mat>(ntheta,ngroup);
  mat betam      = zeros<mat>(size(betall));
  mat betas      = zeros<mat>(size(betall));
  
  unsigned int ngrid = (nint+1)*(nint+1);
  colvec f0xgrid_vec = zeros<colvec>(ngrid);
  colvec fxgridj_vec = zeros<colvec>(ngrid);
  
  cube fxgridm   = zeros<cube>(nint+1,nint+1,ngroup);
  cube fxgrids   = zeros<cube>(nint+1,nint+1,ngroup);
  mat fxgrid_old = zeros<mat>(nint+1,nint+1);
  mat f0xgridm   = zeros<mat>(nint+1,nint+1);
  mat f0xgrids   = zeros<mat>(nint+1,nint+1);
  mat f0Bxgridm  = zeros<mat>(nint+1,nint+1);
  mat f0Bxgrids  = zeros<mat>(nint+1,nint+1);
  colvec fxdatas = zeros<colvec>(ntot);
  
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
  
  
  // local
  colvec yresid = zeros<colvec>(ntot);
  colvec ymean  = zeros<colvec>(ntot);
  mat mat1 = zeros<mat>(1,nint+1);
  mat mat2 = zeros<mat>(nint+1,1);
  colvec yresidj, rj, v, rVec, bn, bj;
  mat rMat, ranMat, vbi, vbi12, vbni, vbni12;
  double sse, ck;
  
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
  double gammaj, gpar, nz;
  uvec zc, zvi, zv, z, zi;
  colvec resid_new, resid_old, fxgrid_new, thvk, fxj, theta0k, fxdata_new,
  fxdata_old, theta_old, thetak_old, gamvec, ths, met_var_new, met_var,
  met_std;
  colvec theta_new  = zeros<colvec>(ntheta);
  
  
  unsigned int kbot, tbot;
  colvec vpar, thvb;
  
  // theta0
  colvec t0var;
  double xi_vn, xi_bn;
  
  colvec resid, resid2, gamk;
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
  
  // main
  
  
  
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
      
      if (bflagHBsigma) {
        // HB error variance
        sigma2_rn = sigma2_alpha + nobs(j);
        sigma2_sn = sigma2_beta  + accu(square(rj));
        s2        = sigma2_sn / (2 * randg(distr_param(sigma2_rn / 2, 1.0)));
        sigma2(j) = s2;
        sigma(j)  = sqrt(s2);
      } // End generate sigma2 in HB model
    } // End Loop over groups to generate beta_j and sigma2_j
    // HB model for sigma if not probit
    if (bflagHBsigma) {
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
        sigma2_met(1)++; // pmet
        sigma2_met.tail(sigma2_met(0)) = met_var_new;
      }
      AdaptMetUpdate(sigma2_met);
    }
    
    yresid = ydata - valpha - wbeta - fxdata;
    
    // Generate sigma2 for homogenous error variance
    if (!bflagHBsigma) {
      //Probit or else
      // Generate homogeneous sigma
      sigma2_sn = sigma2_s0 + sum(square(yresid));
      sigma2    = sigma2_sn/(2 * randg(distr_param(sigma2_rn / 2, 1.0))) * ones<colvec>(ngroup);
      sigma     = sqrt(sigma2);
    } // End generate homogeneous sigma
    
    
    //------------------------------------------
    // Generate Phi (Beta = Zdata*Phi + delta)
    phi_vbni   = kron(lambdai, ztz) + phi_v0i;
    phi_vbni12 = chol(phi_vbni);
    ztb        = zdata.t() * betall;
    for (unsigned int j = 0; j < nparw; j++) {
      if (nparw > 1) {
        a = ztb * lambdai.col(j);
      } else {
        a = ztb * lambdai(0, 0);
      }
      
      phi_bn.subvec(nparz*j, nparz*(j+1)-1) = a;
      
    }
    rMat = randn(phidim, 1);
    phi_vec  = solve(phi_vbni, phi_bn + phi_vbni12.t()*rMat);
    for (unsigned int i = 0; i < nparw; i++) {
      phi.col(i) = phi_vec.subvec(nparz*i, nparz*(i+1)-1);
    }
    
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
    
    //------------------------------------------
    // Do HB nonparametric model
    yresid = ydata - valpha - wbeta;   // Residuals without fxdata
    // Free model: do not need numerical integration
    for (unsigned int j = 0; j < ngroup; j++) { // Loop over groups to generate theta_j
      uvec gidx = iptHB1[j];
      rj    = yresid.elem(gidx);
      phixj = phixdata.rows(gidx);
      phi2j = phi2.slice(j); // Get phi_j'phi_j
      
      thetaj  = thetall.col(j);
      gammaj  = gammall(j);
      gamjvec = exp(-freq_sum*gammaj);
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
        colvec eigval = zeros<colvec>(ntheta);
        colvec eigvec = zeros<colvec>(ntheta,ntheta);
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
      fxgridj_vec        = phixgrid * thetaj;
      fxgridall.slice(j) = reshape(fxgridj_vec, nint+1,nint+1);  // fj computed on xgrid
      phixj              = phixdata.rows(gidx);
      fxj                = phixj * thetaj;     // fj computed at observations
      fxdata.rows(gidx)  = fxj;
    } // End loop over groups to generate theta_j
    
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
      // theta_jk ~ N(theta_0k,tau_k^2*eta_j^2*exp(-k*gamma_j) 
      thvk = tau2all(k)*exp(-freq_sum(k)*gammall);  // Var(theta_jk) for j = 1, ..., J
      // No intercepts
      xi_vn     = 1/(sum(1/thvk) + 1/th0v(k)); // variance
      xi_bn     = xi_vn*accu(thetall.row(k).t()/thvk); // mean
      theta0(k) = xi_bn + sqrt(xi_vn)*randn();
    } // End loop to generate theta_0k
    
    //--------------------------------------------------------------------------
    // Compute upper-level f0
    f0xgrid_vec = phixgrid * theta0;
    f0xgrid     = reshape(f0xgrid_vec, nint+1,nint+1);
    
    //-----------------------------------------------------------------
    // End nonparametric model
    //-----------------------------------------------------
    
    // Generate tau2all for HB Smoothing Prior
    for (unsigned int k=0; k<ntheta; k++) {
      resid_tau  = thetall.row(k).t() - theta0(k);
      resid2_tau = square(resid_tau);
      gamk       = exp(-freq_sum(k)*gammall);
      tau2_sn    = tau2_s0 + accu(resid2_tau/gamk);
      tau2all(k) = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
    }
    // Get std dev
    tauall = sqrt(tau2all);
    
    
    //--------------------------------------------------------------------
    // Generate eta0^2 where theta_0k ~ N(0,eta0^2*(1+kall/gamma_beta)^(-gamma_alpha))
    resid2_eta = square(theta0)/gam0vec;
    sse        = accu(resid2_eta);
    eta02_sn   = eta02_s0 + sse;
    eta02      = eta02_sn/(2*randg(distr_param(eta02_rn/2, 1.0)));
    eta0       = sqrt(eta02);
    th0v       = eta02*gam0vec;
    
    
    //--------------------------------------------------------------------
    // Generate gamma for each population
    tbot = 0;
    for (unsigned int j=0; j<ngroup; j++) {  // Loop over groups to generate gammaj with slice sampling
      gammaj    = gammall(j);
      gamvec    = exp(-freq_sum.subvec(tbot,ntheta-1)*gammaj);
      resid     = thetall(span(tbot,ntheta-1),j) - theta0.subvec(tbot,ntheta-1);
      resid2    = square(resid);
      colvec ck = resid2/(2*tau2all.subvec(tbot,ntheta-1));
      
      // Worry about ck = 0
      z         = find(ck == 0);   // Find index where ck == 0
      zi        = find(ck != 0);   // Find index where ck == 0
      nz        = z.n_elem;
      if (nz == ntheta) {  // all theta's are zeros!
        // Reset gamma_j 
        gammaj     = 1;
        gammall(j) = gammaj;
      } else {
        
        if(nz>0) ck(z) = ones<colvec>(nz);   // Set zeros to 1 for the time being
        u1   = randu(ntheta, 1);
        bn   = gammaj + (log(ck - log(u1)%gamvec) - log(ck))/freq_sum.subvec(tbot,ntheta-1);
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
    kbot             = 0;
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
    gam0vec_new    = pow(1+freq_sum/gamma_beta_new, -gamma_alpha_new);
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
      thetallg.slice(isave)    = thetall;
      theta0g.row(isave)       = theta0.t();
      fxdatam                 += fxdata;
      fxdatas                 += square(fxdata);
      fxgridm                 += fxgridall;
      fxgrids                 += square(fxgridall);
      f0xgridm                += f0xgrid;
      f0xgrids                += square(f0xgrid);
      phig.row(isave)          = phi_vec.t();
      lambdag.row(isave)       = vectorise(lambda).t();
      betam                   += betall;
      betas                   += square(betall);
      betallg.slice(isave)     = betall;
      
      isave++;
    } // End save MCMC
  } // end of mcmc loop
  
  cout << "MCMC is done!" << endl;
  
  
  // Compute summary statistics
  sigma2_met(1) /= (nskip*smcmc);
  gamma_met(1)  /= (nskip*smcmc);
  
  
  Rcpp::List metg = Rcpp::List::create(Rcpp::Named("sigma2_met")     = sigma2_met,
                                       Rcpp::Named("gamma_met")      = gamma_met);
  
  
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
                            Rcpp::Named("mcmcg")    = mcmcg);         // 12
  
}
