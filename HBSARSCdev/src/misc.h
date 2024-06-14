#include <Rmath.h>
#include <RcppArmadillo.h>

#ifndef MISC
#define MISC

arma::colvec rcpp_CumTrap(arma::colvec f, double delta);
arma::mat rcpp_CumInt2D(arma::mat fmat, arma::colvec xdelta);
double rcpp_trapint(arma::colvec f, double xdelta);
double rcpp_Int2D(arma::mat fmat, arma::colvec xdelta);
arma::colvec rcpp_GetSquish(arma::colvec omega, double psi, arma::colvec xgrid);

double rcpp_rndtna(double mu, double sigma, double xtop);
double rcpp_rndtnb(double mu, double sigma, double xbot);
arma::colvec rcpp_rndtna2(arma::colvec mu, double sigma, double xbot);
arma::colvec rcpp_rndtnb2(arma::colvec mu, double sigma, double xbot);
double rcpp_rndtnab(double mu, double sigma, double xbot, double xtop);
double rcpp_rndtnab(arma::colvec& mu, double sigma, double xbot, double xtop);
arma::colvec rcpp_rndtnab(arma::colvec& mu, arma::mat& sigma, double xbot, double xtop);

arma::colvec rcpp_Getfx(Rcpp::List fpars);
arma::colvec rcpp_GetMonofxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetMonoSamefxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetMonoDiffxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetUfxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetSfxgrid(Rcpp::List fpars, int d3);
arma::colvec rcpp_GetUseg(arma::colvec x, arma::colvec gz, double omega, double xdelta);
arma::colvec rcpp_GetMMfxgrid(Rcpp::List fpars);

arma::colvec rcpp_GetSCfxobs(arma::umat& xinx, arma::vec& xover, arma::colvec& fxgrid);

arma::colvec rcpp_cholSol(arma::colvec& b, arma::mat& U);

// metropolis
arma::colvec GetMetVec(double w, double m0, double alpha, unsigned int ndim);
void SetMetVec(arma::colvec &Met, double w, double m0, double alpha, unsigned int ndim);
arma::colvec AdaptMetVar(arma::colvec x);
void AdaptMetUpdate(arma::colvec &x);
arma::colvec AdaptMetUpdate2(arma::colvec x);
void UpdateMet(arma::colvec &x, unsigned int nb, double fct);
arma::colvec UpdateMet2(arma::colvec x, unsigned int nb, double fct);

arma::colvec rcpp_getsmooth(double tau2, double gampar, arma::colvec& egamkall, arma::colvec theta2, 
                            double wk, double tau2_s0, double tau2_rn, arma::colvec freq_sum);
Rcpp::List rcpp_getsmoothSC(double tau2, double gampar, arma::colvec theta, 
                            arma::colvec theta2, arma::colvec freqs, 
                            arma::colvec egamkall, double wk, double tau2_s0, 
                            double tau2_rn, bool bflagZ);

arma::colvec rcpp_fdata_apx(arma::mat fmat, arma::mat f12apx_weights, arma::umat xinx, arma::mat xover, unsigned int nint);

void rcpp_get_phixgrid_2D(arma::mat& phixgrid, arma::mat& xgrid, arma::colvec& xrange,
                          arma::colvec& xmin, arma::mat& freq_all, double phix0,
                          double phiy0, unsigned int ntheta, unsigned int nint);

double GetSmoothGamma(double gampar,arma::colvec theta, arma::colvec theta0,
                      arma::colvec tau2, double gamma_alpha, double gamma_beta, 
                      double wk, arma::colvec kall, double gmax);

Rcpp::List GetGammaHyper(double gamma_alpha, double gamma_beta, double gamma_mu, 
                         double gamma_sigma, double gamma_sigma2, arma::colvec gammall, 
                         arma::colvec theta0, double eta02, double gamma_sum, 
                         double lgamma_sum, arma::colvec kall, arma::colvec gamma_met,
                         double mu_bot, double mu_top, double sigma_bot, double sigma_top,
                         double gamma_mu_m0, double gamma_mu_v02, double gamma_sigma_m0, double gamma_sigma_v02,
                         int ingroup, arma::colvec freq_sum, double gmax);
arma::colvec rcpp_GetMainfxdata(arma::uvec xinx , arma::colvec xover, arma::colvec fxgrid);

arma::mat get_phixgrid2D(arma::mat xgrid, double phix0, double phiy0, arma::mat freq_all, arma::vec xmin,  
                         arma::vec xrange, unsigned int ntheta, unsigned int nint);
#endif
