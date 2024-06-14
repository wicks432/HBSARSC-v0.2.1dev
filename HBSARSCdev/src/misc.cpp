#include "misc.h"
#define pi 3.141592653589793238462643383279502884197

using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// computes cumulative integral of function
// [[Rcpp::export]]
arma::colvec rcpp_CumTrap(arma::colvec f, double delta) {
  unsigned int nr = f.n_elem;
  colvec fcum = arma::zeros<colvec> (nr);
  if (nr > 1) {
    fcum.subvec(1,nr-1) = cumsum(f.subvec(0,nr-2)) + cumsum(f.subvec(1,nr-1));
    fcum *= (delta/2); 
  }
  
  return fcum;
}


// trapidzoid integration
// Int f(x)*dx = xdelta*(f1/2 + f2 + ... + fn-1 + fn/2)
// [[Rcpp::export]]
double rcpp_trapint(arma::colvec f, double xdelta) {
  unsigned int nr = f.n_elem;
  double a = sum(f);
  a -= ((f[0]+f[nr-1])/2);
  return a*xdelta;
}

/* Int2D
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
#*********************************************/
// [[Rcpp::export]]
double rcpp_Int2D(arma::mat fmat, arma::colvec xdelta) {
  unsigned int ny = fmat.n_cols;
  unsigned int nx = fmat.n_rows;
  double fint = 0.0;
  mat f1 = zeros<mat>(nx,  ny-1);
  mat f2 = zeros<mat>(nx-1,ny-1);
  // Compute average of the 4 fmat values at each cell
  f1    = fmat.head_cols(ny-1) + fmat.tail_cols(ny-1);
  f2    = f1.head_rows(nx-1)   + f1.tail_rows(nx-1);
  fint  = accu(f2)/4;
  fint *= (xdelta(0)*xdelta(1));  // Multiply by length and width of grids.
  
  return fint;
}


/*----------------------------------------------------
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
#*****************************************************/
// [[Rcpp::export]]
arma::mat rcpp_CumInt2D(arma::mat fmat, arma::colvec xdelta) {
  unsigned int ny = fmat.n_cols;
  unsigned int nx = fmat.n_rows;
  mat fint = zeros<mat>(nx,  ny);
  mat f1   = zeros<mat>(nx,  ny-1);
  mat f2   = zeros<mat>(nx-1,ny-1);
  // Compute average of the 4 fmat values at each cell
  // Add 0 to first row and first column
  f1    = fmat.head_cols(ny-1) + fmat.tail_cols(ny-1);
  f2    = f1.head_rows(nx-1)   + f1.tail_rows(nx-1);
  fint(span(1,nx-1), span(1,ny-1)) = f2;
  fint /= 4;
  // Cumsum volumns of blocks in rows and columns
  fint  = cumsum(fint.t(), 0).t();// CAUTION  apply & cumsum dumps results in columns
  fint  = cumsum(fint, 0);
  fint *= (xdelta(0)*xdelta(1)); // Multiply by length and width of grids.
  return fint;
}


// return -idx
// [[Rcpp::export]]
arma::uvec get_ex_idx(unsigned int len, arma::uvec idx) {
  if (len <= idx.n_elem) {
    return NULL;
  }
  
  uvec v = zeros<uvec>(len);
  for (unsigned int i=0; i<idx.n_elem; i++) {
    v(idx(i)) = 1;
  }
  
  uvec ex_vec (len - idx.n_elem);
  int j=0;
  for (unsigned int i=0; i<len; i++) {
    if (v(i) == 0) {
      ex_vec(j) = i;
      j++;
    }
  }
  
  return ex_vec;
}

// [[Rcpp::export]]
arma::uvec get_seq(unsigned int a, unsigned int b) {
  unsigned int len = b-a+1;
  uvec seq = zeros<uvec>(len);
  for (unsigned int i=0; i<len; i++) {
    seq(i) = a + i;
  }
  
  return seq;
}


// generates form a truncated above normal distribution
// [[Rcpp::export]]
double rcpp_rndtna(double mu, double sigma, double xtop) {
  double u = randu();
  double fa = 0.0;
  double fb = normcdf(xtop, mu, sigma);
  double p = fa + u*(fb-fa);
  // Worry about p = 0 or p = 1
  if (p < 0.00001) p = 0.00001;
  if (p > 0.99999) p = 0.99999;
  double x = R::qnorm(p, mu, sigma, true, false);
  
  if (x > xtop) x = xtop;
  
  return x;
}


// generates form a blow normal distribution
// [[Rcpp::export]]
arma::colvec rcpp_rndtna2(arma::colvec mu, double sigma, double xtop) {
  unsigned int n  = mu.n_elem;
  colvec fa = zeros<colvec> (n);
  colvec x  = zeros<colvec> (n);
  
  arma::colvec u = arma::randu(n);
  
  for (unsigned int i=0; i<n; i++) {
    double fa = 0.0;
    double fb = normcdf(xtop, mu(i), sigma);
    double p = fa + u(i)*(fb-fa);
    // Worry about p = 0 or p = 1
    if (p < 0.00001) p = 0.00001;
    if (p > 0.99999) p = 0.99999;
    x(i) = R::qnorm(p, mu(i), sigma, true, false);
    
    if (x(i) > xtop) x(i) = xtop;
    
  }
  
  return x;
}



// generates form a blow normal distribution
// [[Rcpp::export]]
double rcpp_rndtnb(double mu, double sigma, double xbot) {
  double u = arma::randu();
  double fa = arma::normcdf(xbot, mu, sigma);
  double fb = 1.0;
  double p = fa + u*(fb-fa);
  // Worry about p = 0 or p = 1
  if (p < 0.00001) p = 0.00001;
  if (p > 0.99999) p = 0.99999;
  double x = R::qnorm(p, mu, sigma, true, false);
  
  if (x < xbot) x = xbot;
  
  return x;
}


// generates form a blow normal distribution
// [[Rcpp::export]]
arma::colvec rcpp_rndtnb2(arma::colvec mu, double sigma, double xbot) {
  unsigned int n  = mu.n_elem;
  colvec fa = zeros<colvec> (n);
  colvec x  = zeros<colvec> (n);
  
  arma::colvec u = arma::randu(n);
  
  for (unsigned int i=0; i<n; i++) {
    double fa = normcdf(xbot, mu(i), sigma);
    double fb = 1.0;
    double p = fa + u(i)*(fb-fa);
    // Worry about p = 0 or p = 1
    if (p < 0.00001) p = 0.00001;
    if (p > 0.99999) p = 0.99999;
    x(i) = R::qnorm(p, mu(i), sigma, true, false);
    
    if (x(i) < xbot) x(i) = xbot;
    
  }
  
  
  return x;
}


// generates form a truncated normal distribution
// [[Rcpp::export]]
double rcpp_rndtnab(double mu, double sigma, double xbot, double xtop) {
  double u = R::runif(0.0,1.0); //randu();
  double fa = normcdf(xbot, mu, sigma);
  double fb = normcdf(xtop, mu, sigma);
  // Need to gaurd against bad falues for fa and fb
  double x0 = (xbot+xtop)/2; // Good guess that stasifies xbot < x < xtop
  double delta = fa < fb ? fb-fa : 0;
  double p = fa + u*delta;
  // Worry about p = 0 or p = 1
  if (p < 0.00001) p = 0.00001;
  if (p > 0.99999) p = 0.99999;
  double x = R::qnorm(p, mu, sigma, true, false);
  if (fa >= fb) x = x0;
  if (x < xbot) x = xbot;
  if (x > xtop) x = xtop;
  
  return x;
}


/***************************************************
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
***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetSquish(arma::colvec omega, double psi, arma::colvec xgrid) {
  unsigned int nint  = xgrid.n_elem-1;
  unsigned int nflex = omega.n_elem;
  double deven = nflex%2 == 0 ? 1.0 : 0.0;
  colvec hfun  = deven*ones<colvec>(nint+1);
  colvec z     = zeros<colvec>(nint+1);
  colvec hj    = zeros<colvec>(nint+1);
  // Loop over inflection points
  for	(unsigned int j=0; j<nflex; j++) {
    z = psi*(xgrid-omega(j));
    uvec idz1 = find(z < -50);
    uvec idz2 = find(z >  50);
    
    if (idz1.n_elem > 0) z.elem(idz1) = -50*ones<colvec>(idz1.n_elem);  
    if (idz2.n_elem > 0) z.elem(idz2) =  50*ones<colvec>(idz2.n_elem);  
    
    z  = exp(z);
    hj = (1-z)/(1+z);
    
    hfun += pow(-1,j)*hj;
  }
  
  return hfun;
}


/***************************************************
# GetMonofxgrid
#   Compute increasing f = int_a^x exp[Z(s)] ds on xgrid
# INPUT
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetMonofxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool bflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec dfx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    dfx = exp(z);
  } else {
    dfx = square(z);
  }
  colvec fx = iflagpn*rcpp_CumTrap(dfx,delta);
  if(bflagCenter) { // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


/***************************************************
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
# ***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetMonoSamefxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  colvec xgrid     = fpars["xgrid"];
  double mslope    = fpars["mslope"];
  double xmin      = fpars["xmin"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool bflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    d2fx = exp(z);
  } else {
    d2fx = square(z);
  }
  colvec dfx = rcpp_CumTrap(d2fx,delta);
  colvec fx  = rcpp_CumTrap(dfx, delta) + (xgrid-xmin)*mslope;
  fx *= iflagpn;
  if (bflagCenter) { // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}



/***************************************************
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
# ***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetMonoDiffxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  colvec xgrid     = fpars["xgrid"];
  double mslope    = fpars["mslope"];
  double xmin      = fpars["xmin"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool bflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    d2fx = exp(z);
  } else {
    d2fx = square(z);
  }
  colvec dfx = reverse(rcpp_CumTrap(reverse(d2fx),delta));
  colvec fx  = rcpp_CumTrap(dfx, delta) + (xgrid-xmin)*mslope;
  fx *= iflagpn;
  if (bflagCenter) { // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


/***************************************************
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
# ***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetSfxgrid(Rcpp::List fpars, int d3) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  colvec hfun      = fpars["hfun"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool bflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec xgrid     = fpars["xgrid"];
  double xmin      = fpars["xmin"];
  colvec d2fx;
  colvec dfx;
  colvec fx;
  double c2;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    d2fx = exp(z);
  } else {
    d2fx = square(z);
  }
  // d3 =  1 -> hfun is "squish down" for "S" 
  // d3 = -1 -> hfun is "squish up" for "rotated S"
  d2fx = d3*hfun%d2fx;
  dfx  = rcpp_CumTrap(d2fx, delta);
  c2   = std::min(0.0, dfx.min());
  fx   = rcpp_CumTrap(dfx,delta);
  fx   = iflagpn*(fx - c2*(xgrid-xmin));
  
  if (bflagCenter) { // center f
    // Compute mean 
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


/***************************************************
# Get convex, U-shaped or concave, inverted U shaped on segement x
# INPUT
#    x     = grid between neighboring inflection points
#    f     = h(x)*g(Z(x)) computed on x
#    omega = stationary point
# ***************************************************/
 // [[Rcpp::export]]
 arma::colvec rcpp_GetUseg(arma::colvec x, arma::colvec gz, double omega, double xdelta) {
   colvec fx;
   colvec zg1;
   double tmp;
   uvec idbw;
   
   fx   = rcpp_CumTrap(gz, xdelta);
   fx   = rcpp_CumTrap(fx, xdelta);
   idbw = find(x <= omega);         // indicies for x <= omega
   zg1  = gz.elem(idbw);            // g(Z) from xmin to omega
   tmp  = rcpp_trapint(zg1,xdelta); // int_a^omega g(Z)
   fx   = (fx - tmp*(x-x(0)));
   
   return fx;
}


/***************************************************
# GetMMfxgrid
#   Multi-modal   functions
#   iflagpn =  1 U, inverted U, ... or min, max, min
#   iflagpn = -1 Inverted U, U, ... or max, min, max
# INPUT
#    Input 
#       fpars is a list
#   Output
#       fx             = function computed on xgrid
# ***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetMMfxgrid(Rcpp::List fpars) {
  mat phix           = fpars["phix"];
  colvec theta       = fpars["theta"];
  colvec hfun        = fpars["hfun"];
  colvec omega       = fpars["omega"];
  int iflagpn        = fpars["iflagpn"];
  double delta       = fpars["delta"];
  bool bflagCenter   = fpars["iflagCenter"];
  double range       = fpars["range"];
  bool bflagZ        = fpars["iflagZ"];
  colvec xgrid       = fpars["xgrid"];
  uvec id_stat       = fpars["id_stat"];
  uvec idb           = fpars["idb"];
  unsigned int nstat = fpars["nstat"];
  
  unsigned int id1 = 1;
  unsigned int id2 = 1;
  
  colvec d2fx;
  colvec dfx;
  uvec idk;
  colvec xk;
  colvec d2fxk;
  colvec fxk;
  colvec fxk_tmp;
  colvec z       = phix * theta;
  colvec omega_s = omega.elem(id_stat);
  if (bflagZ) {
    d2fx  = hfun%exp(z);
  } else {
    d2fx  = hfun%square(z);
  }
  colvec fx = zeros<colvec>(xgrid.n_elem); // Note: f(xmin) = 0.
  
  for	(unsigned int k=0; k<nstat; k++) {
    idk     = get_seq(idb(k),idb(k+1));
    xk      = xgrid.elem(idk);
    d2fxk   = d2fx.elem(idk);
    fxk_tmp = rcpp_GetUseg(xk,d2fxk,omega_s(k),delta);
    id2     = id1 + fxk_tmp.n_elem-2;
    if (fxk_tmp.n_elem>0) {
      if (fxk_tmp.n_elem==1) {
        fx(id1) = fxk_tmp(0);
        id2 = id1;
      } else {
        fxk                = fxk_tmp.tail(fxk_tmp.n_elem-1) + fx(id1-1);
        fx.subvec(id1,id2) = fxk;
        id1 = id2+1;
      }
    }
  } // End get f on segments
  
  fx *= -iflagpn;
  
  if (bflagCenter) { // center f
    // Compute mean 
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


/***************************************************
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
# ***************************************************/
// [[Rcpp::export]]
arma::colvec rcpp_GetUfxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  colvec hfun      = fpars["hfun"];
  colvec xgrid     = fpars["xgrid"];
  colvec omegavec  = fpars["omega"];
  double xmin      = fpars["xmin"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool bflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec dfx;
  colvec fx;
  double c;
  double zintaw;
  double omega = omegavec(0);
  uvec idbw;
  colvec zg1;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    dfx = exp(z);
  } else {
    dfx = square(z);
  }
  
  idbw   = find(xgrid<=omega);      // indiciec for x <= omega
  zg1    = dfx.elem(idbw);          // g(Z) from xmin to omega
  zintaw = rcpp_trapint(zg1,delta); // int_a^omega g(Z)
  fx     = rcpp_CumTrap(dfx,delta);
  fx     = rcpp_CumTrap(fx,delta);
  fx     = -iflagpn*(fx - zintaw*(xgrid-xmin));
  
  if (bflagCenter) { // center f
    // Compute mean 
    c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_Getfx(Rcpp::List fpars) {
  colvec fx;
  
  int fshape = fpars["iflagsc"];
  
  switch (fshape)
  {
  case 0: // free function
  {
    mat phix     = fpars["phix"];
    colvec theta = fpars["theta"];
    fx = phix * theta;
    break;
  }
  case 1: // monotone function
  {  
    fx = rcpp_GetMonofxgrid(fpars);
    break;
  }
  case 2: // Increasing convex or decreasing and concave
  {  
    fx = rcpp_GetMonoSamefxgrid(fpars);
    break;
  }
  case 3: // Increasing and concave or decreasing and convex
  {
    fx = rcpp_GetMonoDiffxgrid(fpars);
    break;
  }
  case 4: // U-shaped function
  {
    fx = rcpp_GetUfxgrid(fpars);
    break;
  }
  case 5: // S-shaped (increasing,convex-to-concave) or (decreasing,concave-to-convex)
  {
    fx = rcpp_GetSfxgrid(fpars,1);
    break;
  }
  case 6: // S-shaped (increasing,concave-to-convex) or (decreasing,convex-to-concave)
  {
    fx = rcpp_GetSfxgrid(fpars,-1);
    break;
  }
  case 7: // multiple extremes
  {
    fx = rcpp_GetMMfxgrid(fpars);
    break;
  }
  default:
    std::cout << "not available function type" << std::endl;
    break;
  }
  
  return fx;
}



// Compute f at observations
// [[Rcpp::export]]
arma::colvec rcpp_GetSCfxobs(arma::umat& xinx, arma::vec& xover, arma::colvec& fxgrid) {
  unsigned int n = xinx.n_elem;
  colvec fxobs   = zeros<colvec>(n);
  
  colvec f0 = fxgrid.elem(xinx.col(0));
  colvec f1 = fxgrid.elem(xinx.col(1));
  
  fxobs = f0 + (f1-f0)%xover;
  
  return fxobs;
}


// [[Rcpp::export]]
arma::colvec rcpp_cholSol (arma::colvec& b, arma::mat& U) {
  return (solve(trimatu(U), solve(trimatl(U.t()), b)));
}



// w=0.5, m0.001, alpha=4.0, ndim=1
// [[Rcpp::export]]
arma::colvec GetMetVec(double w, double m0, double alpha, unsigned int ndim) {
  double beta    = m0*(alpha-1);
  int pmet       = 0;
  int icount     = 0;
  colvec beta_AM = beta*ones<colvec>(ndim);
  colvec m_AM    = m0*ones<colvec>(ndim);
  colvec met_var = m0*ones<colvec>(ndim);
  colvec Met     = ones<colvec>(7+3*ndim);
  
  Met(0) = ndim;
  Met(1) = pmet;
  Met(2) = icount;
  Met(3) = alpha;
  Met(4) = beta;
  Met(5) = w;
  Met(6) = m0;
  Met.subvec(7,          7 +   ndim-1) = beta_AM;
  Met.subvec(7 +   ndim, 7 + 2*ndim-1) = m_AM;
  Met.subvec(7 + 2*ndim, 7 + 3*ndim-1) = met_var;
    
  return Met;
}


// w=0.5, m0.001, alpha=4.0, ndim=1
// [[Rcpp::export]]
void SetMetVec(arma::colvec &Met, double w, double m0, double alpha, unsigned int ndim) {
  double beta    = m0*(alpha-1);
  int pmet       = 0;
  int icount     = 0;
  colvec beta_AM = beta*ones<colvec>(ndim);
  colvec m_AM    = m0*ones<colvec>(ndim);
  colvec met_var = m0*ones<colvec>(ndim);
  
  Met(0) = ndim;
  Met(1) = pmet;
  Met(2) = icount;
  Met(3) = alpha;
  Met(4) = beta;
  Met(5) = w;
  Met(6) = m0;
  Met.subvec(7,          7 +   ndim-1) = beta_AM;
  Met.subvec(7 +   ndim, 7 + 2*ndim-1) = m_AM;
  Met.subvec(7 + 2*ndim, 7 + 3*ndim-1) = met_var;
}



// [[Rcpp::export]]
arma::colvec AdaptMetVar(arma::colvec x) {
  unsigned int ndim = x(0);
  colvec beta_AM    = x.subvec(7, 7+ndim-1);
  double alpha      = x(3);
  //colvec met_var    = beta_AM/randg(ndim, distr_param(alpha, 1.0));
  colvec met_var = zeros<colvec>(ndim);
  for (unsigned int i = 0; i < ndim; i++) {
    met_var(i) = beta_AM(i)/R::rgamma(alpha, 1.0);
  }
      
  return met_var;
}


// [[Rcpp::export]]
void AdaptMetUpdate(arma::colvec &x) {
  unsigned int ndim = x(0);
  int icount        = x(2);
  colvec m_AM       = x.subvec(7 +   ndim, 7 + 2*ndim-1);
  colvec met_var    = x.subvec(7 + 2*ndim, 7 + 3*ndim-1);
  double w          = x(5);
  double alpha      = x(3);
  colvec mnew;
  
  x(2)                        = icount+1;
  x.subvec(7+ndim,7+2*ndim-1) = m_AM + (met_var - m_AM)/x(2); // met_var
  mnew                        = w*met_var + (1-w)*x.subvec(7+ndim,7+2*ndim-1);
  x.subvec(7,7+ndim-1)        = (alpha-1)*mnew; // beta_AM
}


// [[Rcpp::export]]
arma::colvec AdaptMetUpdate2(arma::colvec x) {
  colvec Met = x;
  unsigned int ndim = Met(0);
  int icount        = Met(2);
  colvec m_AM       = Met.subvec(7 +   ndim, 7 + 2*ndim-1);
  colvec met_var    = Met.subvec(7 + 2*ndim, 7 + 3*ndim-1);
  double w          = Met(5);
  double alpha      = Met(3);
  colvec mnew;
  
  Met(2)                        = icount+1;
  Met.subvec(7+ndim,7+2*ndim-1) = m_AM + (met_var - m_AM)/Met(2); // met_var
  mnew                          = w*met_var + (1-w)*Met.subvec(7+ndim,7+2*ndim-1);
  Met.subvec(7,7+ndim-1)        = (alpha-1)*mnew; // beta_AM
  return Met;
}


// fct=10
// [[Rcpp::export]]
void UpdateMet(arma::colvec& x, unsigned int nb, double fct) {
  double       pmet = x(1);
  unsigned int ndim = x(0);
  double       m0   = x(6);
  
  pmet /= nb;
  if (pmet > 0.6) { // Too few acceptances, decrease metm 
    m0 *= fct;
    SetMetVec(x, 0.5, m0, 4.0, ndim);
  } else if (pmet < 0.3) {  // Too few acceptances, decrease metm 
    m0 /= fct;
    if (m0 < 1e-8) m0 = 1e-8;
    SetMetVec(x, 0.5, m0, 4.0, ndim);
  } else {
    x(1) = 0;
  }
}


// fct=10
// [[Rcpp::export]]
arma::colvec UpdateMet2(arma::colvec x, unsigned int nb, double fct) {
  double       pmet = x(1);
  unsigned int ndim = x(0);
  double       m0   = x(6);
  colvec Met = x;

  pmet /= nb;
  if (pmet > 0.6) { // Too few acceptances, decrease metm 
    m0 *= fct;
    SetMetVec(Met, 0.5, m0, 4.0, ndim);
  } else if (pmet < 0.3) {  // Too few acceptances, decrease metm 
    m0 /= fct;
    if (m0 < 1e-8) m0 = 1e-8;
    SetMetVec(Met, 0.5, m0, 4.0, ndim);
  } else {
    Met(1) = 0;
  }
  
  return Met;
}


// [[Rcpp::export]]
arma::colvec rcpp_getsmooth(double tau2, double gampar, arma::colvec& egamkall, arma::colvec theta2, 
                           double wk, double tau2_s0, double tau2_rn, arma::colvec freq_sum) {
  
  unsigned int nb, ntheta0;
  double tau2_sn, tau, bmin, u2;
  colvec ck, u1, bn;
  colvec out = zeros<colvec>(3);
  uvec idtheta0, iz0;
  
  // Full conditional of tau2 is IG(rn/2,sn/2)
  // rgamma(n, shape, rate = 1, scale = 1/rate)
  // The Gamma distribution with parameters shape = a and scale = s has density
  // f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
  nb      = egamkall.n_elem;
  tau2_sn = tau2_s0 + accu(theta2/egamkall); // Posterior rate parameters
  tau2    = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
  tau     = sqrt(tau2);
  
  // Generate gamma with slice sampling
  ck        = theta2/(2*tau2);
  // Worry about theta_k = 0.  Big hassel, but it matters because we take logs
  idtheta0  = find(ck == 0);
  ntheta0   = idtheta0.n_elem;
  if (ntheta0 == nb) {
    gampar  = .1; // all rows = 0. punt!
  } else {
    if (ntheta0>0) ck(idtheta0) = ones<colvec>(ntheta0);  // set component to 1
    u1      = randu(nb,1); //Uniform rvs
    bn 		  = gampar + (log(ck - log(u1)%egamkall) - log(ck))/freq_sum;
    
    bmin    = min(bn);
    u2		  = randu();
    gampar  = bmin + log(u2 + (1-u2)*exp(-wk*bmin))/wk;
  }
  
  egamkall  = exp(-gampar*freq_sum);
  iz0       = find(egamkall < 1e-15);   // Guard against numerical errors
  if (iz0.n_elem>0) egamkall(iz0) = 1e-15*ones<colvec>(iz0.n_elem);
  
  out(0) = tau2;
  out(1) = tau;
  out(2) = gampar;
  
  return out;
}


// [[Rcpp::export]]
Rcpp::List rcpp_getsmoothSC(double tau2, double gampar, arma::colvec theta, 
                            arma::colvec theta2, arma::colvec freqs, 
                            arma::colvec egamkall, double wk, double tau2_s0, 
                            double tau2_rn, bool bflagZ) {
  
  unsigned int nb, ntheta0;
  double tau2_sn, xipar, bmin, u2;
  colvec ck, u1, bn, theta2_drop, egamkall_drop;
  uvec idtheta0, iz0;
  
  // Full conditional of tau2 is IG(rn/2,sn/2)
  // rgamma(n, shape, rate = 1, scale = 1/rate)
  // The Gamma distribution with parameters shape = a and scale = s has density
  // f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
  nb      = egamkall.n_elem;
  if (!bflagZ) {
    xipar = log(theta(0));
    theta2(0) = pow(xipar,2);
  }
  
  tau2_sn = tau2_s0 + accu(theta2/egamkall); // Posterior rate parameters
  tau2    = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
  
  // Generate gamma with slice sampling
  // Drop theta[1], which is the constant
  theta2_drop   = theta2.tail(nb-1);
  egamkall_drop = egamkall.tail(nb-1);
  
  ck        = theta2_drop/(2*tau2);
  // Worry about theta_k = 0.  Big hassel, but it matters because we take logs
  idtheta0  = find(ck == 0);
  ntheta0   = idtheta0.n_elem;
  if (ntheta0 == nb) {
    gampar  = .1; // all rows = 0. punt!
  } else {
    if (ntheta0>0) ck(idtheta0) = ones<colvec>(ntheta0);  // set component to 1
    u1      = randu(nb-1,1); //Uniform rvs
    bn 		  = gampar + (log(ck - log(u1)%egamkall_drop) - log(ck))/freqs.tail(nb-1);
    
    bmin    = min(bn);
    u2  		= randu();
    gampar  = bmin + log(u2 + (1-u2)*exp(-wk*bmin))/wk;
  }
  
  egamkall_drop  = exp(-gampar*freqs);
  iz0 = find(egamkall_drop < 1e-15);   // Guard against numerical errors
  if (iz0.n_elem>0) egamkall_drop(iz0) = 1e-15*ones<colvec>(iz0.n_elem);
  
  
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("tau2")     = tau2,
                                      Rcpp::Named("gampar")   = gampar,
                                      Rcpp::Named("egamkall") = egamkall_drop);
  
  return out;
  
}


// ******************************************************************************
// ******************************************************************************
// fdata_apx 
//   Approximate f at data by using f computed on xgrid
//   Weighted average of corner points for grid cell that data are in
// INPUT
//   fmat        = f computed on xgrid
//   fapx_weight = weights for the four corners based
//   xinx        = identifies indices of grid where data belongs
// RETURN
//   fdata       = approxmation of f at data
// CALL
//   fdata  = fdata_apx(fmat,f12apx_weights,xinx)
// ******************************************************************************
// [[Rcpp::export]]
arma::colvec rcpp_fdata_apx(arma::mat fmat, arma::mat f12apx_weights, arma::umat xinx, arma::mat xover, unsigned int nint) {
  unsigned int nobs = xinx.n_rows;       // number of observations
  unsigned int i1, i2;
  double fi;
  // Approximate fdata_anal by using grid
  colvec fdata = zeros<colvec>(nobs);
  for	(unsigned int i=0; i<nobs; i++) {
    i1 = xinx(i,0); // Grid index for x1
    i2 = xinx(i,1); // Grid index for x2
    fi = 0;         // Weighted average of neighbors
    // Not on right & not on top boundary
    if ((i1 < nint) && (i2 < nint)){
      fi = fmat(i1,i2)*f12apx_weights(i,0) + fmat(i1+1,i2)*f12apx_weights(i,1) + fmat(i1,i2+1)*f12apx_weights(i,2) + fmat(i1+1,i2+1)*f12apx_weights(i,3);
    }
    // Observation on right boundary and not top
    if ((i1 == nint)&&(i2 < nint)) {
      fi = fmat(i1,i2)*(1-xover(i,1)) + fmat(i1,i2+1)*xover(i,1);
    }
    // Observation on top boundary and not right
    if ((i1 < nint)&&(i2 == nint)) {
      fi = fmat(i1,i2)*(1-xover(i,0)) + fmat(i1+1,i2)*xover(i,0);
    } 
    //Observation on right & top boundary: fmat is exact
    if (i1==nint && i2==nint) fi = fmat(i1,i2);
    fdata(i) = fi;
  }
  
  return fdata;
} // End fdata_apx to approximate fdata with f at grid


// [[Rcpp::export]]
void rcpp_get_phixgrid_2D(arma::mat& phixgrid, arma::mat& xgrid, arma::colvec& xrange,
                          arma::colvec& xmin, arma::mat& freq_all, double phix0,
                          double phiy0, unsigned int ntheta, unsigned int nint) {
  
  double xi, yj, phixi, phiyj;
  int k1, k2;
  
  unsigned int icount = 0;
  for	(unsigned int j=0; j<nint+1; j++) {
    yj = xgrid(j,1);
    for	(unsigned int i=0; i<nint+1; i++) {
      xi = xgrid(i,0);
      // Loop over indices where j+k <= kmax
      for	(unsigned int kk=0; kk<ntheta; kk++) {
        k1 = freq_all(kk,0);   // Frequency for xi
        k2 = freq_all(kk,1);   // Frequency for yj
        if (k1==0) {
          phixi  = phix0;
        } else {
          phixi  = sqrt(2/xrange(0))*cos(k1*pi*(xi-xmin(0))/xrange(0));
        }
        if (k2==0) {
          phiyj  = phiy0;
        } else {
          phiyj  = sqrt(2/xrange(1))*cos(k2*pi*(yj-xmin(1))/xrange(1));
        }
        phixgrid(icount,kk) = phixi*phiyj;
      }  // end loop over indices
      icount++;
    }  // End loop over xi
  } // End lover over yj
  
  return;
}


//-------------------------------------------------------------------
// GetSmoothGamma
//    Generate gamma_j smoothing parameter for population j
//    theta_{j,k} ~ N(theta_{0,k},tau_k^2*exp(-k*gamma_j)
//    gamma_j     ~ G(gamma_alpha,gamma_beta)
//    Use slice sampling
//-------------------------------------------------------------------
// [[Rcpp::export]]
double GetSmoothGamma(double gampar, arma::colvec theta, arma::colvec theta0,
                      arma::colvec tau2, double gamma_alpha, double gamma_beta, 
                      double wk, arma::colvec kall, double gmax) {
  
  colvec gamvec, resid2, ck, u1, bn, bg;
  uvec id0, z, zi;
  double u2, bmin, bmax, gpar;
  unsigned int nz;
  
  unsigned int ntheta = theta.n_elem;
  gamvec = exp(-kall*gampar);
  // Check gamvecj for rounding to 0
  id0 = find(gamvec < 1e-10);
  if (id0.n_elem > 0) gamvec(id0) = 1e-10*ones<colvec>(id0.n_elem);
  resid2     = square(theta - theta0);
  ck         = resid2/(2*tau2);
  // Worry about ck = 0
  z          = find(ck == 0);   // Find index where ck == 0
  zi         = find(ck != 0);
  nz         = z.n_elem;
  if (nz == ntheta) {  // all theta's are zeros!
    // Reset gamma_j 
    gampar   = 1.0;
  } else {
    if (nz>0) ck(z).ones(); // Set zeros to 1 for the time being
    u1 = randu(ntheta);
    bn = gampar + (log(ck - log(u1)%gamvec) - log(ck))/kall;
    if (nz>0) bn = bn(zi);   // drop the z's that are 0
    
    // Slice sampling for gamma^(alpha-1)
    u2   = randu();
    bmin = gampar*pow(u2, 1/(gamma_alpha - 1));
    bg = zeros<colvec>(bn.n_elem+1);
    bg.head(bn.n_elem) = bn;
    bg.tail(1) = gmax;
    bmax = min(bg);  // Sometimes gamma wanders off to a large number.  Limit it to gmax
    if (bmin>bmax) {
      gampar = 1.0;
    } else {
      u2     = randu();
      gpar   = wk - gamma_beta;
      gampar = bmax + log(u2 + (1.0-u2)*exp((bmin-bmax)*gpar))/gpar;
    }  
  } // End if ck all zeros
  return gampar;
}


//-----------------------------------------------------------------------
// GetGammaHyper
//   Generate hyper parameters for heterogeneity distribution of gamma
//   gamma_j  ~ G(alpha,beta)
//   theta_0k ~ N(0,eta0^2*(1+k/beta)^(-alpha))
//   Generate gamma_mu = gamma_alpha/gamma_beta and gamma_sigma2 = gamma_mu/gamma_beta
//   Sufficient statistics from gamma_j ~ G(alpha,beta)
//      gamma_sum  = sum(gamma_j)
//      lgamma_sum = sum(log(gamma_j))
//-----------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List GetGammaHyper(double gamma_alpha, double gamma_beta, double gamma_mu, 
                         double gamma_sigma, double gamma_sigma2, arma::colvec gammall, 
                         arma::colvec theta0, double eta02, double gamma_sum, 
                         double lgamma_sum, arma::colvec kall, arma::colvec gamma_met,
                         double mu_bot, double mu_top, double sigma_bot, double sigma_top,
                         double gamma_mu_m0, double gamma_mu_v02, double gamma_sigma_m0, double gamma_sigma_v02,
                         int ingroup, arma::colvec freq_sum, double gmax) {
  
  colvec theta02, gam0vec_old, gam0vec_new, gam0vec, th0v, th0v_old, th0v_new, met_var_new;
  double var_met_g0, std_met_g0, gamma_mu_new, gamma_sigma_new, gamma_sigma2_new, 
  gamma_alpha_new, gamma_beta_new, gamma_prob_new, gamma_prob_old, testpa, 
  pnew, pold;
  
  theta02          = square(theta0);
  met_var_new      = AdaptMetVar(gamma_met); // Variance for Metropolis
  var_met_g0       = 5.66*met_var_new(0);
  std_met_g0       = sqrt(var_met_g0);      // STD DEV for Metropolis
  
  gamma_mu_new     = rcpp_rndtnab(gamma_mu,std_met_g0,mu_bot,mu_top);
  gamma_sigma_new  = rcpp_rndtnab(gamma_sigma,std_met_g0,sigma_bot,sigma_top);
  gamma_sigma2_new = pow(gamma_sigma_new,2);
  gamma_beta_new   = gamma_mu_new/gamma_sigma2_new;
  gamma_alpha_new  = gamma_beta_new*gamma_mu_new;
  gamma_prob_new   = R::pgamma(gmax,gamma_alpha_new, 1/gamma_beta_new, 1, 0);
  gamma_prob_old   = R::pgamma(gmax,gamma_alpha,     1/gamma_beta,     1, 0);
  
  // Likelihood for gamma_j ~ G(alpha,beta)I(gamma_j < gmax)
  testpa = 0;
  testpa = testpa - ingroup*(log(gamma_prob_new) - log(gamma_prob_old));  // Normalizing constant
  testpa = testpa + ingroup*gamma_alpha_new*log(gamma_beta_new);
  testpa = testpa - ingroup*gamma_alpha*log(gamma_beta);
  testpa = testpa - ingroup*lgamma(gamma_alpha_new);
  testpa = testpa + ingroup*lgamma(gamma_alpha);
  testpa = testpa + (gamma_alpha_new - gamma_alpha)*lgamma_sum;
  testpa = testpa - (gamma_beta_new  - gamma_beta)*gamma_sum;
  
  // Likelihood theta_0j ~ N(0,eta0^2*(1+k/gamma_beta)^(-gamma_alpha))
  gam0vec_old    = pow(1+freq_sum/gamma_beta,     -gamma_alpha);
  gam0vec        = gam0vec_old;
  gam0vec_new    = pow(1+freq_sum/gamma_beta_new, -gamma_alpha_new);
  th0v_new       = eta02*gam0vec_new;
  th0v_old       = eta02*gam0vec_old;
  th0v           = th0v_old;
  theta02        = square(theta0);
  testpa         = testpa - accu(theta02/(th0v_new))/2 + accu(theta02/(th0v_old))/2;
  testpa         = testpa - accu(log(th0v_new))/2 + accu(log(th0v_old))/2;
  
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
  pnew = normcdf(sigma_top,gamma_sigma,std_met_g0)     - normcdf(sigma_bot,gamma_sigma,std_met_g0);        // Normalizing constant for gamma_sigma_new
  pold = normcdf(sigma_top,gamma_sigma_new,std_met_g0) - normcdf(sigma_bot,gamma_sigma_new,std_met_g0);    // Normalizing constant for gamma_sigma_old
  testpa = testpa + log(pnew) - log(pold);
  
  if(log(randu())<testpa){
    gamma_alpha  = gamma_alpha_new;
    gamma_beta   = gamma_beta_new;
    gamma_mu     = gamma_mu_new;
    gamma_sigma  = gamma_sigma_new;
    gamma_sigma2 = gamma_sigma2_new;
    gam0vec      = gam0vec_new;
    th0v         = th0v_new;
    //gamma0       = gamma_mu;
    gamma_met(1)++; // pmet
    gamma_met.tail(gamma_met(0)) = met_var_new;
  }  // End generate gamma_alpha and gamma_beta
  gamma_met = AdaptMetUpdate2(gamma_met);
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("gamma_alpha")  = gamma_alpha,
                                      Rcpp::Named("gamma_beta")   = gamma_beta,
                                      Rcpp::Named("gamma_mu")     = gamma_mu,
                                      Rcpp::Named("gamma_sigma")  = gamma_sigma,
                                      Rcpp::Named("gamma_sigma2") = gamma_sigma2,
                                      Rcpp::Named("gam0vec")      = gam0vec,
                                      Rcpp::Named("th0v")         = th0v,
                                      Rcpp::Named("gamma_met")    = gamma_met);
  
  
  return ret;
}



//----------------------------------------------------
//----------------------------------------------------
// GetMainfxdata
//   Compute f at observations for main effects
//   INPUT
//      xgrid   = grid for x
//      xinx    = gives index for xi in xgrid
//                xgrid[xinx[i]] < xdata[i] <= xgrid[xinx[i]+1]
//      xover   = (xdata - xgrid[xinx])/xdelta is excess over boundary
//      idn0    = index of xdata that is not on grid
//      fxgrid  = f computed on xgrid
//                Note: f is increasing on xgrid
//   OUTPUT
//      fxdata   = f computed at observation
// fxdata = GetMainfxdata(xgrid,xdelta,xinx,xover,fxgrid)
//------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec rcpp_GetMainfxdata(arma::uvec xinx , arma::colvec xover, arma::colvec fxgrid) {
  unsigned int n = xinx.n_elem;
  colvec fxdata  = zeros<colvec> (n);
  uvec idn0      = find(xover > 0);
  
  if (idn0.n_elem < n) { // Got data on grid points
    uvec idx = get_ex_idx(n, idn0);
    fxdata.elem(idx) = fxgrid.elem(xinx.elem(idx));
  }
  
  // overshoot
  colvec y0 = fxgrid.elem(xinx.elem(idn0));
  colvec y1 = fxgrid.elem((xinx.elem(idn0)+1));
  fxdata.elem(idn0) = y0 + (y1-y0)%xover.elem(idn0);  // Linearly interpolate  
  
  return(fxdata);
}  // End GetMainfxdata

// [[Rcpp::export]]
arma::mat get_phixgrid2D(arma::mat xgrid, double phix0, double phiy0, arma::mat freq_all, arma::vec xmin,  
                         arma::vec xrange, unsigned int ntheta, unsigned int nint) {
  
  unsigned int ngrid2D = (nint+1)*(nint+1);
  mat phixgrid = zeros<mat> (ngrid2D,ntheta);
  unsigned int icount = 0;
  int k1, k2;
  double xi, yj, phiyj, phixi;
  
  for (unsigned int j = 0; j < nint+1; j++) {
    yj = xgrid(j,1);
    for (unsigned int i = 0; i < nint+1; i++) {
      xi     = xgrid(i,0);
      // Loop over indices where j+k <= kmax
      for (unsigned int kk = 0; kk < ntheta; kk++) {
        k1      = freq_all(kk,0);   // Frequency for xi
        k2      = freq_all(kk,1);   // Frequency for yj
        if (k1==0) {
          phixi = phix0;
        } else {
          phixi = sqrt(2/xrange(0))*cos(k1*pi*(xi-xmin(0))/xrange(0));
        }
        if (k2==0) {
          phiyj = phiy0;
        } else {
          phiyj = sqrt(2/xrange(1))*cos(k2*pi*(yj-xmin(1))/xrange(1));
        }
        phixgrid(icount,kk) = phixi*phiyj;
      }  // end loop over indices
      icount++;
    }  // End loop over xi
  } // End lover over yj
  
  return phixgrid;
}