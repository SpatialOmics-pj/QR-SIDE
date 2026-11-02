// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <boost/math/tools/minima.hpp>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(BH)]]

// cal_xx --- calculate intermediate results
//  ud_xx --- update equations in the algorithm

//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace Rcpp;
using namespace arma;
using namespace std;
using boost::math::tools::brent_find_minima;

const double xi = 1e-8;



// [[Rcpp::export]] 
double sum_nu_lognu0(arma::cube& nu, arma::vec& nu0, arma::mat& rho) {
  int p_m = nu.n_rows, K = nu.n_cols, n = nu.n_slices;
  double tmp = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < p_m; j++)
    {
      for (int k = 0; k < K; k++)
      {
        if((rho(j, k) != 0) && (nu0(k) != 0)) tmp += nu(j, k, i) * log(nu0(k));
      }
    }
  }
  return tmp;
} 

// [[Rcpp::export]] 
double cal_nu_entropy(arma::cube& nu, arma::mat& rho) {
  int p_m = nu.n_rows, K = nu.n_cols, n = nu.n_slices;
  double tmp = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < p_m; j++)
    {
      for (int k = 0; k < K; k++)
      {
        if(rho(j, k) != 0) tmp += -nu(j, k, i) * log(nu(j, k, i) + (nu(j, k, i)==0)) - (1-nu(j, k, i)) * log((1-nu(j, k, i)) + (nu(j, k, i)==1));
      }
    }
  }
  return tmp;
}


 

//call solve.QP from quadprog package
List solve_QP(mat Dmat, vec dvec, mat Amat, vec bvec, int meq = 0){

  // get namespace of package
  Rcpp::Environment pkg = Environment::namespace_env("quadprog");
  // get function from package
  Rcpp::Function f = pkg["solve.QP"];
  
  List out =  f(_["Dmat"]=Dmat, _["dvec"]=dvec, _["Amat"]=Amat, _["bvec"]=bvec, _["meq"]=meq);
  return out; 
}

struct funcdouble
{
  double Xmij, Ni, Lambdaj, pibmu;
  funcdouble(double Xmij, double Ni, double Lambdaj, double pibmu) : Xmij(Xmij), Ni(Ni), Lambdaj(Lambdaj), pibmu(pibmu) {}

  double operator()(double const& x)
  {
    return -Xmij*(x + log(Ni)) + exp(log(Ni) + x) + 0.5*(x*x - 2*x*pibmu)/Lambdaj;
  }
};
  


// [[Rcpp::export]] 
double ud_loglambda_mu(double Xmij, double Ni, double Lambdaj, double pibmu) {
  
  const int double_bits = std::numeric_limits<double>::digits;
  std::pair<double, double> out_x = brent_find_minima(funcdouble(Xmij, Ni, Lambdaj, pibmu), -5.0, 15.0, double_bits);
  // Show all double precision decimal digits and trailing zeros.
  if (false){
  }
  double loglambda_mu = out_x.first;
  return loglambda_mu;
}

// [[Rcpp::export]] 
arma::vec ud_nu0(arma::cube& nu, arma::mat& rho) { 
  int p_m = nu.n_rows, K = nu.n_cols, n = nu.n_slices;
  rowvec rho_nz = sum(rho, 0);
  vec nu0 = ones<vec>(K);
  for (int k = 0; k < K; k++)
  {
    double tmp = 0;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < p_m; j++)
      {
          tmp += rho(j, k) * nu(j, k, i);
      }
    }
    nu0(k) = tmp/(n * rho_nz(k));
  }
  return nu0;
}


// [[Rcpp::export]]
arma::cube cal_mu_m(arma::mat& m, arma::mat& rho, arma::cube& nu, arma::mat& delta){
  cube mu_m(size(nu));
  int n = nu.n_slices;
  for (int i = 0; i < n; i++){      
      mu_m.slice(i) = m + rho % nu.slice(i) % delta;
  }
  return mu_m;
}

// [[Rcpp::export]]
arma::cube cal_bmu(arma::mat& beta, arma::cube& mu_m){
  int n = mu_m.n_slices, p_m = mu_m.n_rows, T = beta.n_rows; 
  cube bmu = cube(T, p_m, n, fill::zeros);
  for (int i = 0; i < n; i++){
      for (int j = 0; j < p_m; j++){
        rowvec mu_m_i_j = mu_m.slice(i).row(j);
        bmu.slice(i).col(j) = sum(beta % repmat(mu_m_i_j, T, 1), 1);
      }
    }
    return bmu;
}

// [[Rcpp::export]]
arma::cube cal_b2mu(arma::mat& beta2, arma::mat& rho, arma::cube& nu, arma::mat& delta2)  {
  int n = nu.n_slices, p_m = nu.n_rows, T = beta2.n_rows; 
  cube b2mu = cube(T, p_m, n, fill::zeros);
  cube mu_diff(size(nu));
  for (int i = 0; i < n; i++){
      mu_diff.slice(i) = rho % nu.slice(i) % delta2 - rho % nu.slice(i) % nu.slice(i) % delta2;
      for (int j = 0; j < p_m; j++){
        rowvec mu_diff_i_j = mu_diff.slice(i).row(j);
        b2mu.slice(i).col(j) = sum(beta2 % repmat(mu_diff_i_j, T, 1), 1);
      }
    }
    return b2mu;
}


// [[Rcpp::export]]
arma::mat  ud_pi(arma::vec& invLambda, arma::mat& loglambda_mu, arma::cube& bmu, arma::cube& b2mu, arma::vec& pi0) {
  int n = loglambda_mu.n_rows, p_m = loglambda_mu.n_cols, T = pi0.n_elem;
  mat pi = ones<mat>(n, T);
  mat phi = zeros<mat>(n, T);
  mat log_relative_pi = zeros<mat>(T, n);
  for (int i = 0; i < n; i++)
  {
    for (int t = 0; t < T; t++)
    {
      double tmp1 = 0, tmp2 = 0;
      for (int j = 0; j < p_m; j++)
      {
        tmp1 += invLambda(j) * pow(loglambda_mu(i, j) - bmu(t, j, i), 2);
        tmp2 += invLambda(j) * b2mu(t, j, i);
      }
      phi(i, t) = -0.5*tmp1 - 0.5*tmp2;   
    }
    //mat log_relative_pi = zeros<mat>(T, n);
    for (int t = 1; t < T; t++)
    {
      log_relative_pi(t, i) = log(pi0(t)/pi0(0)) + phi(i,t) - phi(i,0);   
    }
      

    for (int t = 0; t<T; t++){
      pi(i, t) =1.0 /(accu(exp(log_relative_pi.col(i) - log_relative_pi(t, i))));
    }

    for (int t = 0; t < T; t++){
      if (pi(i,t) < xi) pi(i,t)=xi;
    }
    pi.row(i) = pi.row(i)/sum(pi.row(i));  

  }
  return pi;
}


// [[Rcpp::export]]
arma::vec ud_Lambda(arma::mat& loglambda_mu, arma::mat& loglambda_s2, arma::mat& pi, arma::cube& bmu, arma::cube& b2mu)  {
  int n = pi.n_rows, T = pi.n_cols, p_m = loglambda_mu.n_cols;
  vec Lambda(p_m);
  for (int j = 0; j < p_m; j++)
  {
    double tmp2 = 0, tmp3 = 0;
    for (int i = 0; i < n; i++)
    {
      for (int t = 0; t < T; t++)
      {
        tmp2 += pi(i, t) * pow(loglambda_mu(i, j)- bmu(t, j, i), 2); 
        tmp3 += pi(i, t) * b2mu(t, j, i);
      }
    }
    Lambda(j) = sum(loglambda_s2.col(j))/n + tmp2/n + tmp3/n;
  }
  return Lambda;
}

// [[Rcpp::export]]
void ud_nu(arma::mat& pi, arma::vec& invLambda, arma::mat& loglambda_mu, arma::cube& nu, 
            arma::mat& beta, arma::mat& beta2, arma::mat& m, arma::mat& rho, arma::mat& delta, arma::mat& delta2, arma::vec& nu0) {
  int n = pi.n_rows, T = pi.n_cols, p_m = loglambda_mu.n_cols, K = m.n_cols;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < p_m; j++)
    {
      for (int k = 0; k < K; k++)
      {
        if (rho(j, k) != 0)
        {
          double tmp1 = 0, tmp2 = 0;
          for (int t = 0; t < T; t++)
          {
            double tmp3 = 0;
            for (int k_ = 0; k_ < K; k_++)
            {
              if(k_ != k) tmp3 += beta(t, k_) * (m(j, k_) + rho(j, k_)*nu(j, k_, i)*delta(j, k_)); 
            }
            tmp1 += pi(i, t) * invLambda(j) * (loglambda_mu(i, j) - tmp3 - beta(t, k) * m(j, k)) * beta(t, k) * rho(j, k) * delta(j, k);
            tmp2 += pi(i, t) * invLambda(j) * beta2(t, k) * rho(j, k) * delta2(j, k);
          }
          double phi = tmp1 - 0.5*tmp2 + log(nu0(k) / (1 - nu0(k)));
          nu(j, k, i) = 1.0/(1.0 + exp(-phi));
          //if (nu(j, k, i) > (1.0 - xi)) nu(j, k, i) = 1.0 - xi;
        }
      }
    }
  }
}

// [[Rcpp::export]]
void ud_m(arma::mat& pi, arma::mat& loglambda_mu, arma::cube& nu, 
            arma::mat& beta, arma::mat& beta2, arma::mat& m, arma::mat& rho, arma::mat& delta) {
  int n = pi.n_rows, T = pi.n_cols, p_m = loglambda_mu.n_cols, K = m.n_cols;
  for (int j = 0; j < p_m; j++)
  {
    for (int k = 0; k < K; k++)
    {
      double tmp1 = 0, tmp2 = 0;
      for (int i = 0; i < n; i++)
      {
        for (int t = 0; t < T; t++)
        {
          double tmp3 = 0;
          for (int k_ = 0; k_ < K; k_++)
          {
            if(k_ != k) tmp3 += beta(t, k_) * (m(j, k_) + rho(j, k_)*nu(j, k_, i)*delta(j, k_)); 
          }
          tmp1 += pi(i, t) * (loglambda_mu(i, j) - tmp3 - beta(t, k) * rho(j, k) * nu(j, k, i) * delta(j, k)) * beta(t, k);
          tmp2 += pi(i, t) * beta2(t, k);
        }
      }
      m(j, k) = tmp1/tmp2;  
    }
  }  
}


// [[Rcpp::export]]
void ud_delta(arma::mat& pi, arma::mat& loglambda_mu, arma::cube& nu, 
            arma::mat& beta, arma::mat& beta2, arma::mat& m, arma::mat& rho, arma::mat& delta) {
  int n = pi.n_rows, T = pi.n_cols, p_m = loglambda_mu.n_cols, K = m.n_cols;

  for (int j = 0; j < p_m; j++)
  {
    for (int k = 0; k < K; k++)
    {
      if (rho(j, k) != 0)  {
        double tmp1 = 0, tmp2 = 0;
        for (int i = 0; i < n; i++)
        {
          for (int t = 0; t < T; t++)
          {
            double tmp3 = 0;
            for (int k_ = 0; k_ < K; k_++)
            {
              if(k_ != k) tmp3 += beta(t, k_) * (m(j, k_) + rho(j, k_) * nu(j, k_, i)*delta(j, k_)); 
            }

            tmp1 += pi(i, t) * (loglambda_mu(i, j) - tmp3 - beta(t, k) * m(j, k)) * nu(j, k, i) * beta(t, k);
            tmp2 += pi(i, t) * nu(j, k, i) * beta2(t, k);
          }
        }
        delta(j, k) = (tmp1/tmp2); 
        if (delta(j, k) < 0) delta(j, k) = 0;
       } 
    }
  }  
}


// [[Rcpp::export]]
double cal_Dmat(arma::mat& pi, arma::mat& loglambda_mu, arma::vec& invLambda, arma::mat& m, arma::mat& rho, 
              arma::cube& nu, arma::mat& delta, arma::cube& mu_m, int t, int k, int k_)  {
  int n = pi.n_rows, p_m = m.n_rows;
  double Dmat_kk_ = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < p_m; j++)
    {
      if(k == k_)  Dmat_kk_ += pi(i, t) * invLambda(j) * (m(j, k)*m(j, k) + 2*m(j, k)*rho(j, k)*nu(j, k, i)*delta(j, k) + rho(j, k)*nu(j, k, i)*delta(j, k)*delta(j, k));
      if(k != k_)  Dmat_kk_ += pi(i, t) * invLambda(j) * mu_m(j, k, i) *  mu_m(j, k_, i); 
    }
  }
  return Dmat_kk_;
}

// [[Rcpp::export]]
arma::mat ud_beta(arma::mat& pi, arma::mat& loglambda_mu, arma::vec& invLambda, arma::mat& m, arma::mat& rho, 
            arma::cube& nu, arma::mat& delta, arma::cube& mu_m) {
  int n = pi.n_rows, T = pi.n_cols, p_m = m.n_rows, K = m.n_cols;
  mat beta = zeros<mat>(T, K);
  for (int t = 0; t < T; t++){
    // calculate Dmat and dvec
    mat Dmat = zeros<mat>(K, K);
    rowvec dvec = zeros<rowvec>(K);
    for (int k = 0; k < K; k++)
    {
      for (int k_ = k; k_ < K; k_++)
      {
        Dmat(k, k_) = cal_Dmat(pi, loglambda_mu, invLambda, m, rho, nu, delta, mu_m, t, k, k_);
        if (k_ != k) Dmat(k_, k) = Dmat(k, k_);  
      }
    }
    if(false) {
    }
    for (int k = 0; k < K; k++)
    {
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < p_m; j++)
        {
          dvec(k) += pi(i, t) * invLambda(j) * loglambda_mu(i, j) * mu_m(j, k, i);
        }
      }
    }

    // call solver
    mat Amat = zeros<mat>(K+1, K);
    Amat.row(0) = ones<rowvec>(K);
    Amat.diag(-1) = ones<rowvec>(K);
    vec bvec = zeros<vec>(K+1);
    bvec(0) = 1;
    double sc = norm(Dmat);
    List out = solve_QP(Dmat/sc, trans(dvec)/sc, trans(Amat), bvec, 1);
    vec solution = out["solution"];
    beta.row(t) = trans(solution);  
    for (int k = 0; k < K; k++){
      if (beta(t, k) < xi) beta(t, k) = xi;
    }
    beta.row(t)  = beta.row(t) /sum(beta.row(t));  
  }  
  return beta;
}



// ELBO
// [[Rcpp::export]]
double calELBO(arma::vec& N, arma::mat& X_m, arma::mat& rho, int T, arma::mat& loglambda_mu, arma::mat& loglambda_s2, 
                arma::vec& Lambda, arma::vec& invLambda, arma::mat& pi, arma::vec& pi0, 
                arma::mat& beta, arma::mat& beta2, arma::cube& nu, arma::cube& nu2, arma::vec& nu0, 
                arma::cube& mu_m, arma::mat& delta, arma::mat& delta2){

    
    int n = X_m.n_rows;
    int p_m = X_m.n_cols;
                  
    mat loglambda_mu_logN = loglambda_mu + repmat(log(N), 1, p_m);   
    double possion_term = accu(X_m % loglambda_mu_logN);    
    double possion_term2 = -accu(repmat(N, 1, p_m) % exp(loglambda_mu + 0.5 * loglambda_s2));
          
    vec term = zeros<vec>(8);
    term(0) = -0.5 * n * sum(log(Lambda));
    term(2) = -0.5 * accu(loglambda_s2 % repmat(trans(invLambda), n, 1)); 
    term(4) = accu(pi % repmat(log(trans(pi0)), n, 1));
    term(6) = 0.5 * accu(log(loglambda_s2)) - accu(pi % log(pi));  
    cube tmp1 = 1.0 - nu;
    vec tmp2 = 1.0 - nu0;
    term(5) = sum_nu_lognu0(nu, nu0, rho) + sum_nu_lognu0(tmp1, tmp2, rho); 
    term(7) = cal_nu_entropy(nu, rho); 

    for (int i = 0; i<n; i++){
      mat nu_i = nu.slice(i);
      mat mu_m_i = mu_m.slice(i);
      mat nu2_i = nu2.slice(i);  
      for (int j = 0; j < p_m; j++){
        for (int t = 0; t < T; t++){
          term(1) += -0.5 * pi(i,t) * invLambda(j) * pow(loglambda_mu(i, j) - sum(beta.row(t) % mu_m_i.row(j)), 2);

          term(3) += -0.5 * pi(i,t) * invLambda(j) * sum(beta2.row(t) % (rho.row(j) % nu_i.row(j) % delta2.row(j) - 
                                                                    rho.row(j) % nu2_i.row(j) % delta2.row(j)));                                          
        }                      
      }
    }
    
     
    double ELBO = accu(term) + possion_term + possion_term2;
    return(ELBO);
}



// [[Rcpp::export]]
List DeconvMarker(arma::mat& X_m, arma::mat& rho, int T, arma::mat beta, 
                          arma::mat delta, arma::mat pi, arma::vec pi0, 
                          arma::vec Lambda, arma::mat m, int max_iter = 100, double eps = 1e-5, bool verbose = true){
  
    
  int n = X_m.n_rows;
  int p_m = X_m.n_cols;
  int K = rho.n_cols;
  int iter = 1;
  
  rowvec ELBO = zeros<rowvec>(1, 100000);
  ELBO(0) = -datum::inf;
  
  mat TERM = zeros<mat>(100000, 9);
  TERM.row(0) = TERM.row(0)-datum::inf;

  mat beta2 = beta % beta;
  //mat m = zeros<mat>(p_m, K);
  //mat delta = zeros<mat>(p_m, K);
  mat delta2 = delta % delta;

  //cube nu = cube(p_m, K, n, fill::randu);
  cube nu = cube(p_m, K, n, fill::ones);

  nu = 0.5 * nu;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < p_m; j++)
    {
      for (int k = 0; k < K; k++)
      {
        if (rho(j, k) == 0) nu(j, k, i) = 0;
      }
    }
  }

  cube nu2 = nu % nu;
  vec nu0 = ud_nu0(nu, rho); 
  //cube mu_m = cube(p_m, K, n, fill::zeros);
  cube mu_m = cal_mu_m(m, rho, nu, delta);
  //cube bmu = cube(T, p_m, n, fill::zeros);
  cube bmu = cal_bmu(beta, mu_m);
  //cube b2mu = cube(T, p_m, n, fill::zeros);
  cube b2mu = cal_b2mu(beta2, rho, nu, delta2);

  mat pibmu = ones<mat>(n, p_m);

  vec invLambda = 1.0/Lambda;
  mat loglambda_mu = ones<mat>(n, p_m);
  mat loglambda_s2 = ones<mat>(n, p_m);
  vec N = sum(X_m, 1); 
  N = ones<vec>(n);
        
    
  for (iter = 1; iter < max_iter; iter++){

    ////////////////////////////////// E-step
    for (int i = 0; i < n; i++){
      for (int j = 0; j < p_m; j++){
        // update loglambda_mu 
        rowvec pi_i = pi.row(i);
        pibmu(i, j) = sum(trans(pi_i) % bmu.slice(i).col(j));
        loglambda_mu(i,j) = ud_loglambda_mu(X_m(i,j),  N(i), Lambda(j), pibmu(i, j));

        //update loglambda_s2
        loglambda_s2(i,j) = 1.0/(N(i)*exp(loglambda_mu(i,j)) + invLambda(j));
      }
    }   

    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);  
      //Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
      ELBO(iter) = part;
      TERM(iter, 0) = part  - ELBO(iter-1);
    }
    // pi
    pi = ud_pi(invLambda, loglambda_mu, bmu, b2mu, pi0);   
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 1) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    // nu 
    ud_nu(pi, invLambda, loglambda_mu, nu, beta, beta2, m, rho, delta, delta2, nu0);
    mu_m = cal_mu_m(m, rho, nu, delta);
    nu2 = nu % nu;
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 2) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    ////////////////////////////////// Variational M-step
    // m
    ud_m(pi, loglambda_mu, nu, beta, beta2, m, rho, delta);
    mu_m = cal_mu_m(m, rho, nu, delta);
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 3) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    // delta
    ud_delta(pi, loglambda_mu, nu, beta, beta2, m, rho, delta);
    delta2 = delta % delta;
    mu_m = cal_mu_m(m, rho, nu, delta);
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 4) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    // beta 
    beta = ud_beta(pi, loglambda_mu, invLambda, m, rho, nu, delta, mu_m);
    beta2 = beta % beta; 
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 5) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    // pi0
    for (int t = 0; t<T; t++){
      pi0(t) = sum(pi.col(t))/n;
    }
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 6) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    // nu0
    nu0 = ud_nu0(nu, rho); 
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 7) = part  - ELBO(iter);
      ELBO(iter) = part;
    }

    // Lambda  
    bmu = cal_bmu(beta, mu_m);
    b2mu = cal_b2mu(beta2, rho, nu, delta2);
    Lambda = ud_Lambda(loglambda_mu, loglambda_s2, pi, bmu, b2mu);
    invLambda = 1.0/Lambda;
    if(verbose) {
      double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      TERM(iter, 8) = part  - ELBO(iter);
      ELBO(iter) = part;
    }
    //ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
    
    ELBO(iter) = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
    if( abs(ELBO(iter)  - ELBO(iter-1))/ abs(ELBO(iter-1)) < eps) break; 
    }
  
  
  List resList = List::create(
    Rcpp::Named("loglambda_mu") = loglambda_mu,
    Rcpp::Named("loglambda_s2") = loglambda_s2,
    Rcpp::Named("pi") = pi,
    Rcpp::Named("pi0") = pi0,
    Rcpp::Named("nu") = nu,
    Rcpp::Named("nu0") = nu0,
    Rcpp::Named("bmu") = bmu,
    Rcpp::Named("b2mu") = b2mu,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("delta") = delta,
    Rcpp::Named("m") = m,
    Rcpp::Named("Lambda") = Lambda, 
    Rcpp::Named("TERM") = TERM.rows(1, iter-1), 
    Rcpp::Named("dLogLik") = ELBO(iter) - ELBO(iter-1),
    Rcpp::Named("loglik") = ELBO.subvec(1,iter-1));

  return(resList);
  
}
 
