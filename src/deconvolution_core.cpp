// Core implementation from this repository's bundled SpatialDecon_1.0.tar.gz.
// Keep it in QRSIDE so installation cannot select an unrelated SpatialDecon.
// Removed diagnostic submatrix printing, which failed for dimensions below 5.
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <boost/math/tools/minima.hpp>
#include <math.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

namespace qrside_core {
using namespace Rcpp;
using namespace arma;
using namespace std;
using boost::math::tools::brent_find_minima;

const double log2pi = std::log(2.0 * M_PI);
double X_count_ij, N_i, alpha_i, Wz_mu_ji, Lambda_j;


//call solve.QP from quadprog package
List solve_QP(mat Dmat, vec dvec, mat Amat, vec bvec, int meq = 0){

  // get namespace of package
  Rcpp::Environment pkg = Environment::namespace_env("quadprog");
  // get function from package
  Rcpp::Function f = pkg["solve.QP"];
  
  List out =  f(_["Dmat"]=Dmat, _["dvec"]=dvec, _["Amat"]=Amat, _["bvec"]=bvec, _["meq"]=meq);
  return out; 
}

double calELBO(int T, int K, int q, arma::vec& N, arma::mat& X_count, arma::mat& loglambda_mu, arma::mat& loglambda_s2, arma::vec& invLambda, arma::vec& Lambda, arma::mat& Wz_mu, arma::vec& alpha, arma::mat& z_mu, arma::mat& z_s2, arma::mat& transW_invLambda_W, arma::field<cube>& b, arma::field<cube>& m, arma::mat mu, arma::mat& pi, arma::mat& pi0, arma::mat& epsilon2, arma::mat& sigma2, arma::mat& beta, arma::mat& beta2, bool verbose){
    // ELBO
    
    int n = X_count.n_rows;
    int p = X_count.n_cols;
              
    mat loglambda_mu_logN = loglambda_mu.each_col() + log(N);   
    double possion_term = accu(X_count%loglambda_mu_logN);    
    double possion_term2 = -accu(repmat(N, 1, p)%exp(loglambda_mu + 0.5*loglambda_s2));
      
      
    if(verbose){
    }
    
    double dimreduction_term = 0;
    vec d_invLambda = loglambda_s2 * invLambda;
    vec loglambda_mu_alpha_We = zeros<vec>(p);
    for (int i = 0; i < n; i++){
        loglambda_mu_alpha_We = trans(loglambda_mu.row(i)) - alpha(i) - Wz_mu.col(i);
        dimreduction_term = dimreduction_term  - 0.5*accu(loglambda_mu_alpha_We%invLambda%loglambda_mu_alpha_We) - 0.5*d_invLambda(i) - 0.5*accu(diagvec(transW_invLambda_W)%trans(z_s2.row(i))); 
    }
    dimreduction_term = dimreduction_term - 0.5*n*accu(log(Lambda));
      
      
      
    if(verbose){
    }
    
    mat term = zeros<mat>(4, 1);
    term(0,0) = accu(pi % repmat(log(trans(pi0)), n, 1));
    term(3,0) = - accu(pi%log(pi)) + 0.5*accu(log(loglambda_s2)) + 0.5*accu(log(z_s2));      
    
    for (int i = 0; i<n; i++){
      for (int t = 0; t<T; t++){
        mat b_it = b(i).slice(t);
        mat m_it = m(i).slice(t);
        for (int j = 0; j<q; j++){
          term(1, 0) = term(1, 0) + pi(i,t) * ( -0.5 * log(epsilon2(i,j)) 
                                                  - 0.5* ( pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m_it.row(j)) + z_s2(i,j) )/epsilon2(i,j));
          for (int k = 0; k< K; k++){
            term(2, 0) = term(2, 0) + pi(i,t) * ( -0.5 * log(sigma2(k,j)) 
                                                 - 0.5 * ( pow((b_it(j,k) - mu(k,j)),2) + m_it(j,k) )/sigma2(k,j)  );
          }                      
          for (int k = 0; k< K; k++){
            term(3) = term(3) + pi(i,t) * ( 0.5*log(m_it(j,k)) );
          } 
          
        }
      }
    }
    
    
    if(verbose){
    }
    
    
    double ELBO = accu(term) + possion_term + possion_term2 + dimreduction_term;
    return(ELBO);
}


List CountDeconvolution(arma::mat& X_count, arma::mat Z, arma::mat W, int K, int T, arma::mat beta, arma::mat mu, arma::mat pi, arma::mat pi0, arma::vec Lambda, arma::vec alpha, arma::mat epsilon2, int max_iter = 100, double eps = 1e-5,  bool verbose = true, bool fixed_loglambda_mu = false, bool fixed_loglambda_s2 = false, bool fixed_z_mu = false, bool fixed_z_s2 = false,  bool fixed_m = false, bool fixed_b = false, bool fixed_pi = false, bool fixed_W = false, bool fixed_alpha = false, bool fixed_Lambda = false,  bool fixed_pi0 = false, bool fixed_epsilon2 = false, bool fixed_mu = false, bool fixed_sigma2 = false, bool fixed_beta2 = false, bool alpha_spot_specific = false,  double cutoff = 1){
  
    
  int n = X_count.n_rows;
  int p = X_count.n_cols;
  int iter = 1;
  int q = Z.n_cols;
  
  if(verbose){
  }
    
  rowvec ELBO = zeros<rowvec>(1, 100000);
  ELBO(0) = -datum::inf;
  //mat pi = ones<mat>(n, T)/T;
  //mat pi0 = ones<mat>(T, 1)/T;
  
  if(verbose){
  }
  //mat beta = zeros<mat>(T, K);
  //beta.col(0) = ones<vec>(T);
  mat beta2 = beta%beta;
  
  if(verbose){
  }
  field<cube> b(n);
  field<cube> m(n);
  
  mat m_initial = ones<mat>(q, K)*0.1;
  for (int i=0; i<n; i++){
    b(i) = cube(q, K, T);
    m(i) = cube(q, K, T);
    for (int t = 0; t < T; t++){
      b(i).slice(t) = trans(mu);
      m(i).slice(t) = m_initial;
    }

  }
  cube m2 = cube(q, K, n);
  cube b2 = cube(q, K, n);
  for (int i = 0; i<n; i++){
    b2.slice(i) = trans(mu);
    m2.slice(i) = m_initial;
  }
  
  
  
  if(verbose){
  }
    
  //mat epsilon2 = ones<mat>(n, q);
  mat sigma2 = ones<mat>(K, q)*0.1;
  //mat mu = zeros<mat>(K, P);

  mat term = zeros<mat>(4, 1);
    
  mat diff1 = zeros<mat>(n,T);
  mat diff2 = zeros<mat>(n,T);   
    

  //Lambda = ones<vec>(p); // the default p is 2000, the default q is 15
  vec invLambda = 1/Lambda;
  //mat W = ones<mat>(p, q)*0.0001;
  mat z_mu = Z;
  mat z_s2 = ones<mat>(n, q);
  mat loglambda_mu = ones<mat>(n, p);
  mat loglambda_s2 = ones<mat>(n, p);
  //vec alpha = zeros<vec>(n);
  vec N = sum(X_count, 1); 
  //N = ones<vec>(n);
  mat Wz_mu = W*trans(z_mu); // p*n
  mat log_relative_pi = zeros<mat>(T, n);
    
  mat transW = trans(W);
  mat transW_invLambda = transW.each_row() % trans(invLambda);
  mat transW_invLambda_W = transW_invLambda * W;
    
  //double X_count_ij, N_i, alpha_i, Wz_mu_ji, Lambda_j; 
    
  if(verbose){
  }
    
    
    
  for (iter = 1; iter < max_iter; iter++){

    // update loglambda_mu and loglambda_s2
    for (int i = 0; i < n; i++){
        for (int j = 0; j < p; j++){
            
            X_count_ij = X_count(i,j);
            N_i = N(i);
            alpha_i = alpha(i);
            Wz_mu_ji = Wz_mu(j,i);
            Lambda_j = Lambda(j);
  
            struct funcdouble
            {
                double operator()(double const& x)
                {
                    return -X_count_ij*(x + log(N_i)) + exp(log(N_i) + x) + 0.5*(x - alpha_i - Wz_mu_ji)*(x - alpha_i - Wz_mu_ji)/Lambda_j;
                }
            };
            
            if (fixed_loglambda_mu != true){
              const int double_bits = std::numeric_limits<double>::digits;
              std::pair<double, double> out_x = brent_find_minima(funcdouble(), -50.0, 50.0, double_bits);
              // Show all double precision decimal digits and trailing zeros.
              if (i==0 & j < 10){
              }
              loglambda_mu(i,j) = out_x.first;
            }                 
            
            if (fixed_loglambda_s2 != true){
              loglambda_s2(i,j) = 1.0/(N(i)*exp(loglambda_mu(i,j)) + 1.0/Lambda(j));
            }
        }
    }

    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
    }
                                                
    // update z_mu and z_s2
    for (int i = 0; i < n; i++){
      if (fixed_z_s2 != true){
        z_s2.row(i) = 1.0/(trans(diagvec(transW_invLambda_W)) + 1.0/epsilon2.row(i));
      }
      
      if (fixed_z_mu != true){
        mat b2_i = b2.slice(i); //q-by-K
        vec pi_epsilon2_b2_beta = zeros<vec>(q);
        for (int t = 0; t < T; t++){
          pi_epsilon2_b2_beta = pi_epsilon2_b2_beta + pi(i,t) * diagmat(1/epsilon2.row(i)) * b2_i *  trans(beta.row(t));
        }
        mat X_z_mu =  transW_invLambda_W + diagmat(1/epsilon2.row(i));
        vec Y_z_mu = pi_epsilon2_b2_beta + transW_invLambda * (trans(loglambda_mu.row(i)) - alpha(i));
        z_mu.row(i) = trans(solve(X_z_mu, Y_z_mu));
      }
     }

      
    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    }
      
    diff1 = zeros<mat>(n,T);
    diff2 = zeros<mat>(n,T);
    
    if(verbose){
    }
    
    for (int i = 0; i<n; i++){
      
      // update b2
      mat b2_i = b2.slice(i);
      if (fixed_b != true){
        for (int j = 0; j<q; j++){
          mat X = zeros(K,K);
          vec Y = zeros(K,1);
          for (int t = 0; t<T; t++){
            X = X + pi(i,t)*(trans(beta.row(t))*beta.row(t)/epsilon2(i,j) + diagmat(1.0/sigma2.col(j)));
            Y = Y + pi(i,t)*(mu.col(j)/sigma2.col(j) + z_mu(i,j)/epsilon2(i,j)*trans(beta.row(t)));
          }
          b2_i.row(j) = trans(solve(X, Y));
        }
        b2.slice(i) = b2_i;
      }

      // update m2
      mat m2_i = m2.slice(i);
      if (fixed_m != true){
        for (int j = 0; j<q; j++){
          for (int k = 0; k< K; k++){
            m2_i(j,k) =  1.0/ sum(pi.row(i)%(trans(beta2.col(k)) / epsilon2(i,j) + 1.0/sigma2(k,j)));
          }
        }
        m2.slice(i) = m2_i;
      }
      
      for (int t = 0; t<T; t++){
        b(i).slice(t) = b2.slice(i);
      }
      
      for (int t = 0; t<T; t++){
        m(i).slice(t) = m2.slice(i);
      }    
      
      if (fixed_pi != true){
        for (int t = 0; t<T; t++){
          mat b_it = b(i).slice(t);
          for (int j = 0; j< q; j++){
            diff1(i,t) = diff1(i,t) + (pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m2_i.row(j)) + z_s2(i,j) )/(2*epsilon2(i,j));
          }
          
          if(verbose){
          }
          // update pi

          
          if ( fixed_pi != true){
            for (int t = 1; t<T; t++){
              log_relative_pi(t,i) = log(pi0(t)/pi0(0)) - diff1(i,t) + diff1(i,0);
            }
            for (int t = 0; t<T; t++){
              pi(i, t) =1.0 /(accu(exp(log_relative_pi.col(i) - log_relative_pi(t,i))));
            }
            
            for (int t = 0; t < T; t++){
              if (pi(i,t) < 0.00001){
                pi(i,t)=0.00001;
              }
            }
            pi.row(i) = pi.row(i)/sum(pi.row(i));
            
            if (iter == 5){
            }
          }
        }
      }
      
    }
      
    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    }
    
    if(verbose){
    }
    
    
    // Variational M-step
    // update W
    if (fixed_W != true){
      mat X_W = zeros(q, q);
      mat Y_W = zeros(p, q);
      for (int i = 0; i < n; i++){
        X_W = X_W + trans(z_mu.row(i)) * z_mu.row(i) + diagmat(z_s2.row(i));
        Y_W = Y_W + (trans(loglambda_mu.row(i)) - alpha(i)) * z_mu.row(i);
      }
      W = trans(solve(X_W, trans(Y_W)));
      Wz_mu = W*trans(z_mu); // p*n
    }
    transW = trans(W);
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W; 
      
      
    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    }
      
    // update Lambda
    if (fixed_Lambda != true){
      mat loglambda_mu_alpha_Wz = loglambda_mu - repmat(alpha, 1, p) - trans(Wz_mu);
      mat loglambda_mu_alpha_Wz2_loglambda_s2 = loglambda_mu_alpha_Wz%loglambda_mu_alpha_Wz + loglambda_s2;
      vec diagWtdiagsumfW = sum(W.each_row()%sum(z_s2, 0)%W, 1); //p-by-1
      Lambda = trans(sum(loglambda_mu_alpha_Wz2_loglambda_s2, 0))/n + diagWtdiagsumfW/n;
      invLambda = 1.0/Lambda;
    }
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W; 
  
      
    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    }
    
    if(verbose){
    }    
    
    // update alpha
    if (fixed_alpha != true){
      if (alpha_spot_specific == true){
        for (int i = 0; i < n; i++){
          alpha(i) = sum((trans(loglambda_mu.row(i)) - Wz_mu.col(i))%invLambda)/sum(invLambda);
        }
      }else{
          mat tloglambda_mu_Wz_mu = trans(loglambda_mu) - Wz_mu;
          alpha(0) = accu(tloglambda_mu_Wz_mu.each_col()%invLambda)/(n*sum(invLambda));
          for (int i = 1; i < n; i++){
            alpha(i) = alpha(0);
          }
      }
    }
    
    
    
    
    
    if(verbose){
    }
    
    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    }
    
    if ( fixed_pi0 != true){
      for (int t = 0; t<T; t++){
        pi0(t) = sum(pi.col(t))/n;
      }
    }
    
    if(verbose){
    }
    
    if(verbose){
      ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    }
    
    
      // if (fixed_epsilon2 != true){
      //   for (int i = 0; i<n; i++){
      //     for (int j = 0; j<q; j++){
      //       double Exp1 = 0;
      //       for (int t = 0; t<T; t++){
      //         mat b_it = b(i).slice(t);
      //         mat m_it = m(i).slice(t);
      //         Exp1 = Exp1 +  pi(i,t) * (pow((z_mu(i,j) - sum(beta.row(t)%b_it.row(j))),2) + sum( beta2.row(t)%m_it.row(j)) + z_s2(i,j)); 
      //       }
      //       epsilon2(i,j) =  Exp1  / sum(pi.row(i));
      //     }
      //   }  
      // }
      
      
      if (fixed_epsilon2 != true){
        for (int j = 0; j<q; j++){
          double Exp1 = 0;
          for (int i = 0; i<n; i++){
            for (int t = 0; t<T; t++){
              mat b_it = b(i).slice(t);
              mat m_it = m(i).slice(t);
              Exp1 = Exp1 +  pi(i,t) * (pow((z_mu(i,j) - sum(beta.row(t)%b_it.row(j))),2) + sum( beta2.row(t)%m_it.row(j)) + z_s2(i,j)); 
            }
          }
          epsilon2(0,j) =  Exp1  / accu(pi);
        }  
      }
      if (fixed_epsilon2 != true){
        for (int j = 0; j<q; j++){
          for (int i = 1; i < n; i++){
            epsilon2(i,j) =  epsilon2(0,j); 
          }
        }
      }

      
      if(verbose){
      }
      
      if(verbose){
        ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
        if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
      }
      
      mat pi_b_it = zeros<mat>(q, K);
      
      if (fixed_mu != true){
        for (int i = 0; i<n; i++){
          for (int t = 0; t < T; t++){
            mat b_it = b(i).slice(t);
            pi_b_it = pi_b_it + pi(i,t) * b_it;
          }
        }
        mu = trans(pi_b_it) / accu(pi);
      }
      
      if(verbose){
        ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
        if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
      }
      
      if(verbose){
      }
      
      if (fixed_sigma2 != true){
        for (int j = 0; j< q; j++){
          for (int k = 0; k< K; k++){
            double Exp2 = 0;
            for (int i = 0; i<n; i++){
              for (int t = 0; t < T; t++){
                mat b_it = b(i).slice(t);
                mat m_it = m(i).slice(t);
                Exp2 = Exp2 + sum(pi(i,t)*(pow((b_it(j,k) - mu(k,j)),2) + m_it(j,k)));
              }
            }
            sigma2(k,j) = Exp2/accu(pi);
          }
        }
      }
      
      if(verbose){
        ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
        if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
      }
      
      if(verbose){
      }
      
      if (iter > 4){
        if ((ELBO(iter-2)  - ELBO(iter-3)) < cutoff){
          if (fixed_beta2 != true){
            for (int t = 0; t < T; t++){
              mat Dmat = zeros<mat>(K, K);
              rowvec dvec = zeros<rowvec>(K);
              for (int i = 0; i<n; i++){
                mat b_it = b(i).slice(t);
                mat m_it = m(i).slice(t);
                for (int j = 0; j < q; j++){
                  Dmat = Dmat + pi(i,t)/epsilon2(i,j) *(trans(b_it.row(j))*b_it.row(j) + diagmat(trans(m_it.row(j))));
                  dvec = dvec + z_mu(i,j)*pi(i,t)/epsilon2(i,j) * b_it.row(j);
                }
              }
              //min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.
              
              mat Amat = zeros<mat>(K+1, K);
              Amat.row(0) = ones<rowvec>(K);
              Amat.diag(-1) = ones<rowvec>(K);
              vec bvec = zeros<vec>(K+1);
              bvec(0) = 1;
              double sc = norm(Dmat);
              List out = solve_QP(Dmat/sc, trans(dvec)/sc, trans(Amat), bvec, 1);
              vec solution = out["solution"];
              beta.row(t) = trans(solution);
            }
            beta2 = beta%beta;
          }
        }
      }
      
      
    transW = trans(W);
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W;                    
    
    ELBO(iter) = calELBO(T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b, m, mu, pi, pi0, epsilon2, sigma2, beta, beta2, verbose);
      
    if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    
    if( (ELBO(iter)  - ELBO(iter-1))  < -1e-7){
      if (verbose) Rcpp::Rcout << "Likelihood stopped increasing; returning the current fit.\n";
      break;
    }
      
    if( abs(ELBO(iter)  - ELBO(iter-1))/ abs(ELBO(iter-1)) < eps) break;
    }
  
  
  List resList = List::create(
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("Lambda") = Lambda, 
    Rcpp::Named("W") = W,  
    Rcpp::Named("loglambda_mu") = loglambda_mu,
    Rcpp::Named("loglambda_s2") = loglambda_s2,
    Rcpp::Named("z_mu") = z_mu,
    Rcpp::Named("z_s2") = z_s2,  
    Rcpp::Named("pi") = pi,
    Rcpp::Named("pi0") = pi0,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("b") = b,
    Rcpp::Named("m") = m,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("epsilon2") = epsilon2,
    Rcpp::Named("diff1") = diff1,
    Rcpp::Named("dLogLik") = ELBO(iter) - ELBO(iter-1),
    Rcpp::Named("loglik") = ELBO.subvec(1,iter-1));

  return(resList);
  
}






vec dmvnrm(const mat& x,  
           rowvec mean,  
           mat sigma, 
           bool logd = false,
           int cores = 1) { 
  //omp_set_num_threads(cores);
  int n = x.n_rows;
  int xdim = x.n_cols;
  vec out(n);
  mat rooti = inv(chol(sigma, "lower"));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(xdim/2) * log2pi;
#pragma omp parallel for schedule(static) 
  for (int i=0; i < n; i++) {
    vec z = rooti * trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * sum(z%z) + rootisum;     
  }  
  
  if (logd==false) {
    out=exp(out);
  }
  return(out);
}



sp_mat get_spNbs(ivec y, const sp_mat& Adj) {   
  // row is for pixel.
  //output a sparse matrix, i-th row contains labels of neighbor_i. 
  // Make const iterator
  arma::sp_mat::const_iterator start = Adj.begin();
  //arma::sp_mat::const_iterator end = Adj.end();
  
  // Calculate number of nonzero points
  int n = Adj.n_nonzero;
  
  sp_mat spNbs(y.n_elem, y.n_elem);    //neiborhood state matrix, matched with Adj.
  
  arma::sp_mat::const_iterator it = start; 
  for(int i = 0; i < n; ++i)
  {
    spNbs(it.row(), it.col()) = y(it.col());
    ++it; // increment
  }
  
  return spNbs;
}


arma::mat calXenergy2D_sp(arma::ivec x, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  int n = x.n_rows;
  arma::sp_mat spNbs = get_spNbs(x, Adj);
  arma::sp_mat spNbs_t = spNbs.t();  // transform spNbs to iterate by column.
  arma::mat Ux(n, K);
  int i, k;
  for (k = 0; k < K; k++)
  {
    for (i = 0; i < n; i++)
    {
      arma::sp_mat col(spNbs_t.col(i));
      double n_sameS = 0;
      //double n_nb = 0;
      //int nn = std::distance(col.begin(), col.end());
      int nn = col.n_nonzero;
      for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
        n_sameS += (*j) == (k+1);
        //n_nb++;
      }
      Ux(i, k) = alpha(k) + beta * (nn - n_sameS)/2;
    }
  }

  return Ux;
  
}

double obj_beta(const arma::ivec& y, const arma::mat& gam, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  
  mat Ux = calXenergy2D_sp(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
  mat pxgn = normalise(exp(-Ux), 1, 1); // set all rowSums to be ONE.
  return accu(gam % (log(pxgn)));
}


Rcpp::List runICM_sp(const arma::mat &Uz, arma::ivec y, const arma::sp_mat& Adj,
                      arma::vec beta0_y, double beta_y, arma::mat epsilon2, int maxIter_ICM)	{
  int n = y.n_rows, T = Uz.n_cols;
  int iter;
  
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  mat Uy(n, T);
  mat U(n, T);
  //--------------------------------------------------------------------------------	
  // ICM algrithm
  //--------------------------------------------------------------------------------
  int Iteration = 1;
  

  for (iter = 1; iter < maxIter_ICM; iter++) {
    
    Uy = calXenergy2D_sp(y, Adj, T, beta0_y, beta_y);
    Uy = normalise(exp(-Uy), 1, 1);
    Uy = -log(Uy);
    //Uy = Uy;
    

    U = Uy + Uz;
    vec Umin = min(U, 1);
    uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      break;
    }
  }

  if (iter == maxIter_ICM) {
    Iteration = iter - 1;
  } else {
    Iteration = iter;
  }

  vec energy = Energy.subvec(1, Iteration);

  mat pygn = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE.
  

  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("U") = U,
    Rcpp::Named("Uy") = Uy,
    Rcpp::Named("Uz") = Uz,
    Rcpp::Named("pygn") = pygn,
    Rcpp::Named("energy") = energy);
  
  return output; 
  
}


double calELBO2(arma::ivec y, arma::sp_mat Adj, int T, int K, int q, arma::vec& N, arma::mat& X_count, arma::mat& loglambda_mu, arma::mat& loglambda_s2, arma::vec& invLambda, arma::vec& Lambda, arma::mat& Wz_mu, arma::vec& alpha, arma::mat& z_mu, arma::mat& z_s2, arma::mat& transW_invLambda_W, arma::cube& b2, arma::cube& m2, arma::mat mu, arma::mat& pi, arma::mat& pi0, arma::mat& epsilon2, arma::mat& sigma2, arma::mat& beta, arma::mat& beta2, arma::vec& beta0, double beta_y, bool verbose){
  // ELBO
  
  int n = X_count.n_rows;
  int p = X_count.n_cols;
  
  mat loglambda_mu_logN = loglambda_mu.each_col() + log(N);   
  double possion_term = accu(X_count%loglambda_mu_logN) - accu(repmat(N, 1, p)%exp(loglambda_mu + 0.5*loglambda_s2)); 
  
  if(verbose){
  }
  
  double dimreduction_term = 0;
  vec d_invLambda = loglambda_s2 * invLambda;
  vec loglambda_mu_alpha_We = zeros<vec>(p);
  for (int i = 0; i < n; i++){
    loglambda_mu_alpha_We = trans(loglambda_mu.row(i)) - alpha(i) - Wz_mu.col(i);
    dimreduction_term = dimreduction_term  - 0.5*accu(loglambda_mu_alpha_We%invLambda%loglambda_mu_alpha_We) - 0.5*d_invLambda(i) - 0.5*accu(diagvec(transW_invLambda_W)%trans(z_s2.row(i))); 
  }
  dimreduction_term = dimreduction_term - 0.5*n*accu(log(Lambda));
  
  if(verbose){
  }
  
  mat term = zeros<mat>(4, 1);
  term(0,0) = accu(pi % repmat(log(trans(pi0)), n, 1));
  term(3,0) = -accu(pi%log(pi)) + 0.5*accu(log(loglambda_s2)) + 0.5*accu(log(z_s2));      
  
  for (int i = 0; i<n; i++){
    mat b_it = b2.slice(i);
    mat m_it = m2.slice(i);
    for (int t = 0; t<T; t++){
      for (int j = 0; j<q; j++){
        term(1, 0) = term(1, 0) + pi(i,t) * ( -0.5 * log(epsilon2(i,j)) 
                                                - 0.5* ( pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m_it.row(j)) + z_s2(i,j) )/epsilon2(i,j));
        for (int k = 0; k< K; k++){
          term(2, 0) = term(2, 0) + pi(i,t) * ( -0.5 * log(sigma2(k,j)) 
                                                  - 0.5 * ( pow((b_it(j,k) - mu(k,j)),2) + m_it(j,k) )/sigma2(k,j)  );
        }                      
        for (int k = 0; k< K; k++){
          term(3) = term(3) + pi(i,t) * ( 0.5*log(m_it(j,k)) );
        } 
        
      }
    }
  }
  
  if(verbose){
  }
  
  mat Uy = calXenergy2D_sp(y, Adj, T, beta0, beta_y);
  mat pygn = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE.
    
  double pots_term = accu(pi%(-log(pygn)));
  
  
  double ELBO = accu(term) + possion_term + dimreduction_term + pots_term;
  return(ELBO);
}




List STCountDecon(arma::mat& X_count, arma::ivec& y, arma::sp_mat& Adj, arma::mat Z, arma::mat W, int K, int T, arma::mat beta, arma::mat mu, arma::mat pi, arma::mat pi0, arma::vec Lambda, arma::vec alpha, arma::mat epsilon2, arma::vec beta0_y, arma::vec beta_y_grid, double beta_y = 0, int max_iter = 100, int maxIter_ICM = 50, double eps = 1e-5,  bool verbose = true, bool fixed_loglambda_mu = false, bool fixed_loglambda_s2 = false, bool fixed_z_mu = false, bool fixed_z_s2 = false,  bool fixed_m = false, bool fixed_b = false, bool fixed_pi = false, bool fixed_W = false, bool fixed_alpha = false, bool fixed_Lambda = false,  bool fixed_pi0 = false, bool fixed_epsilon2 = false, bool fixed_mu = false, bool fixed_sigma2 = false, bool fixed_beta2 = false, bool fixed_beta_y = false, bool alpha_spot_specific = false,  double cutoff = 1){
  
  
  int n = X_count.n_rows;
  int p = X_count.n_cols;
  int q = Z.n_cols;
  int iter = 1;

  if(verbose){
  }

  rowvec ELBO = zeros<rowvec>(1, 100000);
  ELBO(0) = -datum::inf;
  mat beta2 = beta%beta;
  
  mat m_initial = ones<mat>(q, K)*0.1;
  cube m2 = cube(q, K, n);
  cube b2 = cube(q, K, n);
  for (int i = 0; i<n; i++){
    b2.slice(i) = trans(mu);
    m2.slice(i) = m_initial;
  }
  
  mat sigma2 = ones<mat>(K, q)*0.1;
  mat term = zeros<mat>(4, 1);
  mat diff1 = zeros<mat>(n,T);
  mat Uz;
  mat U;
  
  vec invLambda = 1.0/Lambda;
  mat z_mu = Z;
  mat z_s2 = ones<mat>(n, q);
  mat loglambda_mu = ones<mat>(n, p);
  mat loglambda_s2 = ones<mat>(n, p);
  vec N = sum(X_count, 1); 
  mat Wz_mu = W*trans(z_mu); // p*n
  mat transW = trans(W);
  mat transW_invLambda = transW.each_row() % trans(invLambda);
  mat transW_invLambda_W = transW_invLambda * W;
  vec LogLik = zeros<vec>(10000);
  mat pzgy = zeros<mat>(n, T);
  mat pygn = zeros<mat>(n, T);
  mat log_relative_pi = zeros<mat>(T, n);
  
  for (iter = 1; iter < max_iter; iter++){
    
    // E-step.
    // update loglambda_mu and loglambda_s2
    for (int i = 0; i < n; i++){
      for (int j = 0; j < p; j++){
        
        X_count_ij = X_count(i,j);
        N_i = N(i);
        alpha_i = alpha(i);
        Wz_mu_ji = Wz_mu(j,i);
        Lambda_j = Lambda(j);
        
        struct funcdouble
        {
          double operator()(double const& x)
          {
            return -X_count_ij*(x + log(N_i)) + exp(log(N_i) + x) + 0.5*(x - alpha_i - Wz_mu_ji)*(x - alpha_i - Wz_mu_ji)/Lambda_j;
          }
        };
        
        if (fixed_loglambda_mu != true){
          const int double_bits = std::numeric_limits<double>::digits;
          std::pair<double, double> out_x = brent_find_minima(funcdouble(), -50.0, 50.0, double_bits);
          // Show all double precision decimal digits and trailing zeros.
          // if (i==0 & j < 10){
          // }
          loglambda_mu(i,j) = out_x.first;
        }                 
        
        if (fixed_loglambda_s2 != true){
          loglambda_s2(i,j) = 1.0/(N(i)*exp(loglambda_mu(i,j)) + invLambda(j));
        }
      }
    }
  
    
    // update z_mu and z_s2
    for (int i = 0; i < n; i++){
      if (fixed_z_s2 != true){
        z_s2.row(i) = 1.0/(trans(diagvec(transW_invLambda_W)) + 1.0/epsilon2.row(i));
      }
      
      if (fixed_z_mu != true){
        mat b2_i = b2.slice(i); //q-by-K
        vec pi_epsilon2_b2_beta = zeros<vec>(q);
        for (int t = 0; t < T; t++){
          pi_epsilon2_b2_beta = pi_epsilon2_b2_beta + pi(i,t) * diagmat(1/epsilon2.row(i)) * b2_i *  trans(beta.row(t));
        }
        mat X_z_mu =  transW_invLambda_W + diagmat(1/epsilon2.row(i));
        vec Y_z_mu = pi_epsilon2_b2_beta + transW_invLambda * (trans(loglambda_mu.row(i)) - alpha(i));
        z_mu.row(i) = trans(solve(X_z_mu, Y_z_mu));
      }
    }

    
    diff1 = zeros<mat>(n,T);
    
    if(verbose){
    }
    
    // ICM Step
    // energy of z
    Uz = zeros<mat>(n,T);
    for (int i = 0; i< n; i++){
      mat m2_i = m2.slice(i);
      for (int t = 0; t<T; t++){
        mat b_it = b2.slice(t);
        for (int j = 0; j< q; j++){
          diff1(i,t) = diff1(i,t) + (pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m2_i.row(j)) + z_s2(i,j) )/(2*epsilon2(i,j));
        }
      }
    }
    vec sum_log_epsilon2 = sum(log2pi + log(epsilon2),1);
    mat diff3 = zeros(n, T);
    for (int t = 0; t<T; t++){
      diff3.col(t) = sum_log_epsilon2;
    }
    Uz = diff1 + 0.5 * diff3;
    Uz = Uz.each_row() - trans(log(pi0));
    
    // update y and pygn
    List fitICM = runICM_sp(Uz, y, Adj, beta0_y, beta_y, epsilon2, maxIter_ICM);
    if (verbose){
    }
    ivec yHat = fitICM["y"]; 
    y = yHat;
    mat UHat = fitICM["U"]; 
    U = UHat;
    mat UzHat = fitICM["Uz"]; 
    Uz = UzHat;
    pzgy = Uz; // save for latter use, will delete later.
    mat pygnHat = fitICM["pygn"];
    pygn = pygnHat;
    vec energy = fitICM["energy"];
    LogLik(iter) = min(energy);
    
    
    diff1 = zeros<mat>(n,T);

    for (int i = 0; i<n; i++){
      
      // update b2
      mat b2_i = b2.slice(i);
      if (fixed_b != true){
        for (int j = 0; j<q; j++){
          mat X = zeros(K,K);
          vec Y = zeros(K,1);
          for (int t = 0; t<T; t++){
            X = X + pi(i,t)*(trans(beta.row(t))*beta.row(t)/epsilon2(i,j) + diagmat(1.0/sigma2.col(j)));
            Y = Y + pi(i,t)*(mu.col(j)/sigma2.col(j) + z_mu(i,j)/epsilon2(i,j)*trans(beta.row(t)));
          }
          b2_i.row(j) = trans(solve(X, Y));
        }
        b2.slice(i) = b2_i;
      }
      
      // update m2
      mat m2_i = m2.slice(i);
      if (fixed_m != true){
        for (int j = 0; j<q; j++){
          for (int k = 0; k< K; k++){
            m2_i(j,k) =  1.0/ sum(pi.row(i)%(trans(beta2.col(k)) / epsilon2(i,j) + 1.0/sigma2(k,j)));
          }
        }
        m2.slice(i) = m2_i;
      }
      
      
      if (fixed_pi != true){
        for (int t = 0; t<T; t++){
          mat b_it = b2.slice(i);
          for (int j = 0; j< q; j++){
            diff1(i,t) = diff1(i,t) + (pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m2_i.row(j)) + z_s2(i,j) )/(2*epsilon2(i,j));
          }
        }
          
        // update pi
        if ( fixed_pi != true){
          for (int t = 1; t<T; t++){
            log_relative_pi(t,i) = log(pi0(t)/pi0(0)) - diff1(i,t) + diff1(i,0) + log(pygn(i,t)) - log(pygn(i,0));
          }
          for (int t = 0; t<T; t++){
            pi(i, t) =1.0 /(accu(exp(log_relative_pi.col(i) - log_relative_pi(t,i))));
          }
          for (int t = 0; t < T; t++){
            if (pi(i,t) < 0.00001){
              pi(i,t)=0.00001;
            }
          }
          pi.row(i) = pi.row(i)/sum(pi.row(i));
        }
      }
    }
    
    
    // Variational M-step
    // update W
    if (fixed_W != true){
      mat X_W = zeros(q, q);
      mat Y_W = zeros(p, q);
      for (int i = 0; i < n; i++){
        X_W = X_W + trans(z_mu.row(i)) * z_mu.row(i) + diagmat(z_s2.row(i));
        Y_W = Y_W + (trans(loglambda_mu.row(i)) - alpha(i)) * z_mu.row(i);
      }
      W = trans(solve(X_W, trans(Y_W)));
      Wz_mu = W*trans(z_mu); // p*n
    }
    if (verbose){
    }
    transW = trans(W);
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W; 
    
    
    // update Lambda
    if (fixed_Lambda != true){
      mat loglambda_mu_alpha_Wz = loglambda_mu - repmat(alpha, 1, p) - trans(Wz_mu);
      mat loglambda_mu_alpha_Wz2_loglambda_s2 = loglambda_mu_alpha_Wz%loglambda_mu_alpha_Wz + loglambda_s2;
      vec diagWtdiagsumfW = sum(W.each_row()%sum(z_s2, 0)%W, 1); //p-by-1
      Lambda = trans(sum(loglambda_mu_alpha_Wz2_loglambda_s2, 0))/n + diagWtdiagsumfW/n;
      invLambda = 1.0/Lambda;
    }
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W; 
    
    
    
    // update alpha
    if (fixed_alpha != true){
      if (alpha_spot_specific == true){
        for (int i = 0; i < n; i++){
          alpha(i) = sum((trans(loglambda_mu.row(i)) - Wz_mu.col(i))%invLambda)/sum(invLambda);
        }
      }else{
        mat tloglambda_mu_Wz_mu = trans(loglambda_mu) - Wz_mu;
        alpha(0) = accu(tloglambda_mu_Wz_mu.each_col()%invLambda)/(n*sum(invLambda));
        for (int i = 1; i < n; i++){
          alpha(i) = alpha(0);
        }
      }
    }
    
    // update pi0
    if ( fixed_pi0 != true){
      for (int t = 0; t<T; t++){
        pi0(t) = sum(pi.col(t))/n;
      }
    }
    

    // update epsilon2
    if (fixed_epsilon2 != true){
      for (int j = 0; j<q; j++){
        double Exp1 = 0;
        for (int i = 0; i<n; i++){
          for (int t = 0; t<T; t++){
            mat b_it = b2.slice(i);
            mat m_it = m2.slice(i);
            Exp1 = Exp1 +  pi(i,t) * (pow((z_mu(i,j) - sum(beta.row(t)%b_it.row(j))),2) + sum( beta2.row(t)%m_it.row(j)) + z_s2(i,j)); 
          }
        }
        epsilon2(0,j) =  Exp1  / accu(pi);
      }  
    }
    if (fixed_epsilon2 != true){
      for (int j = 0; j<q; j++){
        for (int i = 1; i < n; i++){
          epsilon2(i,j) =  epsilon2(0,j); 
        }
      }
    }
    
    // update mu
    mat pi_b_it = zeros<mat>(q, K);
    if (fixed_mu != true){
      for (int i = 0; i<n; i++){
        mat b_i = b2.slice(i);
        for (int t = 0; t < T; t++){
          pi_b_it = pi_b_it + pi(i,t) * b_i;
        }
      }
      mu = trans(pi_b_it) / n;
    }
    
    
    if (fixed_sigma2 != true){
      for (int j = 0; j< q; j++){
        for (int k = 0; k< K; k++){
          double Exp2 = 0;
          for (int i = 0; i<n; i++){
            mat b_i = b2.slice(i);
            mat m_i = m2.slice(i);
            Exp2 = Exp2 + pow((b_i(j,k) - mu(k,j)),2) + m_i(j,k);
          }
          sigma2(k,j) = Exp2/n;
        }
      }
    }
  
    
    if (iter > 4){
      if ((ELBO(iter-2)  - ELBO(iter-3)) < cutoff){
        if (fixed_beta2 != true){
          for (int t = 0; t < T; t++){
            mat Dmat = zeros<mat>(K, K);
            rowvec dvec = zeros<rowvec>(K);
            for (int i = 0; i<n; i++){
              mat b_it = b2.slice(i);
              mat m_it = m2.slice(i);
              for (int j = 0; j < q; j++){
                Dmat = Dmat + pi(i,t)/epsilon2(i,j) *(trans(b_it.row(j))*b_it.row(j) + diagmat(trans(m_it.row(j))));
                dvec = dvec + z_mu(i,j)*pi(i,t)/epsilon2(i,j) * b_it.row(j);
              }
            }
            //min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.
            
            mat Amat = zeros<mat>(K+1, K);
            Amat.row(0) = ones<rowvec>(K);
            Amat.diag(-1) = ones<rowvec>(K);
            vec bvec = zeros<vec>(K+1);
            bvec(0) = 1;
            double sc = norm(Dmat);
            List out = solve_QP(Dmat/sc, trans(dvec)/sc, trans(Amat), bvec, 1);
            vec solution = out["solution"];
            beta.row(t) = trans(solution);
          }
          beta2 = beta%beta;
        }
      }
    }
    if (verbose){
    }

    // update beta_y: grid search.
    if (fixed_beta_y == false){
      int n_beta = beta_y_grid.n_elem;
      vec objBetaVec(n_beta);
      for(int k=0; k < n_beta; ++k){
        objBetaVec(k) = obj_beta(y, pi, Adj, T, beta0_y, beta_y_grid(k));
      }
      beta_y = beta_y_grid(index_max(objBetaVec));
    }
    
    ELBO(iter) = calELBO2(y, Adj, T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b2, m2, mu, pi, pi0, epsilon2, sigma2, beta, beta2, beta0_y, beta_y, verbose);
    
    if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    
    if( (ELBO(iter)  - ELBO(iter-1))  < -1e-7){
      if (verbose) Rcpp::Rcout << "Likelihood stopped increasing; returning the current fit.\n";
      break;
    }
    
    if( abs(ELBO(iter)  - ELBO(iter-1))/ abs(ELBO(iter-1)) < eps) break;
  }
  
  
  List ICMList = List::create(
    Rcpp::Named("pygn") = pygn,
    Rcpp::Named("pzgy") = pzgy,
    Rcpp::Named("U") = U,
    Rcpp::Named("Uz") = Uz,
    Rcpp::Named("log_relative_pi") = log_relative_pi
  );
  
  List intermediate = List::create(
    Rcpp::Named("diff1") = diff1,
    Rcpp::Named("log_relative_pi") = log_relative_pi
  );
  
  List resList = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("Lambda") = Lambda, 
    Rcpp::Named("W") = W,  
    Rcpp::Named("loglambda_mu") = loglambda_mu,
    Rcpp::Named("loglambda_s2") = loglambda_s2,
    Rcpp::Named("z_mu") = z_mu,
    Rcpp::Named("z_s2") = z_s2,  
    Rcpp::Named("pi") = pi,
    Rcpp::Named("pi0") = pi0,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("b") = b2,
    Rcpp::Named("m") = m2,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("epsilon2") = epsilon2,
    Rcpp::Named("beta_y") = beta_y,
    Rcpp::Named("intermediate") = intermediate,
    Rcpp::Named("ICMList") = ICMList,
    Rcpp::Named("ELBO") = ELBO.subvec(1,iter-1));
  
  return(resList);
  
}



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

arma::cube cal_mu_m(arma::mat& m, arma::mat& rho, arma::cube& nu, arma::mat& delta){
  cube mu_m(size(nu));
  int n = nu.n_slices;
  for (int i = 0; i < n; i++){      
    mu_m.slice(i) = m + rho % nu.slice(i) % delta;
  }
  return mu_m;
}

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

struct funcdouble_m
{
  double Xmij, Ni, Lambdaj, pibmu;
  funcdouble_m(double Xmij, double Ni, double Lambdaj, double pibmu) : Xmij(Xmij), Ni(Ni), Lambdaj(Lambdaj), pibmu(pibmu) {}
  
  double operator()(double const& x)
  {
    return -Xmij*(x + log(Ni)) + exp(log(Ni) + x) + 0.5*(x*x - 2*x*pibmu)/Lambdaj;
  }
};


double ud_loglambda_mu(double Xmij, double Ni, double Lambdaj, double pibmu) {
  
  const int double_bits = std::numeric_limits<double>::digits;
  std::pair<double, double> out_x = brent_find_minima(funcdouble_m(Xmij, Ni, Lambdaj, pibmu), -5.0, 15.0, double_bits);
  // Show all double precision decimal digits and trailing zeros.
  if (false){
  }
  double loglambda_mu = out_x.first;
  return loglambda_mu;
}


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


// ELBO
double calELBO3(arma::vec& N, arma::mat& X_m, arma::mat& rho, int T, arma::mat& loglambda_mu, arma::mat& loglambda_s2, 
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
  //term(4) = accu(pi % repmat(log(trans(pi0)), n, 1));
  term(4) = 0;
  term(6) = 0.5 * accu(log(loglambda_s2));// - accu(pi % log(pi));  
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



List JointSTCountDecon(arma::mat& X_m_count, arma::mat& rho, arma::mat& X_count, arma::ivec& y, arma::sp_mat& Adj, arma::mat Z, arma::mat W, int K, int T, arma::mat beta, arma::mat mu, arma::mat pi, arma::vec pi0, arma::vec Lambda, arma::vec alpha, arma::mat epsilon2, arma::vec beta0_y, arma::vec beta_y_grid, arma::vec Lambda_m, arma::mat delta, arma::mat m, double beta_y = 0, int max_iter = 100, int maxIter_ICM = 50, double eps = 1e-5,  bool verbose = true, bool fixed_loglambda_mu = false, bool fixed_loglambda_s2 = false, bool fixed_z_mu = false, bool fixed_z_s2 = false,  bool fixed_m = false, bool fixed_b = false, bool fixed_pi = false, bool fixed_W = false, bool fixed_alpha = false, bool fixed_Lambda = false,  bool fixed_pi0 = false, bool fixed_epsilon2 = false, bool fixed_mu = false, bool fixed_sigma2 = false, bool fixed_beta2 = false, bool fixed_beta_y = false, bool alpha_spot_specific = false,  double cutoff = 1){
  
  int p_m = X_m_count.n_cols;
  
  
  int n = X_count.n_rows;
  int p = X_count.n_cols;
  int q = Z.n_cols;
  int iter = 1;
  
  if(verbose){
  }
  
  
  rowvec ELBO = zeros<rowvec>(1, 100000);
  ELBO(0) = -datum::inf;
  mat beta2 = beta%beta;
  
  mat m_initial = ones<mat>(q, K)*0.1;
  cube m2 = cube(q, K, n);
  cube b2 = cube(q, K, n);
  for (int i = 0; i<n; i++){
    b2.slice(i) = trans(mu);
    m2.slice(i) = m_initial;
  }
  
  mat sigma2 = ones<mat>(K, q)*0.1;
  mat term = zeros<mat>(4, 1);
  mat diff1 = zeros<mat>(n,T);
  mat diff2 = zeros<mat>(n,T);
  mat Uz;
  mat U;
  
  vec invLambda = 1.0/Lambda;
  mat z_mu = Z;
  mat z_s2 = ones<mat>(n, q);
  mat loglambda_mu = ones<mat>(n, p);
  mat loglambda_s2 = ones<mat>(n, p);
  vec N = sum(X_count, 1); 
  mat Wz_mu = W*trans(z_mu); // p*n
  mat transW = trans(W);
  mat transW_invLambda = transW.each_row() % trans(invLambda);
  mat transW_invLambda_W = transW_invLambda * W;
  vec LogLik = zeros<vec>(10000);
  mat pzgy = zeros<mat>(n, T);
  mat pygn = zeros<mat>(n, T);
  mat log_relative_pi = zeros<mat>(T, n);
  
  //mat m = zeros<mat>(p_m, K);
  mat delta2 = delta % delta;
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
  cube mu_m = cal_mu_m(m, rho, nu, delta);
  cube bmu = cal_bmu(beta, mu_m);
  cube b2mu = cal_b2mu(beta2, rho, nu, delta2);
  
  mat pibmu = ones<mat>(n, p_m);
  
  vec invLambda_m = 1.0/Lambda_m;
  mat loglambda_mu_m = ones<mat>(n, p_m);
  mat loglambda_s2_m = ones<mat>(n, p_m);
  vec N_m = sum(X_m_count, 1); 
  N_m = ones<vec>(n);
  
  for (iter = 1; iter < max_iter; iter++){
    
    
    // marker E-step.
    for (int i = 0; i < n; i++){
      for (int j = 0; j < p_m; j++){
        // update loglambda_mu 
        rowvec pi_i = pi.row(i);
        pibmu(i, j) = sum(trans(pi_i) % bmu.slice(i).col(j));
        loglambda_mu_m(i,j) = ud_loglambda_mu(X_m_count(i,j),  N_m(i), Lambda_m(j), pibmu(i, j));
        
        //update loglambda_s2
        loglambda_s2_m(i,j) = 1.0/(N_m(i)*exp(loglambda_mu_m(i,j)) + invLambda_m(j));
      }
    } 
    
    // nu 
    ud_nu(pi, invLambda_m, loglambda_mu_m, nu, beta, beta2, m, rho, delta, delta2, nu0);
    mu_m = cal_mu_m(m, rho, nu, delta);
    nu2 = nu % nu;
    if(verbose) {
      //double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      //TERM(iter, 2) = part  - ELBO(iter);
      //ELBO(iter) = part;
    }
    
    // m
    ud_m(pi, loglambda_mu_m, nu, beta, beta2, m, rho, delta);
    mu_m = cal_mu_m(m, rho, nu, delta);
    if(verbose) {
      //double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      //TERM(iter, 3) = part  - ELBO(iter);
      //ELBO(iter) = part;
    }
    
    // delta
    ud_delta(pi, loglambda_mu_m, nu, beta, beta2, m, rho, delta);
    delta2 = delta % delta;
    mu_m = cal_mu_m(m, rho, nu, delta);
    bmu = cal_bmu(beta, mu_m);
    if(verbose) {
      //double part = calELBO(N, X_m, rho, T, loglambda_mu, loglambda_s2, Lambda, invLambda, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
      //TERM(iter, 4) = part  - ELBO(iter);
      //ELBO(iter) = part;
    }
     
    
    // update loglambda_mu and loglambda_s2
    for (int i = 0; i < n; i++){
      for (int j = 0; j < p; j++){
        
        X_count_ij = X_count(i,j);
        N_i = N(i);
        alpha_i = alpha(i);
        Wz_mu_ji = Wz_mu(j,i);
        Lambda_j = Lambda(j);
        
        struct funcdouble
        {
          double operator()(double const& x)
          {
            return -X_count_ij*(x + log(N_i)) + exp(log(N_i) + x) + 0.5*(x - alpha_i - Wz_mu_ji)*(x - alpha_i - Wz_mu_ji)/Lambda_j;
          }
        };
        
        if (fixed_loglambda_mu != true){
          const int double_bits = std::numeric_limits<double>::digits;
          std::pair<double, double> out_x = brent_find_minima(funcdouble(), -50.0, 50.0, double_bits);
          // Show all double precision decimal digits and trailing zeros.
          // if (i==0 & j < 10){
          // }
          loglambda_mu(i,j) = out_x.first;
        }                 
        
        if (fixed_loglambda_s2 != true){
          loglambda_s2(i,j) = 1.0/(N(i)*exp(loglambda_mu(i,j)) + invLambda(j));
        }
      }
    }
    
    
    // update z_mu and z_s2
    for (int i = 0; i < n; i++){
      if (fixed_z_s2 != true){
        z_s2.row(i) = 1.0/(trans(diagvec(transW_invLambda_W)) + 1.0/epsilon2.row(i));
      }
      
      if (fixed_z_mu != true){
        mat b2_i = b2.slice(i); //q-by-K
        vec pi_epsilon2_b2_beta = zeros<vec>(q);
        for (int t = 0; t < T; t++){
          pi_epsilon2_b2_beta = pi_epsilon2_b2_beta + pi(i,t) * diagmat(1/epsilon2.row(i)) * b2_i *  trans(beta.row(t));
        }
        mat X_z_mu =  transW_invLambda_W + diagmat(1/epsilon2.row(i));
        vec Y_z_mu = pi_epsilon2_b2_beta + transW_invLambda * (trans(loglambda_mu.row(i)) - alpha(i));
        z_mu.row(i) = trans(solve(X_z_mu, Y_z_mu));
      }
    }
    
    
    diff1 = zeros<mat>(n,T);
    diff2 = zeros<mat>(n,T);
    
    if(verbose){
    }
    
    // ICM Step
    // energy of z
    Uz = zeros<mat>(n,T);
    for (int i = 0; i< n; i++){
      mat m2_i = m2.slice(i);
      for (int t = 0; t<T; t++){
        mat b_it = b2.slice(t);
        for (int j = 0; j< q; j++){
          diff1(i,t) = diff1(i,t) + (pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m2_i.row(j)) + z_s2(i,j) )/(2*epsilon2(i,j));
        }
      }
    }
    for (int i = 0; i< n; i++){
      for (int t = 0; t < T; t++)
      {
        double tmp1 = 0, tmp2 = 0;
        for (int j = 0; j < p_m; j++)
        {
          tmp1 += invLambda_m(j) * pow(loglambda_mu_m(i, j) - bmu(t, j, i), 2);
          tmp2 += invLambda_m(j) * b2mu(t, j, i);
        }
        diff2(i, t) = 0.5*tmp1 + 0.5*tmp2;   
      }
    }
    
    vec sum_log_epsilon2 = sum(log2pi + log(epsilon2),1);
    mat diff3 = zeros(n, T);
    for (int t = 0; t<T; t++){
      diff3.col(t) = sum_log_epsilon2;
    }
    Uz = diff1 + 0.5 * diff3 + diff2;
    Uz = Uz.each_row() - trans(log(pi0));
    
    // update y and pygn
    List fitICM = runICM_sp(Uz, y, Adj, beta0_y, beta_y, epsilon2, maxIter_ICM);
    if (verbose){
    }
    ivec yHat = fitICM["y"]; 
    y = yHat;
    mat UHat = fitICM["U"]; 
    U = UHat;
    mat UzHat = fitICM["Uz"]; 
    Uz = UzHat;
    pzgy = Uz; // save for latter use, will delete later.
    mat pygnHat = fitICM["pygn"];
    pygn = pygnHat;
    vec energy = fitICM["energy"];
    LogLik(iter) = min(energy);
    
    
    diff1 = zeros<mat>(n,T);
    diff2 = zeros<mat>(n,T);
    
    for (int i = 0; i<n; i++){
      
      // update b2
      mat b2_i = b2.slice(i);
      if (fixed_b != true){
        for (int j = 0; j<q; j++){
          mat X = zeros(K,K);
          vec Y = zeros(K,1);
          for (int t = 0; t<T; t++){
            X = X + pi(i,t)*(trans(beta.row(t))*beta.row(t)/epsilon2(i,j) + diagmat(1.0/sigma2.col(j)));
            Y = Y + pi(i,t)*(mu.col(j)/sigma2.col(j) + z_mu(i,j)/epsilon2(i,j)*trans(beta.row(t)));
          }
          b2_i.row(j) = trans(solve(X, Y));
        }
        b2.slice(i) = b2_i;
      }
      
      // update m2
      mat m2_i = m2.slice(i);
      if (fixed_m != true){
        for (int j = 0; j<q; j++){
          for (int k = 0; k< K; k++){
            m2_i(j,k) =  1.0/ sum(pi.row(i)%(trans(beta2.col(k)) / epsilon2(i,j) + 1.0/sigma2(k,j)));
          }
        }
        m2.slice(i) = m2_i;
      }
      
      
      if (fixed_pi != true){
        for (int t = 0; t<T; t++){
          mat b_it = b2.slice(i);
          for (int j = 0; j< q; j++){
            diff1(i,t) = diff1(i,t) + (pow((z_mu(i,j) - sum(beta.row(t) % b_it.row(j))),2) + sum(beta2.row(t)%m2_i.row(j)) + z_s2(i,j) )/(2*epsilon2(i,j));
          }
        }
        
        
        for (int t = 0; t < T; t++)
        {
          double tmp1 = 0, tmp2 = 0;
          for (int j = 0; j < p_m; j++)
          {
            tmp1 += invLambda_m(j) * pow(loglambda_mu_m(i, j) - bmu(t, j, i), 2);
            tmp2 += invLambda_m(j) * b2mu(t, j, i);
          }
          diff2(i, t) = 0.5*tmp1 + 0.5*tmp2;   
        }
        
        
        // update pi
        if ( fixed_pi != true){
          for (int t = 1; t<T; t++){
            log_relative_pi(t,i) = log(pi0(t)/pi0(0)) - diff1(i,t) + diff1(i,0)  + log(pygn(i,t)) - log(pygn(i,0)) - diff2(i,t) + diff2(i,0);
          }
          for (int t = 0; t<T; t++){
            pi(i, t) =1.0 /(accu(exp(log_relative_pi.col(i) - log_relative_pi(t,i))));
          }
          for (int t = 0; t < T; t++){
            if (pi(i,t) < 0.00001){
              pi(i,t)=0.00001;
            }
          }
          pi.row(i) = pi.row(i)/sum(pi.row(i));
        }
      }
    }
    
    
    // Variational M-step
    // update W
    if (fixed_W != true){
      mat X_W = zeros(q, q);
      mat Y_W = zeros(p, q);
      for (int i = 0; i < n; i++){
        X_W = X_W + trans(z_mu.row(i)) * z_mu.row(i) + diagmat(z_s2.row(i));
        Y_W = Y_W + (trans(loglambda_mu.row(i)) - alpha(i)) * z_mu.row(i);
      }
      W = trans(solve(X_W, trans(Y_W)));
      Wz_mu = W*trans(z_mu); // p*n
    }
    if (verbose){
    }
    transW = trans(W);
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W; 
    
    
    // update Lambda
    if (fixed_Lambda != true){
      mat loglambda_mu_alpha_Wz = loglambda_mu - repmat(alpha, 1, p) - trans(Wz_mu);
      mat loglambda_mu_alpha_Wz2_loglambda_s2 = loglambda_mu_alpha_Wz%loglambda_mu_alpha_Wz + loglambda_s2;
      vec diagWtdiagsumfW = sum(W.each_row()%sum(z_s2, 0)%W, 1); //p-by-1
      Lambda = trans(sum(loglambda_mu_alpha_Wz2_loglambda_s2, 0))/n + diagWtdiagsumfW/n;
      invLambda = 1.0/Lambda;
    }
    transW_invLambda = transW.each_row() % trans(invLambda);
    transW_invLambda_W = transW_invLambda * W; 
    
    
    
    // update alpha
    if (fixed_alpha != true){
      if (alpha_spot_specific == true){
        for (int i = 0; i < n; i++){
          alpha(i) = sum((trans(loglambda_mu.row(i)) - Wz_mu.col(i))%invLambda)/sum(invLambda);
        }
      }else{
        mat tloglambda_mu_Wz_mu = trans(loglambda_mu) - Wz_mu;
        alpha(0) = accu(tloglambda_mu_Wz_mu.each_col()%invLambda)/(n*sum(invLambda));
        for (int i = 1; i < n; i++){
          alpha(i) = alpha(0);
        }
      }
    }
    
    // update pi0
    if ( fixed_pi0 != true){
      for (int t = 0; t<T; t++){
        pi0(t) = sum(pi.col(t))/n;
      }
    }
    
    
    // update epsilon2
    if (fixed_epsilon2 != true){
      for (int j = 0; j<q; j++){
        double Exp1 = 0;
        for (int i = 0; i<n; i++){
          for (int t = 0; t<T; t++){
            mat b_it = b2.slice(i);
            mat m_it = m2.slice(i);
            Exp1 = Exp1 +  pi(i,t) * (pow((z_mu(i,j) - sum(beta.row(t)%b_it.row(j))),2) + sum( beta2.row(t)%m_it.row(j)) + z_s2(i,j)); 
          }
        }
        epsilon2(0,j) =  Exp1  / accu(pi);
      }  
    }
    if (fixed_epsilon2 != true){
      for (int j = 0; j<q; j++){
        for (int i = 1; i < n; i++){
          epsilon2(i,j) =  epsilon2(0,j); 
        }
      }
    }
    
    // update mu
    mat pi_b_it = zeros<mat>(q, K);
    if (fixed_mu != true){
      for (int i = 0; i<n; i++){
        mat b_i = b2.slice(i);
        for (int t = 0; t < T; t++){
          pi_b_it = pi_b_it + pi(i,t) * b_i;
        }
      }
      mu = trans(pi_b_it) / n;
    }
    
    
    if (fixed_sigma2 != true){
      for (int j = 0; j< q; j++){
        for (int k = 0; k< K; k++){
          double Exp2 = 0;
          for (int i = 0; i<n; i++){
            mat b_i = b2.slice(i);
            mat m_i = m2.slice(i);
            Exp2 = Exp2 + pow((b_i(j,k) - mu(k,j)),2) + m_i(j,k);
          }
          sigma2(k,j) = Exp2/n;
        }
      }
    }
    
    
    if (iter > 4){
      if ((ELBO(iter-2)  - ELBO(iter-3)) < cutoff){
        if (fixed_beta2 != true){
          for (int t = 0; t < T; t++){
            mat Dmat = zeros<mat>(K, K);
            rowvec dvec = zeros<rowvec>(K);
            for (int i = 0; i<n; i++){
              mat b_it = b2.slice(i);
              mat m_it = m2.slice(i);
              for (int j = 0; j < q; j++){
                Dmat = Dmat + pi(i,t)/epsilon2(i,j) *(trans(b_it.row(j))*b_it.row(j) + diagmat(trans(m_it.row(j))));
                dvec = dvec + z_mu(i,j)*pi(i,t)/epsilon2(i,j) * b_it.row(j);
              }
            }
            
            mat Dmat_m = zeros<mat>(K, K);
            rowvec dvec_m = zeros<rowvec>(K);
            for (int k = 0; k < K; k++)
            {
              for (int k_ = k; k_ < K; k_++)
              {
                Dmat_m(k, k_) = cal_Dmat(pi, loglambda_mu_m, invLambda_m, m, rho, nu, delta, mu_m, t, k, k_);
                if (k_ != k) Dmat_m(k_, k) = Dmat_m(k, k_);  
              }
            }
            for (int k = 0; k < K; k++)
            {
              for (int i = 0; i < n; i++)
              {
                for (int j = 0; j < p_m; j++)
                {
                  dvec_m(k) += pi(i, t) * invLambda_m(j) * loglambda_mu_m(i, j) * mu_m(j, k, i);
                }
              }
            }
            Dmat = Dmat + Dmat_m;
            dvec = dvec + dvec_m;
            
            //min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.
            
            mat Amat = zeros<mat>(K+1, K);
            Amat.row(0) = ones<rowvec>(K);
            Amat.diag(-1) = ones<rowvec>(K);
            vec bvec = zeros<vec>(K+1);
            bvec(0) = 1;
            double sc = norm(Dmat);
            List out = solve_QP(Dmat/sc, trans(dvec)/sc, trans(Amat), bvec, 1);
            vec solution = out["solution"];
            beta.row(t) = trans(solution);
          }
          beta2 = beta%beta;
        }
      }
    }
    if (verbose){
    }
    
    // update beta_y: grid search.
    if (fixed_beta_y == false){
      int n_beta = beta_y_grid.n_elem;
      vec objBetaVec(n_beta);
      for(int k=0; k < n_beta; ++k){
        objBetaVec(k) = obj_beta(y, pi, Adj, T, beta0_y, beta_y_grid(k));
      }
      beta_y = beta_y_grid(index_max(objBetaVec));
    }
    
    ELBO(iter) = calELBO2(y, Adj, T, K, q, N, X_count, loglambda_mu, loglambda_s2, invLambda, Lambda, Wz_mu, alpha, z_mu, z_s2, transW_invLambda_W, b2, m2, mu, pi, pi0, epsilon2, sigma2, beta, beta2, beta0_y, beta_y, verbose) + 
      calELBO3(N_m, X_m_count, rho, T, loglambda_mu_m, loglambda_s2_m, Lambda_m, invLambda_m, pi, pi0, beta, beta2, nu, nu2, nu0, mu_m, delta, delta2);
    
    if (verbose) Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", iter, ELBO(iter), (ELBO(iter)  - ELBO(iter-1)));
    
    if( (ELBO(iter)  - ELBO(iter-1))  < -1e-7){
      if (verbose) Rcpp::Rcout << "Likelihood stopped increasing; returning the current fit.\n";
      break;
    }
    
    if( abs(ELBO(iter)  - ELBO(iter-1))/ abs(ELBO(iter-1)) < eps) break;
  }
  
  
  List ICMList = List::create(
    Rcpp::Named("pygn") = pygn,
    Rcpp::Named("pzgy") = pzgy,
    Rcpp::Named("U") = U,
    Rcpp::Named("Uz") = Uz,
    Rcpp::Named("log_relative_pi") = log_relative_pi
  );
  
  List intermediate = List::create(
    Rcpp::Named("diff1") = diff1,
    Rcpp::Named("diff2") = diff2,
    Rcpp::Named("log_relative_pi") = log_relative_pi,
    Rcpp::Named("m") = m,
    Rcpp::Named("loglambda_mu_m") = loglambda_mu_m,
    Rcpp::Named("loglambda_s2_m") = loglambda_s2_m,
    Rcpp::Named("nu") = nu,
    Rcpp::Named("nu0") = nu0,
    Rcpp::Named("bmu") = bmu,
    Rcpp::Named("b2mu") = b2mu,
    Rcpp::Named("delta") = delta,
    Rcpp::Named("Lambda_m") = Lambda_m
  );
  
  List resList = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("Lambda") = Lambda, 
    Rcpp::Named("W") = W,  
    Rcpp::Named("loglambda_mu") = loglambda_mu,
    Rcpp::Named("loglambda_s2") = loglambda_s2,
    Rcpp::Named("z_mu") = z_mu,
    Rcpp::Named("z_s2") = z_s2,  
    Rcpp::Named("pi") = pi,
    Rcpp::Named("pi0") = pi0,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("b2") = b2,
    Rcpp::Named("m2") = m2,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("epsilon2") = epsilon2,
    Rcpp::Named("beta_y") = beta_y,
    Rcpp::Named("intermediate") = intermediate,
    Rcpp::Named("ICMList") = ICMList,
    Rcpp::Named("ELBO") = ELBO.subvec(1,iter-1));
  
  return(resList);
  
}

} // namespace qrside_core

// [[Rcpp::export(name = ".count_deconvolution")]]
Rcpp::List qrside_CountDeconvolution(arma::mat& X_count, arma::mat Z, arma::mat W, int K, int T, arma::mat beta, arma::mat mu, arma::mat pi, arma::mat pi0, arma::vec Lambda, arma::vec alpha, arma::mat epsilon2, int max_iter = 100, double eps = 1e-5,  bool verbose = true, bool fixed_loglambda_mu = false, bool fixed_loglambda_s2 = false, bool fixed_z_mu = false, bool fixed_z_s2 = false,  bool fixed_m = false, bool fixed_b = false, bool fixed_pi = false, bool fixed_W = false, bool fixed_alpha = false, bool fixed_Lambda = false,  bool fixed_pi0 = false, bool fixed_epsilon2 = false, bool fixed_mu = false, bool fixed_sigma2 = false, bool fixed_beta2 = false, bool alpha_spot_specific = false,  double cutoff = 1) {
  return qrside_core::CountDeconvolution(X_count, Z, W, K, T, beta, mu, pi, pi0, Lambda, alpha, epsilon2, max_iter, eps, verbose, fixed_loglambda_mu, fixed_loglambda_s2, fixed_z_mu, fixed_z_s2, fixed_m, fixed_b, fixed_pi, fixed_W, fixed_alpha, fixed_Lambda, fixed_pi0, fixed_epsilon2, fixed_mu, fixed_sigma2, fixed_beta2, alpha_spot_specific, cutoff);
}

// [[Rcpp::export(name = ".joint_st_count_decon")]]
Rcpp::List qrside_JointSTCountDecon(arma::mat& X_m_count, arma::mat& rho, arma::mat& X_count, arma::ivec& y, arma::sp_mat& Adj, arma::mat Z, arma::mat W, int K, int T, arma::mat beta, arma::mat mu, arma::mat pi, arma::vec pi0, arma::vec Lambda, arma::vec alpha, arma::mat epsilon2, arma::vec beta0_y, arma::vec beta_y_grid, arma::vec Lambda_m, arma::mat delta, arma::mat m, double beta_y = 0, int max_iter = 100, int maxIter_ICM = 50, double eps = 1e-5,  bool verbose = true, bool fixed_loglambda_mu = false, bool fixed_loglambda_s2 = false, bool fixed_z_mu = false, bool fixed_z_s2 = false,  bool fixed_m = false, bool fixed_b = false, bool fixed_pi = false, bool fixed_W = false, bool fixed_alpha = false, bool fixed_Lambda = false,  bool fixed_pi0 = false, bool fixed_epsilon2 = false, bool fixed_mu = false, bool fixed_sigma2 = false, bool fixed_beta2 = false, bool fixed_beta_y = false, bool alpha_spot_specific = false,  double cutoff = 1) {
  return qrside_core::JointSTCountDecon(X_m_count, rho, X_count, y, Adj, Z, W, K, T, beta, mu, pi, pi0, Lambda, alpha, epsilon2, beta0_y, beta_y_grid, Lambda_m, delta, m, beta_y, max_iter, maxIter_ICM, eps, verbose, fixed_loglambda_mu, fixed_loglambda_s2, fixed_z_mu, fixed_z_s2, fixed_m, fixed_b, fixed_pi, fixed_W, fixed_alpha, fixed_Lambda, fixed_pi0, fixed_epsilon2, fixed_mu, fixed_sigma2, fixed_beta2, fixed_beta_y, alpha_spot_specific, cutoff);
}
