// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


arma::vec compute_dist(const arma::mat &pts){
  
  arma::uword n_pts = pts.n_rows;
  
  arma::vec dist ((n_pts*n_pts - n_pts)/2);
  
  for(arma::uword i = 0; i < n_pts - 1; i++){
      for(arma::uword j = i + 1; j < n_pts; j++){
        dist(n_pts*i - (i + 1)*i/2 + j - i - 1) = std::sqrt( arma::sum(arma::square(pts.row(i) - pts.row(j))) );
    }
  }
  dist = dist*(dist.size()/arma::sum(dist));
  return dist;
}

// [[Rcpp::export]]
NumericVector dist_Cpp(SEXP embeddings, int idx, int npcs, double pow){
  
  XPtr<BigMatrix> xpEmb (embeddings);
  MatrixAccessor<double> acess (*xpEmb);
  int n_emb = xpEmb->ncol();
  
  idx = idx - 1;
  
  arma::mat emb1(acess[idx],  xpEmb->nrow()/npcs, npcs, false);
  
  arma::vec d1 = compute_dist(emb1);
  
  NumericVector out (n_emb);
  for(int i = idx + 1; i < n_emb; i++){
    arma::mat emb2(acess[i], xpEmb->nrow()/npcs, npcs, false);
    arma::vec d2 = compute_dist(emb2);
    out[i] = std::pow(arma::sum(arma::pow(d1 - d2, pow)/d1.size()), 1/pow);
  }
  return out;
}