#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


std::vector<double> compute_dist(double* emb, int n_pts, int npcs){

  std::vector<double> dist ((n_pts*n_pts - n_pts)/2, 0.0);

  for(int k = 0; k < npcs; k++){
    for(int i = 0; i < n_pts - 1; i++){
      for(int j = i + 1; j < n_pts; j++){
        dist[n_pts*i - (i + 1)*i/2 + j - i - 1] +=
          (emb[k*n_pts + i] - emb[k*n_pts + j])*(emb[k*n_pts + i] - emb[k*n_pts + j]);
      }
    }
  }
  for(int i = 0; i < dist.size(); i++) dist[i] = std::sqrt(dist[i]);

  return dist;
}

// [[Rcpp::export]]
NumericVector dist_Cpp(SEXP embeddings, int idx, int npcs){

  XPtr<BigMatrix> xpEmb (embeddings);
  MatrixAccessor<double> acess (*xpEmb);
  int n_emb = xpEmb->ncol();
  NumericVector out (n_emb);

  idx = idx - 1;

  std::vector<double> d1 = compute_dist(acess[idx],  xpEmb->nrow()/npcs, npcs);

  for(int i = idx + 1; i < n_emb; i++){
    double val = 0.0;
    std::vector<double> d2 = compute_dist(acess[i],  xpEmb->nrow()/npcs, npcs);
    
    for(int j = 0; j < d1.size(); j++) val += std::abs(d1[j] - d2[j]);

    out[i] = val/d1.size();
  }
  return out;
}