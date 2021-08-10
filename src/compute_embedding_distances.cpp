#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


std::vector<std::vector<int>> line_to_triang(int k1, int k2, int n_cells){
  
  std::vector<std::vector<int>> out (2);
  
  if(k1 > k2) stop("Error: k1 > k2");
  
  if((k1 < 0) || (k1 >= (n_cells*(n_cells - 1))/2) || (k2 < 0) || (k2 >= (n_cells*(n_cells - 1))/2)){
    stop("Error: k1 or k2 out of bounds");
  }
 
  int pos_ud = 0;
  int i = 0;
  while(pos_ud <= k2 && i < n_cells -1){
    if( (i == n_cells -1) || ( (pos_ud <= k1) && (k1 < pos_ud + n_cells - i - 1))){
      out[0] = {i, k1 - pos_ud + i + 1};
    }
    if((i == n_cells -1) || (k2 < pos_ud + n_cells - i - 1)){
      out[1] = {i, k2 - pos_ud + i + 1};
    }
    pos_ud += n_cells - i - 1;
    i++;
  }
  
  if((k1 != n_cells*out[0][0] - (out[0][0] + 1)*out[0][0]/2 + out[0][1] - out[0][0] - 1) ||
      (k2 != n_cells*out[1][0] - (out[1][0] + 1)*out[1][0]/2 + out[1][1] - out[1][0] - 1))
  {
    stop("Error: something is wrong");
  }
  return out;
}


std::vector<double> compute_dist(double* emb, std::vector<std::vector<int>> ends, int n_cells, int npcs){

  
  int i1 = ends[0][0];
  int j1 = ends[0][1];
  
  int i2 = ends[1][0];
  int j2 = ends[1][1];
  
  int k1 = n_cells*i1 - (i1 + 1)*i1/2 + j1 - i1 - 1;
  int k2 = n_cells*i2 - (i2 + 1)*i2/2 + j2 - i2 - 1;
  std::vector<double> dist (k2 - k1 + 1, 0.0);
  
  
  for(int col = 0; col < npcs; col++){
    for(int i = i1; i <= i2; i++){
      int left_bound = (i == i1) ? j1 : (i + 1);
      int right_bound = (i == i2) ? j2 : (n_cells -1);
      for(int j = left_bound; j <= right_bound; j++){
        dist[n_cells*i - (i + 1)*i/2 + j - i - 1 - k1] +=
          (emb[col*n_cells + i] - emb[col*n_cells + j])*(emb[col*n_cells + i] - emb[col*n_cells + j]);
      }
    }
  }
  for(int i = 0; i < dist.size(); i++) dist[i] = std::sqrt(dist[i]);

  return dist;
}

// [[Rcpp::export]]
NumericVector partial_dist_Cpp(SEXP embeddings, int k1, int k2, int npcs, int stride){
  XPtr<BigMatrix> xpEmb (embeddings);
  MatrixAccessor<double> acess (*xpEmb);
  int n_emb = xpEmb->ncol();
  int n_cells = xpEmb->nrow()/npcs;
  NumericVector out ((n_emb*(n_emb - 1))/2);
  
  int n_strides = (k2 - k1 + 1) / stride;
  
  if(n_strides*stride < (k2 - k1 + 1)) n_strides++;
  
  // passing to 0-index
  k1--;
  k2--;
  
  for(int s = 0; s < n_strides; s++){
    std::vector<std::vector<double>> distance (n_emb);
    int stride_size = std::min(k2 - s*stride - k1 + 1, stride);
    
    std::vector<std::vector<int>> ends = \
      line_to_triang(k1 + s*stride, k1 + s*stride + stride_size - 1, n_cells);
    
    for(int n = 0; n < n_emb; n++){
      distance[n] = compute_dist(acess[n], ends, n_cells, npcs);
    }
    
    for(int i = 0; i < n_emb; i ++){
      for(int j = i + 1; j < n_emb; j ++){
        for(int k = 0; k < stride_size; k++){
          out[n_emb*i - (i + 1)*i/2 + j - i - 1] += std::abs(distance[i][k] - distance[j][k]);
        }
      }
    }
  }
  return out;
}