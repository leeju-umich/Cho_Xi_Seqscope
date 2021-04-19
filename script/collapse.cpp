#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
typedef std::vector<double> stdvec;

// arma::sp_mat col_sp(const arma::sp_mat& x,const arma::uvec& index) {
//   int n_cols = index.n_elem;
//   arma::sp_mat x_subset(x.n_rows,index.n_elem);
//   for(int i=0; i<n_cols; i++){
//     x_subset.col(i) = x.col(index(i));
//   }
//   return x_subset;}


//' Rcpp function to collapse HDMIs
//' @param m_tile a sparse matrix of counts
//' @param HDMIind a vector of indicator of HDMIs to m_tile
//' @param interv a vector of interval length
// [[Rcpp::export]]
arma::sp_mat collapse(arma::sp_mat& m_tile, NumericVector HDMIind, NumericVector interv)
{

  //put sanity check???
  int nrow = m_tile.n_rows;
  int ncol = interv.size();
  //NumericMatrix collapseM(nrow,ncol);
  arma::sp_mat collapseM(nrow,ncol-1);

  for (int i=0; i<ncol-1; i++)
    {
      int start =interv[i] ;
      int end = interv[i+1]-1;


      NumericVector ind = HDMIind[Rcpp::Range(start,end)];
      //std::cout<<ind.length()<<std::endl;
      arma::uvec inds = as<arma::uvec>(ind);

      arma::uvec onesvec = ones<uvec>(ind.length());
      arma::sp_mat tempMat = m_tile.cols(inds-onesvec );
      //arma::sp_mat tempMat = m_tile.cols(inds );
      arma::mat sumRows = tempMat * ones(tempMat.n_cols, 1);
      collapseM.col(i) = sumRows.col(0);
      //arma::sp_mat rowSums = sum(tempMat, 1);

    }
  return collapseM;

}
