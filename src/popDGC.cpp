#include <Rcpp.h>

//[[Rcpp::export]]
namespace Rcpp {
class dgCMatrix {
public:
  IntegerVector i, p, Dim;
  NumericVector x;
  List Dimnames;

  // constructor
  dgCMatrix(S4 mat) {
    i = mat.slot("i");
    p = mat.slot("p");
    x = mat.slot("x");
    Dim = mat.slot("Dim");
    Dimnames = mat.slot("Dimnames");
  };

  // column iterator
  class col_iterator {
  public:
    int index;
    col_iterator(dgCMatrix& g, int ind) : parent(g) { index = ind; }
    bool operator!=(col_iterator x) { return index != x.index; };
    col_iterator& operator++(int) { ++index; return (*this); };
    int row() { return parent.i[index]; };
    int col() { return column; };
    double& value() { return parent.x[index]; };
  private:
    dgCMatrix& parent;
    int column;
  };
  col_iterator begin_col(int j) { return col_iterator(*this, p[j]); };
  col_iterator end_col(int j) { return col_iterator(*this, p[j + 1]); };

};

template <> dgCMatrix as(SEXP mat) { return dgCMatrix(mat); }

template <> SEXP wrap(const dgCMatrix& sm) {
  S4 s(std::string("dgCMatrix"));
  s.slot("i") = sm.i;
  s.slot("p") = sm.p;
  s.slot("x") = sm.x;
  s.slot("Dim") = sm.Dim;
  s.slot("Dimnames") = sm.Dimnames;
  return s;
}
}

//[[Rcpp::export]]
Rcpp::dgCMatrix R_to_Cpp_to_R(Rcpp::dgCMatrix& mat, Rcpp::NumericVector r, Rcpp::NumericVector c, Rcpp::NumericVector n, int ncol){
  for (int i = 0; i < c.length(); ++i){
    mat.x.push_back(n[i]);
    mat.i.push_back(r[i]);
    for (int j = c[i]; j < ncol; ++j){
      mat.p[j]++;
    }
  }
  return mat;
}

//[[Rcpp::export]]
Rcpp::dgCMatrix r_to_Cpp_to_R(Rcpp::dgCMatrix& mat, Rcpp::NumericVector r, Rcpp::NumericVector c, Rcpp::NumericVector n, int ncol){
  for (int i = 0; i < c.length(); ++i){
    if(i % 10000 == 0)
      std::cout << i << "\n";
    mat.x.push_back(n[i]);
    mat.i.push_back(r[i]);
    for (int j = c[i]+1; j < ncol+1; ++j){
      mat.p[j]++;
    }
  }
  return mat;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix getNumMat(Rcpp::NumericVector r, Rcpp::NumericVector c, Rcpp::NumericVector n, double nrow, double ncol){
  Rcpp::NumericMatrix mat(nrow, ncol);
  for (int i = 0; i < c.length(); ++i){
    mat[r[i], c[i]] = n[i];
  }
  return mat;
}
