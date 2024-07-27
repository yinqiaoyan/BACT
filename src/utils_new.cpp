#include <RcppArmadillo.h>
#include <stdio.h>
#include <cmath>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


bool IsInVec(int x, arma::vec& v) {
  bool res=false;
  int i, l=v.size();

  for (i=0; i<l; i++) {
    if (x == v(i)) {
      res = true;
      break;
    }
  }

  return res;
}



// [[Rcpp::export]]
arma::uvec FindNeighbors_Rcpp(int i, arma::vec& X_loc, arma::vec& Y_loc, std::string platform) {  // platform = "ST" or "Visium"
  int x_tmp=X_loc(i), y_tmp=Y_loc(i), num_count=0;
  bool bool_tmp;
  arma::vec inds_vec;
  arma::vec tmpY;
  arma::uvec res;

  if (platform == "ST")
  {
    inds_vec = -1 * ones<arma::vec>(4);

    tmpY = Y_loc( find(x_tmp == X_loc));
    bool_tmp = IsInVec(y_tmp - 1, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find(x_tmp == X_loc && (y_tmp - 1) == Y_loc));
      num_count += 1;
    }

    // tmpY = Y_loc( find(x_tmp == X_loc));
    bool_tmp = IsInVec(y_tmp + 1, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find(x_tmp == X_loc && (y_tmp + 1) == Y_loc));
      num_count += 1;
    }

    tmpY = Y_loc( find((x_tmp - 1) == X_loc));
    bool_tmp = IsInVec(y_tmp, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find((x_tmp - 1) == X_loc && y_tmp == Y_loc));
      num_count += 1;
    }

    tmpY = Y_loc( find((x_tmp + 1) == X_loc));
    bool_tmp = IsInVec(y_tmp, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find((x_tmp + 1) == X_loc && y_tmp == Y_loc));
      // num_count += 1;
    }
  }
  else if (platform == "Visium")
  {
    inds_vec = -1 * ones<arma::vec>(6);

    tmpY = Y_loc( find(x_tmp == X_loc));
    bool_tmp = IsInVec(y_tmp + 2, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find(x_tmp == X_loc && (y_tmp + 2) == Y_loc));
      num_count += 1;
    }

    // tmpY = Y_loc( find(x_tmp == X_loc));
    bool_tmp = IsInVec(y_tmp - 2, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find(x_tmp == X_loc && (y_tmp - 2) == Y_loc));
      num_count += 1;
    }

    tmpY = Y_loc( find((x_tmp - 1) == X_loc));
    bool_tmp = IsInVec(y_tmp + 1, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find((x_tmp - 1) == X_loc && (y_tmp + 1) == Y_loc));
      num_count += 1;
    }

    // tmpY = Y_loc( find((x_tmp - 1) == X_loc));
    bool_tmp = IsInVec(y_tmp - 1, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find((x_tmp - 1) == X_loc && (y_tmp - 1) == Y_loc));
      num_count += 1;
    }

    tmpY = Y_loc( find((x_tmp + 1) == X_loc));
    bool_tmp = IsInVec(y_tmp + 1, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find((x_tmp + 1) == X_loc && (y_tmp + 1) == Y_loc));
      num_count += 1;
    }

    // tmpY = Y_loc( find((x_tmp + 1) == X_loc));
    bool_tmp = IsInVec(y_tmp - 1, tmpY);
    if (bool_tmp) {
      // inds_vec.insert_rows(num_count, 1);
      inds_vec(num_count) = conv_to<int>::from(find((x_tmp + 1) == X_loc && (y_tmp - 1) == Y_loc));
    }
  }

  res = conv_to<arma::uvec>::from(inds_vec( find(inds_vec > -1) ));

  return res;
}




List FindFiniteSet_c_Rcpp(double u, double dpAlpha, arma::vec dpXi){   // new changed!!! (remove M0) date:20221116
  // Rcpp::Environment env2("package:stats");
  // Rcpp::Function Rrbeta = env2["rbeta"];
  int m = 0, num_count = 0, dpXi_size = dpXi.size();
  double mass = 1, weight, rest;
  // double tmpBeta;
  arma::vec finiteSet;

  while (true) {
    // if (m > M0 - 1) {   // new changed!!! date:20221116
    if (m > dpXi_size - 1) {
      dpXi.insert_rows(m, 1);
      dpXi(m) = R::rbeta(1, dpAlpha);
    }

    weight = mass * dpXi(m);
    if (weight >= u) {
      finiteSet.insert_rows(num_count, 1);
      finiteSet(num_count) = m;
      num_count += 1;
    }

    rest = mass * (1 - dpXi(m));
    if (rest < u) {
      break;
    }
    else {
      mass = rest;
      m += 1;
    }
  }

  return Rcpp::List::create(Rcpp::Named("finiteSet") = finiteSet + 1,
                            Rcpp::Named("dpXi") = dpXi);
}




double ComputeLogQ_c_Rcpp(int i, int m, double tmpBeta,
                          double psi1, double psi2, double geoq,
                          arma::uvec& c_vec, arma::uvec& L1_Ids, arma::uvec& L2_Ids,
                          arma::uvec& tmpNei, arma::vec& X_loc, arma::vec& Y_loc) {

  double sum_nei_m=0, p1_tmp, p2_tmp, denom;
  double res;

  sum_nei_m = sum(c_vec(tmpNei) == m) * tmpBeta;
  denom = geoq / (1 - pow(1 - geoq, m - 1) * geoq);
  if (L1_Ids(i) == m)
    p1_tmp = psi1;
  else
    p1_tmp = (1 - psi1) * pow(1 - geoq, L1_Ids(i) - 1) * denom;
  if (L2_Ids(i) == m)
    p2_tmp = psi2;
  else
    p2_tmp = (1 - psi2) * pow(1 - geoq, L2_Ids(i) - 1) * denom;

  res = log(p1_tmp) + log(p2_tmp) + sum_nei_m;

  return res;
}


// [[Rcpp::export]]
double ComputePottsDist_Rcpp(double tmpBeta, arma::vec c_vec,
                             arma::vec& X_loc, arma::vec& Y_loc,
                             std::string platform) {
  // Rcpp::Environment env = Rcpp::Environment::global_env();
  // Rcpp::Function FindNeighbors = env["FindNeighbors"];
  int i, numOfData = c_vec.size();
  double sum_nei_m=0;
  arma::uvec tmpNei, tmpNei2;

  for (i=0; i<numOfData; i++) {
    tmpNei = FindNeighbors_Rcpp(i, X_loc, Y_loc, platform);
    sum_nei_m += sum(c_vec(tmpNei) == c_vec(i));
  }

  return sum_nei_m / 2 * tmpBeta;
}


// [[Rcpp::export]]
double ComputePottsDist_sc_Rcpp(double tmpBeta, arma::vec c_vec,
                                arma::vec& X_loc, arma::vec& Y_loc,
                                arma::mat& neighbor_indices) {
  // Rcpp::Environment env = Rcpp::Environment::global_env();
  // Rcpp::Function FindNeighbors = env["FindNeighbors"];
  int i, numOfData = c_vec.size();
  double sum_nei_m=0;
  arma::uvec tmpNei, tmpNei2;

  for (i=0; i<numOfData; i++) {
    tmpNei = arma::conv_to<arma::uvec>::from(neighbor_indices.row(i) - 1);
    sum_nei_m += sum(c_vec(tmpNei) == c_vec(i));
  }

  return sum_nei_m / 2 * tmpBeta;
}





// [[Rcpp::export]]
List UpdateClIds(double dpAlpha, double pottsBeta,     // new changed!!! date:20221116 (remove M0)
                 double psi1, double psi2, double geoq,
                 arma::vec Pi, arma::vec dpXi,
                 arma::uvec clIds, arma::uvec& L1_Ids, arma::uvec& L2_Ids,
                 arma::vec& X_loc, arma::vec& Y_loc,
                 std::string platform) {
  // Rcpp::Environment env = Rcpp::Environment::global_env();
  // Rcpp::Environment env2("package:base");
  // Rcpp::Function FindFiniteSet_c = env["FindFiniteSet_c"];
  // Rcpp::Function ComputeLogQ_c_v1 = env["ComputeLogQ_c_v1"];
  // Rcpp::Function Rsample = env2["sample"];
  // bool dmh = false;
  int i, m, numOfData = clIds.size(), set_size;
  double u_i;
  arma::vec finiteSet_c, q_c;
  arma::uvec vec_nei;
  List res;
  //
  for (i = 0; i < numOfData; i++) {
    // i = 0;
    u_i = R::runif(0, Pi(clIds(i) - 1));   // randu() * Pi(clIds(i))

    // res = FindFiniteSet_c_Rcpp(u_i, dpAlpha, M0, dpXi);
    res = FindFiniteSet_c_Rcpp(u_i, dpAlpha, dpXi);   // new changed!!! date:20221116
    finiteSet_c = as<arma::vec>(res["finiteSet"]);
    dpXi = as<arma::vec>(res["dpXi"]);

    // cout << "q_c" << endl;
    // cout << finiteSet_c << endl;
    set_size = finiteSet_c.size();
    q_c = zeros<arma::vec>(set_size);
    vec_nei = FindNeighbors_Rcpp(i, X_loc, Y_loc, platform);

    for (m = 0; m < set_size; m++) {
      q_c(m) = ComputeLogQ_c_Rcpp(i, finiteSet_c(m), pottsBeta, psi1, psi2,
          geoq, clIds, L1_Ids, L2_Ids, vec_nei, X_loc, Y_loc);
    }
    // // cout << q_c << endl;
    // //
    q_c = q_c - max(q_c);
    q_c = exp(q_c) / sum(exp(q_c));
    // // cout << q_c << endl;
    if (set_size > 1) {
      clIds(i) = Rcpp::RcppArmadillo::sample(finiteSet_c, 1, false, q_c)(0);
    }
    else {
      clIds(i) = finiteSet_c(0);
    }

  }
  return Rcpp::List::create(Rcpp::Named("clIds") = clIds,
                            Rcpp::Named("dpXi") = dpXi);
  // return u_vec;
}



// [[Rcpp::export]]
List UpdateClIds_sc(double dpAlpha, double pottsBeta,     // new changed!!! date:20221116 (remove M0)
                    double psi1, double psi2, double geoq,
                    arma::vec Pi, arma::vec dpXi,
                    arma::uvec clIds, arma::uvec& L1_Ids, arma::uvec& L2_Ids,
                    arma::vec& X_loc, arma::vec& Y_loc,
                    arma::mat& neighbor_indices) {

  int i, m, numOfData = clIds.size(), set_size;
  double u_i;
  arma::vec finiteSet_c, q_c;
  arma::uvec vec_nei;
  List res;

  for (i = 0; i < numOfData; i++) {
    // i = 0;
    u_i = R::runif(0, Pi(clIds(i) - 1));   // randu() * Pi(clIds(i))

    // res = FindFiniteSet_c_Rcpp(u_i, dpAlpha, M0, dpXi);
    res = FindFiniteSet_c_Rcpp(u_i, dpAlpha, dpXi);   // new changed!!! date:20221116
    finiteSet_c = as<arma::vec>(res["finiteSet"]);
    dpXi = as<arma::vec>(res["dpXi"]);

    // cout << "q_c" << endl;
    // cout << finiteSet_c << endl;
    set_size = finiteSet_c.size();
    q_c = zeros<arma::vec>(set_size);
    vec_nei = arma::conv_to<arma::uvec>::from(neighbor_indices.row(i) - 1);

    for (m = 0; m < set_size; m++) {
      q_c(m) = ComputeLogQ_c_Rcpp(i, finiteSet_c(m), pottsBeta, psi1, psi2,
          geoq, clIds, L1_Ids, L2_Ids, vec_nei, X_loc, Y_loc);
    }
    // // cout << q_c << endl;
    // //
    q_c = q_c - max(q_c);
    q_c = exp(q_c) / sum(exp(q_c));
    // // cout << q_c << endl;
    if (set_size > 1) {
      clIds(i) = Rcpp::RcppArmadillo::sample(finiteSet_c, 1, false, q_c)(0);
    }
    else {
      clIds(i) = finiteSet_c(0);
    }

  }

  return Rcpp::List::create(Rcpp::Named("clIds") = clIds,
                            Rcpp::Named("dpXi") = dpXi);
  // return u_vec;
}





//// L1
static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

arma::vec dmvnrm_arma_mc(arma::mat const &x,
                         arma::rowvec const &mean,
                         arma::mat const &sigma,
                         bool const logd = false) {  // int const cores = 1
  using arma::uword;
  // omp_set_num_threads(cores);
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  // #pragma omp parallel for schedule(static) private(z)
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}



