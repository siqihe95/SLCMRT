#include <RcppArmadillo.h>
#include "mcmc_samples.h"
#include "generate_atable.h"
#include "model_check.h"
using namespace Rcpp;

//' @title Generate Random Delta Matrix
//' @description Generate a random Delta matrix with specified dimensions and parameters.
//'
//' This function generates a random Delta matrix with specified dimensions and parameters,
//' satisfying identifiability constraints.
//'
//' @param J Number of items.
//' @param K Number of attributes.
//' @param w Probability parameter specifying sparsity.
//' @return A random Delta matrix.
//' @export
//' @examples
//' J <- 10
//' K <- 3
//' w <- 0.5
// [[Rcpp::export]]
arma::mat random_D(unsigned int J,unsigned int K, double w){
  unsigned int P = pow(2,K);
  //Generate two identity matrices
  arma::mat I_K(K,K,arma::fill::eye);
  arma::mat Two_I_K = arma::join_cols(I_K,I_K);
  arma::mat zeros_2K = arma::zeros<arma::mat>((2*K),(P-K-1));
  arma::mat Q0 = arma::join_rows(Two_I_K ,zeros_2K);

  //generate Q1
  unsigned int R_Q1 = J-2*K;
  arma::mat U1 = arma::randu<arma::mat>(R_Q1,K);
  arma::mat Q1 = arma::zeros<arma::mat>(R_Q1,K);

  //fix elements so columns are nonzero
  arma::vec row_ks = arma::randi<arma::vec>(K,arma::distr_param(0,R_Q1-1));
  for(unsigned int k=0;k<K;k++){
    Q1(row_ks(k),k) = 1;
  }

  Q1.elem(arma::find(U1 < w)).fill(1.0);

  //Generating the remaining elements of Q in Q2
  arma::mat U2 = arma::randu<arma::mat>(R_Q1,(P-K-1));
  arma::mat Q2 = arma::zeros<arma::mat>(R_Q1,(P-K-1));
  Q2.elem(arma::find(U2 < w)).fill(1.0);

  arma::mat Q = arma::join_rows(Q1,Q2);
  arma::vec r_sum = arma::sum(Q,1);
  arma::uvec r_ind = find(r_sum==0);
  unsigned int r_num = r_ind.size();

  //fix elements so rows are nonzero
  if(r_num>0){
    arma::vec col_ks = arma::randi<arma::vec>(r_num,arma::distr_param(0,(P-2)));
    for(unsigned int i=0;i<r_num;i++){
      Q(r_ind(i),col_ks(i)) = 1;
    }
  }
  Q = arma::join_cols(Q0,Q);
  arma::vec one_J = arma::ones<arma::vec>(J);

  //Q
  arma::uvec p = arma::linspace<arma::uvec>(0,(J-1),J);
  for(unsigned int j=0;j<J;j++){
    p(j)=j;
  }
  p = arma::shuffle(p);
  return Q.rows(p);
}


//' @title Generate Random Delta Matrix with Strict Identifiability
//' @description Generate a random Delta matrix with strict identifiability condition and specified dimensions and parameters.
//'
//' This function generates a random Delta matrix with strict identifiability condition.
//'
//' @param J Number of items.
//' @param K Number of attributes.
//' @param M Number of levels per attribute.
//' @param t Probability parameter specifying sparsity.
//' @param DtoQtable Table indicating the number of attributes per item.
//' @return A random Delta matrix satisfying the strict identifiability condition.
//' @export
// [[Rcpp::export]]
arma::mat random_D_strict(int J, int K, int M, double t, arma::mat DtoQtable) {
  int nClass = pow(M, K);

  // construct three simple structure D
  arma::mat I_K = arma::zeros<arma::mat>(K, M-1);
  I_K.row(0).fill(1);
  for(int i = 1; i < K; i++) {
    arma::mat temp = arma::zeros<arma::mat>(K, M-1);
    temp.row(i).fill(1);
    I_K = join_rows(I_K, temp);
  }
  arma::mat two_I_K = arma::join_cols(arma::join_cols(I_K, I_K), I_K);
  arma::mat zero2_2K = arma::zeros<arma::mat>(3*K, nClass-(M-1)*K-1);
  arma::mat D0 = arma::join_rows(two_I_K, zero2_2K);

  // construct the remainder with random filled 0/1
  int R_D1 = J-3*K;
  arma::mat U1 = arma::randu<arma::mat>(R_D1, nClass-1);
  arma::mat D1 = arma::zeros<arma::mat>(R_D1, nClass-1);
  D1.elem(find(U1 < t)).ones();

  arma::mat D = arma::join_cols(D0, D1);
  arma::mat intercept = arma::ones<arma::mat>(J, 1);
  D = arma::join_rows(intercept, D);

  arma::uvec main = arma::find(arma::sum(DtoQtable, 0) == 1);
  arma::uvec inter = arma::find(arma::sum(DtoQtable, 0) > 1 );
  arma::uvec attorder = arma::join_cols(arma::uvec({0}), main, inter);

  // permute cols to the Atable main*interaction terms order
  D = D.cols(attorder);

  // shuffle rows
  // int n = D.n_rows;
  // arma::uvec indices = arma::shuffle(arma::linspace<arma::uvec>(0, n-1, n));
  // D = D.rows(indices);

  return D;
}



//' @title Permute Atable Indices
//' @description Permute indices of the attribute profile table based on specified parameters.
//'
//' This function permutes the indices of the attribute profile table (\code{Atable}) based on the provided parameters including
//' the number of classes (\code{nClass}), the number of attributes (\code{K}), the number of levels of attribute (\code{M}),
//' the vector of values (\code{vv}) and the permutation vector (\code{perm}).
//'
//' @param nClass Number of classes.
//' @param K Number of attributes.
//' @param M Levels of attribute.
//' @param vv Vector of values.
//' @param perm Permutation vector.
//' @return A vector representing the permuted indices of the attribute profile table.
//' @export
// [[Rcpp::export]]
arma::vec permuteAtableIndices(unsigned int nClass,unsigned int K,unsigned int M,
                              const arma::vec& vv, const arma::vec& perm){
 unsigned int order = K;
 unsigned int level = M-1;
 arma::vec vvperm(K);
 arma::vec fullorigperm=arma::linspace(0,nClass-1,nClass);
 for(unsigned int k=0;k<K;k++){
   vvperm(k)=vv(perm(k));
 }
 arma::vec attrnum(nClass);
 arma::vec attrlevel(nClass);
 arma::vec fullpermindices(nClass);
 for(unsigned int cr=0;cr<nClass;cr++){
   arma::vec alpha_r = inv_gen_bijectionvector(K,M,cr);
   double nof0s=0.;
   for(unsigned int k=0;k<K;k++){
     nof0s+=1.*(alpha_r(k)==0);
   }
   attrnum(cr)=1.*(nof0s>double(K-order)-1);
   attrlevel(cr)=1.*(arma::max(alpha_r) <= double(level) || double(K-nof0s) < 2);
   arma::vec alpha_perm(K);
   fullpermindices(cr)=arma::accu(alpha_r%vvperm);
 }
 arma::uvec finalcols = find(attrnum == 1  && attrlevel==1);
 arma::vec origperm=fullorigperm(finalcols);
 arma::vec reducedpermindices=fullpermindices(finalcols);
 arma::vec permindices(origperm.n_elem);
 for(unsigned int p=0;p<origperm.n_elem;p++){
   double origval=origperm(p);
   for(unsigned int pp=0;pp<origperm.n_elem;pp++){
     if(origval==reducedpermindices(pp)){
       permindices(p)=pp;
     }
   }
 }
 return permindices;
}


//' @title Permute Attribute Order
//' @description Permute the order of attributes in the attribute profile table.
//'
//' This function permutes the order of attributes in the attribute profile table (\code{perm}) based on the provided parameters including
//' the number of classes (\code{nClass}), the number of attributes (\code{K}), the number of levels of attribute (\code{M}).
//'
//' @param nClass Number of classes.
//' @param K Number of attributes.
//' @param M Levels of attribute.
//' @param perm Matrix representing the attribute profile table, permat = permuteGeneral(c(0:(K-1))).
//' @return A matrix representing the permuted indices of the attribute profile table.
//' @export
// [[Rcpp::export]]
arma::mat permute_attr_order(unsigned int nClass, unsigned int K, unsigned int M, arma::mat perm){
 arma::vec vv(K);
 vv(0) = 1;
 for(int i=1; i<K; i++) {
   vv(i) = vv(i-1) * M;
 }
 arma::mat perindex = arma::zeros<arma::mat>(nClass, perm.n_rows);
 for (int i = 0; i < perm.n_rows; i++) {
   arma::vec permute = arma::conv_to<arma::vec>::from(perm.row(i));
   arma::vec indices = permuteAtableIndices(nClass, K, M, vv, permute);
   perindex.col(i) = arma::conv_to<arma::colvec>::from(indices) + 1;
 }
 return perindex;
}

//' @title Generate Y Matrix
//' @description Generate a matrix of categorical responses based on probability cubes and class labels.
//'
//' This function generates a matrix of categorical responses (\code{Y}) based on the probability cubes (\code{PY_a}) and class labels (\code{CLASS}).
//'
//' @param PY_a A 3-dimensional cube containing probabilities for each class.
//' @param CLASS A numeric vector containing class labels for each respondent.
//' @param N Number of respondents.
//' @param J Number of response categories.
//' @param P Number of intervals for each response category.
//' @return A matrix representing the generated categorical responses.
//' @export
// [[Rcpp::export]]
arma::mat gen_Y(arma::cube& PY_a, arma::vec& CLASS,unsigned int N,unsigned int J,unsigned int P){
 arma::cube prob(J,N,P);
 arma::mat Y(J,N);
 for(unsigned int j=0;j<J;j++){
   for(unsigned int i=0;i<N;i++){
     unsigned int Class0 = CLASS(i);
     arma::mat PY_acc=PY_a.subcube(0,Class0,0,J-1,Class0,P);
     for(unsigned int p=0;p<P;p++){
       prob(j,i,p) = PY_a(j,Class0,p+1)-PY_a(j,Class0,p);
     }
     Y(j,i) = rmultinomial(prob.subcube(j,i,0,j,i,P-1));
   }
 }
 return Y;
}


//' @title Simulate Y Matrix
//' @description Simulate a matrix of categorical responses based on provided parameters.
//'
//' This function simulates a matrix of categorical responses (\code{Y}) based on the provided parameters including the number of respondents (\code{N}),
//' the number of items (\code{J}), the number of levels of attribute (\code{M}), the number of classes (\code{nClass}), class labels (\code{CLASS}),
//' attribute profile table (\code{Atable}), item response coefficient matrix (\code{BETA}), and item response threshold matrix (\code{Kappa}).
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param M Levels of attribute.
//' @param nClass Number of classes.
//' @param CLASS A numeric vector containing class labels for each respondent.
//' @param Atable A matrix representing the attribute profile table.
//' @param BETA A matrix representing the item response coefficient matrix.
//' @param Kappa A matrix representing the item response threshold matrix.
//' @return A matrix of (\code{N}) by (\code{J}) representing the simulated categorical responses.
//' @export
// [[Rcpp::export]]
arma::mat simY(unsigned int N,unsigned int J,unsigned int M,
              unsigned int nClass,
              const arma::vec& CLASS,const arma::mat& Atable,
              const arma::mat& BETA,const arma::mat& Kappa){
 arma::cube PY_a(J,nClass,(M+1));
 PY_a.slice(0)=arma::zeros<arma::mat>(J,nClass);
 PY_a.slice(M)=arma::ones<arma::mat>(J,nClass);
 for(unsigned int cc=0;cc<nClass;cc++){
   arma::rowvec a_alpha=Atable.row(cc);
   for(unsigned int j=0;j<J;j++){
     double aBj = arma::accu(a_alpha%BETA.row(j));
     for(unsigned int m=0;m<(M-1);m++){
       PY_a(j,cc,m+1)=R::pnorm(Kappa(j,m),aBj,1.,1,0);
     }
   }
 }
 arma::mat Y(N,J);
 for(unsigned int i=0;i<N;i++){
   double class_i = CLASS(i) - 1.0;
   for(unsigned int j=0;j<J;j++){
     arma::vec cumulativepsij=PY_a.tube(j,class_i);
     arma::vec psij(M);
     for(unsigned int m=0;m<M;m++){
       psij(m)=cumulativepsij(m+1)-cumulativepsij(m);
     }
     Y(i,j)=rmultinomial(psij);
   }
 }
 return Y;
}



//' Compute the log-likelihood ratio test statistic for response time.
//'
//' This function computes the log-likelihood ratio test statistic.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param nClass Number of latent classes.
//' @param logt A matrix of log-transformed response time.
//' @param Eta A vector of item response time parameters.
//' @param Tau A vector of latent speed coefficients.
//' @param Atable A matrix of item attribute profiles.
//' @param H A matrix of item-attribute response time coefficients.
//' @param sigma2_t A vector of item response time variances.
//' @param pis Vector of latent class proportion parameters.
//' @return Log-likelihood ratio test statistic.
//' @export
// [[Rcpp::export]]
double RTLL(const int& N, const int& J, const int& nClass, const arma::mat& logt,
            const arma::vec& Eta, const arma::vec& Tau, const arma::mat& Atable,
            const arma::mat& H, const arma::vec& sigma2_t, const arma::vec& pis) {
  double m2ll = 0.0;
  for (int i = 0; i < N; ++i) {
    double Li = 0.0;
    for (int c = 0; c < nClass; ++c) {
      double temp = 1.0;
      for (int j = 0; j < J; ++j) {
        double mu = Eta(j) * Tau(i) + sum(Atable.row(c) % H.row(j));
        temp *= R::dnorm(logt(i,j), mu, sqrt(sigma2_t(j)),0);
      }
      Li += pis(c) * temp;
    }
    m2ll += log(Li);
  }
  return -2.*m2ll;
}


//' Compute the -2*log-likelihood for responses.
//'
//' This function computes the -2*log-likelihood for responses.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param P Number of response categories.
//' @param nClass Number of latent classes.
//' @param Y A matrix of observations with respondents by items.
//' @param pis Vector of latent class proportion parameters.
//' @param PY_a Cube of the joint probability distribution.
//' @return -2*log-likelihood.
//' @export
// [[Rcpp::export]]
double m2LL(unsigned int N,unsigned int J,unsigned int P,
            unsigned int nClass,const arma::mat& Y,
            const arma::vec& pis,const arma::cube& PY_a){

  double m2ll=0.;
  for(unsigned int i=0;i<N;i++){
    arma::vec Yi=Y.col(i);
    double Li=0.;
    for(unsigned int cc=0;cc<nClass;cc++){
      double py_a=1.;
      arma::mat PY_acc=PY_a.subcube(0,cc,0,J-1,cc,P);
      for(unsigned int j=0;j<J;j++){
        double Yij=Yi(j);
        py_a*=(PY_acc(j,Yij+1)-PY_acc(j,Yij));
      }
      Li+=py_a*pis(cc);
    }
    m2ll+=log(Li);
  }
  return -2.*m2ll;
}
