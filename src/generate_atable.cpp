#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "generate_atable.h"
using namespace Rcpp;

//' @title Get Mode
//' @description Get the mode of a vector.
//'
//' This function calculates the mode of a numeric vector.
//'
//' @param v A numeric vector.
//' @return The mode of the input vector as an integer.
//' @export
// [[Rcpp::export]]
int getmode(NumericVector v) {
 NumericVector uniqv = unique(v);
 int maxIndex = 0;
 int maxValue = 0;

 for (int i = 0; i < uniqv.size(); i++) {
   int count = std::count(v.begin(), v.end(), uniqv[i]);
   if (count > maxValue) {
     maxIndex = i;
     maxValue = count;
   }
 }

 return uniqv[maxIndex];
}

//' @title Generate Random Normal Samples
//' @description Generate random samples from a normal distribution with given mean and standard deviation.
//'
//' This function generates random samples from a normal distribution with given mean (\code{mu}) and standard deviation (\code{sigma}).
//'
//' @param mu A numeric vector specifying the mean of the normal distribution.
//' @param sigma A numeric vector specifying the standard deviation of the normal distribution.
//' @return A numeric vector of random samples.
//' @export
// [[Rcpp::export]]
arma::vec  myrnorm(arma::vec mu, arma::vec sigma) {
 unsigned int n = mu.size();
 arma::vec a = as<arma::vec>(Rcpp::rnorm(n))%sigma + mu;
 return (a);
}


//' @title Generate Bernoulli Random Sample
//' @description Generate a Bernoulli random sample with given probability of success.
//'
//' This function generates a single Bernoulli random sample with the specified probability (\code{prob}) of success.
//'
//' @param prob A numeric value specifying the probability of success.
//' @return An integer indicating the outcome of the Bernoulli trial (0 or 1).
//' @export
// [[Rcpp::export]]
int rbernoulli_one (double prob){
 double u;
 u = R::runif(0,1);
 int z;
 if(u < prob){
   z = 1;}
 else{z = 0;}

 return z;
}



//' @title Generate Truncated Normal Sample
//' @description Generate a truncated normal sample with specified mean, standard deviation, truncation direction, and truncation bound.
//'
//' This function generates a single truncated normal sample with the specified mean (\code{mu}), standard deviation (\code{sigma}),
//' truncation direction (\code{y}), and truncation bound (\code{c}).
//'
//' @param mu A numeric value specifying the mean of the normal distribution.
//' @param sigma A numeric value specifying the standard deviation of the normal distribution.
//' @param y A numeric value specifying the truncation direction (0 for upper bound, 1 for lower bound).
//' @param c A numeric value specifying the truncation bound.
//' @return A numeric value representing the truncated normal sample.
//' @export
// [[Rcpp::export]]
double rTruncNorm_1b(double mu,double sigma, double y,double c){
 Rcpp::NumericVector Y = wrap(y);
 Rcpp::NumericVector u = Rcpp::runif(1);
 Rcpp::NumericVector sd_mu = wrap((c-mu)/sigma);
 Rcpp::NumericVector r = Rcpp::pnorm(sd_mu);
 Rcpp::NumericVector z = Rcpp::qnorm((r+u*(1-r))*Y+u*r*(1-Y));
 double Z = z(0);
 Z = Z * sigma + mu;
 if(Z>999){Z = c;}
 if(Z<-999){Z = c;}
 return Z;}



//' @title Generate Dirichlet Sample
//' @description Generate a random sample from a Dirichlet distribution with given parameters.
//'
//' This function generates a single random sample from a Dirichlet distribution with parameters specified by the vector \code{deltas}.
//'
//' @param deltas A numeric vector specifying the parameters of the Dirichlet distribution.
//' @return A numeric vector representing the generated Dirichlet sample.
//' @export
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec& deltas){
 unsigned int C = deltas.n_elem;
 arma::vec Xgamma(C);
 for(unsigned int c = 0;c < C;c++){
   Xgamma(c) = R::rgamma(deltas(c),1.0);
 }
 return Xgamma/sum(Xgamma);
}


//' @title Generate Truncated Normal Sample
//' @description Generate a truncated normal sample with specified mean, standard deviation, multinomial result, and cumulative probabilities
//' given the multinomial distribution.
//'
//' This function generates a single truncated normal sample with the specified mean (\code{mean}), standard deviation (\code{sd}),
//' response category (\code{w}), and cumulative probabilities (\code{ps}).
//'
//' @param mean A numeric value specifying the mean of the normal distribution.
//' @param sd A numeric value specifying the standard deviation of the normal distribution.
//' @param w A numeric value specifying the response category.
//' @param ps A numeric vector of cumulative probabilities.
//' @return A numeric value representing the truncated normal sample.
//' @export
// [[Rcpp::export]]
double rTruncNorm(double mean,double sd, double w,
                 const arma::vec& ps){
 double uZ = R::runif(0,1);
 double p0 = ps(w);
 double p1 = ps(w+1);
 double pz = p0 + uZ*(p1-p0);
 double Z = R::qnorm(pz, mean, sd, 1, 0);
 return Z;
}

//' @title Generate Multinomial Sample
//' @description Generate a multinomial sample with given probabilities.
//'
//' This function generates a single multinomial sample with the specified probabilities (\code{ps}).
//'
//' @param ps A numeric vector specifying the probabilities of each category.
//' @return An unsigned integer representing the generated multinomial sample.
//' @export
// [[Rcpp::export]]
 unsigned int rmultinomial(const arma::vec& ps){
   unsigned int C = ps.n_elem;
   arma::vec pps = ps/sum(ps);
   double u = R::runif(0,1);
   arma::vec cps = cumsum(pps);
   arma::vec Ips = arma::zeros<arma::vec>(C);
   Ips.elem(arma::find(cps < u) ).fill(1.0);
   return sum(Ips);
 }



//' @title Get Bijection Vector
//' @description Get the bijection vector for a given number of attributes and level of attribute.
//'
//' This function calculates the bijection vector for a specified number of attributes (K) and level of attribute (M).
//'
//' @param K Number of attributes.
//' @param M Level of attribute.
//' @return The bijection vector.
//' @examples
//' gen_bijectionvector(3, 2)
//' @export
// [[Rcpp::export]]
arma::vec gen_bijectionvector(unsigned int K,unsigned int M) {
 arma::vec vv(K);
 for(unsigned int k=0;k<K;k++){
   vv(k) = pow(M,K-k-1);
 }
 return vv;
}

//' @title Inverse Generate Bijection Vector
//' @description Inverse function of generating the bijection vector for a given number of attributes and level of attribute.
//'
//' This function calculates the inverse of the function that generates the bijection vector for a specified number of attributes (K),
//' level of attribute (M), and cumulative value (CL).
//'
//' @param K Number of attributes.
//' @param M Level of attribute.
//' @param CL Cumulative value.
//' @return A numeric vector representing the inverse bijection vector.
//' @export
// [[Rcpp::export]]
arma::vec inv_gen_bijectionvector(unsigned int K,unsigned int M,double CL){
arma::vec alpha(K);
for(unsigned int k=0;k<K;k++){
 double Mpow = pow(M,K-k-1);
 double ak=0.;
 while(( (ak+1.)*Mpow<=CL) & (ak<M)){
   ak +=1.;
 }
 alpha(k) = ak;
 CL = CL - Mpow*alpha(k);
}
return alpha;
}



//' @title Generate Inverse Bijection Table
//' @description Generate a table of inverse bijection vectors for given parameters.
//'
//' This function generates a table of inverse bijection vectors for a specified number of attributes (K), level of attribute (M), and number of classes (nClass).
//'
//' @param K Number of attributes.
//' @param M Level of attribute.
//' @param nClass Number of classes.
//' @return A matrix representing the generated inverse bijection table.
//' @export
// [[Rcpp::export]]
arma::mat CL_gen_invbijection_table(unsigned int K,unsigned int M,unsigned int nClass){
arma::mat CLtable(K,nClass);
for(unsigned int cc=0;cc<nClass;cc++){
 CLtable.col(cc) = inv_gen_bijectionvector(K,M,cc);
}
return CLtable;
}


//' @title Generate Attribute Table
//' @description Generate attribute-related tables based on specified parameters.
//'
//' This function generates attribute-related tables including the attribute table (\code{Atable}), lower bound table (\code{LBtable}),
//' index table (\code{finalcols}), D-to-Q table (\code{DtoQtable}), adjacency table (\code{adjtable}), and item index table (\code{itemindexnum})
//' based on the provided parameters including the number of classes (\code{nClass}), the number of attributes (\code{K}),
//' the number of levels of attribute (\code{M}).
//'
//' @param nClass Number of classes.
//' @param K Number of attributes.
//' @param M Levels of attribute.
//' @return A list containing the generated tables.
//' @export
// [[Rcpp::export]]
Rcpp::List GenerateAtable(unsigned int nClass,unsigned int K,unsigned int M){
unsigned int order = K;
unsigned int level = M-1;
arma::mat FullAtable(nClass,nClass);
arma::mat FullLBtable(nClass,nClass);
arma::mat FullDtoQtable(K,nClass);
arma::mat Fulladjtable(nClass,nClass);
arma::vec attrnum(nClass);
arma::vec attrlevel(nClass);
arma::vec itemindex(nClass);
for(unsigned int cr=0;cr<nClass;cr++){
 arma::vec alpha_r = inv_gen_bijectionvector(K,M,cr);
 double nof0s=0.;
 for(unsigned int k=0;k<K;k++){
   nof0s+=1.*(alpha_r(k)==0);
   FullDtoQtable(k,cr)=1.*(alpha_r(k)>0);
 }
 attrnum(cr)=1.*(double(K-nof0s-order) < 1);
 attrlevel(cr)=1.*(arma::max(alpha_r) <= double(level) || double(K-nof0s) < 2);
 double isZeroTwo=0.;
 for(unsigned int j=0;j<K;j++){
   isZeroTwo+=1.*(alpha_r[j]==0);
   isZeroTwo+=1.*(alpha_r[j]==(M-1));
 }
 itemindex(cr)=1.*(isZeroTwo==K);
 for(unsigned int cc=0;cc<nClass;cc++){
   arma::vec alpha_c = inv_gen_bijectionvector(K,M,cc);
   double mindiff = arma::min(alpha_r-alpha_c);
   FullAtable(cr,cc) = 1.*(mindiff>-1);
   double maxdiff = arma::accu(abs(alpha_c-alpha_r));
   FullLBtable(cr,cc) = 1.*(maxdiff==1)*(mindiff<0)+2.*(mindiff>-1)*(maxdiff!=0);
   Fulladjtable(cr,cc)=1.*(maxdiff==1);
 }
}
arma::uvec finalcols = find(attrnum == 1 && attrlevel==1);
arma::uvec itemindexnum = find(itemindex==1);
arma::mat Atable=FullAtable.cols(finalcols);
arma::mat LBtable=FullLBtable.cols(finalcols);
arma::mat DtoQtable=FullDtoQtable.cols(finalcols);
arma::mat adjtable=Fulladjtable.submat(finalcols,finalcols);
//return FullAtable.cols(finalcols);
return Rcpp::List::create(Rcpp::Named("Atable",Atable),
                         Rcpp::Named("LBtable",LBtable),
                         Rcpp::Named("finalcols",finalcols),
                         Rcpp::Named("DtoQtable",DtoQtable),
                         Rcpp::Named("adjtable",adjtable),
                         Rcpp::Named("itemindex",itemindexnum)
);
}
