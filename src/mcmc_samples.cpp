#include <RcppArmadillo.h>
#include "generate_atable.h"
#include "mcmc_samples.h"
using namespace Rcpp;

//' @title Update the latent speed coefficients Tau
//' @description Update the latent speed coefficients Tau using MCMC sampling.
//'
//' This function updates the latent speed coefficients Tau using MCMC sampling.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param Eta A numeric vector of item response time parameters.
//' @param sigma2_t A numeric vector of item response time variances.
//' @param logt A numeric matrix of log-transformed response time.
//' @param alpha A numeric matrix of attribute profiles by respondents.
//' @param H A numeric matrix of item-attribute response time coefficients.
//' @return A numeric vector.
//' @export
// [[Rcpp::export]]
arma::vec Update_Tau(unsigned int N, unsigned int J, const arma::vec& Eta, const arma::vec& sigma2_t,
                    const arma::mat& logt, const arma::mat& alpha,
                    const arma::mat& H) {

 arma::vec Tau(N, arma::fill::zeros);

 double sigma2_tau_ast = 1.0 / (1.0 + arma::accu(Eta % Eta/ sigma2_t) );

 for (unsigned int i = 0; i < N; i++) {
   double temp = 0.0;
   for (unsigned int j = 0; j < J; j++) {
     double sigma2_tj = sigma2_t(j);
     temp += Eta[j] * (logt(i,j) - arma::accu(alpha.row(i) % H.row(j)))/sigma2_tj;
   }
   double mu_tau_ast = sigma2_tau_ast * temp;
   Tau(i) = R::rnorm(mu_tau_ast, sqrt(sigma2_tau_ast));
 }

 return Tau;
}


//' @title Update Sigma2_t
//' @description Update Sigma2_t using MCMC sampling.
//'
//' This function updates Sigma2_t using MCMC sampling.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param Eta A numeric vector of item response time parameters.
//' @param Tau A numeric vector of the latent speed Tau values.
//' @param logt A numeric matrix of log-transformed response time.
//' @param alpha A numeric matrix of attribute profiles by respondents.
//' @param H A numeric matrix of item-attribute response time coefficients.
//' @param a_t Shape parameter for the gamma distribution.
//' @param b_t Rate parameter for the gamma distribution.
//' @return A numeric vector of updated Sigma2_t values.
//' @export
// [[Rcpp::export]]
arma::vec Update_sigma2_t(unsigned int N, unsigned int J, const arma::vec& Eta, const arma::vec& Tau,
                         const arma::mat& logt, const arma::mat& alpha, const arma::mat& H,
                         double a_t, double b_t) {
 arma::vec sigma2_t(J, arma::fill::zeros);
 for (unsigned int j = 0; j < J; j++) {
   double temp = 0.0;
   for (unsigned int i = 0; i < N; i++) {
     temp += pow(logt(i,j) - Eta(j) * Tau(i) - arma::accu(alpha.row(i) % H.row(j)), 2.0);
   }
   double precision_t = R::rgamma(a_t + N/ 2.0, 1.0 / (b_t + temp / 2.0));
   double sigma2_tj = 1.0 / precision_t;
   sigma2_t(j) = sigma2_tj;
 }
 return sigma2_t;
}


//' @title Update Eta
//' @description Update Eta using MCMC sampling.
//'
//' This function updates Eta using MCMC sampling.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param Tau A numeric vector of latent speed coefficients.
//' @param logt A numeric matrix of log-transformed response time.
//' @param alpha A numeric matrix of attribute profiles by respondents.
//' @param H A numeric matrix of item-attribute response time coefficients.
//' @param sigma2_t A numeric vector of item response time variances.
//' @param sigma2_eta Variance parameter for the Eta sampling.
//' @return A numeric vector of updated Eta values.
//' @export
// [[Rcpp::export]]
arma::vec Update_Eta(unsigned int N, unsigned int J, const arma::vec& Tau, const arma::mat& logt,
                    const arma::mat& alpha, const arma::mat& H, const arma::vec& sigma2_t, double sigma2_eta) {

 arma::vec Eta(J, arma::fill::zeros);
 for (unsigned int j = 0; j < J; j++) {
   double sigma2_tj = sigma2_t(j);
   double sigma2_eta_ast = 1.0 / (1.0 / sigma2_eta + arma::accu(pow(Tau, 2.0)) / sigma2_tj);
   double temp = 0.0;
   for (unsigned int i = 0; i < N; i++) {
     temp += Tau(i) * (logt(i,j) - arma::accu(alpha.row(i) % H.row(j)));
   }
   double mu_eta_ast = sigma2_eta_ast * temp / sigma2_tj;
   Eta(j) = R::rnorm(mu_eta_ast, sqrt(sigma2_eta_ast));
 }

 return Eta;
}


//' @title Update Delta
//' @description Sample delta.
//'
//' This function samples from the DJP (Double Joint Probability) distribution.
//'
//' @param Bplb Lower bound for Bp.
//' @param omega Omega parameter.
//' @param sigma_B Standard deviation for B.
//' @param Bp Bp parameter.
//' @param ApApp ApApp parameter.
//' @param ApAb ApAb parameter.
//' @param ApZjp ApZjp parameter.
//' @param term_h Additional term related to H.
//' @return A double value.
//' @export
// [[Rcpp::export]]
double djpsamp(double Bplb, double omega, double sigma_B,double Bp,double ApApp,
              double ApAb, double ApZjp, double term_h){

 double term1=R::pnorm(Bplb/sigma_B,0.,1.,0,1);
 double u=R::runif(0,1);

 double var_p = 1/(ApApp+1/pow(sigma_B,2));
 double sigma_p = sqrt(var_p);
 double mu_p = var_p*(ApZjp-ApAb);

 double term2 =R::pnorm((mu_p-Bplb)/sigma_p,0.,1.,1,1);
 double lnexp =  .5*(mu_p*mu_p/var_p);
 double lnB = log(1-omega)-log(omega) + term1 - log(sigma_p) + log(sigma_B) - term2 - lnexp - term_h;
 double djp = 1.*(log(1.-u)-log(u)>lnB);
 return djp;
}


//' @title Compute Lower Bound
//' @description Compute the lower bound for a given parameter.
//'
//' This function computes the lower bound for a given parameter based on provided tables and coefficients.
//'
//' @param nClass Number of classes.
//' @param LBtable Matrix representing the lower bound constraints.
//' @param Atable Matrix representing the attribute profile table.
//' @param p Index of the parameter.
//' @param betaj Row vector of item response coefficient matrix B.
//' @param Bp Bp parameter.
//' @param Bindices Vector of latent class indices.
//' @return The computed lower bound.
//' @export
// [[Rcpp::export]]
double computeLB(unsigned int nClass,const arma::mat& LBtable,const arma::mat& Atable,unsigned int p,
                const arma::rowvec& betaj,double Bp,const arma::uvec& Bindices){
 double Bplb=-100.;
 if(p<nClass-1){
   arma::vec lbmax(2);
   arma::vec LBp=LBtable.col(p);
   arma::uvec lbinds=find(LBp==1);//find lower classes
   arma::rowvec ap=Atable.row(Bindices(p));
   arma::rowvec betajpeq0=betaj;
   betajpeq0(p)=0.;
   double gamp=arma::accu(ap%betajpeq0);
   //double gamp=arma::accu(ap%betaj)-Bp;
   arma::vec gams=(Atable.rows(lbinds))*betaj.t();
   lbmax(0)=arma::max(gams)-gamp;
   arma::uvec lbinds2=find(LBp==2);//find greater or equal classes
   //arma::vec gams2=(Atable.rows(lbinds2))*betaj.t()-Bp;
   arma::vec gamdiff=gamp-(Atable.rows(lbinds2))*betajpeq0.t();
   lbmax(1)=arma::max(gamdiff);
   Bplb=arma::max(lbmax);
 }
 if(p==nClass-1){//need this because there is no 2 in LBtable for last coefficient
   arma::vec LBp=LBtable.col(p);
   arma::uvec lbinds=find(LBp==1);//find lower classes
   arma::rowvec ap=Atable.row(Bindices(p));
   arma::rowvec betajpeq0=betaj;
   betajpeq0(p)=0.;
   double gamp=arma::accu(ap%betajpeq0);
   arma::vec gams=(Atable.rows(lbinds))*betaj.t();
   Bplb=arma::max(gams)-gamp;
 }
 return Bplb;
}


//' @title Update Alpha
//' @description Update alpha values using MCMC sampling.
//'
//' This function updates alpha values based on provided parameters using MCMC sampling.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param nClass Number of classes.
//' @param Y Matrix of categorical responses.
//' @param CLASS Vector of class labels.
//' @param PY_a Cube of conditional probabilities.
//' @param pis Vector of latent class proportion parameters.
//' @param Atable Matrix representing the attribute profile table.
//' @param alpha Matrix of attribute profiles by respondents.
//' @param d0 Dirichlet distribution parameter.
//' @param logt Matrix of log-transformed response times.
//' @param Eta Vector of item response time parameters.
//' @param Tau Vector of latent speed coefficients.
//' @param H Matrix of item-attribute response time coefficients.
//' @param sigma2_t Vector of item response time variances.
//' @return A list containing updated values.
//' @export
// [[Rcpp::export]]
Rcpp::List Update_alpha_SLCMRT(unsigned int N, unsigned int J, unsigned int nClass, const arma::mat& Y,
                              arma::vec& CLASS, arma::cube& PY_a, arma::vec& pis, const arma::mat& Atable, arma::mat& alpha,
                              const arma::vec& d0, const arma::mat& logt,const arma::vec& Eta,const arma::vec& Tau,
                              const arma::mat& H, arma::vec& sigma2_t){

 arma::vec dtilde=arma::zeros<arma::vec>(nClass);

 for(unsigned int i=0;i<N;i++){
   arma::colvec Yi = Y.col(i);
   double denominator=0.;
   arma::vec numerator(nClass);
   for(unsigned int cc=0;cc<nClass;cc++){
     double picc=1.;
     double temp = 1.;
     for(unsigned int j=0;j<J;j++){
       double Yij=Yi(j);
       double sigma2_tj = sigma2_t(j);
       temp *= R::dnorm(logt(i,j), Eta(j)*Tau(i)+sum(Atable.row(cc)%H.row(j)), std::sqrt(sigma2_tj),0);
       picc *=(PY_a(j,cc,Yij+1.)-PY_a(j,cc,Yij));
     }
     numerator(cc)=picc*temp*pis(cc);
     denominator+=picc*temp*pis(cc);
   }
   arma::vec pai=numerator/denominator;
   double class_i = rmultinomial(pai);
   CLASS(i) = class_i;
   arma::rowvec a_alpha=Atable.row(class_i);
   alpha.row(i)=a_alpha;

   dtilde(class_i)+=1.;
   pis = rDirichlet(dtilde+d0);
 }
 return Rcpp::List::create(
   // Rcpp::Named("PY_a")=PY_a,
   Rcpp::Named("pis")=pis,
   Rcpp::Named("CLASS")=CLASS.t(),
   Rcpp::Named("alpha")=alpha);
}


//' @title Simulate Y star.
//' @description Simulate a matrix of categorical responses.
//'
//' This function simulates a matrix of categorical responses using the provided parameters and an adjusted sampling technique.
//'
//' @param N Number of respondents.
//' @param J Number of items.
//' @param Y Matrix of categorical responses.
//' @param CLASS Vector of class labels.
//' @param PY_a Cube of conditional probabilities.
//' @param ABETA Matrix of adjusted item response coefficients.
//' @return A matrix of simulated categorical responses.
//' @export
// [[Rcpp::export]]
arma::mat simYast(unsigned int N, unsigned int J, const arma::mat& Y,
                 arma::vec& CLASS, arma::cube& PY_a, arma::mat& ABETA){
 arma::mat Yast=arma::zeros<arma::mat>(J,N);
 for(unsigned int i=0;i<N;i++){
   arma::vec Yasti=arma::zeros<arma::vec>(J);
   arma::colvec Yi = Y.col(i);
   double class_i=CLASS(i);
   for(unsigned int j=0;j<J;j++){
     double Yij = Yi(j);
     double aiBj=ABETA(j,class_i);
     arma::vec pYij= PY_a.tube(j,class_i);
     double Zji= rTruncNorm(aiBj,1.,Yij,pYij);
     Yasti(j) = Zji;
   }
   Yast.col(i)=Yasti;
 }
 return Yast;
}

//' @title Compute Conditional Probabilities
//' @description Compute conditional probabilities based on provided parameters.
//'
//' This function computes conditional probabilities based on the provided attribute profile table, item response coefficient matrix, and item response threshold matrix.
//'
//' @param J Number of items.
//' @param P Number of levels for each item.
//' @param nClass Number of classes.
//' @param Atable Matrix representing the attribute profile table.
//' @param BETA Matrix representing the item response coefficient matrix.
//' @param Kappa Matrix representing the item response threshold matrix.
//' @return A list containing the adjusted item response coefficient matrix (ABETA) and the cube of conditional probabilities (PY_a).
//' @export
// [[Rcpp::export]]
Rcpp::List computePYa(unsigned int J, unsigned int P, unsigned int nClass, arma::mat & Atable,
                     arma::mat & BETA, arma::mat & Kappa){
 arma::cube PY_a(J,nClass,(P+1));
 PY_a.slice(0) =  arma::zeros<arma::mat>(J,nClass);
 PY_a.slice(P) = arma::ones<arma::mat>(J,nClass);
 arma:: mat ABETA(J,nClass);
 for(unsigned int c=0;c<nClass;c++){
   arma::rowvec a_alpha=Atable.row(c);
   for(unsigned int j=0;j<J;j++){
     double aBj = arma::accu(a_alpha%BETA.row(j));
     ABETA(j,c) = aBj;
     for(unsigned int m=0;m<(P-1);m++){
       PY_a(j,c,m+1)=R::pnorm(Kappa(j,m),aBj,1.,1,0);
     }
   }
 }
 return Rcpp::List::create(
   Rcpp::Named("ABETA")=ABETA,
   Rcpp::Named("PY_a")=PY_a
 );
}


//' @title Check Strict Identifiability Condition
//' @description Check if the Delta matrix satisfies the strict identifiability condition.
//'
//' This function checks if the provided Delta matrix satisfies the strict identifiability condition based on the given parameters.
//'
//' @param DELTA Matrix representing the Delta matrix.
//' @param completematrix Matrix used to check block diagonal matrix.
//' @param M Levels of attribute.
//' @param K Number of attribute.
//' @param ncomp Number of parameters.
//' @param mainorder A vector representing the main order of the parameters.
//' @return A boolean indicating whether the strict identifiability condition is satisfied.
//' @export
// [[Rcpp::export]]
double identify_strictcheck(const arma::mat &DELTA, const arma::mat &completematrix,
                           unsigned int M, unsigned int K,
                           unsigned int ncomp,  const arma::uvec& mainorder){

 arma::mat D = DELTA.cols(mainorder);
 D.shed_col(0);

 // the matrix to check block diagonal matrix
 arma::mat ones_zero_on_diag = completematrix;

 // find the rows with row sum 0 in the interation terms
 // submatrix of these rows are used to find two block diagonal matrices
 arma::vec row_sum = (arma::sum(D.cols(K*(M-1),ncomp-2),1));
 arma::mat Dsub = D.rows(find(row_sum==0));

 // check if 2 block diagonal matrices exist (G1)
 arma::mat I_check = Dsub.cols(0,K*(M-1)-1)*ones_zero_on_diag.t();
 arma::mat I_count = I_check;
 I_count.elem(arma::find(I_check == (M-1)) ).fill(-1.0);
 arma::vec row_sumsum = (arma::sum(I_count,1));
 arma::mat Dsubsub = I_count.rows(find(row_sumsum == -1));
 arma::vec n_ek = (arma::sum(Dsubsub,0)).t();

 // check row sum > 0
 arma::vec r_sum = arma::sum(D,1);

 // check condition (G2)
 arma::mat I_check2 = D.cols(0,K*(M-1)-1)*ones_zero_on_diag.t();
 arma::mat I_count2 = I_check2;
 I_count2.elem(arma::find(I_check2 == (M-1))).fill(1.0);
 I_count2.elem(arma::find(I_check2 < (M-1))).fill(0);
 arma::vec n_ek2 = (arma::sum(I_count2,0)).t();

 double min_ek = (arma::min(n_ek)<-1);
 double min_r = (arma::min(r_sum)>0);
 double min_c = (arma::min(n_ek2)>2);

 return (min_c+min_r+min_ek>2);
}

//' @title Update H, D and B matrix
//' @description Update the HDB matrix with strict identifiability condition.
//'
//' This function updates the HDB matrix with strict identifiability condition using the provided parameters.
//' The H coefficients are negative in default.
//'
//' @param J Number of items.
//' @param N Number of respondents.
//' @param M Levels of attribute.
//' @param K Number of attributes.
//' @param nClass Number of classes.
//' @param sigma2_h Variance parameter for H.
//' @param alpha Matrix of attribute profiles by respondents.
//' @param sigma2_t Vector of variance parameters for the response time.
//' @param logt Matrix of log-transformed response time.
//' @param Atable Matrix representing the attribute profile table.
//' @param H Matrix of item-attribute response time coefficients.
//' @param Eta Vector of item response time coefficients.
//' @param Tau Vector of latent speed coefficients.
//' @param DELTA Matrix representing the item-attribute structure.
//' @param omega Hyperparameter.
//' @param Yast Matrix representing simulated responses.
//' @param BETA Matrix representing the item response coefficient matrix.
//' @param ABETA Matrix representing the item response coefficient matrix.
//' @param sigma_B Variance parameter for B.
//' @param main Vector representing the column index of the main-effect parameters.
//' @param mainorder Vector representing the index of the beta coefficients.
//' @param completematrix Matrix used to check the identifiablity constraints.
//' @param LBtable Matrix representing the lower bound table.
//' @param Bindices Vector representing the latent class indices.
//' @return None. The function updates the HDB matrix in place.
//' @export
// [[Rcpp::export]]
void Update_HDB_strict(unsigned int J, unsigned int N, unsigned int M, unsigned int K,
                      unsigned int nClass, double sigma2_h,
                      const arma::mat& alpha, const arma::vec& sigma2_t, const arma::mat& logt,
                      const arma::mat& Atable, arma::mat& H, const arma::vec& Eta, const arma::vec& Tau,
                      arma::mat& DELTA, const double omega, const arma::mat& Yast,
                      arma::mat& BETA,arma::mat& ABETA, double sigma_B,
                      const arma::uvec& main, const arma::uvec& mainorder,
                      const arma::mat& completematrix, const arma::mat& LBtable,
                      const arma::uvec& Bindices){

 double term_h;
 arma::mat ApA=arma::zeros<arma::mat>(nClass,nClass);
 arma::mat ApZ=arma::zeros<arma::mat>(J,nClass);
 ApA =  alpha.t()*alpha;
 ApZ = Yast*alpha;

 for (unsigned int j = 0; j < J; j++) {
   arma::rowvec ApZj=ApZ.row(j);
   arma::rowvec betaj=BETA.row(j);
   arma::rowvec deltaj=DELTA.row(j);
   arma::rowvec ABETAj=ABETA.row(j);
   double sigma2_tj = sigma2_t(j);

   for (unsigned int c = 0; c < nClass; c++) {
     double sigma2_h_ast = 1.0 / (1.0 / sigma2_h + arma::accu(alpha.col(c)%alpha.col(c)) / sigma2_tj);

     double temp = 0.0;
     for (unsigned int i = 0; i < N; i++) {
       temp += alpha(i, c) * (logt(i, j) - arma::accu(alpha.row(i) % H.row(j)) +
         alpha(i, c) * H(j, c) - Eta(j) * Tau(i));
     }
     double mu_h_ast = (sigma2_h_ast / sigma2_tj) * temp;

     double ApApp=ApA(c,c);
     double Bp = betaj(c);
     double ApAb=arma::accu(ApA.row(c)*betaj.t())-ApApp*Bp;
     double ApZjp = ApZj(c);
     double var_p=1./(ApApp+1./pow(sigma_B,2.));
     double mu_p = var_p*(ApZjp-ApAb);
     double sigma_p = sqrt(var_p);

     double Bplb_h = 0.0;
     // term_h = log(sqrt(sigma2_h_ast / sigma2_h)) + (0.5 * pow(mu_h_ast, 2) / sigma2_h_ast);
     // term_h = log(sqrt(sigma2_h_ast / sigma2_h)) + (0.5 * pow(mu_h_ast, 2) / sigma2_h_ast) - R::pnorm(-Bplb_h/sqrt(sigma2_h),0.,1.,1,1) + R::pnorm((mu_h_ast-Bplb_h)/sigma2_h_ast,0.,1.,1,1);
     term_h = log(sqrt(sigma2_h_ast / sigma2_h)) + (0.5 * pow(mu_h_ast, 2) / sigma2_h_ast) - R::pnorm(Bplb_h/sqrt(sigma2_h),0.,1.,1,1) + R::pnorm((Bplb_h-mu_h_ast)/sigma2_h_ast,0.,1.,1,1);


     if(c<1){
       deltaj(c)=1.;
       betaj(c) = R::rnorm(mu_p,sigma_p);
       H(j, c) = R::rnorm(mu_h_ast, sqrt(sigma2_h_ast));
     }
     else{
       arma::mat D0 = DELTA;
       D0(j,c)= 1.0 - DELTA(j,c);
       double flag = identify_strictcheck(D0, completematrix, M, K, nClass, mainorder);
       double Bplb=-10000.0;
       if(any(c == main)){
         Bplb=0.;
       }
       else{
         Bplb=computeLB(nClass,LBtable,Atable,c,betaj,Bp,Bindices);
       }
       if(flag > 0){
         double  djp = djpsamp(Bplb,omega,sigma_B,Bp,ApApp,ApAb,ApZjp,term_h);
         deltaj(c)=djp;
       }
       if(deltaj(c) > 0.){
         // betaj(c)=R::rnorm(mu_p,sigma_p);
         betaj(c)=rTruncNorm_1b(mu_p,sigma_p,1,Bplb);
         // H(j, c) = R::rnorm(mu_h_ast, sqrt(sigma2_h_ast));
         // H(j, c) = rTruncNorm_1b(mu_h_ast, sqrt(sigma2_h_ast),1,0.0);
         H(j, c) = rTruncNorm_1b(mu_h_ast, sqrt(sigma2_h_ast),0,0.0);
       }
       else{
         betaj(c)=0.;
         H(j, c) = 0.;
       }
     }

   }
   ABETAj=(Atable*betaj.t()).t();
   ABETA.row(j) = ABETAj;
   BETA.row(j) = betaj;
   DELTA.row(j) = deltaj;
 }
}

//' Update_alpha
//'
//' Update the alpha parameters in the SLCM model.
//'
//' @param N An integer specifying the number of observations.
//' @param J An integer specifying the number of items.
//' @param nClass An integer specifying the number of latent classes.
//' @param Y Matrix of categorical responses.
//' @param CLASS Vector of class labels.
//' @param PY_a A cube containing the conditional probabilities for each observation-feature combination and class.
//' @param pis Vector of latent class proportion parameters.
//' @param Atable Matrix representing the attribute profile table.
//' @param alpha Matrix of attribute profiles by respondents.
//' @param d0 A vector of hyper parameters.
//'
//' @return A list containing the updated parameters:
//' \itemize{
//'   \item pis: Updated class proportions.
//'   \item CLASS: Updated class assignments.
//'   \item alpha: Updated alpha parameters.
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List Update_alpha(unsigned int N, unsigned int J, unsigned int nClass, const arma::mat& Y,
                       arma::vec& CLASS, const arma::cube& PY_a, arma::vec& pis, const arma::mat& Atable, arma::mat& alpha,
                       const arma::vec& d0){

 arma::vec dtilde=arma::zeros<arma::vec>(nClass);

 for(unsigned int i=0;i<N;i++){
   arma::colvec Yi = Y.col(i);
   double denominator=0.;
   arma::vec numerator(nClass);
   for(unsigned int cc=0;cc<nClass;cc++){
     double picc=1.;
     for(unsigned int j=0;j<J;j++){
       double Yij=Yi(j);
       picc *=(PY_a(j,cc,Yij+1.)-PY_a(j,cc,Yij));
     }
     numerator(cc)=picc*pis(cc);
     denominator+=picc*pis(cc);
   }
   arma::vec pai=numerator/denominator;
   double class_i = rmultinomial(pai);
   CLASS(i) = class_i;
   arma::rowvec a_alpha=Atable.row(class_i);
   alpha.row(i)=a_alpha;

   dtilde(class_i)+=1.;
   pis = rDirichlet(dtilde+d0);
 }
 return Rcpp::List::create(
   Rcpp::Named("pis")=pis,
   Rcpp::Named("CLASS")=CLASS.t(),
   Rcpp::Named("alpha")=alpha);
}


//' Update_DB
//'
//' Update BETA and DELTA parameters.
//'
//' @param J An integer specifying the number of rows.
//' @param M An integer specifying the number of components.
//' @param K An integer specifying the number of factors.
//' @param alpha A matrix containing alpha values.
//' @param nClass An integer specifying the number of classes.
//' @param Yast A matrix containing Yast values.
//' @param BETA A matrix containing BETA values.
//' @param DELTA A matrix containing DELTA values.
//' @param ABETA A matrix containing ABETA values.
//' @param Atable A matrix containing Atable values.
//' @param sigma_B A double specifying the parameter sigma_B.
//' @param omega A double specifying the parameter omega.
//'
//' @return A list containing updated BETA, DELTA, and ABETA values.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List Update_DB(unsigned int J, unsigned int M,unsigned int K,const arma::mat& alpha,
                    unsigned int nClass, const arma::mat& Yast,
                    arma::mat& BETA, arma::mat& DELTA, arma::mat& ABETA, const arma::mat& Atable,
                    double sigma_B, double omega){

 arma::mat ApA=arma::zeros<arma::mat>(nClass,nClass);
 arma::mat ApZ=arma::zeros<arma::mat>(J,nClass);
 ApA =  alpha.t()*alpha;
 ApZ = Yast*alpha;

 for(unsigned int j=0;j<J;j++){
   arma::rowvec ApZj=ApZ.row(j);
   arma::rowvec betaj=BETA.row(j);
   arma::rowvec deltaj=DELTA.row(j);
   arma::rowvec ABETAj=ABETA.row(j);

   for(unsigned int c=0;c<nClass;c++){
     double ApApp=ApA(c,c);
     double Bp = betaj(c);
     double ApAb=arma::accu(ApA.row(c)*betaj.t())-ApApp*Bp;
     double ApZjp = ApZj(c);

     double var_p=1./(ApApp+1./pow(sigma_B,2.));
     double mu_p = var_p*(ApZjp-ApAb);
     double sigma_p = sqrt(var_p);
     if(c<1){
       deltaj(c)=1.;
       betaj(c) = R::rnorm(mu_p,sigma_p);
     }
     else{
       double Bplb=-1000.0;
       double djp=djpsamp(Bplb,omega,sigma_B,Bp,ApApp,ApAb,ApZjp,0.0);
       deltaj(c)=djp;
       if(djp>0){
         betaj(c)=R::rnorm(mu_p,sigma_p);
       }
       else{
         betaj(c)=0.;
       }
     }
   }

   ABETAj=(Atable*betaj.t()).t();
   ABETA.row(j) = ABETAj;
   BETA.row(j) = betaj;
   DELTA.row(j) = deltaj;
 }
 return Rcpp::List::create(
   Rcpp::Named("BETA")=BETA,
   Rcpp::Named("DELTA")=DELTA,
   Rcpp::Named("ABETA")=ABETA);
}
