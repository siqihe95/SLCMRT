#ifndef mcmc_samples_h
#define mcmc_samples_h

arma::vec Update_Tau(unsigned int N, unsigned int J, const arma::vec& Eta, const arma::vec& sigma2_t,
                     const arma::mat& logt, const arma::mat& alpha,
                     const arma::mat& H);

arma::vec Update_sigma2_t(unsigned int N, unsigned int J, const arma::vec& Eta, const arma::vec& Tau,
                          const arma::mat& logt, const arma::mat& alpha, const arma::mat& H,
                          double a_t, double b_t);

arma::vec Update_Eta(unsigned int N, unsigned int J, const arma::vec& Tau, const arma::mat& logt,
                     const arma::mat& alpha, const arma::mat& H, const arma::vec& sigma2_t, double sigma2_eta);

double djpsamp(double Bplb, double omega, double sigma_B,double Bp,double ApApp,
               double ApAb, double ApZjp, double term_h);

double computeLB(unsigned int nClass,const arma::mat& LBtable,const arma::mat& Atable,unsigned int p,
                 const arma::rowvec& betaj,double Bp,const arma::uvec& Bindices);


Rcpp::List Update_alpha_SLCMRT(unsigned int N, unsigned int J, unsigned int nClass, const arma::mat& Y,
                               arma::vec& CLASS, arma::cube& PY_a, arma::vec& pis, const arma::mat& Atable, arma::mat& alpha,
                               const arma::vec& d0, const arma::mat& logt,const arma::vec& Eta,const arma::vec& Tau,
                               const arma::mat& H, arma::vec& sigma2_t);

arma::mat simYast(unsigned int N, unsigned int J, const arma::mat& Y,
                  arma::vec& CLASS, arma::cube& PY_a, arma::mat& ABETA);


Rcpp::List computePYa(unsigned int J, unsigned int P, unsigned int nClass, arma::mat & Atable,
                      arma::mat & BETA, arma::mat & Kappa);

double identify_strictcheck(const arma::mat &DELTA, const arma::mat &completematrix,
                            unsigned int M, unsigned int K,
                            unsigned int ncomp,  const arma::uvec& mainorder);

void Update_HDB_strict(unsigned int J, unsigned int N, unsigned int M, unsigned int K,
                       unsigned int nClass, double sigma2_h,
                       const arma::mat& alpha, const arma::vec& sigma2_t, const arma::mat& logt,
                       const arma::mat& Atable, arma::mat& H, const arma::vec& Eta, const arma::vec& Tau,
                       arma::mat& DELTA, const double omega, const arma::mat& Yast,
                       arma::mat& BETA,arma::mat& ABETA, double sigma_B,
                       const arma::uvec& main, const arma::uvec& mainorder,
                       const arma::mat& completematrix, const arma::mat& LBtable,
                       const arma::uvec& Bindices);

Rcpp::List Update_alpha(unsigned int N, unsigned int J, unsigned int nClass, const arma::mat& Y,
                        arma::vec& CLASS, const arma::cube& PY_a, arma::vec& pis, const arma::mat& Atable, arma::mat& alpha,
                        const arma::vec& d0);

Rcpp::List Update_DB(unsigned int J, unsigned int M,unsigned int K,const arma::mat& alpha,
                     unsigned int nClass, const arma::mat& Yast,
                     arma::mat& BETA, arma::mat& DELTA, arma::mat& ABETA, const arma::mat& Atable,
                     double sigma_B, double omega);

#endif
