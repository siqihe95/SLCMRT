#ifndef model_check_h
#define model_check_h

arma::mat random_D(unsigned int J,unsigned int K, double w);

arma::mat random_D_strict(int J, int K, int M, double t, arma::mat DtoQtable);

arma::vec permuteAtableIndices(unsigned int nClass,unsigned int K,unsigned int M,
                               const arma::vec& vv, const arma::vec& perm);

arma::mat permute_attr_order(unsigned int nClass, unsigned int K, unsigned int M, arma::mat perm);

arma::mat gen_Y(arma::cube& PY_a, arma::vec& CLASS,unsigned int N,unsigned int J,unsigned int P);

arma::mat simY(unsigned int N,unsigned int J,unsigned int M,
               unsigned int nClass, const arma::vec& CLASS,const arma::mat& Atable,
               const arma::mat& BETA,const arma::mat& Kappa);

double RTLL(const int& N, const int& J, const int& nClass, const arma::mat& logt,
            const arma::vec& Eta, const arma::vec& Tau, const arma::mat& Atable,
            const arma::mat& H, const arma::vec& sigma2_t, const arma::vec& pis);

double m2LL(unsigned int N,unsigned int J,unsigned int P,
            unsigned int nClass,const arma::mat& Y,
            const arma::vec& pis,const arma::cube& PY_a);

#endif
