#ifndef generate_atable_h
#define generate_atable_h

int getmode(Rcpp::NumericVector v);

arma::vec  myrnorm(arma::vec mu, arma::vec sigma);

int rbernoulli_one (double prob);

double rTruncNorm_1b(double mu,double sigma, double y,double c);

arma::vec rDirichlet(const arma::vec& deltas);

double rTruncNorm(double mean,double sd, double w,
                  const arma::vec& ps);

unsigned int rmultinomial(const arma::vec& ps);

arma::vec gen_bijectionvector(unsigned int K,unsigned int M);

arma::vec inv_gen_bijectionvector(unsigned int K,unsigned int M,double CL);

arma::mat CL_gen_invbijection_table(unsigned int K,unsigned int M,unsigned int nClass);

Rcpp::List GenerateAtable(unsigned int nClass,unsigned int K,unsigned int M);

#endif
