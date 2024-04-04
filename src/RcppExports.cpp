// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getmode
int getmode(NumericVector v);
RcppExport SEXP _SLCMRT_getmode(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(getmode(v));
    return rcpp_result_gen;
END_RCPP
}
// myrnorm
arma::vec myrnorm(arma::vec mu, arma::vec sigma);
RcppExport SEXP _SLCMRT_myrnorm(SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(myrnorm(mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rbernoulli_one
int rbernoulli_one(double prob);
RcppExport SEXP _SLCMRT_rbernoulli_one(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(rbernoulli_one(prob));
    return rcpp_result_gen;
END_RCPP
}
// rTruncNorm_1b
double rTruncNorm_1b(double mu, double sigma, double y, double c);
RcppExport SEXP _SLCMRT_rTruncNorm_1b(SEXP muSEXP, SEXP sigmaSEXP, SEXP ySEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rTruncNorm_1b(mu, sigma, y, c));
    return rcpp_result_gen;
END_RCPP
}
// rDirichlet
arma::vec rDirichlet(const arma::vec& deltas);
RcppExport SEXP _SLCMRT_rDirichlet(SEXP deltasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type deltas(deltasSEXP);
    rcpp_result_gen = Rcpp::wrap(rDirichlet(deltas));
    return rcpp_result_gen;
END_RCPP
}
// rTruncNorm
double rTruncNorm(double mean, double sd, double w, const arma::vec& ps);
RcppExport SEXP _SLCMRT_rTruncNorm(SEXP meanSEXP, SEXP sdSEXP, SEXP wSEXP, SEXP psSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ps(psSEXP);
    rcpp_result_gen = Rcpp::wrap(rTruncNorm(mean, sd, w, ps));
    return rcpp_result_gen;
END_RCPP
}
// rmultinomial
unsigned int rmultinomial(const arma::vec& ps);
RcppExport SEXP _SLCMRT_rmultinomial(SEXP psSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type ps(psSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultinomial(ps));
    return rcpp_result_gen;
END_RCPP
}
// gen_bijectionvector
arma::vec gen_bijectionvector(unsigned int K, unsigned int M);
RcppExport SEXP _SLCMRT_gen_bijectionvector(SEXP KSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_bijectionvector(K, M));
    return rcpp_result_gen;
END_RCPP
}
// inv_gen_bijectionvector
arma::vec inv_gen_bijectionvector(unsigned int K, unsigned int M, double CL);
RcppExport SEXP _SLCMRT_inv_gen_bijectionvector(SEXP KSEXP, SEXP MSEXP, SEXP CLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type CL(CLSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_gen_bijectionvector(K, M, CL));
    return rcpp_result_gen;
END_RCPP
}
// CL_gen_invbijection_table
arma::mat CL_gen_invbijection_table(unsigned int K, unsigned int M, unsigned int nClass);
RcppExport SEXP _SLCMRT_CL_gen_invbijection_table(SEXP KSEXP, SEXP MSEXP, SEXP nClassSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    rcpp_result_gen = Rcpp::wrap(CL_gen_invbijection_table(K, M, nClass));
    return rcpp_result_gen;
END_RCPP
}
// GenerateAtable
Rcpp::List GenerateAtable(unsigned int nClass, unsigned int K, unsigned int M);
RcppExport SEXP _SLCMRT_GenerateAtable(SEXP nClassSEXP, SEXP KSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateAtable(nClass, K, M));
    return rcpp_result_gen;
END_RCPP
}
// Update_Tau
arma::vec Update_Tau(unsigned int N, unsigned int J, const arma::vec& Eta, const arma::vec& sigma2_t, const arma::mat& logt, const arma::mat& alpha, const arma::mat& H);
RcppExport SEXP _SLCMRT_Update_Tau(SEXP NSEXP, SEXP JSEXP, SEXP EtaSEXP, SEXP sigma2_tSEXP, SEXP logtSEXP, SEXP alphaSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2_t(sigma2_tSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_Tau(N, J, Eta, sigma2_t, logt, alpha, H));
    return rcpp_result_gen;
END_RCPP
}
// Update_sigma2_t
arma::vec Update_sigma2_t(unsigned int N, unsigned int J, const arma::vec& Eta, const arma::vec& Tau, const arma::mat& logt, const arma::mat& alpha, const arma::mat& H, double a_t, double b_t);
RcppExport SEXP _SLCMRT_Update_sigma2_t(SEXP NSEXP, SEXP JSEXP, SEXP EtaSEXP, SEXP TauSEXP, SEXP logtSEXP, SEXP alphaSEXP, SEXP HSEXP, SEXP a_tSEXP, SEXP b_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type a_t(a_tSEXP);
    Rcpp::traits::input_parameter< double >::type b_t(b_tSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_sigma2_t(N, J, Eta, Tau, logt, alpha, H, a_t, b_t));
    return rcpp_result_gen;
END_RCPP
}
// Update_Eta
arma::vec Update_Eta(unsigned int N, unsigned int J, const arma::vec& Tau, const arma::mat& logt, const arma::mat& alpha, const arma::mat& H, const arma::vec& sigma2_t, double sigma2_eta);
RcppExport SEXP _SLCMRT_Update_Eta(SEXP NSEXP, SEXP JSEXP, SEXP TauSEXP, SEXP logtSEXP, SEXP alphaSEXP, SEXP HSEXP, SEXP sigma2_tSEXP, SEXP sigma2_etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2_t(sigma2_tSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_eta(sigma2_etaSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_Eta(N, J, Tau, logt, alpha, H, sigma2_t, sigma2_eta));
    return rcpp_result_gen;
END_RCPP
}
// djpsamp
double djpsamp(double Bplb, double omega, double sigma_B, double Bp, double ApApp, double ApAb, double ApZjp, double term_h);
RcppExport SEXP _SLCMRT_djpsamp(SEXP BplbSEXP, SEXP omegaSEXP, SEXP sigma_BSEXP, SEXP BpSEXP, SEXP ApAppSEXP, SEXP ApAbSEXP, SEXP ApZjpSEXP, SEXP term_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Bplb(BplbSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_B(sigma_BSEXP);
    Rcpp::traits::input_parameter< double >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< double >::type ApApp(ApAppSEXP);
    Rcpp::traits::input_parameter< double >::type ApAb(ApAbSEXP);
    Rcpp::traits::input_parameter< double >::type ApZjp(ApZjpSEXP);
    Rcpp::traits::input_parameter< double >::type term_h(term_hSEXP);
    rcpp_result_gen = Rcpp::wrap(djpsamp(Bplb, omega, sigma_B, Bp, ApApp, ApAb, ApZjp, term_h));
    return rcpp_result_gen;
END_RCPP
}
// computeLB
double computeLB(unsigned int nClass, const arma::mat& LBtable, const arma::mat& Atable, unsigned int p, const arma::rowvec& betaj, double Bp, const arma::uvec& Bindices);
RcppExport SEXP _SLCMRT_computeLB(SEXP nClassSEXP, SEXP LBtableSEXP, SEXP AtableSEXP, SEXP pSEXP, SEXP betajSEXP, SEXP BpSEXP, SEXP BindicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LBtable(LBtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type betaj(betajSEXP);
    Rcpp::traits::input_parameter< double >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Bindices(BindicesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLB(nClass, LBtable, Atable, p, betaj, Bp, Bindices));
    return rcpp_result_gen;
END_RCPP
}
// Update_alpha_SLCMRT
Rcpp::List Update_alpha_SLCMRT(unsigned int N, unsigned int J, unsigned int nClass, const arma::mat& Y, arma::vec& CLASS, arma::cube& PY_a, arma::vec& pis, const arma::mat& Atable, arma::mat& alpha, const arma::vec& d0, const arma::mat& logt, const arma::vec& Eta, const arma::vec& Tau, const arma::mat& H, arma::vec& sigma2_t);
RcppExport SEXP _SLCMRT_Update_alpha_SLCMRT(SEXP NSEXP, SEXP JSEXP, SEXP nClassSEXP, SEXP YSEXP, SEXP CLASSSEXP, SEXP PY_aSEXP, SEXP pisSEXP, SEXP AtableSEXP, SEXP alphaSEXP, SEXP d0SEXP, SEXP logtSEXP, SEXP EtaSEXP, SEXP TauSEXP, SEXP HSEXP, SEXP sigma2_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma2_t(sigma2_tSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_alpha_SLCMRT(N, J, nClass, Y, CLASS, PY_a, pis, Atable, alpha, d0, logt, Eta, Tau, H, sigma2_t));
    return rcpp_result_gen;
END_RCPP
}
// simYast
arma::mat simYast(unsigned int N, unsigned int J, const arma::mat& Y, arma::vec& CLASS, arma::cube& PY_a, arma::mat& ABETA);
RcppExport SEXP _SLCMRT_simYast(SEXP NSEXP, SEXP JSEXP, SEXP YSEXP, SEXP CLASSSEXP, SEXP PY_aSEXP, SEXP ABETASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ABETA(ABETASEXP);
    rcpp_result_gen = Rcpp::wrap(simYast(N, J, Y, CLASS, PY_a, ABETA));
    return rcpp_result_gen;
END_RCPP
}
// computePYa
Rcpp::List computePYa(unsigned int J, unsigned int P, unsigned int nClass, arma::mat& Atable, arma::mat& BETA, arma::mat& Kappa);
RcppExport SEXP _SLCMRT_computePYa(SEXP JSEXP, SEXP PSEXP, SEXP nClassSEXP, SEXP AtableSEXP, SEXP BETASEXP, SEXP KappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type P(PSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Kappa(KappaSEXP);
    rcpp_result_gen = Rcpp::wrap(computePYa(J, P, nClass, Atable, BETA, Kappa));
    return rcpp_result_gen;
END_RCPP
}
// identify_strictcheck
double identify_strictcheck(const arma::mat& DELTA, const arma::mat& completematrix, unsigned int M, unsigned int K, unsigned int ncomp, const arma::uvec& mainorder);
RcppExport SEXP _SLCMRT_identify_strictcheck(SEXP DELTASEXP, SEXP completematrixSEXP, SEXP MSEXP, SEXP KSEXP, SEXP ncompSEXP, SEXP mainorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type completematrix(completematrixSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncomp(ncompSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mainorder(mainorderSEXP);
    rcpp_result_gen = Rcpp::wrap(identify_strictcheck(DELTA, completematrix, M, K, ncomp, mainorder));
    return rcpp_result_gen;
END_RCPP
}
// Update_HDB_strict
void Update_HDB_strict(unsigned int J, unsigned int N, unsigned int M, unsigned int K, unsigned int nClass, double sigma2_h, const arma::mat& alpha, const arma::vec& sigma2_t, const arma::mat& logt, const arma::mat& Atable, arma::mat& H, const arma::vec& Eta, const arma::vec& Tau, arma::mat& DELTA, const double omega, const arma::mat& Yast, arma::mat& BETA, arma::mat& ABETA, double sigma_B, const arma::uvec& main, const arma::uvec& mainorder, const arma::mat& completematrix, const arma::mat& LBtable, const arma::uvec& Bindices);
RcppExport SEXP _SLCMRT_Update_HDB_strict(SEXP JSEXP, SEXP NSEXP, SEXP MSEXP, SEXP KSEXP, SEXP nClassSEXP, SEXP sigma2_hSEXP, SEXP alphaSEXP, SEXP sigma2_tSEXP, SEXP logtSEXP, SEXP AtableSEXP, SEXP HSEXP, SEXP EtaSEXP, SEXP TauSEXP, SEXP DELTASEXP, SEXP omegaSEXP, SEXP YastSEXP, SEXP BETASEXP, SEXP ABETASEXP, SEXP sigma_BSEXP, SEXP mainSEXP, SEXP mainorderSEXP, SEXP completematrixSEXP, SEXP LBtableSEXP, SEXP BindicesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_h(sigma2_hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2_t(sigma2_tSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Yast(YastSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ABETA(ABETASEXP);
    Rcpp::traits::input_parameter< double >::type sigma_B(sigma_BSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type main(mainSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type mainorder(mainorderSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type completematrix(completematrixSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LBtable(LBtableSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Bindices(BindicesSEXP);
    Update_HDB_strict(J, N, M, K, nClass, sigma2_h, alpha, sigma2_t, logt, Atable, H, Eta, Tau, DELTA, omega, Yast, BETA, ABETA, sigma_B, main, mainorder, completematrix, LBtable, Bindices);
    return R_NilValue;
END_RCPP
}
// Update_alpha
Rcpp::List Update_alpha(unsigned int N, unsigned int J, unsigned int nClass, const arma::mat& Y, arma::vec& CLASS, const arma::cube& PY_a, arma::vec& pis, const arma::mat& Atable, arma::mat& alpha, const arma::vec& d0);
RcppExport SEXP _SLCMRT_Update_alpha(SEXP NSEXP, SEXP JSEXP, SEXP nClassSEXP, SEXP YSEXP, SEXP CLASSSEXP, SEXP PY_aSEXP, SEXP pisSEXP, SEXP AtableSEXP, SEXP alphaSEXP, SEXP d0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d0(d0SEXP);
    rcpp_result_gen = Rcpp::wrap(Update_alpha(N, J, nClass, Y, CLASS, PY_a, pis, Atable, alpha, d0));
    return rcpp_result_gen;
END_RCPP
}
// Update_DB
Rcpp::List Update_DB(unsigned int J, unsigned int M, unsigned int K, const arma::mat& alpha, unsigned int nClass, const arma::mat& Yast, arma::mat& BETA, arma::mat& DELTA, arma::mat& ABETA, const arma::mat& Atable, double sigma_B, double omega);
RcppExport SEXP _SLCMRT_Update_DB(SEXP JSEXP, SEXP MSEXP, SEXP KSEXP, SEXP alphaSEXP, SEXP nClassSEXP, SEXP YastSEXP, SEXP BETASEXP, SEXP DELTASEXP, SEXP ABETASEXP, SEXP AtableSEXP, SEXP sigma_BSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Yast(YastSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ABETA(ABETASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_B(sigma_BSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_DB(J, M, K, alpha, nClass, Yast, BETA, DELTA, ABETA, Atable, sigma_B, omega));
    return rcpp_result_gen;
END_RCPP
}
// random_D
arma::mat random_D(unsigned int J, unsigned int K, double w);
RcppExport SEXP _SLCMRT_random_D(SEXP JSEXP, SEXP KSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(random_D(J, K, w));
    return rcpp_result_gen;
END_RCPP
}
// random_D_strict
arma::mat random_D_strict(int J, int K, int M, double t, arma::mat DtoQtable);
RcppExport SEXP _SLCMRT_random_D_strict(SEXP JSEXP, SEXP KSEXP, SEXP MSEXP, SEXP tSEXP, SEXP DtoQtableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type DtoQtable(DtoQtableSEXP);
    rcpp_result_gen = Rcpp::wrap(random_D_strict(J, K, M, t, DtoQtable));
    return rcpp_result_gen;
END_RCPP
}
// permuteAtableIndices
arma::vec permuteAtableIndices(unsigned int nClass, unsigned int K, unsigned int M, const arma::vec& vv, const arma::vec& perm);
RcppExport SEXP _SLCMRT_permuteAtableIndices(SEXP nClassSEXP, SEXP KSEXP, SEXP MSEXP, SEXP vvSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vv(vvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(permuteAtableIndices(nClass, K, M, vv, perm));
    return rcpp_result_gen;
END_RCPP
}
// permute_attr_order
arma::mat permute_attr_order(unsigned int nClass, unsigned int K, unsigned int M, arma::mat perm);
RcppExport SEXP _SLCMRT_permute_attr_order(SEXP nClassSEXP, SEXP KSEXP, SEXP MSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(permute_attr_order(nClass, K, M, perm));
    return rcpp_result_gen;
END_RCPP
}
// gen_Y
arma::mat gen_Y(arma::cube& PY_a, arma::vec& CLASS, unsigned int N, unsigned int J, unsigned int P);
RcppExport SEXP _SLCMRT_gen_Y(SEXP PY_aSEXP, SEXP CLASSSEXP, SEXP NSEXP, SEXP JSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_Y(PY_a, CLASS, N, J, P));
    return rcpp_result_gen;
END_RCPP
}
// simY
arma::mat simY(unsigned int N, unsigned int J, unsigned int M, unsigned int nClass, const arma::vec& CLASS, const arma::mat& Atable, const arma::mat& BETA, const arma::mat& Kappa);
RcppExport SEXP _SLCMRT_simY(SEXP NSEXP, SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP CLASSSEXP, SEXP AtableSEXP, SEXP BETASEXP, SEXP KappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Kappa(KappaSEXP);
    rcpp_result_gen = Rcpp::wrap(simY(N, J, M, nClass, CLASS, Atable, BETA, Kappa));
    return rcpp_result_gen;
END_RCPP
}
// RTLL
double RTLL(const int& N, const int& J, const int& nClass, const arma::mat& logt, const arma::vec& Eta, const arma::vec& Tau, const arma::mat& Atable, const arma::mat& H, const arma::vec& sigma2_t, const arma::vec& pis);
RcppExport SEXP _SLCMRT_RTLL(SEXP NSEXP, SEXP JSEXP, SEXP nClassSEXP, SEXP logtSEXP, SEXP EtaSEXP, SEXP TauSEXP, SEXP AtableSEXP, SEXP HSEXP, SEXP sigma2_tSEXP, SEXP pisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const int& >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logt(logtSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Eta(EtaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2_t(sigma2_tSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pis(pisSEXP);
    rcpp_result_gen = Rcpp::wrap(RTLL(N, J, nClass, logt, Eta, Tau, Atable, H, sigma2_t, pis));
    return rcpp_result_gen;
END_RCPP
}
// m2LL
double m2LL(unsigned int N, unsigned int J, unsigned int P, unsigned int nClass, const arma::mat& Y, const arma::vec& pis, const arma::cube& PY_a);
RcppExport SEXP _SLCMRT_m2LL(SEXP NSEXP, SEXP JSEXP, SEXP PSEXP, SEXP nClassSEXP, SEXP YSEXP, SEXP pisSEXP, SEXP PY_aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type P(PSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type PY_a(PY_aSEXP);
    rcpp_result_gen = Rcpp::wrap(m2LL(N, J, P, nClass, Y, pis, PY_a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SLCMRT_getmode", (DL_FUNC) &_SLCMRT_getmode, 1},
    {"_SLCMRT_myrnorm", (DL_FUNC) &_SLCMRT_myrnorm, 2},
    {"_SLCMRT_rbernoulli_one", (DL_FUNC) &_SLCMRT_rbernoulli_one, 1},
    {"_SLCMRT_rTruncNorm_1b", (DL_FUNC) &_SLCMRT_rTruncNorm_1b, 4},
    {"_SLCMRT_rDirichlet", (DL_FUNC) &_SLCMRT_rDirichlet, 1},
    {"_SLCMRT_rTruncNorm", (DL_FUNC) &_SLCMRT_rTruncNorm, 4},
    {"_SLCMRT_rmultinomial", (DL_FUNC) &_SLCMRT_rmultinomial, 1},
    {"_SLCMRT_gen_bijectionvector", (DL_FUNC) &_SLCMRT_gen_bijectionvector, 2},
    {"_SLCMRT_inv_gen_bijectionvector", (DL_FUNC) &_SLCMRT_inv_gen_bijectionvector, 3},
    {"_SLCMRT_CL_gen_invbijection_table", (DL_FUNC) &_SLCMRT_CL_gen_invbijection_table, 3},
    {"_SLCMRT_GenerateAtable", (DL_FUNC) &_SLCMRT_GenerateAtable, 3},
    {"_SLCMRT_Update_Tau", (DL_FUNC) &_SLCMRT_Update_Tau, 7},
    {"_SLCMRT_Update_sigma2_t", (DL_FUNC) &_SLCMRT_Update_sigma2_t, 9},
    {"_SLCMRT_Update_Eta", (DL_FUNC) &_SLCMRT_Update_Eta, 8},
    {"_SLCMRT_djpsamp", (DL_FUNC) &_SLCMRT_djpsamp, 8},
    {"_SLCMRT_computeLB", (DL_FUNC) &_SLCMRT_computeLB, 7},
    {"_SLCMRT_Update_alpha_SLCMRT", (DL_FUNC) &_SLCMRT_Update_alpha_SLCMRT, 15},
    {"_SLCMRT_simYast", (DL_FUNC) &_SLCMRT_simYast, 6},
    {"_SLCMRT_computePYa", (DL_FUNC) &_SLCMRT_computePYa, 6},
    {"_SLCMRT_identify_strictcheck", (DL_FUNC) &_SLCMRT_identify_strictcheck, 6},
    {"_SLCMRT_Update_HDB_strict", (DL_FUNC) &_SLCMRT_Update_HDB_strict, 24},
    {"_SLCMRT_Update_alpha", (DL_FUNC) &_SLCMRT_Update_alpha, 10},
    {"_SLCMRT_Update_DB", (DL_FUNC) &_SLCMRT_Update_DB, 12},
    {"_SLCMRT_random_D", (DL_FUNC) &_SLCMRT_random_D, 3},
    {"_SLCMRT_random_D_strict", (DL_FUNC) &_SLCMRT_random_D_strict, 5},
    {"_SLCMRT_permuteAtableIndices", (DL_FUNC) &_SLCMRT_permuteAtableIndices, 5},
    {"_SLCMRT_permute_attr_order", (DL_FUNC) &_SLCMRT_permute_attr_order, 4},
    {"_SLCMRT_gen_Y", (DL_FUNC) &_SLCMRT_gen_Y, 5},
    {"_SLCMRT_simY", (DL_FUNC) &_SLCMRT_simY, 8},
    {"_SLCMRT_RTLL", (DL_FUNC) &_SLCMRT_RTLL, 10},
    {"_SLCMRT_m2LL", (DL_FUNC) &_SLCMRT_m2LL, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SLCMRT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}