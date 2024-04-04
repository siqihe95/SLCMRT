#' @title SLCMRT_Gibbs
#'
#' @description Run the Gibbs sampling algorithm for fitting the SLCMRT model.
#'
#' @param Y A matrix of observations with respondents by items.
#' @param logt A matrix of log-transformed response times with respondents by items.
#' @param K An integer specifying the number of attributes.
#' @param M An integer specifying the number of levels of attributes.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param chain_length An integer specifying the total number of iterations.
#' @param rep_each An integer specifying the number of repetitions.
#'
#' @return A list containing estimated parameters including:
#' \itemize{
#'   \item est_Eta: Estimated item response time coefficients.
#'   \item est_Tau: Estimated person speed coefficients.
#'   \item est_CL: Estimated latent classes.
#'   \item est_pis: Estimated class proportions.
#'   \item est_H: Estimated item-attribute response time coefficient matrix.
#'   \item sum_D: Sum of D matrices over iterations.
#'   \item est_sigma2_eta: Estimated hyperparameter of variance.
#'   \item est_sigma2_t: Estimated hyperparameters of variance.
#'   \item est_D: Estimated item-attribute structure matrix.
#'   \item est_Beta: Estimated item response coefficients.
#'   \item est_ABETA: Estimated Alpha times Beta coefficient matrix.
#' }
#'
#' @export
#' @importFrom stats rgamma rnorm
slcmrt_Gibbs = function(Y, logt, K, M, burnin, chain_length, rep_each){

  # extract the model parameters
  J = dim(Y)[1]
  N = dim(Y)[2]
  P = nlevels(factor(Y))
  nClass = M^K
  CL = nClass-1
  Atable = GenerateAtable(nClass,K,M)$Atable
  Bindices = t(GenerateAtable(nClass,K,M)$finalcols)
  DtoQtable = GenerateAtable(nClass,K,M)$DtoQtable
  nClass = length(Bindices)
  main = which(colSums(DtoQtable)==1)
  main_idx = main - 1
  inter = which(colSums(DtoQtable)!=1 & colSums(DtoQtable)!=0)
  attorder = c(1,main,inter)-1
  LBtable <- GenerateAtable(nClass,K,M)$LBtable

  completematrix <- matrix(0,K,K*(M-1))
  xindex <- sort(rep(c(1:K),(M-1)))
  yindex <- c(1:(K*(M-1)))
  for(k in 1:length(xindex)){
    completematrix[xindex[k],yindex[k]] <- 1
  }

  LogL = 9999999999999999

  for (s in 1:rep_each){
    chain_m_burn = chain_length-burnin
    class_dist = matrix(0, N, chain_m_burn)
    pi_dist = matrix(0, nClass, chain_m_burn)
    D_dist = array(0, dim = c(J, nClass, chain_m_burn))

    Eta_dist = matrix(0, J, chain_m_burn)
    Tau_dist = matrix(0, N, chain_m_burn)
    H_dist = array(0, dim = c(J, nClass, chain_m_burn))
    sigma2_eta_dist = rep(0, chain_m_burn)
    sigma2_t_dist = matrix(0, J, chain_m_burn)

    BETAS = array(0, dim = c(J, nClass, chain_m_burn))
    CLs = matrix(0, N, chain_m_burn)
    PIs = matrix(0, nClass, chain_m_burn)

    kappa = seq(0,2*(P-2),by=2)
    Kappa = matrix(rep(kappa,J),nrow=J,byrow=TRUE)
    Kappa = Kappa/4

    #hyperparameter
    a_t = b_t = 1
    a_H = b_H = 1
    a_B = b_B = 1
    w_0 = w_1 = 1
    a_eta = b_eta = 1
    sigma2_h = 1
    sigma_B = 1
    sigma2_eta = 1
    omega = 0.5

    # Gibbs sampling initialization
    DELTA = random_D_strict(J,K,M,0.2,DtoQtable)
    H = matrix(0,J,nClass)
    BETA = matrix(0,J,nClass)
    non_zero_h0 = NA
    non_zero_beta0 = NA
    for(k in 1:sum(DELTA)){
      non_zero_h0[k] = rTruncNorm_1b(0,1,1,0)
      non_zero_beta0[k] = rTruncNorm_1b(0,1,1,0)
    }
    H[DELTA==1] = non_zero_h0
    BETA[DELTA==1] = non_zero_beta0

    inputs = computePYa(J, P, nClass, Atable, BETA, Kappa)
    ABETA = inputs$ABETA
    PY_a = inputs$PY_a

    Tau = rnorm(N)
    Eta = rnorm(J,mean = 0, sd = sigma2_eta)
    sigma2_eta = 1
    sigma2_t = rep(1,J)

    d0 = rep(1,nClass)
    pis = rDirichlet(d0)
    CLASS = c(1:nClass, sample(1:nClass, N-nClass, replace=TRUE, prob=pis))
    alpha =  Atable[CLASS,]
    CLASS = CLASS-1

    for(t in 1:chain_length){

      # update Yast
      Yast = simYast(N, J, Y, CLASS, PY_a, ABETA)

      # Update pis and alpha and CLASS
      NewA = Update_alpha_SLCMRT(N, J, nClass, Y, CLASS, PY_a, pis, Atable, alpha,
                          d0, logt, Eta, Tau, H, sigma2_t)

      # update H
      NewHDB = Update_HDB_strict(J, N, M, K, nClass, sigma2_h, alpha, sigma2_t, logt, Atable,
                                 H, Eta, Tau, DELTA, omega, Yast, BETA, ABETA, sigma_B, main_idx,
                                 attorder, completematrix, LBtable, Bindices)
      # update PY_a
      newPY = computePYa(J, P, nClass, Atable, BETA, Kappa)
      PY_a <- newPY$PY_a

      # update omega
      sumD <- sum(DELTA)-J
      omega = rDirichlet(c(1,1)+c(sumD,J*(nClass-1)-sumD))[1]

      # update Tau
      Tau = Update_Tau(N,J,Eta,sigma2_t,logt,alpha,H)

      # update sigma2_t
      sigma2_t = Update_sigma2_t(N,J,Eta,Tau,logt,alpha,H,a_t,b_t)

      # update Eta
      Eta = Update_Eta(N,J,Tau,logt,alpha,H,sigma2_t,sigma2_eta)

      # update sigma2_eta
      precision_eta = rgamma(1, shape = a_eta + J/2, rate = b_eta + sum(Eta^2)/2)
      sigma2_eta = 1/precision_eta

      # update sigma_B
      precision_B = rgamma(1, shape = a_B + sum(DELTA)/2, rate = b_B + sum(BETA^2)/2)
      sigma2_B = 1/precision_B
      sigma_B = sqrt(sigma2_B)

      if(t > burnin){
        tmburn = t - burnin
        Eta_dist[ ,tmburn] = Eta
        Tau_dist[ ,tmburn] = Tau
        class_dist[ ,tmburn] = CLASS
        pi_dist[ ,tmburn] = pis
        H_dist[,,tmburn] = H
        D_dist[,,tmburn] = DELTA
        BETAS[,,tmburn] =  BETA
        CLs[, tmburn] = CLASS
        PIs[, tmburn] =  pis
        sigma2_eta_dist[tmburn] = sigma2_eta
        sigma2_t_dist[ ,tmburn] = sigma2_t

      }
    }
    est_Eta = rowMeans(Eta_dist)
    est_Tau = rowMeans(Tau_dist)
    est_CL = apply(class_dist,1,getmode) + 1
    est_pis = rowMeans(pi_dist)
    est_H = apply(simplify2array(H_dist), 1:2, mean)
    est_sigma2_eta = mean(sigma2_eta_dist)
    est_sigma2_t = rowMeans(sigma2_t_dist)
    est_Beta = apply(simplify2array(BETAS), 1:2, mean)
    sum_D = apply(simplify2array(D_dist), 1:2, sum)
    est_D = (sum_D/chain_m_burn>0.5)+0
    newLogL <- RTLL(N, J, nClass, logt, Eta, Tau,
                    Atable, H, sigma2_t, pis)+m2LL(N,J,P,nClass,Y,est_pis,PY_a)
  }
  if(newLogL < LogL){
    LogL = newLogL
    resultlist = list(est_Eta = est_Eta,
                      est_Tau = est_Tau,
                      est_CL = est_CL,
                      est_pis  = est_pis,
                      est_H = est_H,
                      sum_D = sum_D,
                      est_sigma2_eta=est_sigma2_eta,
                      est_sigma2_t=est_sigma2_t,
                      est_D=est_D,
                      est_Beta = est_Beta,
                      est_ABETA = ABETA)
  }
  return(resultlist)
}
