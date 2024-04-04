# Simulation Study Description
#
# This file provides a simulation study based on pre-saved model parameters.
# The simulation is conducted for an item size of 12, with the number of attributes set to 3,
# and each attribute having 2 levels.

library("RcppAlgos")

# 1. Loading Data and initial preprocessing
Beta0 <- readRDS("inst/extdata/Beta0_K3M2.rds")
H0 <- readRDS("inst/extdata/H0_K3M2.rds")
Delta0 = ifelse(H0==0,0,1)
Beta0 = as.matrix(Beta0)
H0 = as.matrix(H0)
J = 12
K = 3
M = 2
N = 5000
P = 5
CL = M^K-1
nClass = M^K
Atable = GenerateAtable(nClass,K,M)$Atable
Bindices = t(GenerateAtable(nClass,K,M)$finalcols)
DtoQtable = GenerateAtable(nClass,K,M)$DtoQtable
nClass = length(Bindices)
main = which(colSums(DtoQtable)==1)
main_idx = main - 1
inter = which(colSums(DtoQtable)!=1 & colSums(DtoQtable)!=0)
attorder = c(1,main,inter)-1
LBtable <- GenerateAtable(nClass,K,M)$LBtable


# 2. Setting model parameters
set.seed(1)
# person speed par@ tau, item speed par@ Eta
Tau0 = rnorm(N)
sigma2_eta0 = 1
Eta0 = rnorm(J,mean = 0, sd = sqrt(sigma2_eta0))

# threshold par@ kappa
kappa = seq(0,2*(P-2),by=2)
Kappa = matrix(rep(kappa,J),nrow=J,byrow=TRUE)
Kappa = Kappa/4

# class parameter@ Class, prob vec@ pi, profile mat@ alpha
d0 = rep(1,nClass)
pi0 = rep(1/nClass,nClass)
CLASS0 = c(1:nClass, sample( 1:nClass, N-nClass, replace=TRUE, prob=pi0))
Alpha0 = Atable[CLASS0,]

completematrix <- matrix(0,K,K*(M-1))
xindex <- sort(rep(c(1:K),(M-1)))
yindex <- c(1:(K*(M-1)))
for(k in 1:length(xindex)){
  completematrix[xindex[k],yindex[k]] <- 1
}


# 3. Simulations of categorical response and response time
sigma2_t = rep(1,J)
logt = matrix(NA,N,J)
for (i in 1:N){
  for (j in 1:J){
    mu = Eta0[j]*Tau0[i] + sum(Alpha0[i,]*H0[j,])
    logt[i,j] =  rnorm(1,mu,sqrt(sigma2_t))
  }
}

# generate item response @Y
PY_a0 = computePYa(J, P, nClass, Atable, Beta0, Kappa)$PY_a
ABETA0 =  computePYa(J, P, nClass, Atable, Beta0, Kappa)$ABETA
Y = gen_Y(PY_a0,CLASS0-1,N,J,P)


# 4.Model estimation via Gibbs sampling
out = slcmrt_Gibbs(Y, logt, K=3, M=2, burnin=1000, chain_length=2000, rep_each = 1)

# 5. Attribute-order rotation and model evaluation
permat = permuteGeneral(c(0:(K-1)))
permat = permute_attr_order(nClass, K, M, permat)
DeltaPro = 0
for(r in 1:ncol(permat)){
  temp = mean(Delta0==out$est_D[,permat[,r]])
  if(temp > DeltaPro){
    DeltaPro = temp
    trueidx = permat[,r]
  }
}

est_Delta = out$est_D[,trueidx]
EAR_Delta <- mean(((est_Delta-Delta0)==0))
TPR_Delta <- mean(((est_Delta-Delta0)==0)[Delta0==1])
FPR_Delta <- 1 - mean(((est_Delta-Delta0)==0)[Delta0==0])

# Absolute bias error
est_Beta = out$est_Beta[,trueidx]
est_H = out$est_H[,trueidx]

abe_Beta = mean(abs(est_Beta-Beta0))
abe_Beta_Delta1 = mean(abs(est_Beta-Beta0)[Delta0==1])
abe_Beta_Delta0 = mean(abs(est_Beta-Beta0)[Delta0==0])

abe_H = mean(abs(est_H-H0))
abe_H_Delta1 = mean(abs(est_H-H0)[Delta0==1])
abe_H_Delta0 = mean(abs(est_H-H0)[Delta0==0])

abe_pis = mean(abs(pi0-out$est_pis))

# Root mean squared error

rmse_Beta = mean((est_Beta-Beta0)^2)
rmse_Beta_Delta1 = mean(((est_Beta-Beta0)^2)[Delta0==1])
rmse_Beta_Delta0 = mean(((est_Beta-Beta0)^2)[Delta0==0])

rmse_H = mean((est_H-H0)^2)
rmse_H_Delta1 = mean(((est_H-H0)^2)[Delta0==1])
rmse_H_Delta0 = mean(((est_H-H0)^2)[Delta0==0])

rmse_pis = mean(abs(pi0-out$est_pis)^2)

sim_error = c(EAR_Delta,TPR_Delta,FPR_Delta,
             abe_Beta,abe_Beta_Delta1,abe_Beta_Delta0,
             abe_H,abe_H_Delta1,abe_H_Delta0,
             rmse_Beta,rmse_Beta_Delta1,rmse_Beta_Delta0,
             rmse_H,rmse_H_Delta1,rmse_H_Delta0,
             abe_pis, rmse_pis)

names(sim_error) = c("EAR_Delta","TPR_Delta","FPR_Delta",
                    "abe_Beta","abe_Beta_Delta1","abe_Beta_Delta0",
                    "abe_H","abe_H_Delta1","abe_H_Delta0",
                    "rmse_Beta","rmse_Beta_Delta1","rmse_Beta_Delta0",
                    "rmse_H","rmse_H_Delta1","rmse_H_Delta0",
                    "abe_pis", "rmse_pis")

round(sim_error,2)
