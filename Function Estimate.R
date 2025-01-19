#' Estimate estimand with delay : 
#'
#' This function performs bootstrap simulations to estimate estimands under different structural models 
#' (`XY`, `XMY`, `XLMY`).
#'
#' @param Boots Integer. The number of bootstrap iterations for the main estimation process.
#' @param B Integer. The number of bootstrap iterations for marker-specific perturbations.
#' @param tmax Numeric. The maximum time for the simulation.
#' @param delta Numeric. The time step interval.
#' @param type Character. The type of structural model to use. Possible values are:
#'   - `"XY"`: Direct relationship between predictors and outcomes.
#'   - `"XMY"`: Indirect effec through the mediator.
#'   - `"XLMY"`: Indirect effec through the mediator and the time varying confounder..
#' @param seed Integer. Seed for random number generation (default = 1).
#' 
#' @return A numeric vector containing the aggregated expected values for 
#'         each bootstrap iteration.
#'
#' @details 
#' The function initializes model parameters, performs bootstrapping to estimate 
#' the variance-covariance matrix and marker-specific means, and adjusts these 
#' estimates through perturbations for each bootstrap sample.
#'
#'
#' @import MASS
#' @import mvtnorm
#' @import lcmm
#' @import dplyr
#' @import parallel
#' @import CInLPN2
#' @import CInLPN
#' @export


################################################################################
# This R script is designed to estimate causal effects in complex longitudinal
# models using a bootstrap approach. It leverages structural and measurement 
# models to evaluate relationships between variables over time. 
# Key components include:
#  
#       Input Parameters: The script accepts inputs such as the number of 
#       bootstrap iterations, simulation steps, time intervals (delta), and 
#       the type of causal relationship (e.g., direct effect, mediation, or 
#       more complex pathways).
#
#       Data Simulation: Synthetic datasets are created with temporal structures
#       and predefined relationships between variables (e.g., X, Y, M, and L).
#
#       Model Estimation: The CInLPN2_estimand function is used to fit the data
#       to structural models. Parameter uncertainty is accounted for by sampling
#       from multivariate normal distributions based on estimated means and 
#       covariance matrices.
#
#       Bootstrap Procedure: At each iteration, conditional expectations of 
#       outcomes are calculated by simulating latent and observed variables.
#       Differences between these expectations are used to compute causal effects.
#
#       Output: The script produces estimated causal effects for various 
#       scenarios, saved as text files for further analysis.
#
################################################################################

estim_X_Y_anti<- function(Boots,B, tmax, delta, type, seed=1){
  if(type=="XY"){ # changer les écritures du type
    a=1
    b=0
    c=0
    d=0
    e=0
    f=0
  }
  if(type=="XMY"){
    a=1
    b=1
    c=0
    d=1
    e=0
    f=0
  }
  if(type=="XLMY"){
    a=1
    b=1
    c=1
    d=1
    e=1
    f=0
  } 
  
  t=seq(0,tmax,by=delta)
  data1 <- data.frame(t=seq(0,tmax,by=delta), X_Y = a, X_M = b, X_L=c, "NUM"=rep(1, length(t)),"L"=0,"M"=0,"Y"=0)
  E_X_Y <- NULL
  
  for (i in 1: Boots){
    set.seed(seed*i)
    thetaB<- rmvnorm(n=1, mean=CI_app$best, sigma = ma)
    THETA <- rep(NA,length(CI_app$best)+length(which(CI_app$posfix==1)))
    THETA[-which(CI_app$posfix==1)] <- thetaB
    THETA[which(CI_app$posfix==1)] <- CI_app$coefficients[which(CI_app$posfix==1)]
    indexparaFixeUser <- which(CI_app$posfix==1)
    paraFixeUser <- THETA[indexparaFixeUser] 
    
    CI_estim <- CInLPN2_estimand(structural.model = list(fixed.LP0 = ~1+X_L|1+X_M|1+X_Y, fixed.DeltaLP = L|M|Y~1+X_L|1+X_M|1+X_Y,
                                                         random.DeltaLP = ~1|1|1, trans.matrix=~1, delta.time= delta),
                                 measurement.model = list(link.functions = list(links = c(NULL,NULL,NULL),knots = list(NULL,NULL,NULL))),
                                 parameters = list(paras.ini = THETA, Fixed.para.index = indexparaFixeUser,
                                                   Fixed.para.values = paraFixeUser),
                                 option = list(parallel=F,nproc=6, print.info=F, makepred=F),
                                 Time="t",subject="NUM",data= data1, TimeDiscretization = F, cholesky = T)
    
    VC <-CI_estim$VC
    MU <- CI_estim$Mu
    
    
    mu_L <- MU[sequence(length(t),by=3),1]
    mu_M <- MU[sequence(length(t),by=3)+1,1] 
    mu_Y_t <- MU[length(t)*3]
    mu_LM <- c(mu_L,mu_M) 
    
    VC_L <- VC[sequence(length(t),by=3),sequence(length(t),by=3)]
    VC_M <- VC[sequence(length(t),by=3)+1,sequence(length(t),by=3)+1]
    VC_Y_t <- VC[dim(VC)[2],dim(VC)[2]]
    
    covLM <- VC[sequence(length(t),by=3)+1,sequence(length(t),by=3)]
    covLY <-VC[sequence(length(t),by=3)+2,sequence(length(t),by=3)]
    covMY <-VC[sequence(length(t),by=3)+2,sequence(length(t),by=3)+1]
    covMY_t <- covMY[dim(covMY)[1],]
    covLY_t <- covLY[dim(covLY)[1],]
    
    VC_LM <- matrix(NA, nrow=dim(covLM)[2]*2, ncol=dim(covLM)[2]*2)
    VC_LM[1:dim(covLM)[2],1:dim(covLM)[2]] <- VC_L 
    VC_LM[(dim(covLM)[2]+1):(2*dim(covLM)[2]),(dim(covLM)[2]+1):(2*dim(covLM)[2])] <- VC_M 
    VC_LM[1:dim(covLM)[2],(dim(covLM)[2]+1):(2*dim(covLM)[2])] <- covLM 
    VC_LM[(dim(covLM)[2]+1):(2*dim(covLM)[2]),1:dim(covLM)[2]] <- covLM 
    
    cov_LM_Y_t <- matrix(c(covLY_t,covMY_t),ncol=2)
    
    
    E_f <- NULL 
    
    
    
    for (i in 1:B){
      set.seed(seed*i)
      ml_tb0 <- mvrnorm(1,rep(0,length(mu_LM)),VC_LM)
      ml_tb <- mu_LM + ml_tb0
      ml_tbneg <-  mu_LM - ml_tb0
      mat_temp <- t(as.matrix(c(covLY_t,covMY_t))) %*% solve(VC_LM)
      E <- mu_Y_t +      mat_temp %*% (as.matrix(c(ml_tb-mu_LM)))
      Eneg <- mu_Y_t + mat_temp %*% (as.matrix(c(ml_tbneg-mu_LM)))
      E_f <- append(E_f, c(E,Eneg)) 
    }
    Estimand <- mean(E_f)
    
    
    data2 <- data.frame(t=seq(0,tmax,by=delta), X_M = e , X_Y = d, X_L=f,"NUM"=rep(1, length(t)),"L"=0,"M"=0,"Y"=0)
    
    
    CI_estim <- CInLPN2_estimand(structural.model = list(fixed.LP0 = ~1+X_L|1+X_M|1+X_Y, fixed.DeltaLP = L|M|Y~1+X_L|1+X_M|1+X_Y,
                                                         random.DeltaLP = ~1|1|1, trans.matrix=~1, delta.time= delta),
                                 measurement.model = list(link.functions = list(links = c(NULL,NULL,NULL),knots = list(NULL,NULL,NULL))),
                                 parameters = list(paras.ini = THETA, Fixed.para.index = indexparaFixeUser,
                                                   Fixed.para.values = paraFixeUser),
                                 option = list(parallel=F,nproc=6, print.info=F, makepred=F),
                                 Time="t",subject="NUM",data= data2, TimeDiscretization = F, cholesky = T)
    
    
    
    VC <-CI_estim$VC
    MU <- CI_estim$Mu
    
    mu_L <- MU[sequence(length(t),by=3),1]
    mu_M <- MU[sequence(length(t),by=3)+1,1] 
    mu_Y_t <- MU[length(t)*3]
    mu_LM <- c(mu_L,mu_M) 
    
    VC_L <- VC[sequence(length(t),by=3),sequence(length(t),by=3)]
    VC_M <- VC[sequence(length(t),by=3)+1,sequence(length(t),by=3)+1]
    VC_Y_t <- VC[dim(VC)[2],dim(VC)[2]]
    
    covLM <- VC[sequence(length(t),by=3)+1,sequence(length(t),by=3)]
    covLY <-VC[sequence(length(t),by=3)+2,sequence(length(t),by=3)]
    covMY <-VC[sequence(length(t),by=3)+2,sequence(length(t),by=3)+1]
    covMY_t <- covMY[dim(covMY)[1],]
    covLY_t <- covLY[dim(covLY)[1],]
    
    VC_LM <- matrix(NA, nrow=dim(covLM)[2]*2, ncol=dim(covLM)[2]*2)
    VC_LM[1:dim(covLM)[2],1:dim(covLM)[2]] <- VC_L 
    VC_LM[(dim(covLM)[2]+1):(2*dim(covLM)[2]),(dim(covLM)[2]+1):(2*dim(covLM)[2])] <- VC_M 
    VC_LM[1:dim(covLM)[2],(dim(covLM)[2]+1):(2*dim(covLM)[2])] <- covLM 
    VC_LM[(dim(covLM)[2]+1):(2*dim(covLM)[2]),1:dim(covLM)[2]] <- covLM 
    
    cov_LM_Y_t <- matrix(c(covLY_t,covMY_t),ncol=2)
    
    E_f <- NULL 
    
    
    
    for (i in 1:B){
      set.seed(seed*i)
      ml_tb0 <- mvrnorm(1,rep(0,length(mu_LM)),VC_LM)
      ml_tb <- mu_LM + ml_tb0
      ml_tbneg <-  mu_LM - ml_tb0
      mat_temp <- t(as.matrix(c(covLY_t,covMY_t))) %*% solve(VC_LM)
      E <- mu_Y_t +      mat_temp %*% (as.matrix(c(ml_tb-mu_LM)))
      Eneg <- mu_Y_t + mat_temp %*% (as.matrix(c(ml_tbneg-mu_LM)))
      E_f <- append(E_f, c(E,Eneg)) 
    }
    Estimand1 <- mean(E_f)
    
    E_X_Y0 <-Estimand-Estimand1
    E_X_Y <- append(E_X_Y0, E_X_Y) 
  }
  return(E_X_Y)
} 
