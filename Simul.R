#' Une_Simul
#' 
#' This function simulates longitudinal data when the time is in delay based on mixed-effects models and covariates. 
#' It generates time-series data for multiple individuals, including random effects, fixed 
#' and random covariates, and introduces missing observations and visits based on specified rates.
#' 
#' @param seed Integer, seed for reproducibility.
#' @param txMissOutC Numeric, rate of missing observations in the output data.
#' @param txMissVisit Numeric, rate of missing visits.
#' @param K Integer, number of dimensions for longitudinal responses (e.g., number of markers).
#' @param I Integer, number of simulated individuals.
#' @param TempsFin Numeric, total duration of the simulation.
#' @param DeltaT Numeric, time interval between consecutive points.
#' @param DeltaTestim Numeric, time interval for selecting final observations.
#' @param fixed_X0.models formula, models for fixed covariates at initial time.
#' @param fixed_DeltaX.models formula, models for fixed covariates on time differences.
#' @param randoms_X0.models formula, models for random effects at initial time.
#' @param randoms_DeltaX.models  formula, models for random effects on time differences.
#' @param mod_trans.model Model for the transition matrix (modA_mat).
#' @param predictorsA Parameters for defining the predictors used in the transition matrix.
#' @param para_mu0 Vector, parameters for fixed covariates at the initial time.
#' @param para_mu Vector, parameters for fixed covariates on time.
#' @param matD Matrix, covariance matrix for random effects.
#' @param q Vector, dimensions of random effects for each covariate.
#' @param vec_alpha_ij Vector, coefficients for transition models (Aj).
#' @param Sig Matrix, covariance matrix for measurement errors.
#' @param paraEtha0 Vector, parameters for fixed terms in responses.
#' @param paraEtha1 Matrix, coefficients for covariate-dependent terms in responses.
#' 
#' @return A data.frame containing the simulated data, including covariates, visits, 
#'         and longitudinal responses for each individual.
#'
#' @details 
#' - Covariates are generated at the initial time and evolve over time according to the provided models.
#' - Random effects are sampled from a multivariate normal distribution.
#' - Longitudinal responses are generated for each individual at each visit, based on covariates, 
#'   random effects, and measurement errors.
#' - Missing data are introduced based on the specified rates for visits (`txMissVisit`) 
#'   and observations (`txMissOutC`).
#' 
#' @examples
#' set.seed(123)
#' data <- Une_Simul(seed = 123, txMissOutC = 0.1, txMissVisit = 0.2, K = 3, I = 100, 
#'                   TempsFin = 10, DeltaT = 1, DeltaTestim = 1, 
#'                   fixed_X0.models = NULL, fixed_DeltaX.models = NULL, 
#'                   randoms_X0.models = NULL, randoms_DeltaX.models = NULL, 
#'                   mod_trans.model = NULL, predictorsA = NULL, para_mu0 = rep(0, 3), 
#'                   para_mu = rep(0, 3), matD = diag(6), q = c(3, 3), 
#'                   vec_alpha_ij = rep(0, 3), Sig = diag(3), 
#'                   paraEtha0 = rep(0, 3), paraEtha1 = diag(3))
#'
#' @export


Une_Simul <- function(seed, txMissOutC=0, txMissVisit=0, K, I, TempsFin, DeltaT, DeltaTestim, fixed_X0.models,
                      fixed_DeltaX.models, randoms_X0.models, randoms_DeltaX.models, mod_trans.model, 
                      predictorsA, para_mu0, para_mu, matD, q, vec_alpha_ij, Sig,paraEtha0, paraEtha1){
  
  
  f_na <- function(X,tau_c){
    K <- length(X)
    pas <- (1-tau_c)/K
    r <- runif(1)
    X.na <- X
    if(r > tau_c){
      j <-0
      while(j<= K){
        if(r >= (tau_c +j*pas) && r < (tau_c +(j+1)*pas)){
          X.na[j+1] <- NA
        }
        j <- j+1
      }
    }
    
    return(X.na)
  }
  set.seed(seed=seed)
  Time <- seq(from = 0, to = TempsFin, by = DeltaT)
  Visit <- seq(0,(length(Time)-1))
  data <- NULL
  matEtha1 <- paraEtha1*diag(K)
  #   p.i <- 0
  for(i in 1:I){
    Y_ij <-NULL
    X <- NULL 
    X0 <- 0
    X1 <- 1
    X <- rep(as.numeric(rbinom(n = 1,size = 1,prob = 0.6)))
    data_covariate_i <- as.data.frame(cbind(id=rep(i, length(Time)),X,X0,X1, Time, Visit))
    cols_data_covariate_i <- colnames(data_covariate_i)
    out <- create_x0_x_z0_z_modA_mat(data=data_covariate_i, fixed_X0.models=fixed_X0.models, fixed_DeltaX.models=fixed_DeltaX.models,
                                     randoms_X0.models=randoms_X0.models, randoms_DeltaX.models=randoms_DeltaX.models, 
                                     mod_trans.model=mod_trans.model, subject="id", Time="Time")
   
    x0i <- as.matrix(out$x0) 
    xi <- as.matrix(out$x) 
 
    z0i <- as.matrix(out$z0) 
    zi <- as.matrix(out$z) 
    
    modA_mat_i <- as.matrix(out$modA_mat)
 
    Re <- mvrnorm(n = 1, mu = rep(0,(sum(q)+K)), Sigma = matD)
    
    
    wi <- Re[1:K]
    ui <- Re[(K+1):(sum(q)+K)]
    
    
    X_ij <- x0i%*%para_mu0 + z0i%*%wi
    Y_ij <- rbind(Y_ij, as.numeric(paraEtha0 + matEtha1%*%(X_ij + mvrnorm(n = 1, mu = rep(0,K), Sigma = Sig))))
    
    X_ij_1 <- X_ij
    Aj_1 <- ConstrA(K=K, t=0, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
    for(j in 1:(length(Time)-1)){ 
      
      X_ij <- DeltaT*xi[(j*K+1):((j+1)*K),]%*%para_mu + DeltaT*zi[(j*K+1):((j+1)*K),]%*%ui + Aj_1%*% X_ij_1
      Y_ij <- rbind(Y_ij, as.numeric(paraEtha0 + matEtha1%*%(X_ij + mvrnorm(n = 1, mu = rep(0,K), Sigma = Sig))))
     
      X_ij_1 <- X_ij
      Aj_1 <- ConstrA(K=K, t=j, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
    }
    Y_ij <- as.data.frame(Y_ij)
    colnames(Y_ij)<- paste("Y",1:K, sep="")
    
    data_i_cplt <- as.data.frame(cbind(data_covariate_i,Y_ij))
    data_i_cplt <- data_i_cplt[which(Time%%DeltaTestim==0),]
    data <- rbind(data,data_i_cplt)
  }
 
  missVisit <- sample(0:dim(data)[1], round(txMissVisit*dim(data)[1]), replace=FALSE)
  if(txMissVisit!=0){
    data <- data[-c(missVisit),]
  }
  
  for( i in 1: dim(data)[1]){
    data[i,-c(1:5)] <- f_na(data[i,-c(1:5)],(1-txMissOutC))
  }
  return(data)
}

