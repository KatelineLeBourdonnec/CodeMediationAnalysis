Une_Simul_age_unif <- function(seed, txMissOutC=0, txMissVisit=0, K, I, DeltaT, DeltaTestim, fixed_X0.models,
                               fixed_DeltaX.models, randoms_X0.models, randoms_DeltaX.models, mod_trans.model, 
                               predictorsA, para_mu0, para_mu, matD, q, vec_alpha_ij, Sig,paraEtha0, paraEtha1){
  set.seed(seed=seed)
  
  
  f_na <- function(X,tau_c){
    #tau_c : tau de couverture
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
  AgeMin <- 65
  AgeMax <- 75
  followup <- 20
  TimeSeq <- seq(AgeMin, AgeMax+followup,by=DeltaT) - AgeMin
  Visit <- seq(0,(length(TimeSeq)-1))
  
  data <- NULL
  matEtha1 <- paraEtha1*diag(K)
  #   p.i <- 0
  for(i in 1:I){
    
    # générer en fonction du temps
    
    AgeDebut <- runif(1,AgeMin,AgeMax) 
    AgeFin <- AgeDebut + runif(1,DeltaT,followup) 
    Time_i <- seq(from = AgeDebut, to = AgeFin, by = 2) 
    Time_i <- Time_i + rnorm(length(Time_i), mean = 0, sd = 0.5)
    Time_i <- Time_i - AgeMin 
    Visit_i <- seq(0,(length(Time_i)-1))
    Vis_seq_i <-sapply(1:length(Time_i), function(x){which.min(abs(TimeSeq - Time_i[x]))}) 
    # Visit <- Time
    Y_ij_seq <-NULL
    # Simulation des covariables
    X <- NULL 
    
    X0 <- 0
    X1 <- 1
    X <- rep(as.numeric(rbinom(n = 1,size = 1,prob = 0.6)))
    
    C <- NULL
    C <- rep(as.numeric(rbinom(1,1,0.4)))
    data_covariate_i <- as.data.frame(cbind(id=rep(i, length(TimeSeq)),X,X0,X1,C, TimeSeq, Visit))
    cols_data_covariate_i <- colnames(data_covariate_i)
    out <- create_x0_x_z0_z_modA_mat(data=data_covariate_i, fixed_X0.models=fixed_X0.models, fixed_DeltaX.models=fixed_DeltaX.models,
                                     randoms_X0.models=randoms_X0.models, randoms_DeltaX.models=randoms_DeltaX.models, 
                                     mod_trans.model=mod_trans.model, subject="id", Time="TimeSeq")
    #xi pour l'individu i
    x0i <- as.matrix(out$x0)
    x0i <- x0i[,-1]# covariable sur valeur ini inntercept ? 
    xi <- as.matrix(out$x)
    xi <- xi[,-1]
    # covariable sur la pente ? 
    #zi pour l'individu i
    z0i <- as.matrix(out$z0) #  sur EA intercept ? 
    zi <- as.matrix(out$z) 
    # modA_mat_i pour l'individu i
    modA_mat_i <- as.matrix(out$modA_mat)
    # générations des effets aléatoires
    Re <- mvrnorm(n = 1, mu = rep(0,(sum(q)+K)), Sigma = matD)
    
    
    wi <- Re[1:K]
    ui <- Re[(K+1):(sum(q)+K)]
    
    # calcul de X0
    X_ij <- x0i%*%para_mu0 + z0i%*%wi
    Y_ij_seq <- rbind(Y_ij_seq, as.numeric(paraEtha0 + matEtha1%*%(X_ij + mvrnorm(n = 1, mu = rep(0,K), Sigma = Sig))))
    ## boucle pour la recurrence
    X_ij_1 <- X_ij
    Aj_1 <- ConstrA(K=K, t=0, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
    for(j in 1:(length(TimeSeq)-1)){ 
      # il faut sauter xi[1:K,] qui correspond a xi0
      X_ij <- DeltaT*xi[(j*K+1):((j+1)*K),]%*%para_mu + DeltaT*zi[(j*K+1):((j+1)*K),]%*%ui + Aj_1%*% X_ij_1
      Y_ij_seq <- rbind(Y_ij_seq, as.numeric(paraEtha0 + matEtha1%*%(X_ij + mvrnorm(n = 1, mu = rep(0,K), Sigma = Sig))))
      #maj X_ij_1 et A_j_1
      X_ij_1 <- X_ij
      Aj_1 <- ConstrA(K=K, t=j, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
    }
    #Y_ij_seq <- as.data.frame(Y_ij_seq)
    Y_ij <- as.data.frame(Y_ij_seq)
    Y_ij[-which(Visit %in% Vis_seq_i),] <- NA
    colnames(Y_ij)<- paste("Y",1:K, sep="")
    ### data complet
    data_i_cplt <- as.data.frame(cbind(data_covariate_i,Y_ij))
    #data_i_cplt <- data_i_cplt[which(Time%%DeltaTestim==0),]
    data <- rbind(data,data_i_cplt)
  }
  #   data[,Time] <- rep(seq(from = 0,(length(data$Time[which(data$id==1)])-1),by = 1),I)
  #Visite et marqueurs manquantes
  ### visite manquante : taux de visites manquantes: txMissVisit
  missVisit <- sample(0:dim(data)[1], round(txMissVisit*dim(data)[1]), replace=FALSE)
  if(txMissVisit!=0){
    data <- data[-c(missVisit),]
  }
  ### observation manquantes
  for( i in 1: dim(data)[1]){
    data[i,-c(1:6)] <- f_na(data[i,-c(1:6)],(1-txMissOutC))
  }
  return(data)
}
