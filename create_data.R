create_x0_x_z0_z_modA_mat <- function(data, fixed_X0.models, fixed_DeltaX.models, randoms_DeltaX.models, 
                                      randoms_X0.models, mod_trans.model, subject, Time){
  colnames <- colnames(data)
  tau <- unique(sort(data$Visit)) # il y a forcement un $visit car creer dans une_simul (tau = vecteur du nb de visite)
  nb.outcomes <- length(fixed_DeltaX.models)
  all.Y <- 1:nb.outcomes
  I <- 1
  ## construction de la matrice x0=====================================
  x_cov<- NULL
  for(k in 1:nb.outcomes){
    indY_x <- rep(all.Y[k], dim(data)[1])
    data_x_cov_i <- cbind(data, indY_x)
    x_cov <- rbind(x_cov, data_x_cov_i)
  }
  x_cov <- x_cov[order(x_cov[,subject],x_cov[,Time] ),]
  
  ##only for x0 #####
  x0 <- NULL
  nb_x0_k <- NULL
  col_k<-list()
  x0_cov <- x_cov[which(x_cov$Visit==0),]
  #
  ##
  indY_x0 <- x0_cov$indY_x     # pour l intercept, indY_X0 <- 1,1,1,... (x le nombre d'individu)
  for(k in 1:nb.outcomes){
    r<-as.formula(paste(subject, fixed_X0.models[k], sep="~-1+"))
    x0k<-model.matrix(r,data=x0_cov)
    nb_x0_k <- c(nb_x0_k,ncol(x0k))
    if(length(x0k)==0){
      col <- paste(k,"zero",sep="")
      x0k<-matrix(assign(col,rep(0,dim(x0_cov)[1])))
      nb_x0_k <- c(nb_x0_k,ncol(x0k))
    }
    
    colnames<-colnames(x0k)
    colnames<-paste(all.Y[k],colnames,sep="")
    colnames(x0k) <-colnames
    col_k <-c(col_k,list(colnames))
    x0<-cbind(x0,as.matrix(x0k))
  }
  x0 <-cbind(indY_x0,x0)
  ### remplissage avec les zeros
  tous_col_x0 <-unlist(col_k)
  for(i in 1:nrow(x0)){
    col_i <- unlist(col_k[[x0[i,"indY_x0"]]])
    col_0<-tous_col_x0[which(!(tous_col_x0 %in% col_i))]
    x0[i,col_0]<-0 #  passer pour optimisation
  }
  x0 <- as.matrix(x0)
  colnames <- colnames(x0)
  #x0 <- as.matrix(x0[,-c(1)])
  #colnames(x0) <- colnames[-c(1)]
  
  ##only for x #####
  x <- NULL
  nb_x_k <- NULL
  col_k<-list()
  indY_x <- x_cov$indY_x
  for(k in 1:nb.outcomes){
    r<-as.formula(paste(subject, fixed_DeltaX.models[k], sep="~-1+"))
    xk<-model.matrix(r,data=x_cov)
    nb_x_k <- c(nb_x_k,ncol(xk))
    if(length(xk)==0){
      col <- paste(k,"zero",sep="")
      xk<-matrix(assign(col,rep(0,dim(x_cov)[1])))
      nb_x_k <- c(nb_x_k,ncol(xk))
    }
    colnames<-colnames(xk)
    colnames<-paste(all.Y[k],colnames,sep="")
    colnames(xk) <-colnames
    col_k <-c(col_k,list(colnames))
    x<-cbind(x,as.matrix(xk))
  }
  x <-cbind(indY_x,x)
  ### remplissage avec les zeros
  tous_col_x <-unlist(col_k)
  for(i in 1:nrow(x)){
    col_i <- unlist(col_k[[x[i,"indY_x"]]])
    col_0<-tous_col_x[which(!(tous_col_x %in% col_i))]
    x[i,col_0]<-0 # z  passer pour optimisation
  }
  
  x <- as.matrix(x)
  colnames <- colnames(x)
  # x <- as.matrix(x[,-c(1)])
  #colnames(x) <- colnames[-c(1)]
  
  
  #===================================================================
  #     # construction des  matrice z0 et z========================
  data_z_cov <- data[, c(subject,Time,"Visit")]
  z_cov<- NULL
  for(k in 1:nb.outcomes){
    indY_z <- rep(all.Y[k], dim(data_z_cov)[1])
    data_z_cov_i <- cbind(data_z_cov, indY_z)
    z_cov <- rbind(z_cov, data_z_cov_i)
  }
  z_cov <- z_cov[order(z_cov[,subject],z_cov$Visit ),]
  #   z_cov[,Time] <- z_cov[,Time]*DeltaT ######
  
  #### only for z0 ####
  z0_cov <- z_cov[which(z_cov$Visit==0),]
  indY_z0 <- z0_cov$indY_z
  z0 <- NULL
  col_k<-list()
  q0 <- NULL
  nb_paraDw <- 0
  for(k in 1:nb.outcomes){
    r<-as.formula(paste(subject,randoms_X0.models[k], sep="~-1+"))
    z0k<-model.matrix(r,data=z0_cov)
    if(length(z0k)==0){
      col <- paste(k,"zero",sep="")
      z0k<-matrix(assign(col,rep(0,dim(z0_cov)[1])))
    }
    colnames<-colnames(z0k)
    colnames<-paste(all.Y[k],colnames,sep="")
    colnames(z0k) <-colnames
    col_k <-c(col_k,list(colnames))
    z0<-cbind(z0,z0k)
    q0 <- c(q0,ncol(z0k))
    #     nb_paraDw <- nb_paraDw + ncol(z0k)*(ncol(z0k)+1)/2
  }
  z0 <-cbind(indY_z0,z0)
  ### remplissage avec les zeros
  tous_col_z0 <-unlist(col_k)
  for(i in 1:nrow(z0)){
    col_i <- unlist(col_k[[z0[i,"indY_z0"]]])
    col_0<-tous_col_z0[which(!(tous_col_z0 %in% col_i))]
    z0[i,col_0]<-0 # z  passer pour optimisation
  }
  z0 <- z0[,-c(1)]
  z0 <- as.matrix(z0)
  
  
  #### only for z ####
  indY_z <- z_cov$indY_z
  z <- NULL
  col_k<-list()
  q <- NULL
  nb_paraDu <- 0
  for(k in 1:nb.outcomes){
    r<-as.formula(paste(subject,randoms_DeltaX.models[k], sep="~-1+"))
    zk<-model.matrix(r,data=z_cov)
    if(length(zk)==0){
      col <- paste(k,"zero",sep="")
      zk<-matrix(assign(col,rep(0,dim(z_cov)[1])))
    }
    colnames<-colnames(zk)
    colnames<-paste(all.Y[k],colnames,sep="")
    colnames(zk) <-colnames
    col_k <-c(col_k,list(colnames))
    z<-cbind(z,zk)
    q <- c(q,ncol(zk))
    #     nb_paraDu <- nb_paraDu + ncol(zk)*(ncol(zk)+1)/2
  }
  z <-cbind(indY_z,z)
  ### remplissage avec les zeros
  tous_col_z <-unlist(col_k)
  for(i in 1:nrow(z)){
    col_i <- unlist(col_k[[z[i,"indY_z"]]])
    col_0<-tous_col_z[which(!(tous_col_z %in% col_i))]
    z[i,col_0]<-0 # z  passer pour optimisation
  }
  z <- z[,-c(1)]
  z <- as.matrix(z)
  #============================================================
  # calcul de la matrice de design pour le mod?le de transition
  f<-as.formula(paste(subject,mod_trans.model, sep="~-1+"))# subject pour juste avoir un premier membre pour la formule
  modA_mat<-model.matrix(f,data=data)
  dim(modA_mat)
  return(list(x0= x0, x= x, z0= as.matrix(z0), z= as.matrix(z), modA_mat= as.matrix(modA_mat)))
}