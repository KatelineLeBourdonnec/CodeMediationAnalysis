################################################################################
#
#
#
#                                     SCENARIO 1D
#
#
#
#
#  This script is designed to generate the results associated with Scenario 1D.
#  
#  In this scenario, we have defined the parameters by setting delta to 0.1,
#  the dropout rate to 10 and the number of subject to 250
#
#
# Each scenario is replicated 250 times to ensure robustness and reliability 
# in the results.
#
#
################################################################################

Sys.setenv(R_LIBS_USER = "~/R/x86_64-pc-linux-gnu-library/4.2")
options(rgl.useNULL=TRUE)
args <- commandArgs(TRUE)
job <- as.numeric(args[1])
rep <- as.numeric(args[2])

#### LOAD PACKAGES :####
require(MASS)
require(splines)
require(mvtnorm)
require(lcmm)
require(dplyr)
require(parallel)
require(CInLPN)

#### LOAD FUNCTION

source(create_data.R)
source(Simul.R)

#### PARAMETERS :  ####
K <- 3 # number of marker
I <- 250 # number of subject

seed <- rep
TempsFin <- 5
DeltaT <- 0.1
DeltaTestim <- 0.1 

para_mu0 <- c(0.5, 1.80, 0.2, 0.90, 0.6, 1.5) 
para_mu <- c(0.10, 0.20, 0.2, 0.80, 0.8, 0.40) 
para_RE <- rep(0,21)

A <- matrix(c(1,0.2,0.1,0,0,0,
              0.2,2,0.1,0,0,0,
              0.1,0.1,3,0,0,0,
              0,0,0,1,0,0,
              0,0,0,0,2,0,
              0,0,0,0,0,3), nrow = 6) 
L <- t(chol(A))

lower_triangular <- L[lower.tri(L, diag = TRUE)]
para_RE <- lower_triangular
vec_alpha_ij <-c(0, 0, 0,0.3,0,0,0.4,0.5, 0) # Transition matrix

para_B <- c(0,0,0)
para_Sig <- c(0.3, 0.6, 0.2) # SD 
paraEtha0 <- c(0, 0,0)
paraEtha1 <- c(1, 1,1)       

###### 
#Construction matD, matB et  Sig
q <- c(1,1,1) # nb of RE
matD <- DparChol(q = (sum(q)+K), qvector= para_RE)
Sig <- (para_Sig)^2*diag(K)
predictorsA <- NULL

## trend models
fixed_X0.models <- c("1 + X", "1 + X","1 + X")
fixed_DeltaX.models <- c("1 + X", "1 + X","1 + X")

randoms_X0.models <- c("1", "1","1")
randoms_DeltaX.models <- c("1", "1","1")
mod_trans.model = "1"

###### call of function une_simul ########
subject="id"
Time="Time"
all.preds=c("X") 

data_r <- Une_Simul(seed=(seed*3333), txMissOutC=0.2, txMissVisit=0.1, K=K, I=I, TempsFin=TempsFin, DeltaT=DeltaT, DeltaTestim=DeltaTestim,
                    fixed_X0.models = fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                    randoms_X0.models = randoms_X0.models, randoms_DeltaX.models = randoms_DeltaX.models, 
                    mod_trans.model = mod_trans.model, predictorsA = predictorsA, para_mu0=para_mu0, 
                    para_mu=para_mu, matD, q=q, vec_alpha_ij=vec_alpha_ij, Sig=Sig, paraEtha0, paraEtha1)

colnames(data_r)<-c(colnames(data_r)[1:6],"L","M","Y")

epsa <- 0.0001
epsb <- 0.0001
epsd <- 0.0001


MULT_L <-  hlme(L ~ (X)*Time, random =~1+Time, subject='id', data = data_r)
MULT_M <-  hlme(M ~ (X)*Time, random =~1+Time, subject='id', data = data_r)
MULT_Y <-  hlme(Y ~  (X)*Time, random =~1+Time, subject='id', data = data_r)
coefL <- MULT_L$best
coefM <- MULT_M$best
coefY <- MULT_Y$best
# 

indexparaFixeUser <- c(16:18,21:23,25:27,29,30,32,34:36,38:39,42,46:51)
coefMU0 <- c(coefL[1:2],coefM[1:2],coefY[1:2])
coefMU <- c(coefL[3:4],coefM[3:4],coefY[3:4])
para_ini_fixe <- c(coefMU0,coefMU)
para_ini_2MARQ <- c(para_ini_fixe,lower_triangular,rep(0,9),coefL[8],coefM[8],coefY[8],0,1,0,1,0,1)
paraFixeUser <- para_ini_2MARQ[indexparaFixeUser]


res <-
  CInLPN(structural.model = list(fixed.LP0 = ~ 1 + X| 1 + X| 1 + X,
                                 fixed.DeltaLP =L| M | Y  ~ 1 +  X | 1 + X| 1 + X,
                                 random.DeltaLP = ~ 1|1|1,
                                 trans.matrix = ~ 1 ,
                                 delta.time = DeltaTestim),
         measurement.model = list(link.functions = list(links = c(NULL,NULL,NULL),
                                                        knots = list(NULL, NULL,NULL))),
         
         parameters = list(paras.ini = para_ini_2MARQ, Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser),
         option = list(parallel = T, nproc = 1, print.info = F, epsa= epsa, epsb = epsb, epsd = epsd,
                       makepred = F),
         Time = "Time",
         subject = "id",
         data = data_r,
         TimeDiscretization = F, cholesky = T)


CI_app <- res

save(CI_app,file=paste("Review_delaisD_", rep,".Rdata",sep=""))