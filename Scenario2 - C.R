################################################################################
#
#
#
#                                     SCENARIO 2C
#
#
#
#
#  This script is designed to generate the results associated with Scenario 1A.
#  
#  In this scenario, we have defined the parameters by setting delta to 1,
#  and the dropout rate to 10 with MAR dropout and a baseline confounder.
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
source(Simul_age.R)

#### PARAMETERS :  ####
K <- 3 # number of marker
I <- 500 # number of subject

DeltaT <- 1 #discretization step
DeltaTestim <- 1 
Delta <- 1

seed <- rep

para_mu0 <- c(-.33, -0.36, 0.4, 0.3, 0.35, 0.2, -0.9, 0.01,0.9) # initial level (int, Bx, intb, Bxb)
para_mu <- c(0.28, 0.06,0.02, -0.36, -0.02,0.3, 1.16, 0.34,0.5) # change over time
para_RE <- rep(0,21)

L <- matrix(c(-1,-0.42,-0.33,0,0,0,
              0.01,-0.31,0,0,0,0,
              0.02,0,0,0,0,0,
              0,0,0,-0.3,0,0,
              0,0,0,0,0.3,0,
              0,0,0,0,0,-0.76), nrow = 6) # matrix

lower_triangular <- L[lower.tri(L, diag = TRUE)]
para_RE <- lower_triangular

vec_alpha_ij <-c(0, 0, 0,-0.4,0,0,0.28,-1.71, 0) # transition matrix

para_B <- c(0,0,0) # error AR
para_Sig <- c(1, -0.4, 1.2) # SD 
paraEtha0 <- c(0, 0,0)
paraEtha1 <- c(1, 1,1)        

#Construction of matD, matB and  Sig
q <- c(1,1,1) # number of RE on DeltaX 
matD <- DparChol(q = (sum(q)+K), qvector= para_RE)
Sig <- (para_Sig)^2*diag(K)
predictorsA <- NULL 

## trend models
fixed_X0.models <- c("1 + X+C", "1 + X+C","1 + X+C")
fixed_DeltaX.models <- c("1 + X+C", "1 + X+C","1 + X+C")

randoms_X0.models <- c("1", "1","1")
randoms_DeltaX.models <- c("1", "1","1")
mod_trans.model = "1"

###### CALL OF FUNCTION Une_Simul_age_unif ########
subject="id"
Time="Time"
all.preds=c("X", "C")

data_r <- Une_Simul_age_unif(seed=(seed*3333), txMissOutC=0, txMissVisit=0, K=K, I=I, DeltaT=DeltaT, DeltaTestim=DeltaTestim,
                             fixed_X0.models = fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                             randoms_X0.models = randoms_X0.models, randoms_DeltaX.models = randoms_DeltaX.models, 
                             mod_trans.model = mod_trans.model, predictorsA = predictorsA, para_mu0=para_mu0, 
                             para_mu=para_mu, matD, q=q, vec_alpha_ij=vec_alpha_ij, Sig=Sig, paraEtha0, paraEtha1)

data_rY <- data_r %>%
  filter(!is.na(Y3)|!is.na(Y2)|!is.na(Y1))
table(data_rY$id)


colnames(data_r)<-c(colnames(data_r)[1:7],"L","M","Y")

quantile(data_r$Y, probs = c(0.6, 0.7,0.8,0.85, 0.9, 0.95,1), na.rm = T)


data_r1 <- data_r %>% filter(!is.na(Y))


data_r1 <- data_r1 %>%
  group_by(id) %>%
  filter(row_number() <= first(which(Y > 350), default = n())) %>%
  ungroup()



epsa <- 0.0001
epsb <- 0.0001
epsd <- 0.0001


coef <-c(-0.407963144, -0.175653645 , 0.207046863,  0.291790972,  0.402947144 , 0.075850815 ,-0.888157762,  0.045105665,
         0.435478430 , 0.373228090 ,-0.008527491,  0.014264427, -0.297428751, -0.155828467,  0.273820038,  1.084452335,
         0.417328477 , 0.599799181 ,-0.818419348, -0.504694611, -1.304953666,  0.061890452, -0.296546538,  0.554132690,
         0.025267308  ,0.199115339 , 0.763642656 ,-0.304796869 , 0.252215493 ,-0.276056601, -0.398405299,  0.269221196,
         -1.712525874  ,1.009956995, -0.392884210  ,1.204037468)
indexparaFixeUser <- c(22:23,27,29,31:32,35:36,37:39,40:42,44:45,48,52:57)
nb <- seq(1,57)
nb_est <- nb[-indexparaFixeUser]

para_ini_2MARQ <- rep(NA,57)

para_ini_2MARQ[nb_est] <- coef
para_ini_2MARQ[indexparaFixeUser] <- c(rep(0,17),0,1,0,1,0,1)
paraFixeUser <- para_ini_2MARQ[indexparaFixeUser]



Delta <- 1
DeltaTestim <- 1

data_r1 <- as.data.frame(data_r1)


res <-
  CInLPN(structural.model = list(fixed.LP0 = ~ 1 + X +C| 1 + X+C| 1 + X+C,
                                 fixed.DeltaLP =L| M | Y  ~ 1 +  X+C | 1 + X+C| 1 + X+C,
                                 random.DeltaLP = ~ 1|1|1,
                                 trans.matrix = ~ 1 ,
                                 delta.time = DeltaTestim),
         measurement.model = list(link.functions = list(links = c(NULL,NULL,NULL),
                                                        knots = list(NULL, NULL,NULL))),
         
         parameters = list(paras.ini = para_ini_2MARQ, Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser),
         option = list(parallel = T, nproc = 10, print.info = T, epsa= epsa, epsb = epsb, epsd = epsd,
                       makepred = F),
         Time = "TimeSeq",
         subject = "id",
         data = data_r1,
         TimeDiscretization = F, cholesky = T)

CI_app <- res


save(CI_app,file=paste("Review_MAR2", rep,".Rdata",sep=""))
