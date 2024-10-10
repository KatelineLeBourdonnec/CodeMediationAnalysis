###### OUVRIR UNE SIMUL.R

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

#### LOAD FUNCTION : ####
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

data_r <- Une_Simul_age_unif(seed=(seed*3333), txMissOutC=0.1, txMissVisit=0.2, K=K, I=I, DeltaT=DeltaT, DeltaTestim=DeltaTestim,
                             fixed_X0.models = fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                             randoms_X0.models = randoms_X0.models, randoms_DeltaX.models = randoms_DeltaX.models, 
                             mod_trans.model = mod_trans.model, predictorsA = predictorsA, para_mu0=para_mu0, 
                             para_mu=para_mu, matD, q=q, vec_alpha_ij=vec_alpha_ij, Sig=Sig, paraEtha0, paraEtha1)

colnames(data_r)<-c(colnames(data_r)[1:7],"L","M","Y")

# Random selection of 3 visits by subject
data_t <- data_r %>%
  group_by(id) %>%
  filter(!is.na(L) & !is.na(M) & !is.na(Y)) %>%
  filter(n() >= 3) %>%
  sample_n(3)  

data_t2 <- data_r %>%
  group_by(id) %>%
  filter(!is.na(L) & !is.na(M) & !is.na(Y)) %>%
  filter(n()<=2)

data_t <- rbind(data_t,data_t2)

res_fusion <- left_join(data_r, data_t, by = c("id", "TimeSeq"), suffix = c("_old", "_new"), relationship="many-to-many")

data_r1 <- data_r 

data_r1$M <- res_fusion$M_new
data_r1$L <- res_fusion$L_new


epsa <- 0.0001
epsb <- 0.0001
epsd <- 0.0001

MULT_L <-  hlme(L ~ (X+C)*TimeSeq, random =~1+TimeSeq, subject='id', data = data_r)
MULT_M <-  hlme(M ~ (X+C)*TimeSeq, random =~1+TimeSeq, subject='id', data = data_r)
MULT_Y <-  hlme(Y ~  (X+C)*TimeSeq, random =~1+TimeSeq, subject='id', data = data_r)
coefL <- MULT_L$best
coefM <- MULT_M$best
coefY <- MULT_Y$best


indexparaFixeUser <- c(16:18,21:23,25:27,29,30,32,34:36,38:39,42,46:51)+6
coefMU0 <- c(coefL[1:3],coefM[1:3],coefY[1:3])
coefMU <- c(coefL[4:6],coefM[4:6],coefY[4:6])
para_ini_fixe <- c(coefMU0,coefMU)
para_ini_2MARQ <- c(para_ini_fixe,lower_triangular,rep(0,9),coefL[10],coefM[10],coefY[10],0,1,0,1,0,1)
paraFixeUser <- para_ini_2MARQ[indexparaFixeUser]

data_r1 <- data_r1[-which(is.na(data_r1$L)&is.na(data_r1$M)&is.na(data_r1$Y)), ]


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
save(CI_app,file=paste("Review_AgeA_INCOMP_COV", rep,".Rdata",sep=""))
