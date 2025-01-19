
########### Example Code for Aggregating Results ###############################
#
#
# This repository contains an R script that demonstrates how to aggregate 
# simulation results from multiple files
#
## How to Use
#
# 1. Download or clone the repository.
# 2. Place your input files (e.g., `path_REF_*.txt`) in the appropriate directory.
# 3. Adjust the file paths in the script to match your directory structure.
# 4. Run the script in an R environment.
#
## Requirements
#
#- R (version 4.0 or higher)
##- Optional: Customize the script for specific simulations or data processing.
#########################################################################################

# here enter your path files 
setwd("C:/Users/PATH")

# here enter your names files
pattern <- "Path_1A" 

files <- list.files(pattern=pattern)
newpath <- NULL
for(f1 in files){
  path <- read.table(f1)
  newpath <- rbind(newpath,path)
}

newpath$ID <- sort(rep(1:500,length(files)))

meanbyreplique <- aggregate(x= newpath[,-16],     
                            
                            by = list(newpath$ID),      
                            
                            FUN = mean)
SDbyreplique <- aggregate(x= newpath[,-16],     
                          
                          by = list(newpath$ID),      
                          
                          FUN = sd)


# Import true value

pathTrue <- read.table("paths_true.txt")

trueVal <- pathTrue
timeT <- seq(0,5,0.1)

meanbyreplique <- meanbyreplique[,-1]
SDbyreplique <- SDbyreplique[,-1]
moy_XY <- meanbyreplique[,seq(from=1,to=max(timeT)*3,by=3)]
sd_XY <-  SDbyreplique[,seq(from=1,to=max(timeT)*3,by=3)]

moy_XMY <- meanbyreplique[,seq(from=2,to=max(timeT)*3,by=3)]
sd_XMY <-  SDbyreplique[,seq(from=2,to=max(timeT)*3,by=3)]

moy_XLMY <- meanbyreplique[,seq(from=3,to=max(timeT)*3,by=3)]
sd_XLMY <-  SDbyreplique[,seq(from=3,to=max(timeT)*3,by=3)]

M_P_XY <- apply(moy_XY,2,mean)
M_P_XMY <- apply(moy_XMY,2,mean)   
M_P_XLMY <- apply(moy_XLMY,2,mean)

Total <- M_P_XY + M_P_XMY + M_P_XLMY
Joint <- M_P_XMY + M_P_XLMY

trueVal <- rbind(trueVal,timeT)
time <- c(1:5)
T_XY <- trueVal[1,which(trueVal[4,]%in%time)]
T_XMY <- trueVal[2,which(trueVal[4,]%in%time)]  
T_XLMY <- trueVal[3,which(trueVal[4,]%in%time)]  
T_XY+T_XMY+T_XLMY

biaisXY <- (M_P_XY-T_XY)/T_XY *100
biaisXMY <- (M_P_XMY-T_XMY)/T_XMY *100
biaisXLMY <- (M_P_XLMY-T_XLMY)/T_XLMY *100

moyXY <- cbind("Time"=c(rep(1,500),rep(2,500),rep(3,500),rep(4,500),rep(5,500)),"Effect"=c(moy_XY$XY1,moy_XY$XY2,moy_XY$XY3,moy_XY$XY4,moy_XY$XY5))
moyXMY <- cbind("Time"=c(rep(1,500),rep(2,500),rep(3,500),rep(4,500),rep(5,500)),"Effect"=c(moy_XMY$XMY1,moy_XMY$XMY2,moy_XMY$XMY3,moy_XMY$XMY4,moy_XMY$XMY5))
moyXLMY <- cbind("Time"=c(rep(1,500),rep(2,500),rep(3,500),rep(4,500),rep(5,500)),"Effect"=c(moy_XLMY$XLMY1,moy_XLMY$XLMY2,moy_XLMY$XLMY3,moy_XLMY$XLMY4,moy_XLMY$XLMY5))

moyXY <- as.data.frame(moyXY)
moyXMY <- as.data.frame(moyXMY)
moyXLMY <- as.data.frame(moyXLMY)

moyXY$Time <- as.character(moyXY$Time)
moyXMY$Time <- as.character(moyXMY$Time)
moyXLMY$Time <- as.character(moyXLMY$Time)

TXY <- as.numeric(c(rep(T_XY[1],500),rep(T_XY[2],500),rep(T_XY[3],500),rep(T_XY[4],500),rep(T_XY[5],500)))
TXMY <- as.numeric(c(rep(T_XMY[1],500),rep(T_XMY[2],500),rep(T_XMY[3],500),rep(T_XMY[4],500),rep(T_XMY[5],500)))
TXLMY <- as.numeric(c(rep(T_XLMY[1],500),rep(T_XLMY[2],500),rep(T_XLMY[3],500),rep(T_XLMY[4],500),rep(T_XLMY[5],500)))

