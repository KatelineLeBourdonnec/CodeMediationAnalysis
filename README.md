---
  title: "README"
output:
  html_document: default
pdf_document: default
date: "2024-11-02"
---
  
  
  
  ## Tables
  - [Information](#Information)
    - [Author](#author)
      - [Documentation](#Documentation)
        - [Packages](#packages)
          - [How to use scripts](#How to use the scripts)
            1) For application
          2) To simulate 
        3) To combine simulation results
      - [Session Info](#Session info)
        
        ## Information
        
        The code was written in R with the following software versions : 
          R version 4.3.3 (2024-02-29 ucrt) -- "Angel Food Cake"
        Copyright (C) 2024 The R Foundation for Statistical Computing
        Platform: x86_64-w64-mingw32/x64 (64-bit)
        
        ## Author 
        
        Corresponding author : Kateline Le Bourdonnec,
        kateline.le-bourdonnec@u-bordeaux.fr
        
        ## Documentation
        
        The "CodeMediationAnalysis.Proj" project contains the R codes used to write the paper "Continuous-time mediation analysis for repeatedly measured mediators and outcomes " to find the same 
        results. 
        
        ## packages 
        
        load packages :
          ```{r, message=F, warning=F, results = 'hide'}
        options(repos = "https://cran.r-project.org")
       
        install.packages("MASS")
        install.packages("mvtnorm")
        install.packages("ggplot2")
        install.packages("dplyr")
        install.packages("parallel")
        install.packages("CInLPN")
        install.packages("lcmm")
        install.packages("devtools")
        install_github("https://github.com/KatelineLeBourdonnec/CInLPN2")
        
      
        ```
        
        ## How to use scripts  
        

        2. For application : 
        As the data used for our application is confidential and requires investigator agreements to be shared, we have not provided them.
       
          
          2. To simulate : 
          
          To run simulation you need to geneData function (in [geneData](geneData.R) file) and bootstrap_param function (in [boostrap_param](boostrap_param.R)
                                                                                                                         
         Simulation (file [Function_for_Cont_Log-Sub](CODE_BJ_Function_Simu_Cont_Log-Sub.R),
                     [Function_for_Log-Res](CODE_BJ_Function_Simu_Res.R), [Function_for_Linear-Sub](CODE_BJ_Function_Simu_Log-Linear.R), [Function_for_Log-Sub](CODE_BJ_Function_Simu_Log-Sub.R))
         were run on a parallel computing server to reduce calculation time. 

These files contain 4 scripts for function : log/res, log/sub, linear/sub, log/sub (when X continuous).

We have provided the simulation results in the following files: 

- [cont_log_sub](cont_log_sub)

- [linear_sub](linear_sub)

- [log_res](log_res) 

- [log_sub](log_sub)

In these file, each .Rdata return a list with respectively in this order :

* CoefTrue,CoefTrueInter,sdTrue, sdTrueInter <- coefficients and standard error of the true regression model (except for log_res)

* CoefNaifInter0,CoefNaif0,sdNaif0,sdNaifInter0 <- coefficients and standard error of the naive regression model

* Coef2stageM,Coef2stageInter,sd2stageM,sd2stageInter <- coefficients and standard error of the IV regression model

* r2_1,F_1 <- R² statistic and Fisher statistic

* var_corrig,se_corr <- corrected variance and corrected standard error

* TC,TCint <- Coverage rate for coef2stageM and Coef2stageInt (slope)

* VarIntra,VarInter,VarTot,VarIntraInt,VarInterInt,VarTotInt <- variance intra inter subject and total variance

* resu[[1]],resu[[2]] <- bootstrap coefficients


To obtain identical results, simply identify the seed (or job) and replica numbers indicated at the start of each function script.


<font color="red"> All these statistical indicators were used for the analyses conducted, but only the naive and two-stages coefficients over-time are employed for the reproducibility of Figure 2 in script [Script : simulation data](BJ_script_combine_simu.R). </font>

3. To combine simulation : [Script : simulation data](BJ_script_combine_simu.R)

Code for loading simulations that have run on a compute server
and retrieving coefficients over time for naive and two-stage models. 
Different scenarios are represented :
            * Binary / Continuous
            
            * N = 2000, 6000, 20K
            
            * alpha = 2, 3, 4
            
            * 4 methods : Naive, Log/Res, Line/Subs, Log/Subs

## Session info

