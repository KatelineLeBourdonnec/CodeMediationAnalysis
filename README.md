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
        - [How to use scripts](#How-to-use-the-scripts)
          1) To simulate
          2) For application
          3) To combine simulation results
      - [Session Info](#Session-info)

## Information
The code was written in R using the following software versions:
- **R version 4.3.3** (2024-02-29 ucrt) -- "Angel Food Cake"
- Copyright (C) 2024 The R Foundation for Statistical Computing
- Platform: x86_64-w64-mingw32/x64 (64-bit)

## Author
Corresponding author: **Kateline Le Bourdonnec**,  
Email: [kateline.le-bourdonnec@u-bordeaux.fr](mailto:kateline.le-bourdonnec@u-bordeaux.fr)

## Documentation
The "CodeMediationAnalysis.Proj" project contains the R scripts used to replicate the results from the paper "Continuous-time mediation analysis for repeatedly measured mediators and outcomes."

## Packages

To install the required packages, run the following code:

```r
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


1. To simulate:
To run the simulation, you need to use the [simul](Simul.R) function to generate data, the CInLPN function from the CInLPN package to estimate the model, and the [estimate function](Function Estimate.R), or respectively [Simul_age](Simul_age.R) and [Estimate_function_age](Function Estimate Age.R) if you are working with age data rather than delays.

The following simulation files were run on a parallel computing server to reduce computation time:

[Scenario1A](Scenario1 - A.R)
[Scenario1B](Scenario1 - B.R)
[Scenario1C](Scenario1 - C.R)
[Scenario1D](Scenario1 - D.R)
[Scenario1E](Scenario1 - E.R)
[Scenario2A](Scenario2 - A.R)
[Scenario2B](Scenario2 - B.R)
[Scenario2C](Scenario2 - C.R)
[Scenario2D](Scenario2 - D.R)

2. For application:
The data used for our application are confidential and require investigator agreements for sharing. Therefore, we have not provided them.

3.To combine simulation results:

Use the [Script]("Concat_result.R") to load the simulation results that were run on the compute server and retrieve mean estimands, bias, ..

## Session :

R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8   
[3] LC_MONETARY=French_France.utf8 LC_NUMERIC=C                  
[5] LC_TIME=French_France.utf8    

time zone: Europe/Paris
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.3.3    tools_4.3.3       rstudioapi_0.16.0
