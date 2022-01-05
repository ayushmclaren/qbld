<p align="center">
    <a href="https://cran.r-project.org/web/packages/qbld/index.html">
    <img src="https://img.shields.io/cran/v/qbld?style=flat-square"
         alt="CRAN version">
    <a href="https://github.com/ayushmclaren/qbld/actions">
    <img src="https://github.com/ayushmclaren/qbld/workflows/R-CMD-check/badge.svg"    
         alt="RCMD check()">    
    <a href="https://travis-ci.com/github/ayushmclaren/qbld">
    <img src="https://travis-ci.com/ayushmclaren/qbld.svg?branch=master"
         alt="Build status (Travis CI)">
    <a href="https://www.codacy.com/manual/ayushmclaren/qbld?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ayushmclaren/qbld&amp;utm_campaign=Badge_Grade">
    <img src="https://app.codacy.com/project/badge/Grade/3ee1436280ad4736bf4f6f909bf881fd"    
         alt="Codacy Badge"> 
    <a href="https://summerofcode.withgoogle.com/projects/#6628115486343168">
    <img src="https://img.shields.io/badge/Google-Funded-success?style=flat&logo=Google"
         alt="GSoC project">
      <a href="https://www.r-project.org/">
    <img src="https://img.shields.io/badge/100%25--blue?style=flat&logo=R"
         alt="The R project for statistical computing"> 
    <a href="https://github.com/ayushmclaren/qbld/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/ayushmclaren/qbld"
         alt="GitHub License"> </a>   
     <a href="https://cran.r-project.org/web/packages/qbld/index.html">
    <img src="https://cranlogs.r-pkg.org/badges/grand-total/qbld"
         alt="Downloads"> </a>   
</p> 




# qbld : Quantile Regression for Binary Longitudinal Data
The R package `qbld` implements the Bayesian quantile regression model for binary longitudinal data (QBLD) developed in [Rahman and Vossmeyer (2019)](https://www.emerald.com/insight/content/doi/10.1108/S0731-90532019000040B009/full/html). The model handles both fixed and random effects and implements both a blocked and an unblocked Gibbs sampler for posterior inference.

Author: [Ayush Agarwal](https://www.linkedin.com/in/ayushmclaren/)\[aut, cre\], [Dootika Vats](http://home.iitk.ac.in/~dootika/)\[ctb\]

## Installation
This R package is on CRAN, and its preferred URL is <https://CRAN.R-project.org/package=qbld>.
```{r}
#Install from CRAN
install.packages("qbld")

#load the package
library(qbld)
```

To download this development repo,  through the `devtools` package:

```{r}
# install.packages("devtools")
library(devtools)
devtools::install_github("ayushmclaren/qbld")
```
We recommend updating to the latest R version, with appropriate compilation tools.

## Help using qbld
*  Manual:- [qbld-manual](https://cran.r-project.org/web/packages/qbld/qbld.pdf)  
*  Vignette:- [Using qbld](https://cran.r-project.org/web/packages/qbld/vignettes/qbld.pdf)  

## Citation
Please run `citation("qbld")` after loading the package for citation details.

This package was supported by [Google Summer of Code, 2020](https://summerofcode.withgoogle.com/projects/#6628115486343168)

Special thanks to the mentors for all the support throughout the summer:

*  Prof. Dootika Vats <dootika@iitk.ac.in>,     
*  Dr. Adam Maidman <abmaidman@gmail.com>.    
