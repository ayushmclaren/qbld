# qbld : Quantile Regression for Binary Longitudinal Data

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a3b8e6d79f2c4d76bc2ddaffd21917d7)](https://app.codacy.com/manual/ayushmclaren/qbld?utm_source=github.com&utm_medium=referral&utm_content=ayushmclaren/qbld&utm_campaign=Badge_Grade_Dashboard)

The R package qbld implements the Bayesian quantile regression model for binary longitudinal data (QBLD) developed in [Rahman and Vossmeyer (2019)](https://arxiv.org/abs/1909.05560). The model handles both fixed and random effects and implements both a blocked and an unblocked Gibbs sampler for posterior inference.

Author: [Ayush Agarwal](https://www.linkedin.com/in/ayushmclaren/)\[aut, cre\], [Dootika Vats](http://home.iitk.ac.in/~dootika/)\[ctb\]

## Installation
To download this development repo,  through the `devtools` package:

```{r}
# install.packages("devtools")
library(devtools)
devtools::install_github("ayushmclaren/qbld")
```
We recommend updating to the latest R version, with appropriate compilation tools.

## Help using qbld
* Manual:- [qbld-manual](https://github.com/ayushmclaren/ExplainIt/blob/master/qbld-manual.pdf)  
* Vignette:- [Using qbld](https://github.com/ayushmclaren/ExplainIt/blob/master/Using%20qbld.pdf)  

## Citation
Please run `citation("qbld")` after loading the package for citation details.

This package was supported by [Google Summer of Code, 2020](https://summerofcode.withgoogle.com/projects/#6628115486343168)

Special thanks to the mentors for all the support throughout the summer:

* Prof. Dootika Vats(dootika@iitk.ac.in),     
* Dr. Adam Maidman(abmaidman@gmail.com)  
