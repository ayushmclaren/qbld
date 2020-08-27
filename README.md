# qbld : Quantile Regression for Binary Longitudinal Data
An R package follows [Rahman and Vossmeyer (2019)](https://arxiv.org/abs/1909.05560) as its motivating literature, and contributes by extending the various methodologies in quantile framework, to a hierarchical Bayesian quantile regression model for binary longitudinal data (QBLD) and proposing a Markov chain Monte Carlo (MCMC) algorithm to estimate the model. 
The model handles both common (fixed) and individual-specific (random) parameters (commonly referred to as mixed effects in statistics). The algorithm implements a blocking, and an unblocking procedure that is computationally efficient and the distributions involved allow for straightforward calculations of covariate effects.

Author: [Ayush Agarwal](https://www.linkedin.com/in/ayushmclaren/)\[aut, cre\], [Dootika Vats](http://home.iitk.ac.in/~dootika/)\[ctb\]

# Installation
To download this development repo,  through the the `devtools` package:

```{r}
# install.packages("devtools")
library(devtools)
devtools::install_github("ayushmclaren/qbld")
```
I recommend updating to the latest R version, with appropriate compilation tools.

# Help using qbld
* Manual:- [qbld-manual](https://github.com/ayushmclaren/ExplainIt/blob/master/qbld-manual.pdf)
* Vignette:- [Using qbld](https://github.com/ayushmclaren/ExplainIt/blob/master/Using%20qbld.pdf)

# Citation
Please run `citation("qbld")` after loading the package for citation details.
This project was developed through [Google Summer of Code, 2020](https://summerofcode.withgoogle.com/projects/#6628115486343168)
