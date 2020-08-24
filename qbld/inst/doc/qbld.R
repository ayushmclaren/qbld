## ----aldmix--------------------------------------------------------------
library(qbld)
set.seed(10)

#generate 1e4 samples
ald.sample <- raldmix(n = 1e4, mu = 0, sigma = 1, p = 0.5)
plot(density(ald.sample), main="AL(0,1,0.5)")

## additional functions
ald.density <- daldmix(c(4,5),mu = 0,sigma = 1,p = 0.5)
ald.cdf <- paldmix(c(1,4),mu = 0,sigma = 1,p = 0.5,lower.tail=TRUE)
ald.quantile <- qaldmix(0.5,mu = 0,sigma = 1,p = 0.5,lower.tail=TRUE)

## ----gig-----------------------------------------------------------------

# random generation
gig.sample <- rgig(n = 1e4, lambda = 0.5, a = 1, b = 2)
plot(density(gig.sample),main="GIG(1,2,0.5)")

# density
gig.density <- dgig(x = 1, a = 1, b = 2, p = 0.5, log_density = FALSE)

## ----data----------------------------------------------------------------
data(airpollution)
str(airpollution)

## ----Block,results='hide'------------------------------------------------

##modelling the output :- Blocked 
#no burn, no verbose, no summary

output.block <- model.qbld(fixed_formula = wheeze~smoking+I(age^2)+age, 
                           data = airpollution, id="id", 
                           random_formula = ~counts+1, p=0.25, 
                           nsim=5000, method="block", burn=0, 
                           summarize=FALSE, verbose=FALSE) 

## ----qbldclass-----------------------------------------------------------
str(output.block)

## ----Unblock-------------------------------------------------------------

##modelling the output :- Unblocked 
#Using burn, no verbose, and summary

output.unblock <- model.qbld(fixed_formula = wheeze~smoking+I(age^2)+age, 
                           data = airpollution, id="id", 
                           random_formula = ~counts+1, p=0.25, 
                           nsim=5000, method="Unblock", burn=0.2, 
                           summarize=TRUE, verbose=FALSE) 

## ----qbldsummaryclass----------------------------------------------------
summary.unblock = summary(output.unblock, 
                          quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                          epsilon=0.10)
str(summary.unblock)

## ----plot----------------------------------------------------------------
plot(output.block, trace = TRUE, density = TRUE, 
     auto.layout = TRUE, ask = NULL)

