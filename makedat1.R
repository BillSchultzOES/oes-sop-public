## The file to make the dat1 simulated dataset
library(randomizr)
library(mvtnorm)
library(dplyr)
set.seed(20180101)

## This seems very clear but would be happy to do this using DeclareDesign too
## Define size of the experiment
N <- 100
## Create a pre-existing characteristic of people that predicts their outcome
covmeans <- rep(0,20)
varcov <- matrix(2,nrow=20,ncol=20)
diag(varcov) <- 5
covmat <- rmvnorm(N,mean=covmeans,sigma=varcov)
summary(as.vector(cor(covmat)))

dat1 <- data.frame(covmat)
names(dat1) <- paste("cov",1:20,sep="")

## Define the baseline outcome as kind of skewed
dat1 <- dat1 %>% mutate(y0noise = round(rchisq(n=nrow(dat1),df=1)),
                        y0 = ifelse(y0noise==0,0, y0noise+ 2*sd(y0noise)*cov2),
                        y0 = y0*(y0>0))

## Test the covariate to baseline outcome relationship
summary(lm(y0~cov2,data=dat1))$r.squared
##var(fitted(lm1))/var(y0)
## Define a simple treatment effect
tau <- 1.5*sd(dat1$y0)
## Make the counterfactual outcome in the treatment condition
## it has less variation than the outcome in the baseline/control condition.
dat1$y1 <- with(dat1,mean(y0) + tau + (y0 - mean(y0))/2 +  sd(y0noise) * cov2)

with(dat1,boxplot(list(y0=y0,y1=y1)))

## Randomly assign treatment by drawing from an urn --- so that each treatment
## arm has the same size across possible randomizations
dat1$Z <- complete_ra(N,m=25)
## We see y1 when Z=1 and y0 when Z=0 (the experimental assignment reveals to
## us one of the two potential outcomes for each person), and each person is
## independent of each other.
dat1$Y <- with(dat1,Z*y1 + (1-Z)*y0)
dat1$id <- 1:N
head(dat1)

write.csv(dat1,file="dat1.csv")
