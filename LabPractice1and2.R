#/*******************************************************************************

#This code relates to a grant application to seek funding for replicating the 
# sleep deprivation study in Belenkey et al. 2003, which studied the effect of
# sleep deprivation on average reaction time over time.
#This code performs a simulation of data from Belenkey et al. 2003, 
#sampling from a linear mixed effect model to 1) calualte power for a replication study, 
#and 2) investigate the anticipated distribution of the effect estimate beta1. 

# Lines 31 to 42 set up the work space (load packages, load data, get original study results)

#Problem 1: Power calculation
# Lines x to x address Problem 1, where I am calculating power under the following settings:
# Collect data for 8 nights (lines 60/64)
# From 40 subject (line 50)
# With a B1 = 3, since we anticipate a weaker association between Day and Reaction (line 77)
# With a larger error due to heterogeneity of the new sample. eij ~ N(0,30^2) (line 64)
# With random intercept and random slope anticipated to be independent.  cor(b0,i,b1,i)=0 (line 68)
# Power for different sample sizes, including 40, is plotted on line 96
# Power for N=40 is specifically extracted on line 100

#Problem 2: Anticipate distribution of B1
# Lines 105 to 151 address Problem 2, where I calculate the anticipated distribution of B1,
# assuming the same settings in lab problem 1, and using the same linear mixed model. 
# On line 108 I set the sample size to 40, and on line 111 I set the sampling to repeat 200 times. 
# B1 estimates are extracted on line 139
# The distribution of B1 is plotted on line 150.

  
#/*******************************************************************************
rm(list=ls())
library(MASS)
library(lmer)
library(lme4)
library(lmerTest)
library(dplyr)

#load data
data(sleepstudy)

#Get estimates from original study
lmer(Reaction ~ Days + (Days | Subject), sleepstudy)



### LAB 1 Probelm

# We vary the sample size from 5 to 40. When we later plot this we will see what the power is at N = 40
p.matrix = c()
for (N in seq(5,40,by=5)) { # On this line I setting my sample size. Here I look at sample sizes 5-40 (by 5 person increments). I could set this just to 40 if needed. 
  p.temp.vector = c()
  # We repeat simulation for 500 times for each N (Step 3)
  for (i in 1:500) {
    # Set seed
    set.seed(i)
    # Step 1: generate data
    # Generate id
    dat = data.frame(id = paste0(615,1:N))
    # Long format: each subject was followed by 8 days
    dat = dat %>% slice(rep(1:n(), each=8)) # Set days to 8
    # Make Day variable
    dat = dat %>% group_by(id) %>% mutate(Days = 1:n()) %>% ungroup()
    # Simulate random error
    dat = dat %>% mutate(err = rnorm(8*N,mean=0,sd=30)) #Set days to 8 and error to 30
    
    # Simulate (correlated) subject-level random effects for intercepts and slopes
    ## Covariance matrix
    S1 = diag(c(24.741, 5.922)) %*% matrix(c(1, 0, 0, 1), nrow=2) %*% diag(c(24.741, 5.922)) #Make intercept and slope independent
    ## Simulate realization of random intercept and slope
    U1 = mvrnorm(N, mu=c(0,0), Sigma=S1) %>% as.data.frame()
    ## Add identifier (subject id)
    U1 = U1 %>% bind_cols(id = paste0(615,1:N))
    ## Merge subject-level random effects back to data
    dat = dat %>% left_join(U1,by="id")
    
    # Simulate the outcome: Reaction_ij
    dat = dat %>% mutate(Reaction = (251.405 + V1) + (3 +V2)*Days + err) #Reduce effect of B1 to 3. 
    
    # Step 2: test the null hypothesis
    mod = lmer(Reaction ~ Days + (Days | id), dat)
    p.value = summary(mod)$coef["Days","Pr(>|t|)"]
    # Save p value
    p.temp.vector = c(p.temp.vector,p.value)
  }
  # Save p value vector for each N
  p.matrix = cbind(p.matrix,p.temp.vector)
}
# Matrix => data.frame
p.matrix = p.matrix %>% as.data.frame()
# Add column names
names(p.matrix) = seq(5,40,by=5)
# Step 4: calculate power
power = p.matrix %>% summarise_all(function(x) mean(x<0.05))

#plot the distribution of power for sample sized 5 to 40
plot(seq(5,40,by=5),power,xlab = "Sample size (N)",ylab="Power",type = "b", pch = 19)
abline(h=0.8,lty=2)

#Specifically pull out the power for N = 40
power_N40 <- power$`40`
#[1] 0.736



# Lab Problem 2

b.matrix = c()
for (N in 40) { # Set sample size to 40
  b.temp.vector = c()
  # We repeat simulation for 200 times (Step 3)
  for (i in 1:200) {
    # Set seed
    set.seed(i)
    # Step 1: generate data
    # Generate id
    dat = data.frame(id = paste0(615,1:N))
    # Long format: each subject was followed by 8 days
    dat = dat %>% slice(rep(1:n(), each=8))
    # Make Day variable
    dat = dat %>% group_by(id) %>% mutate(Days = 1:n()) %>% ungroup()
    # Simulate random error
    dat = dat %>% mutate(err = rnorm(8*N,mean=0,sd=30))
    
    # Simulate (correlated) subject-level random effects for intercepts and slopes
    ## Covariance matrix
    S1 = diag(c(24.741, 5.922)) %*% matrix(c(1, 0, 0, 1), nrow=2) %*% diag(c(24.741, 5.922))
    ## Simulate realization of random intercept and slope
    U1 = mvrnorm(N, mu=c(0,0), Sigma=S1) %>% as.data.frame()
    ## Add identifier (subject id)
    U1 = U1 %>% bind_cols(id = paste0(615,1:N))
    ## Merge subject-level random effects back to data
    dat = dat %>% left_join(U1,by="id")
    
    # Simulate the outcome: Reaction_ij
    dat = dat %>% mutate(Reaction = (251.405 + V1) + (3 +V2)*Days + err)
    
    # Step 2: test the null hypothesis
    mod = lmer(Reaction ~ Days + (Days | id), dat)
    b.value = summary(mod)$coef["Days","Estimate"] #Extract B1 for Days
    # Save p value
    b.temp.vector = c(b.temp.vector,b.value)
  }
  # Save p value vector for each N
  b.matrix = cbind(b.matrix,b.temp.vector)
}
# Matrix => data.frame
b.matrix = b.matrix %>% as.data.frame()

#Plot distribution of B1 estimates
hist(b.matrix$b.temp.vector)
abline(v=3,lwd=3)
