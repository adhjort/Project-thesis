### ABOUT
# Author:   Anders Hjort
# Date:     24.10.18
# About:    Code to simulate a bayesian normal multivariate distribution of beta coefficients 
#           in a linear model y = beta0 + x1*beta1

### PACKAGES
library(MASS)
library(ggplot2)
library(mvtnorm)
library(ellipse)
library(gridExtra)

#### PRIOR DISTIRBUTION ###
n = 30 #obsevations
p = 1 #covariates (in addition to intercept)
mu = c(0,0) 
cov = 0.8
sigma = matrix(c(1,cov,cov,1), nrow = 2)

#### LIKELIHOOD / OBSERVATIONS ###
# Generate design matrix of n observations with p covariates + intercept column
X = matrix(rnorm(n*p), ncol = p)
X = cbind(rep(1,n),X)

# Create the true beta-coefficients used to simulate data y = beta0 + x1*beta1+x2*beta2
b = c(1,0.5)

# Noise
eps_var = 1
eps = rnorm(n = n, mean = 0, sd = eps_var)

# "Observe"/generate some values (based on the true betas)
y = X%*%b + eps 

### POSTERIOR ###
# Use posterior mean and posterior covariance mtx calculations from project thesis
post_sigma_inv = solve(sigma) + (1/eps_var)*t(X)%*%X
post_sigma = solve(post_sigma_inv)
post_mu = t(post_sigma)%*%t(t(mu)%*%solve(sigma)+(1/eps_var)*t(y)%*%X)

# Make a grid centered in the prior mu
gridsize = 4 #number of standard deviances from the mean 
xygrid = expand.grid(x = seq(from = mu[1]-gridsize*sigma[1,1], to = mu[1]+gridsize*sigma[1,1], length.out = 100),
                     y = seq(from = mu[2]-gridsize*sigma[2,2], to = mu[2]+gridsize*sigma[2,2], length.out = 100))

# The actual calculations of the prior and posterior
prior = dmvnorm(x = xygrid, mean = mu, sigma = sigma)
posterior = dmvnorm(x = xygrid, mean = post_mu, sigma = post_sigma)
df = data.frame(cbind(xygrid, prior,posterior))

### PLOTTING ###
# Contour plots
contPlot = ggplot(df, aes(x=x, y=y)) + 
  geom_contour(aes(z = prior, col = "prior")) +
  geom_contour(aes(z = posterior, col = "posterior"))+
  coord_fixed(xlim = c(mu[1]-3*sigma[1,1],mu[1]+3*sigma[1,1]),
                               ylim = c(mu[2]-3*sigma[2,2], mu[2]+3*sigma[2,2]), ratio = 1) + 
  labs(x ="beta0",y="beta1", title = "Prior vs. posterior") + 
  scale_color_manual(name="", values = c("prior"="steelblue", "posterior" = "red"))

# Regression plot
fit1 = lm(y ~ X[,-1]) #Remove intercept
summary(fit1)

xx = seq(from=min(X[,2]), to = max(X[,2]), length.out = n)
y_pred = coef(fit1)[1]+coef(fit1)[2]*xx

reg_df = data.frame(x = X[,2],y_obs = y)

regPlot = ggplot(data = reg_df) + 
  geom_point(aes(x=x, y=y_obs, col = "Observed values"))+
  labs(x ="x",y="y", title = "Observations and regression") +
  geom_line(aes(x=xx, y = y_pred, col = "Linear prediction"))+
  scale_color_manual(name="", values = c("Observed values"="steelblue", "Linear prediction" = "red"))

contPlot
regPlot

