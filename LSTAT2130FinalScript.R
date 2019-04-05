setwd("C:/Users/Nicephore/Documents/Q2/BayesianProject")
absences <- 
  read.table("C:/Users/Nicephore/Documents/Q2/BayesianProject/absences.txt", 
             header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)

absences_whole_set <- absences 

#Need to keep the whole data set for use in the last question

absences <- na.omit(absences) #Removal of rows with unavailable values

x <- absences[,3]
y <- absences[,1]
z <- absences[,2]



#2.

#This is the logarithm of the likelihood

loglikelihood <- function(V){
  alpha <- V[1]
  beta_0 <- V[2]
  beta_1 <- V[3]
  tau <- V[4]
  
  q <- exp(beta_0+alpha*z+(beta_1+tau*z)*(x-39)) 
  result <- sum((-q)+y*log(q)-log(factorial(y)))
  return(result)
}

#Trying to calculate the MLEs 
#argmax f(x) = argmin ( -f(x) )

minus_loglikelihood <- function(V){
  result <- - loglikelihood(V)
  return(result)
}

nlm(minus_loglikelihood, c(0,0,0,0), gradtol = 1e-6, fscale = 1, stepmax = 10, iterlim = 50)

 


#Let's work with "uninformative" priors.

prior <- function(V){
  A <- V[1]
  B <- V[2]
  C <- V[3]
  D <- V[4]
  cst <- 10^6
  result <- pnorm(A,0,cst)*pnorm(B,0,cst)*pnorm(C,0,cst)*pnorm(D,0,cst)
  return(result)
}

#3.

logposterior <- function(V){
  
  alpha <- V[1]
  beta_0 <- V[2]
  beta_1 <- V[3] 
  tau <- V[4]
  result <- loglikelihood(V) + log(prior(V))
  return(result)
}

#If one does not want to call another function within function logposterior then:

####logposterior <- function(V){
  
####  alpha <- V[1]
####  beta_0 <- V[2]
####  beta_1 <- V[3] 
####  tau <- V[4]
  
####  q <- exp(beta_0+alpha*z+(beta_1+tau*z)*(x-39)) 
####  loglikelihood <- sum((-q)+y*log(q)-log(factorial(y)))
####  cst <- 10^6
####  prior <- pnorm(alpha,0,cst)*pnorm(beta_0,0,cst)*pnorm(beta_1,0,cst)*pnorm(tau,0,cst)
####  result <- loglikelihood + log(prior)
####  return(result)
####}

#4. (a)

#Implementation inspired by an old assignment from a different module


#M: number of iterations
#burnin: number of initial elements of the samples which will be discarded
#in order to minimise the influence of the initial values fed to the Metropolis function
#std_proposal: 4-vector of standard deviations of the normal 'disruptions' 

Metropolis <- function( M, 
                        burnin, 
                        alphainit,
                        beta_0init,
                        beta_1init,
                        tau_init,
                        std_proposal
)
{
  set.seed(2018)
  result <- matrix(0, nrow = M, ncol = 4)
  accept <- matrix(0, ncol = 4, nrow = M)
  lcur <- logposterior( c(alphainit, beta_0init, beta_1init, tau_init))
  paramProp <- rep(0, 4)
  paramCur  <- c(alphainit, beta_0init, beta_1init, tau_init)
  
  for(cpt1 in 1:M){ # iteration loop
    for(cpt2 in 1:4) { # loop for number of parameters
      paramProp <- paramCur + rnorm(1,0, std_proposal[cpt2])*diag(1,4)[cpt2,] 
      lprop <- logposterior(c(paramProp[1],paramProp[2], paramProp[3], paramProp[4]))  
      
      if(runif(1) < min(1, exp(lprop - lcur))) {
        paramCur <- paramProp 
        lcur <- lprop 
        accept[cpt1, cpt2] <- 1
      }
    }
    result[cpt1, ] <- paramCur
  }
  
  return(list(alpha = result[(burnin+1):M,1], 
              beta_0 = result[(burnin+1):M,2], 
              beta_1 = result[(burnin+1):M,3],
              tau = result[(burnin+1):M,4],
              accept_Rate = list(alpha = mean(accept[(burnin+1):M,1]),
                                 beta_0 = mean(accept[(burnin+1):M,2]),
                                 beta_1 = mean(accept[(burnin+1):M,3])),
                                 tau = mean(accept[(burnin+1):M,4])))
}

M <- Metropolis(M=2.5*10^4, burnin = 1.5*10^4, alphainit = 0.15, beta_0init = 0.15, beta_1init = 0.15, tau_init = 0.15, std_proposal = c(0.22, 0.14, 0.01, 0.016))
M

#All acceptance rates are between 41% and 43%

#Now, we have got to use two different criteria to check the convergence of the chains.

#(b).

library(coda)

MCMC_alpha <- mcmc(data = M$alpha)
MCMC_beta0  <- mcmc(data = M$beta_0)
MCMC_beta1 <- mcmc(data = M$beta_1)
MCMC_tau <- mcmc(data = M$tau)

geweke.diag(MCMC_alpha, frac1=0.1, frac2=0.5)
geweke.diag(MCMC_beta0, frac1=0.1, frac2=0.5)
geweke.diag(MCMC_beta1, frac1=0.1, frac2=0.5)
geweke.diag(MCMC_tau, frac1=0.1, frac2=0.5)


#The z-scores for all parameters  are located between -1.96 and 1.96, so the null hypothesis
#that the chains converged cannot be rejected according to Geweke diagnostic. (See p.15/26 in Chapter 3)

heidel.diag(MCMC_alpha)
heidel.diag(MCMC_beta0)
heidel.diag(MCMC_beta1)
heidel.diag(MCMC_tau)

#Traceplots

par(mfrow = c(2,2))
traceplot(MCMC_alpha, main = "alpha")
traceplot(MCMC_beta0, main = "beta_0")
traceplot(MCMC_beta1, main = "beta_1")
traceplot(MCMC_tau, main = "tau")

#All p-values are considerably larger 5%, so the null hypothesis that the chains converged cannot be rejected according 
#to Heidelberger-Welch diagnostic.

#(c).

#Giving point estimates for each parameter

median(MCMC_alpha) #Estimator of the posterior mean of alpha
median(MCMC_beta0) #Estimator of the posterior mean of beta0
median(MCMC_beta1) #Estimator of the posterior mean of beta1
median(MCMC_tau) #Estimator of the posterior mean of tau

mean(MCMC_alpha) #Estimator of the posterior mean of alpha
mean(MCMC_beta0) #Estimator of the posterior mean of beta0
mean(MCMC_beta1) #Estimator of the posterior mean of beta1
mean(MCMC_tau) #Estimator of the posterior mean of tau

HPDinterval(MCMC_alpha, prob = 0.95) #HPD interval for the posterior distribution of alpha
HPDinterval(MCMC_beta0, prob = 0.95) #HPD interval for the posterior distribution of beta0
HPDinterval(MCMC_beta1, prob = 0.95) #HPD interval for the posterior distribution of beta1
HPDinterval(MCMC_tau, prob = 0.95) #HPD interval for the posterior distribution of tau


#5.

library(rjags)

M=length(x) 
alpha <- rep(0,M)
beta_0 <- rep(0,M)
beta_1 <- rep(0,M)
tau <- rep(0,M)

mu = c()
beta_0
beta_1
tau

for (i in 1:length(x)){
  mu[i] =  exp(beta_0+alpha*z[i]+(beta_1+tau*z[i])*(x[i]-39))
}

mu


# initialization of the parameters


mod= "model{
for (i in 1:length(x)){
y[i] ~ dpois(mu[i])
mu[i] = exp(beta_0+alpha*z[i]+(beta_1+tau*z[i])*(x[i]-39))
}


# specify the priors (uninformative here)
alpha ~ dnorm(0,0.000001)
beta_0 ~ dnorm(0,0.000001)
beta_1 ~ dnorm(0,0.000001)
tau ~ dnorm(0,0.000001)
}"

alpha
beta_0
beta_1
tau

set.seed(1348)
colnames(absences)= c('y','z','x')
data = as.list(absences)
param=c("alpha","beta_0","beta_1","tau")

jags <- jags.model(textConnection(mod),
                   data = data,
                   n.chains = 4, # number of parallele chains
                   n.adapt = 1000 # burn-in
)
update(jags, 10000)


samples <- coda.samples(jags,
                        param,
                        10000)
#plot(samples)

samples
samples[[1]]
apply(samples[[1]],2,median)
apply(samples[[1]],2,mean)



S_alpha <- samples[[1]][,1]
S_beta0 <- samples[[1]][,2]
S_beta1 <- samples[[1]][,3]
S_tau <- samples[[1]][,4]

S_alpha <- mcmc(data = S_alpha)
S_beta0  <- mcmc(data = S_beta0)
S_beta1 <- mcmc(data = S_beta1)
S_tau <- mcmc(data = S_tau)

geweke.diag(S_alpha, frac1=0.1, frac2=0.5)
geweke.diag(S_beta0, frac1=0.1, frac2=0.5)
geweke.diag(S_beta1, frac1=0.1, frac2=0.5)
geweke.diag(S_tau, frac1=0.1, frac2=0.5)

heidel.diag(S_alpha)
heidel.diag(S_beta0)
heidel.diag(S_beta1)
heidel.diag(S_tau)

par(mfrow=c(2,2))
traceplot(MCMC_alpha, main = "jags_alpha")
traceplot(MCMC_beta0, main = "jags_beta_0")
traceplot(MCMC_beta1, main = "jags_beta_1")
traceplot(MCMC_tau, main = "jags_tau")

median(S_alpha)
median(S_beta0)
median(S_beta1)
median(S_tau)

mean(S_alpha)
mean(S_beta0)
mean(S_beta1)
mean(S_tau)

HPDinterval(S_alpha, prob = 0.95) #HPD interval for the posterior distribution of alpha
HPDinterval(S_beta0, prob = 0.95) #HPD interval for the posterior distribution of beta0
HPDinterval(S_beta1, prob = 0.95) #HPD interval for the posterior distribution of beta1
HPDinterval(S_tau, prob = 0.95) #HPD interval for the posterior distribution of tau

#6.

manabsence <- function(age){ 
  result <- 2.48*(1.06^(age-39))
  return(result)
}

womanabsence <- function(age){ 
  result <- 2.48*1.47*(1.043^(age-39))
  return(result)
}

ratio <- function(age){
  result <- womanabsence(age)/manabsence(age)
  return(result)
}

par(mfrow = c(1,3))

plot(manabsence, xlab = "age", xlim = c(15,100))
plot(womanabsence, xlab = "age", xlim = c(15,100))
plot(ratio, xlab = "age", xlim = c(15,100), main = "Ratio womanabsence/manabsence")
abline(a = 1, b = 0, col = "red")

#7.(a)

x <- absences_whole_set[,3]
y <- absences_whole_set[,1]
z <- absences_whole_set[,2]

colnames(absences_whole_set)= c('y','z','x')

M=length(x) 
alpha <- rep(0,M)
beta_0 <- rep(0,M)
beta_1 <- rep(0,M)

tau <- rep(0,M)


mu = c()

for (i in 1:length(x)){
  mu[i] =  exp(beta_0+alpha*z[i]+(beta_1+tau*z[i])*(x[i]-39))
}

mu


# initialization of the parameters


mod= "model{
for (i in 1:length(x)){
y[i] ~ dpois(mu[i])
mu[i] = exp(beta_0+alpha*z[i]+(beta_1+tau*z[i])*(x[i]-39))
}


# specify the priors (uninformative here)
alpha ~ dnorm(0,0.000001)
beta_0 ~ dnorm(0,0.000001)
beta_1 ~ dnorm(0,0.000001)
tau ~ dnorm(0,0.000001)
z[93]  ~ dbern(0.5)
z[94]  ~ dbern(0.5)
z[95]  ~ dbern(0.5)
z[96]  ~ dbern(0.5)
x[97]  ~ dunif(18,65)
x[98]  ~ dunif(18,65)
x[99]  ~ dunif(18,65)
x[100]  ~ dunif(18,65)
}"

alpha
beta_0
beta_1
tau
z[93]
z[94]
z[95]  
z[96]  
x[97]  
x[98] 
x[99]  
x[100] 

set.seed(1348)
data = as.list(absences_whole_set)
param=c("alpha","beta_0","beta_1","tau", "z[93]", "z[94]", "z[95]", "z[96]", "x[97]", "x[98]", "x[99]", "x[100]")

jags <- jags.model(textConnection(mod),
                   data = data,
                   n.chains = 4, # number of parallel chains
                   n.adapt = 1000 # burn-in
)
update(jags, 10000)


samples <- coda.samples(jags,
                        param,
                        10000)
#plot(samples)
#Samples are ordered as follows: alpha, beta_0, beta_1, tau, x[100], x[97], x[98], x[99], z[93], z[94], z[95], z[96]

samples

median(samples[[1]][,1]) #alpha
median(samples[[1]][,2]) #beta_0
median(samples[[1]][,3]) #beta_1
median(samples[[1]][,4]) #tau

median(samples[[1]][,5]) #x[100]
median(samples[[1]][,6]) #x[97]
median(samples[[1]][,7]) #x[98]
median(samples[[1]][,8]) #x[99]

mean(samples[[1]][,9]) #z[93]
mean(samples[[1]][,10]) #z[94]
mean(samples[[1]][,11]) #z[95]
mean(samples[[1]][,12]) #z[96]

