# Restricted dynamic occupancy JS model for FKW 2000-2015
# Pseudospatial covariate 

## run model 
library(tidyverse)
library(coda)
library(shinystan)
library(postpack)
library(R2jags)
library(truncnorm)
library(mvtnorm)



rootdir<-"~/Desktop/JS_investigation/pseudospatialCR"
Date<-Sys.Date()
outdir<-file.path(rootdir,"output", Date) 
dir.create(outdir) #Create dated directory for output 


years<- 2000:2015
T<-length(years)

#read in data
data.orig<-read.csv(file.path(rootdir,"FKW_sightings.csv"),as.is=TRUE)%>%
  dplyr::filter(Distinctiveness >=3)

source(file.path(rootdir,"BuildCapHist.R"))
fkw.all<-captureHistories(data.orig,1,2)%>%
  mutate(Cluster = as.factor(Cluster))
fkw.all<-fkw.all[ , order(names(fkw.all))]



CH<-as.matrix(fkw.all[,c(1:T)])
clst<-as.integer(as.character(fkw.all$Cluster))
n.clust<-max(clst)

ol<-read.table(file.path(rootdir, "overlap.txt"))%>%as.matrix()
 

# Augment data
nz <- 200
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

CH.z <- CH.aug
CH.z[CH.z == 0] <- NA

mod.init.z <- function(ch, notseen){
  z.known <- ch # get known states
  z.known[z.known==notseen] <- NA
  z.init <- as.data.frame(t(z.known)) # transpose
  z.init <- tidyr::fill(z.init, names(z.init)) # fill all subsequent unknown states with previous known state
  z.init <- t(z.init) # transpose back
  z.init[!is.na(z.known)] <- NA # replace known values with NA
  return(z.init)
}


jags.data <- list(y = CH.aug, 
                  n.occasions = dim(CH.aug)[2],
                  M = dim(CH.aug)[1],
                  N.obs = dim(CH)[1],
                  z = CH.z,
                  ol= ol,
                  n.clust = n.clust,
                  clst = c(clst, rep(NA, nz)))
# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), 
                         mean.p = runif(dim(CH)[2], 0, 1), 
                         z = mod.init.z(CH.aug,0),
                         clst = c(rep(NA, dim(CH)[1]),sample(1:3, nz, replace=T)))}

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "b", 
                "Nsuper", "N", "B", "lli",
                "delta")

# MCMC settings
ni <- 40000 
nt <- 3
nb <- 10000 
nc <- 3


jags.model.txt <- function(){
  #--------------------------------------
  # Parameters:
  # phi: survival probability
  # gamma: removal entry probability
  # p: capture probability
  #--------------------------------------


  for (i in 1:M){
    for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
    } #t
    for (t in 1:n.occasions){
      logit(p[i,t]) <- mean.lp[t] + delta*ol[t, clst[i]]
    } #t 
  } #i
  mean.phi ~ dbeta(8, 1)      # Informative prior for phi (high survival rate ~ 0.94)
  delta ~ dnorm(0,0.001)
  
  for(t in 2:n.occasions){
    mean.p[t] ~ dunif(0,1)
    mean.lp[t] <- log(mean.p[t]/(1-mean.p[t]))
  }
  
  mean.p[1] ~ dunif(0,1)        # Constrain p (confounding with gamma)
  mean.lp[1] <- log(mean.p[1]/(1-mean.p[1]))
  
  
  
  for (t in 1:n.occasions){ 
    gamma[t] ~ dunif(0,1)
  } #t
  
  # Clusters
  for(i in 1:M){
    clst[i] ~ dcat(rep((1/n.clust), n.clust))
  }
  
  
  # Likelihood
  for (i in 1:M){
    z[i,1] ~ dbern(gamma[1]) 
    mu1[i] <- z[i,1] * p[i,1] # Observation process 
    y[i,1] ~ dbern(mu1[i])
    
    # Subsequent occasions
    for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
    } #t
  } #i
  
  
  # # Calculate derived population parameters
  
  for (t in 1:n.occasions){ 
    qgamma[t] <- 1-gamma[t] 
  }
  
  cprob[1] <- gamma[1]
  
  for (t in 2:n.occasions){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)]) 
  } #t
  
  psi <- sum(cprob[1:n.occasions])     # Inclusion probability
  
  for (t in 1:n.occasions){ 
    b[t] <- cprob[t] / psi         # Entry probability 
  } #t
  
  for (i in 1:M){
    recruit[i,1] <- z[i,1] 
    for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
    } #t 
  } #i
  for (t in 1:n.occasions){
    N[t] <- sum(z[1:M,t])      # Actual population size
    B[t] <- sum(recruit[1:M,t])         # Number of entries
  } #t
  for (i in 1:M){
    Nind[i] <- sum(z[i,1:n.occasions])
    Nalive[i] <- 1-equals(Nind[i], 0)
  } #i
  Nsuper <- sum(Nalive[1:M])         # Superpopulation size
  
  
  # loglike: for WAIC estimation 
  for(i in 1:N.obs){
    for(t in 2:(n.occasions-1)){
      llit[i,t] <- logdensity.bern(y[i,t],mu3[i,t])
    } # t primary periods
    lli[i] <- sum(llit[i,2:(n.occasions-1)])
  } # individuals N.obs
  
}


m <- jags(data  = jags.data,
          inits = inits,
          parameters.to.save = parameters,
          model.file = jags.model.txt,
          n.chains = nc,
          n.iter = ni,
          n.burnin = nb)
m

post<-mcmcplots::as.mcmc.rjags(m)

waic1.jags <- function(post,n.obs,loglike.grep.txt = "lli"){
  logsumexp <- function(x){ max.x <- max(x); max.x - log(length(x)) +log(sum(exp(x-max.x)))} # this is equivalient to log( 1/S * sum_s[exp{ x }]))    
  # get the columns indices that have 'lli' in their name (for loglike-inidividual)
  ll.col <-grep(loglike.grep.txt,colnames(post[[1]]))
  # get the mcmc values for the log values
  if(length(post)>1){ lli.samp <- as.matrix(do.call("rbind",post)[,ll.col]) } else { lli.samp <- as.matrix(post[[1]][,ll.col])}
  # get the complexity penality (WAIC1; from Gelman et al 2014)
  logElike <- apply(lli.samp,2,logsumexp)
  Eloglike <- apply(lli.samp,2,mean)
  p_waic1 <- 2*sum(logElike-Eloglike)
  # get the lppd log pointwaise predictive density
  # use the log-sum-exp trick to handle underflow issues with exp(loglike) ~= 0
  logsumf.i <- apply(lli.samp,2,logsumexp)
  lppd <- sum(logsumf.i)
  waic1 <- -2*(lppd-p_waic1)
  return(waic1)
}

waic1.jags(post, dim(CH)[1])

summ<- MCMCvis::MCMCsummary(m,params=c("N"))
summ$Year<-years
colnames(summ)<-c("mean","sd","lower","mid","upper","rhat","n.eff","Year")

library(ggplot2)
r<-ggplot(summ, aes(x=Year, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  labs(y=paste("Estimated Number of Distinctive Indivdiuals", sep=""))+
  theme_classic()
r

summ<- MCMCvis::MCMCsummary(m,params=c("mean.p","mean.phi","b","delta"))
summ


save.image(file.path(outdir, "ps_out.RData"))


summ<- MCMCvis::MCMCsummary(m,params=c("mean.p"))
summ$Year<-years
colnames(summ)<-c("mean","sd","lower","mid","upper","rhat","n.eff","Year")

library(ggplot2)
p<-ggplot(summ, aes(x=Year, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  labs(y=paste("Detection probability", sep=""))+
  theme_classic()
p
