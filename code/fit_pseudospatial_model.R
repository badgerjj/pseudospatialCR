# Jolly-Seber model for population analysis of false killer whales

dir<-"~/pseudospatialCR"
source(file.path(dir, "code","helper.R"))

outdir<-file.path(dir,"output")
#dir.create(outdir) #Create directory for output 

Date<-Sys.Date()
outdir<-file.path(outdir, Date) 
dir.create(outdir) #Create dated directory for output 



## Read in capture history data

load(file.path(dir, "data","fkw_CH_2000_2015.RData"))

years<-2000:2015
n.occasions<-length(years)
CH<-as.matrix(fkw.CH[,c(1:n.occasions)])
nind = dim(CH)[1]

clst<-as.integer(as.character(fkw.CH$Cluster))
n.clust<- length(unique(clst))

## Read in overlap information 

overlap<-read.table(file.path(dir,"data","processed_data","overlap_3clust.txt"))

overlap<-reshape2::melt(as.data.frame(overlap))
colnames(overlap)<-c("cluster","overlap")
overlap$year<-rep(1:n.occasions, length.out=length(overlap$cluster))
o<-ggplot(data=overlap, aes(y=overlap, x=year))+
  geom_line(aes(color=cluster))+
  theme_classic()

o

overlap$overlap.s<-(overlap$overlap-mean(overlap$overlap))/sd(overlap$overlap) #standardize

ol<-matrix(NA, nrow=n.occasions, ncol=n.clust)
ol[,1]<-overlap$overlap.s[overlap$cluster=="V1"]
ol[,2]<-overlap$overlap.s[overlap$cluster=="V2"]
ol[,3]<-overlap$overlap.s[overlap$cluster=="V3"]


## Prep data for JAGS model 

nz <-1.5*nind       # Number of pseudoindividuals
y.aug <-rbind(CH,matrix(0,ncol=n.occasions,nrow=nz))      # augmented dataset

y.aug<-cbind(rep(0, nz+nind),y.aug)        #dummy column 
ol<-rbind(rep(0,3), ol) 

y.aug[y.aug==0]<-2         # seen=1, notseen = 2

y.jags=y.aug

jags.model.txt <- function(){
  #--------------------------------------
  # Parameters:
  # phi: survival probability
  # gamma: removal entry probability
  # p: capture probability
  #--------------------------------------
  # States (S):
  # 1 not yet entered
  # 2 alive
  # 3 dead
  # Observations (O):
  # 1 seen
  # 2 not seen
  #--------------------------------------
  
  # Priors and constraints
  
  # priors on capture (NOTE: constraint ON p[1])
  for(t in 1:(n.occasions-1)){ 
    for(i in 1:M){
      logit(p.rv[i,t])<-mean.lp[t] + delta*ol[t, clst[i]]
      
      # set deterministic nodes to enforce constraint
      p[i,(t+1)] <- p.rv[i,t]
    } 
  }
  
  
  p[1:M,1] <- p[1:M,2]
  
  for(t in 1:n.occasions){
    mean.p[t] ~ dunif(0, 1)
    mean.lp[t] <- log(mean.p[t]/(1-mean.p[t]))
  }
  
  for(t in 1:(n.occasions-1)){
    gamma[t] ~ dunif(0,1)
  }
  
  for(i in 1:M){
    clst[i] ~ dcat(rep((1/n.clust), n.clust))
  }
  
  # priors on survival (NOTE: constraint on final phi)
  for(t in c(1:(n.occasions-1))){
    phi[t] <-mean.phi
  }
  
  mean.phi ~ dbeta(8,2)
  
  delta ~ dnorm(0,0.001)
  
  # Define state-transition and observation matrices
  # Define probabilities of state S(t+1) given S(t)
  for (t in 1:(n.occasions-1)){
    ps[1,t,1] <- 1 - gamma[t]
    ps[1,t,2] <- gamma[t]
    ps[1,t,3] <- 0
    ps[2,t,1] <- 0
    ps[2,t,2] <- phi[t]
    ps[2,t,3] <- 1 - phi[t]
    ps[3,t,1] <- 0
    ps[3,t,2] <- 0
    ps[3,t,3] <- 1
    for(i in 1:M){
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 1
      po[2,i,t,1] <- p[i,t]
      po[2,i,t,2] <- 1 - p[i,t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
    }
  }
  # Likelihood
  for (i in 1:M){
    # Define latent state at first occasion
    z[i,1] <- 1 # Make sure that all M individuals are in state 1 at t=1
    
    for (t in 2:n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], t-1, 1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t],i,t-1,1:2])
    } #t
    
  } #i
  
  # # Calculate derived population parameters
  for (t in 1:(n.occasions-1)){
    qgamma[t] <- 1 - gamma[t]
  }
  cprob[1] <- gamma[1]
  for (t in 2:(n.occasions-1)){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } #t
  psi <- sum(cprob[]) # Inclusion probability
  for (t in 1:(n.occasions-1)){
    b[t] <- cprob[t] / psi # Entry probability
  } #t
  for (i in 1:M){
    for (t in 2:n.occasions){
      al[i,t-1] <- equals(z[i,t], 2)
    } #t
    for (t in 1:(n.occasions-1)){
      d[i,t] <- equals(z[i,t] - al[i,t],0)
    } #t
    alive[i] <- sum(al[i,])
  } #i
  for (t in 1:(n.occasions-1)){
    N[t] <- sum(al[,t]) # Actual population size
    B[t] <- sum(d[,t]) # Number of entries
    
  } #t
  for (i in 1:M){
    w[i] <- 1 - equals(alive[i],0)
  } #i
  Nsuper <- sum(w[]) # Superpopulation size
  
  
  
  
  # loglik: for WAIC estimation 
  for(i in 1:N.obs){
    for(t in 1:(n.occasions-1)){
      llit[i,t] <- logdensity.cat(y[i,t], po[z[i,t],i,t,1:2])
    } # t primary periods
    lli[i] <- sum(llit[i,1:(n.occasions-1)])
  } # individuals N.obs
}

jags.data <- list(y = y.jags, 
                  n.occasions = dim(y.jags)[2], 
                  M = dim(y.jags)[1],
                  ol=ol,
                  clst = c(clst, rep(NA, nz)),
                  n.clust=n.clust,
                  N.obs = dim(fkw.CH)[1])
# Initial values
zinit <- cbind(rep(NA, dim(y.jags)[1]), y.jags[,-1])
zinit[zinit==1] <- 2
inits <- function(){list(delta=0.1, 
                         z = zinit, 
                         clst = c(rep(NA, dim(CH)[1]), rcat(nz, prob =rep((1/n.clust), n.clust))))}
# Parameters monitored
parameters <- c("mean.p", "N", "delta", "b", "B", "gamma", "phi")

# MCMC settings
ni <- 10000
nb <- 2500
nc <- 3
m <- jags(data  = jags.data,
          inits = inits,
          parameters.to.save = parameters,
          model.file = jags.model.txt,
          n.chains = nc,
          n.iter = ni,
          n.burnin = nb)
#m

k<-mcmcplots::as.mcmc.rjags(m)%>%as.shinystan()
summ<- MCMCvis::MCMCsummary(m,params=c("N"))
summ$Year<-2000:2015

colnames(summ)<-c("mean","sd","lower","mid","upper","rhat","n.eff","Year")

r<-ggplot(summ, aes(x=Year, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  labs(y="Estimated abundance")+
  theme_classic()
r







