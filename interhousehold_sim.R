traj <-readRDS("intrahouse_traj.rds") #input possible chagas trajectories
# "intrahouse_traj.rds" is produced by running intrahousehold_trajectories.R
#==============================================================================
# Default Parameter Values - 
#
# simlength: number of days in a simulation
# bbfeedingrate: feeding rate of bed bugs (no. of feeding events/day)
# bbprevalence: proportion of households with a bed bug infestation
# bbduration: length of each bed bug infestation
# bugtohumI: probability that an infectious bite infects a susceptible human
# crate: contact rate (amt of time /day each household spends visiting others)
# connect: amount of connectivity in the neighborhood network
#==============================================================================

simlength <- 1000 
bbduration <- 365 
bbfeedingrate <- 1/7 
bbprevalence <- .1 
bugtohumI <- 0.00058 
crate <- 1/30 
connect<-0.02 

#===============================================================================
# Function: set up simulation
#	  -Assign index cases for bed bug infestations
#   -Assign one chagas case
#===============================================================================

Setup<-function() {
  # Build a contact network with all "L" households in the neighborhood 
  NEIGH <- matrix(rbinom(L^2, 1, connect),L) 
  for( i in 1:L)
  {
    NEIGH[i,] <- NEIGH[,i]
  }  
  diag(NEIGH)<-0
  
  # Set up bed bug infestation status of the households in the neighborhood
  S <<- rep(1, L)
  E <<- rep(0, L)
  I <<- rep(0, L)
  recoverday <<-rep(0, L)  # the day that an infested household "recovers"
  
  bb.index <<- sample(L, bbprevalence*L) 
  I[bb.index] <<- 1
  S[bb.index] <<- 0
  recoverday[bb.index] <<- sample(1:bbduration, length(bb.index), replace=TRUE)
  bb.cases <<- bb.index
  
  BBExpose()
  
  # Set up Chagas status in the households:
  chagasI <<- rep(0,L)  # whether or not there is >=1 person with Chagas
  chagas.inf.day <<- rep(NA,L)  # first day that a household member gets Chagas
  chagas.traj<<-rep(NA,L)  # the trajectory drawn from traj array
  
  chagas.index<<-sample(1:L,1)  # choose a single household to have Chagas
  chagasI[chagas.index]<<-1
  chagas.inf.day[chagas.index]<<-0
  chagas.traj[chagas.index]<<-sample(1:dim(traj)[1],1)
  chagas.cases <<- chagas.index
}


#==============================================================================
#  Functions: interhousehold bed bug and Chagas transmission events
#==============================================================================

BBInfest <- function(node, day) {
  # Change the bed bug infestation status of selected households to I.
  #
  # Args:
  #   node: the household to be infested
  #   day: the day in the simulation that this function is called
  I[node] <<- 1
  E[node] <<- 0
  S[node] <<- 0
  bb.cases <<- which(I == 1)
  recoverday[node] <<- day + bbduration
}

BBExpose<-function() {
  # Expose households in contact with infested houses 
  if(length(bb.cases) == 1) {
    index <- NEIGH[, bb.cases] == 1
  }
  if(length(cases) > 1) {
    index <- rowSums(NEIGH[, bb.cases]) >= 1
  }
  if(length(cases) == 0){
    index <- integer(0)
  }
  #index<- index & !I # recovered are not really exposed
  E[index] <<- 1
  #E[cases] <<- 0
  S[index] <<- 0
}

BBRecover<-function(node) {
  # Recover infected households selected by "node"
  I[node] <<- 0
  E[node] <<- 0
  S[node] <<- 1  # recovered households are susceptible again
  bb.cases <<- which(I==1)
}

ChagasInfect <- function(chagashouse, day) {
  # Determine if Chagas free, bed bug infested houses get infected with Chagas
  # from connected households that have both Chagas and bed bugs.
  #
  # Args:
  #   chagashouse: houses with Chagas that may pose a risk to their contacts
  #   day: the day in the simulation that this function is called
  for (house in chagashouse) {
    if (I[house] == 1) {
      contacts <- which(NEIGH[,house] == 1)
      index <- contacts[which(I[contacts] == 1)] 
      random <- runif(length(index), 0, 1)
      test <- random < crate/(sum(NEIGH[, house])) * traj[chagas.traj[house],
                                              day - chagas.inf.day[house],
                                              2] * bbfeedingrate * bugtohumI
      # ^ contact rate/ (number of contacts) * bb prevalence * ...
      # this is the probability that a day of visitation leads to Chagas being
      # transmitted from an infected bug to a susceptible visitor
      chagasI[index][test] <<- 1
      chagas.traj[index][test] <<- sample(1:dim(traj)[1], 
                                          length(chagas.traj[index][test]))
      chagas.inf.day[index][test] <<- day
    }
  }
}

#===================================================================================
# Function: running the simulation
# b: probability that exposure to bed bugs leads to infestation, solved in each
# case so that bbprevalence can stay (approximately) constant throughout the 
# simulation
#===================================================================================

RunSim <- function(bbprevalence, connect, crate) {
  # Runs one entire interhousehold simulation with both bed bug and bed bug-
  # mediated Chagas transmission. 
  #
  # Args:
  #   bbprevalence, connect, and crate as defined in Default Parameter Values
  #   
  # Returns:
  #   Number of Chagas infected households at the end of the simulation
  bbprevalence <<- bbprevalence
  b <<- exp(-7.274+4.408*bbprevalence) #prob that bb exposure leads to infestation
  connect <<- connect
  crate <<- crate
  Setup()
  for(i in 2:simlength) {
    random <- runif(L)
    risk <- E*b		#probability that exposed households become infested
    BBInfest(random < risk, i)
    BBRecover(which(recoverday == i))
    BBExpose()
    ChagasInfect(which(chagasI==1),i)
  }
  return(sum(chagasI))
}

#===============================================================================
# Execution: connectivity sensitivity analysis
#===============================================================================

conrange<-c(seq(0.005,0.05,.005))
constore0.02<-rep(0,length(conrange))
k=1
for (con in conrange){
  jstore<-rep(0,200)
  for(j in 1:200){
    jstore[j]<-RunSim(0.1,con,0.02)
  }
  constore0.02[k]<-sum(jstore>1)
  k=k+1
}

conrange<-c(seq(0.005,0.05,.005))
constore0.01<-rep(0,length(conrange))
k=1
for (con in conrange){
  jstore<-rep(0,200)
  for(j in 1:200){
    jstore[j]<-RunSim(0.1,con,0.01)
  }
  constore0.01[k]<-sum(jstore>1)
  k=k+1
}
conrange<-c(seq(0.005,0.05,.005))
constore0.03<-rep(0,length(conrange))
k=1
for (con in conrange){
  jstore<-rep(0,200)
  for(j in 1:200){
    jstore[j]<-RunSim(0.1,con,0.03)
  }
  constore0.03[k]<-sum(jstore>1)
  k=k+1
}

data1=data.frame(conrange,constore0.01/200)
data2=data.frame(conrange,constore0.02/200)
data3=data.frame(conrange,constore0.03/200)

mod1 <- nls(constore0.01/200 ~ exp(a + b * conrange), data = data1, 
            start = list(a = 0, b = 0))
mod2 <- nls(constore0.02/200 ~ exp(a + b * conrange), data = data2, 
            start = list(a = 0, b = 0))
mod3 <- nls(constore0.03/200 ~ exp(a + b * conrange), data = data3, 
            start = list(a = 0, b = 0))


plot(conrange,constore0.01/200,xlab="Connectivity (c)",ylab="Probability",
     ylim=c(0,max(constore0.03/200)),col="slategray4")
lines(conrange, predict(mod1, list(x = conrange)),lwd=2,col="slategray4")
points(conrange,constore0.02/200,col="firebrick4")
lines(conrange, predict(mod2, list(x = conrange)),lwd=2,col="firebrick4")
points(conrange,constore0.03/200,col="skyblue4")
lines(conrange, predict(mod3, list(x = conrange)),lwd=2,col="skyblue4")
legend(x=.03,y=.17,legend=c("p=0.10","k=0.03","k=0.02","k=0.01"),
       col=c("white","skyblue4","firebrick4","slategray4"),bty="n",lty=c(1,1,1),
       lwd=2)
dev.print(device=pdf, "connect.pdf", width=5, height=5)

#===============================================================================
# Execution: contact rate sensitivity analysis
#===============================================================================

crange<-seq(.02,.2,.02)
cstore.05<-rep(0,length(crange))
k=1
for (c in crange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-RunSim(0.05,.02,c)
  }
  cstore.05[k]<-sum(jstore>1)
  k=k+1
}
cstore.1<-rep(0,length(crange))
k=1
for (c in crange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-RunSim(0.1,.02,c)
  }
  cstore.1[k]<-sum(jstore>1)
  k=k+1
}
crange<-seq(.02,.2,.02)
cstore.2<-rep(0,length(crange))
k=1
for (c in crange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-RunSim(0.2,.02,c)
  }
  cstore.2[k]<-sum(jstore>1)
  k=k+1
}
data1=data.frame(crange,cstore.05/200)
data2=data.frame(crange,cstore.1/200)
data3=data.frame(crange,cstore.2/200)

mod1 <- nls(cstore.05/200 ~ a + b * crange, data = data1, 
            start = list(a = 0, b = 0))
mod2 <- nls(cstore.1/200 ~ a + b * crange, data = data2, 
            start = list(a = 0, b = 0))
mod3 <- nls(cstore.2/200 ~ a + b * crange, data = data3, 
            start = list(a = 0, b = 0))

plot(crange,cstore.05/200,xlab="Contact rate (k)",ylab="Probability",
     ylim=c(0,max(cstore.2/200+.015)),col="slategray4")
lines(crange, predict(mod1, list(x = conrange)),lwd=2,col="slategray4")
points(crange,cstore.1/200,col="firebrick4")
lines(crange, predict(mod2, list(x = conrange)),lwd=2,col="firebrick4")
points(crange,cstore.2/200,col="skyblue4")
lines(crange, predict(mod3, list(x = conrange)),lwd=2,col="skyblue4")
legend(x=.02,y=.42,legend=c("c=0.02","p=0.20","p=0.10","p=0.05"),
       col=c("white","skyblue4","firebrick4","slategray4"),bty="n",lty=c(1,1,1),lwd=2)
dev.print(device=pdf, "contact.pdf", width=5, height=5)

#===============================================================================
# Execution: bed bug prevalence rate sensitivity analysis
#===============================================================================

bbprevrange<-seq(0.02,0.2,0.02)
bstore.02<-rep(0,length(bbprevrange))
k=1
for (bb in bbprevrange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-RunSim(bb,.02,0.02)
  }
  bstore.02[k]<-sum(jstore>1)
  k=k+1
}
bbprevrange<-seq(0.02,0.2,0.02)
bstore.05<-rep(0,length(bbprevrange))
k=1
for (bb in bbprevrange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-RunSim(bb,.02,0.05)
  }
  bstore.05[k]<-sum(jstore>1)
  k=k+1
}
bbprevrange<-seq(0.02,0.2,0.02)
bbstore.10<-rep(0,length(bbprevrange))
k=1
for (bb in bbprevrange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-RunSim(bb,.02,0.1)
  }
  bbstore.10[k]<-sum(jstore>1)
  k=k+1
}
data1=data.frame(bbprevrange,b.02)
data2=data.frame(bbprevrange,b.05)
data3=data.frame(bbprevrange,b.1)

mod1 <- nls(b.02 ~ a + b * bbprevrange, data = data1, start = list(a = 0, b = 0))
mod2 <- nls(b.05 ~ a + b * bbprevrange, data = data2, start = list(a = 0, b = 0))
mod3 <- nls(b.1 ~ a + b * bbprevrange, data = data3, start = list(a = 0, b = 0))


plot(bbprevrange,b.02,xlab="Bed bug prevalence (p)",ylab="Probability",ylim=c(0,max(b.1)),col="slategray4")
lines(bbprevrange, predict(mod1, list(x = bbprevrange)),lwd=2,col="slategray4")
points(bbprevrange,b.05,col="firebrick4")
lines(bbprevrange, predict(mod2, list(x = bbprevrange)),lwd=2,col="firebrick4")
points(bbprevrange,b.1,col="skyblue4")
lines(bbprevrange, predict(mod3, list(x = bbprevrange)),lwd=2,col="skyblue4")
legend(x=.02,y=.32,legend=c("c=0.02","k=0.02","k=0.05","k=0.10"),col=c("white","skyblue4","firebrick4","slategray4"),bty="n",lty=c(1,1,1),lwd=2)
dev.print(device=pdf, "prev.pdf", width=5, height=5)