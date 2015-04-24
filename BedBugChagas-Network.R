
#install.packages("network")
#library(network)
#help("network-package")
#install.packages("animation")
#library(animation)


traj <-readRDS("array.rds") #input possible chagas trajectories
#===================================================================================
#  Parameter Values
#===================================================================================

days<-1000 #length of each simulation
b <-.0015 #probaility of bb infection given exposure per day (ie rate)
duration<-365 #length of each bed bug infestation
bbfeedingrate<-1/7 #bed bugs feed once a week
BtoHIProb<-0.00058 #probability a bite from an infected bug leads to infection

#default values for parameters undergoing sensitivity analysis
bbprevalence=.1 #proportion of households with a bed bug infestation
connect<-0.02 #amount of connectivity in the neighborhood network
crate <- 1/30 #contact rate (amount of time per day spent visiting other households)



#===================================================================================
#  Simulation Setup
#  -setup()  set S to 1, I, E, R to 0
#	-bbprev: an empty vector to keep track of bed bug prevalence at each timestep
#	-Assign index cases for bed bug infestations
# -Assign one chagas case
#===================================================================================


setup<-function()
{
  #build a contact network with all "L" households in the neighborhood 
  NEIGH<-matrix(rbinom(L^2,1,connect),L) 
  for( i in 1:L)
  {
    NEIGH[i,]<-NEIGH[,i]
  }  
  diag(NEIGH)<-0
  
  #Setting up for bed bugs:
  S<<-rep(1,L)
  E<<-rep(0,L)
  I<<-rep(0,L)
  recoverday<<-rep(0,L)
  
  i<-1
  bbprev<<-rep(0,days)
  
  index_case<<-sample(L,bbprevalence*L) #choose the houses to be initially infested with bed bugs
  I[index_case]<<-1
  S[index_case]<<-0
  cases<<-index_case
  #so they don't all recover on the same day:
  recoverday[index_case]<<-sample(1:duration,length(index_case),replace=TRUE)
  expose()
  bbprev[i]<<-sum(I)
  
  #Setting up for Chagas: 
  chagascases<<-rep(0,L) #vector storing whether there is a chagas infection in each household
  chagasdayinf<<-rep(NA,L) #vector storing the day a household gets infected with chagas
  chagastraj<<-rep(NA,L) #the trajectory a Chagas infected household follows (after one person gets infected)
  chagas_index<<-sample(1:L,1) #choose a single household to have Chagas
  chagascases[chagas_index]<<-1
  chagasdayinf[chagas_index]<<-0
  chagastraj[chagas_index]<<-sample(1:dim(traj)[1],1)
}


#===================================================================================
#  Key functions of stochastic simulation
#===================================================================================

infect<-function(node,day)
{
  I[node]<<-1
  E[node]<<-0
  S[node]<<-0
  cases<<-which(I==1)
  recoverday[node]<<-day+duration
}

expose<-function()
{
  
  if(length(cases)==1)
  {
    index<-NEIGH[,cases]==1
  }
  if(length(cases)>1)
  {
    index<-rowSums(NEIGH[,cases])>=1
  }
  if(length(cases)==0){
    index<-integer(0)
  }
  #index<- index & !I # recovered are not really exposed
  E[index]<<-1
  E[cases]<<-0
  S[index]<<-0
}

recover<-function(node)
{
  I[node]<<-0
  E[node]<<-0
  S[node]<<-1 #recovered households are susceptible again
  cases<<-which(I==1)
}

chagas<-function(chagashouse,day)
{
  for (house in chagashouse){
    if (I[house]==1){
      contacts<-which(NEIGH[,house]==1)
      index<-contacts[which(I[contacts]==1)]
      random<-runif(length(index),0,1)
      test<-random<crate/(connect*L)*traj[chagastraj[house],day-chagasdayinf[house],2]*bbfeedingrate*BtoHIProb
      chagascases[index][test]<<-1
      chagastraj[index][test]<<-sample(1:dim(traj)[1],length(chagastraj[index][test]))
      chagasdayinf[index][test]<<-day
    }
  }
}

#===================================================================================
#  Simulation loop
#	-reps: define how many time steps
#	-outputs a movie
#	-after movie plots prevalence over time
#===================================================================================

network<-function(bbprevalence,connect,crate){
  bbprevalence<<-bbprevalence
  b<<-exp(-7.274+4.408*bbprevalence)
  connect<<-connect
  crate<<-crate
  setup()
  for(i in 2:days)
  {
    random<-runif(L)
    risk<-E*b		#multiplying the risk by whether or not expsoed
    infect(random<risk,i)
    recover(which(recoverday==i))
    expose()
    chagas(which(chagascases==1),i)
    bbprev[i] <<- sum(I)

  }
  return(sum(chagascases))
}

####effect of connectivity###
conrange<-c(seq(0.005,0.05,.005))
constore0.02<-rep(0,length(conrange))
k=1
for (con in conrange){
  jstore<-rep(0,200)
  for(j in 1:200){
    jstore[j]<-network(0.1,con,0.02)
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
    jstore[j]<-network(0.1,con,0.01)
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
    jstore[j]<-network(0.1,con,0.03)
  }
  constore0.03[k]<-sum(jstore>1)
  k=k+1
}

data1=data.frame(conrange,constore0.01/200)
data2=data.frame(conrange,constore0.02/200)
data3=data.frame(conrange,constore0.03/200)

mod1 <- nls(constore0.01/200 ~ exp(a + b * conrange), data = data1, start = list(a = 0, b = 0))
mod2 <- nls(constore0.02/200 ~ exp(a + b * conrange), data = data2, start = list(a = 0, b = 0))
mod3 <- nls(constore0.03/200 ~ exp(a + b * conrange), data = data3, start = list(a = 0, b = 0))


plot(conrange,constore0.01/200,xlab="Connectivity (c)",ylab="Probability",ylim=c(0,max(constore0.03/200)),col="slategray4")
lines(conrange, predict(mod1, list(x = conrange)),lwd=2,col="slategray4")
points(conrange,constore0.02/200,col="firebrick4")
lines(conrange, predict(mod2, list(x = conrange)),lwd=2,col="firebrick4")
points(conrange,constore0.03/200,col="skyblue4")
lines(conrange, predict(mod3, list(x = conrange)),lwd=2,col="skyblue4")
legend(x=.03,y=.17,legend=c("p=0.10","k=0.03","k=0.02","k=0.01"),col=c("white","skyblue4","firebrick4","slategray4"),bty="n",lty=c(1,1,1),lwd=2)
dev.print(device=pdf, "connect.pdf", width=5, height=5)
####effect of contact rate###
crange<-seq(.02,.2,.02)
cstore.05<-rep(0,length(crange))
k=1
for (c in crange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-network(0.05,.02,c)
  }
  cstore.05[k]<-sum(jstore>1)
  k=k+1
}
cstore.1<-rep(0,length(crange))
k=1
for (c in crange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-network(0.1,.02,c)
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
    jstore[j]<-network(0.2,.02,c)
  }
  cstore.2[k]<-sum(jstore>1)
  k=k+1
}
data1=data.frame(crange,cstore.05/200)
data2=data.frame(crange,cstore.1/200)
data3=data.frame(crange,cstore.2/200)

mod1 <- nls(cstore.05/200 ~ a + b * crange, data = data1, start = list(a = 0, b = 0))
mod2 <- nls(cstore.1/200 ~ a + b * crange, data = data2, start = list(a = 0, b = 0))
mod3 <- nls(cstore.2/200 ~ a + b * crange, data = data3, start = list(a = 0, b = 0))


plot(crange,cstore.05/200,xlab="Contact rate (k)",ylab="Probability",ylim=c(0,max(cstore.2/200+.015)),col="slategray4")
lines(crange, predict(mod1, list(x = conrange)),lwd=2,col="slategray4")
points(crange,cstore.1/200,col="firebrick4")
lines(crange, predict(mod2, list(x = conrange)),lwd=2,col="firebrick4")
points(crange,cstore.2/200,col="skyblue4")
lines(crange, predict(mod3, list(x = conrange)),lwd=2,col="skyblue4")
legend(x=.02,y=.42,legend=c("c=0.02","p=0.20","p=0.10","p=0.05"),col=c("white","skyblue4","firebrick4","slategray4"),bty="n",lty=c(1,1,1),lwd=2)
dev.print(device=pdf, "contact.pdf", width=5, height=5)
####effect of bbprevalence###
bbprevrange<-seq(0.02,0.2,0.02)
bstore.02<-rep(0,length(bbprevrange))
k=1
for (bb in bbprevrange){
  jstore<-rep(0,200)
  for (j in 1:200){
    jstore[j]<-network(bb,.02,0.02)
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
    jstore[j]<-network(bb,.02,0.05)
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
    jstore[j]<-network(bb,.02,0.1)
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