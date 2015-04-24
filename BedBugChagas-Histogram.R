
simlength = 180 #Number of days in the simulation
reps = 100 #Number of times to repeat the simulation

#What is the makeup of the household (all hosts and vectors)?
numHum = 5 
numBug = 1000

#Initial number infected amongst host and vector
HumI0 = 1
BugI0 = 0

#Parameters
feedInterval = 7 #Average number of days between feeding
BugLifetime = 200 #Average bug lifespan in days
BtoHIProb = 0.00058 #The probability of an infectious bite infecting a human
HtoBIProb = 0.042 #The probability of an infected human infecting a bug from a bite; First three years of infection 
HtoBIProb2 = 0.0048 #After first three years of infection

#Information about distributions to draw feeding days and bug lifetimes
sdFeed <- 2 # sd of days between feeding for vectors
sdBL <- 30 # sd of days for bug lifetime

#The transmission function: one simulation of length "reps", given a certain HtoBP (determined by whether patient has a chronic or acute infectionn)
Transmission <- function(HtoBP, reps){
  HtoBIProb <- HtoBP  
  #vector to store a count of the total number of humans infected at the end of the simulation
  SumHumI <- as.numeric()
  for (rep in 1:reps){
    #Initial values
    ##Vectors of infection status
    HumI <- c(rep(1,HumI0),rep(0,numHum - HumI0))
    BugI <- c(rep(1,BugI0),rep(0,numBug - BugI0))
    ##Vector of bug feeding days
    BugFD <- ceiling(runif(numBug,0,feedInterval))
    ##Death day of all bugs
    BugDD <- ceiling(runif(numBug,0,BugLifetime))
    
  for (day in 1:simlength){
    for (bug in 1:numBug){
      #determine if bug is feeding that day
      if (BugFD[bug]==day){#begin feeding loop
        #determine which human the bug will feed on
        H = sample(1:numHum,1)
        #determine if uninfected bug is feeding on infected human
        if ((HumI[H]==1) & (BugI[bug]==0)){
          if (runif(1)<HtoBIProb){
            BugI[bug] <- 1 #Bug gets infected with probability HtoBIProb
          }
        }
      
        #determine if infected bug is feeding on uninfected human
        if ((HumI[H]==0) & (BugI[bug]==1)){
          if (runif(1)<BtoHIProb){
            HumI[H] <- 1 #Human gets infected with probability BtoHIProb
          }
        } 
      
        #reset feeding day 
        BugFD[bug] = ceiling(rnorm(1, mean=feedInterval+day, sd=sdFeed)) 
      
      }#end feeding loop
    
      if (BugDD[bug] == day){
        BugI[bug] <- 0 #reset infection status
        BugDD[bug] <- ceiling(rnorm(1, mean=BugLifetime+day, sd=sdBL))
      }
    }#end bug loop
  }#end day loop
  SumHumI[rep] = sum(HumI)
  }#end rep loop
  return (SumHumI)
}#end function loop


Trans = Transmission (HtoBIProb,reps)
length(Trans[which(Trans>0)])


NewInf = Trans - 1
saveRDS(NewInf,file="NewInf.rds")
#hist(NewInf,breaks = c(0,1,2,3,4), xlab = "Are", ylab = "Frequency", xlim=c(0,4))
hist(NewInf,breaks=c(-.1,.9,1.9,2.9,3.9),freq=FALSE,xlab="Number of New Human Infections",ylab="Probability Density",main=NA,ylim=c(0,0.5))
dev.print(device=pdf, "hist.pdf", width=5, height=5)
#Assumptions/omissions:
  #No incubation period 
  #Bug is infected for life
  #No seasonal/temperature dependent cycling
#Ricardo Gartland 
