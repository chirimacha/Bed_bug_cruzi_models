# This code saves a number (default value: 10,000) of possible trajectories (number
# of infected humans and bed bugs at each time-step) as a file entitled 
# "intrahouse_traj.rds"

simlength <- 1000 # number of days in the simulation
reps <- 1e4  # number of times to repeat the simulation

# Number of humans and bed bugs in a hosehold:
num.hum <- 5 
num.bug <- 1000

# Initial number of infected humans and bugs:
humI0 <- 1
bugI0 <- 0

# Parameter values:
feed.interval <- 7   # average number of days between feeding
sd.feed <- 2  # sd of days between feeding for vectors
bug.lifetime <- 200  # average bug lifespan in days
sd.buglife <- 30 # sd of days for bug lifetime
bugtohumI <- 0.00058  # probability of an infectious bite infecting a human
humtobugI <- 0.042  # probability of an infected human infecting a bug from a 
  #bite during the first three years of infection 
humtobugI2 <- 0.0048 #After first three years of infection

Simulate <- function(){
  # Runs a single intra-household simulation 
  #
  # Returns:
  #   A data table of the number of the number of humans and bed bugs infected
  #   at each time-step (day) 
  
  # Vectors to store total number of infecteds at the end of the simulation:
  hum.prev <- as.numeric()
  bug.prev <- as.numeric()
    
  # Vectors of infection status:
  humI <- c(rep(1,humI0),rep(0,num.hum - humI0))
  bugI <- c(rep(1,bugI0),rep(0,num.bug - bugI0))
  
  # Vectors of bug feeding and death days:
  bug.feed <- ceiling(runif(num.bug, 0, feed.interval))
  bug.death <- ceiling(runif(num.bug, 0, bug.lifetime))
    
  for (day in 1:simlength) {
    for (bug in 1:num.bug) {
      #determine if bug is feeding that day
      if (bug.feed[bug] == day) {#begin feeding loop
        #determine which human the bug will feed on
        hum <- sample(1:num.hum, 1)
        
        #determine if uninfected bug is feeding on infected human
        if ((humI[hum] == 1) & (bugI[bug] == 0)){
          if (runif(1) < humtobugI){
            bugI[bug] <- 1 #Bug gets infected with probability humtobugI
          }
        }
      
        #determine if infected bug is feeding on uninfected human
        if ((humI[hum] == 0) & (bugI[bug] == 1)){
          if (runif(1) < bugtohumI){
            humI[hum] <- 1 #Human gets infected with probability bugtohumI
          }
        } 
      
        #reset feeding day 
        bug.feed[bug] = ceiling(rnorm(1, mean=feed.interval+day, sd=sdFeed)) 
      
      }#end feeding loop
    
      if (bug.death[bug] == day){
        bugI[bug] <- 0 #reset infection status
        bug.death[bug] <- ceiling(rnorm(1, mean=bug.lifetime+day, sd=sd.buglife))
      }
    }#end bug loop
    hum.prev[day] <- sum(humI)
    bug.prev[day] <- sum(bugI)
  }#end day loop
  data <- cbind(hum.prev, bug.prev)
  return (data)
}#end function loop


output <- array(0,dim=c(reps,simlength,2)) 
for(rep in 1:reps){
  output[rep,,] <- Simulate(humtobugI)
}
# note: takes ~ 24h to get 10,000 reps = 20MIL values
saveRDS(output,file="intrahouse_traj.rds")

# intrahouse_traj.rds is a 3-D array of dimension reps x simlength x 2 
# it stores the number of infected humans and bed bugs at each time step for a
# number of simulations (specified by "reps")