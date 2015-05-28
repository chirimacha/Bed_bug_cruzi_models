# This code makes a histogram of the number of new human infections that result
# from intrahousehold feeding of bed bugs on humans.

simlength <- 180  # number of days per simulation
reps <- 100  # number of times to repeat the simulation

# Makeup of the household:
num.hum <- 5  # number of humans within the house
num.bug <- 1000  # number of bed bugs within the house

# Initial number of infected humans and bed bug:
hum.I0 <- 1
bug.I0 <- 0

# Parameter values:
feed.interval <- 7  # Average number of days between feeding
bug.lifetime <- 200  # Average bug lifespan in days
bugtohum.I <- 0.00058  # The probability of an infectious bite infecting a human
humtobug.I <- 0.042  # The probability of an infected human infecting a bug from a 
  # bite in the first three years of infection
humtobug.I2 <- 0.0048 #After first three years of infection

# Designating distributions to draw bug feeding days and lifetimes:
sd.feed <- 2  # sd of days between feeding for vectors
sd.buglife <- 30  # sd of days for bug lifetime

SimulateReps <- function(humtobug.I, reps){
  # SimulateReps runs the intra-household simulation for a specified number of
  # "reps"
  #
  # Args:
  #   humtobug.I: probability of an infected human infecting a susceptible bug
  #   reps: the number of times to repeat the simulation
  #
  # Returns:
  #   A vector of the number of humans infected at the end of each simulation 
  humtobug.I <- humtobug.I  
  total.hum.I <- as.numeric()
  for (rep in 1:reps) {
    
    # Vectors of infection status:
    hum.I <- c(rep(1,hum.I0),rep(0,num.hum - hum.I0))
    bug.I <- c(rep(1,bug.I0),rep(0,num.bug - bug.I0))
    
    # Vector of bug feeding days:
    bug.feed.days <- ceiling(runif(num.bug,0,feed.interval))
    
    # Death day of all bugs:
    bug.death.days <- ceiling(runif(num.bug,0,bug.lifetime))
    
  for (day in 1:simlength) {
    
    for (bug in 1:num.bug) {
      
      # determine if bug is feeding that day
      if (bug.feed.days[bug] == day) {  # begin feeding loop
        
        # determine which human the bug will feed on
        hum = sample(1:num.hum,1)
        
        # determine if uninfected bug is feeding on infected human
        if ((hum.I[hum]==1) & (bug.I[bug]==0)) {
          
          if (runif(1) < humtobug.I) {
            bug.I[bug] <- 1  # Bug gets infected with probability humtobug.I
          }
        }
      
        # determine if infected bug is feeding on uninfected human
        if ((hum.I[hum] == 0) & (bug.I[bug] == 1)) {
          if (runif(1)<bugtohum.I) {
            hum.I[hum] <- 1  # Human gets infected with probability bugtohum.I
          }
        } 
      
        # reset feeding day 
        bug.feed.days[bug] <- ceiling(rnorm(1, mean = feed.interval + day, 
                                           sd = sd.feed)) 
      
      }  # end feeding loop
    
      if (bug.death.days[bug] == day){
        bug.I[bug] <- 0  # reset infection status
        bug.death.days[bug] <- ceiling(rnorm(1, mean = bug.lifetime + day, 
                                             sd = sd.buglife))
      }
    }  # end bug loop
  }  # end day loop
  
  total.hum.I[rep] = sum(hum.I)
  
  }  # end rep loop
  
  return (total.hum.I)
}  # end function loop


hum.I.store <- SimulateReps (humtobug.I,reps)

new.hum.infections <- hum.I.store - 1 
#saveRDS(NewInf,file="newhuminfections.rds")

# Plotting a histogram of the new human infections:
hist(new.hum.infections, breaks=0:4, freq=FALSE, las = 1,
     xlab = "Number of New Human Infections", ylab = "Probability Density", 
     right=F,main=NA)
#dev.print(device=pdf, "hist.pdf", width=5, height=5)

