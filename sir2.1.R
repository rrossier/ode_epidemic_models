## Define as a function
epi203 <- function(pars, initialConditions){
  
  ## Show parameters
  print(pars)
  
  ## Additional parameters
  times <- seq(from = 0, to = 600, by = 1)              # we want to run the model for 3000 time steps
  yinit <- c(Susc = initialConditions$susc, Infected = initialConditions$infected, Recovered = initialConditions$recovered) # this parameter sets the initial conditions
  
  ## below is the code for the actual model including the equations that you should recognize
  SIR_model <- function(times, yinit, pars){
    
    with(as.list(c(yinit,pars)), {
      
      dSusc      <- birth - beta*Infected*Susc                     - death*Susc
      dInfected  <-         beta*Infected*Susc - recovery*Infected - death*Infected
      dRecovered <-                              recovery*Infected - death*Recovered
      
      return(list(c(dSusc, dInfected, dRecovered)))})
  }
  
  ## run the ode solver for the function specified (function defined above is used)
  ## return the value of each compartment (Susc, Infected, Recovered) for each time step.
  results <- ode(func = SIR_model, times = times, y = yinit, parms = pars)
  results <- as.data.frame(results)
  
  ## Return result
  return(results)
}

##############################################################################

test.pars <- c(beta = 0.1, recovery = 0.005, death = 0.001, birth = 0.001)
initialConditions = list(susc=0.99, infected=0.01, recovered=0)
results   <- epi203(test.pars, initialConditions)

##############################################################################
## Plotting
matplot(results[, 1], results[, 2:4], type="l", lty=1);legend("topright", col=1:3, legend=c("S", "I", "R"), lwd=1)

