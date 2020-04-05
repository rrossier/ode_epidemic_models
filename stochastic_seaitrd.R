library(odin)
library(rjson)

## function to make colors transparent, stolen from adegenet::transp
transp <- function (col, alpha = 0.5) {
  res <- apply(col2rgb(col), 2,
               function(c) rgb(c[1]/255, c[2]/255, 
                               c[3]/255, alpha))
  return(res)
}

## x: instance of odin model
## t: time steps
## n: number of replicates
run_model <- function(x, t = 0:100, n = 1, ...) {
  res <- x$run(t, replicate = n, ...)
  res <- x$transform_variables(res)
  res <- cbind.data.frame(t = res[[1]], res[-1])
  details = fromJSON(x[["ir"]])
  names = do.call(c, lapply(details[["data"]][["variable"]][["contents"]], function(i){return(i$name)}))
  attr(res, "n_compartments") <- length(names) #length(x$names) - 1
  attr(res, "n_replicates") <- n
  attr(res, "compartments") <- names
  class(res) <- c("pretty_odin", class(res))
  res
}


## plot function
plot.pretty_odin <- function(x, palette = NULL, transparency=FALSE,...) {
  ## handle colors
  n_compartments <- attr(x, "n_compartments")
  n_replicates <- attr(x, "n_replicates")
  if(!is.null(palette)){
    pal = colorRampPalette(hcl.colors(n_compartments,palette=palette))
  }else{
    defaults_col <- c("#587CFF", "#BCDA61", "#BA78A7", "#E2473F", "#FDE476", "#85CDE3", "#FDA31A", "#000000")
    pal <- colorRampPalette(defaults_col)
  }
  col_leg <- pal(n_compartments)
  alpha <- max(10 / n_replicates, 0.05)
  if(transparency){
    col <- rep(transp(col_leg, alpha), each = n_replicates)
  }else{
    col <- rep(col_leg, each = n_replicates)
  }
  ## make plot
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 3)
  matplot(x[, 1], x[, -1], xlab = "Time", ylab = "Number of individuals",
          type = "l", col = col, lty = 1, lwd=2, ...)
  legend("topright", lwd = 2, col = col_leg, bty = "n",
         legend = attr(x, "compartments"))
}


seaitrd_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_S_naturalD - n_SE + n_RS
  update(E) <- E + n_SE - n_EAI - n_E_naturalD + n_import
  update(A) <- A + n_EA - n_AR - n_A_naturalD + n_import
  update(I) <- I + n_EI - n_IT - n_IR - n_ID - n_I_naturalD + n_import
  update(Tr) <- Tr + n_IT - n_TR - n_TD - n_T_naturalD
  update(R) <- R + n_AR + n_IR + n_TR - n_RS - n_R_naturalD + n_import
  update(D) <- D + n_ID + n_TD
  update(naturalD) <- naturalD + n_S_naturalD + n_E_naturalD + n_A_naturalD + n_I_naturalD + n_T_naturalD + n_R_naturalD
  
  ## Individual probabilities of transition:
  p_SE <- 1 - exp(-beta * AIT / N)
  p_RS <- 1 - exp(-omega) # R to S
  p_EAI <- 1 - exp(-kappa) # E to A and I
  #p_EA <- 1 - exp(-kappa*(1-rho)) # E to A
  #p_EI <-  1 - exp(-kappa*rho) # E to I
  p_AR <-  1 - exp(-gamma1) # A to R
  p_IT <-  1 - exp(-alpha) # I to T
  p_IR <-  1 - exp(-gamma1) # I to R
  p_ID <-  1 - exp(-delta2) # I to D
  p_TR <-  1 - exp(-gamma2) # T to R
  p_TD <-  1 - exp(-delta1) # T to D
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_RS <- rbinom(R, p_RS)
  n_EAI <- rbinom(E, p_EAI)
  #n_EA <- rbinom(E, p_EA)
  #n_EI <- rbinom(E, p_EI)
  n_AR <- rbinom(A, p_AR)
  n_IT <- rbinom(I, p_IT)
  n_IR <- rbinom(I, p_IR)
  n_ID <- rbinom(I, p_ID)
  n_TR <- rbinom(Tr, p_TR)
  n_TD <- rbinom(Tr, p_TD)
  
  n_S_naturalD <- rbinom(S, mu)
  n_E_naturalD <- rbinom(E, mu)
  n_A_naturalD <- rbinom(A, mu)
  n_I_naturalD <- rbinom(I, mu)
  n_T_naturalD <- rbinom(Tr, mu)
  n_R_naturalD <- rbinom(R, mu)
  
  n_EAI_2d[] <- rmultinom(n_EAI, p)
  p[1] <- 1 - rho
  p[2] <- rho
  dim(p) <- 2
  dim(n_EAI_2d) <- 2
  n_EI <- n_EAI_2d[1]
  n_EA <- n_EAI_2d[2]

  n_import <- rpois(epsilon)
  
  ## Total population size, and number of infecteds
  AIT <- A + I + Tr
  N <- S + E + AIT + R + D
  
  ## Initial states
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(A) <- 0
  initial(I) <- 0
  initial(Tr) <- 0
  initial(R) <- 0
  initial(D) <- 0
  initial(naturalD) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1e6) # susceptibles
  E_ini <- user(1) # exposed
  beta <- user(5e-2) # exposition rate
  delta1 <- user(8e-4) # death rate during treatment
  delta2 <- user(3e-3) # death rate without treatment
  alpha <- user(1e-2) # treating rate
  gamma1 <- user(5e-3) # recovery rate without treatment
  gamma2 <- user(1e-2) # recovery rate during treatment
  kappa <- user(5e-1) # infection rate
  rho <- user(5e-1) # asymptomatic rate
  mu <- user(8e-6) # natural death rate
  omega <- user(5e-3) # rate of waning immunity
  epsilon <- user(1e-1) # import case rate
  
}, verbose = FALSE)
#devtools::find_rtools() # run this line if the return is "I don't see gcc on the PATH"

#set.seed(1)
seaitrd <- seaitrd_generator()
model_results <- run_model(seaitrd, t = 0:2000, n = 100)
plot.pretty_odin(model_results)
sum(model_results[2001,-1])
