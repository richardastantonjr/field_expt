## simulate data to see how the robust design occupancy model data are structured before fitting in
## unmarked. Modified from the Kery and Chandler's unmarked vignette, chapter on the colext 
## function

M <- 24 # Number of sites
J <- 2 # num secondary sample periods
T <- 5 # num primary sample periods
psi <- rep(NA, T) # Occupancy probability
muZ <- z <- array(dim = c(M, T)) # Expected and realized occurrence
y <- array(NA, dim = c(M, J, T)) # Detection histories
set.seed(13973)
psi[1] <- 0.4 # Initial occupancy probability
p <- c(0.3,0.4,0.5,0.5,0.1,0.3,0.5,0.5,0.6,0.2)
phi <- runif(n=T-1, min=0.6, max=0.8) # Survival probability (1-epsilon)
gamma <- runif(n=T-1, min=0.1, max=0.2) # Colonization probability
# Generate latent states of occurrence
 # First year
 z[,1] <- rbinom(M, 1, psi[1]) # Initial occupancy state
# Later years
  for(i in 1:M){ # Loop over sites
    for(k in 2:T){ # Loop over years
      muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1]
      z[i,k] <- rbinom(1, 1, muZ[k])
    }
  }
 # Generate detection/non-detection data
   for(i in 1:M){
    for(k in 1:T){
      prob <- z[i,k] * p[k]
      for(j in 1:J){
        y[i,j,k] <- rbinom(1, 1, prob)
      }
    }
  }
# Compute annual population occupancy
  for (k in 2:T){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
  }
 
 ## y is the array of sites*visits*periods of interest. Analysis using colext function in unmarked
 ## requires the modifications below
 yy <- matrix(y, M, J*T)  ## Note this is easier to arrive at directly when massaging the real data
                          ## a site by visit matrix of detection-nondetection- expanded into an
                          ## array with nspec slices.
 year <- matrix(c("01","02","03","04","05"),nrow(yy), T, byrow=TRUE)
 
 
 ##covariates can also be supplied using the arguments siteCovs, yearlySiteCovs, and obsCovs. 
 ## yearlySiteCovs must have M rows and T columns.
 