

## functions to compute lifespan inequality measures from a life table
## source code derived from Riffe et al. (2023) LifeIneq R package
## available at: https://github.com/alysonvanraalte/LifeIneq


## function to compute e-dagger
edag_fun <- function(age, mx, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  ## truncate life table if missing mx (cohort tables)
  m <- m_trunc <- length(age)
  if (any(is.na(mx))){
    age <- age[!is.na(mx)]
    dx <- dx[!is.na(mx)]
    lx <- lx[!is.na(mx)]
    ex <- ex[!is.na(mx)]
    ax <- ax[!is.na(mx)]
    mx <- mx[!is.na(mx)]
    m_trunc <- length(age)
  }
  
  # length of the age interval
  n <- c(diff(age),1)
  explusone <- c(ex[-1],ex[length(age)])
  # the average remaining life expectancy in each age interval 
  # (as opposed to the beginning of the interval)
  # ends up being roughly half of the ex between ages
  ex_average <- ex + ax / n * (explusone - ex)
  
  if (m_trunc==0){
    rep(NA,m)
  }else{
    if (m == m_trunc){
      rev(cumsum(rev(dx * ex_average))) / lx 
    }else{
      c(rev(cumsum(rev(dx * ex_average))) / lx, rep(NA,m-m_trunc)) 
    }
  }
    
}

## function to compute variance in ages at deaths
var_fun <- function(age, mx, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  ## truncate life table if missing mx (cohort tables)
  m <- m_trunc <- length(age)
  if (any(is.na(mx))){
    age <- age[!is.na(mx)]
    dx <- dx[!is.na(mx)]
    lx <- lx[!is.na(mx)]
    ex <- ex[!is.na(mx)]
    ax <- ax[!is.na(mx)]
    mx <- mx[!is.na(mx)]
    m_trunc <- length(age)
  }
  
  age0 <- age - age[1]
  n   <- length(age)
  out <- rep(NA, n)
  
  for (i in 1:n){
    axAge    <- age0[1:(n+1-i)] + ax[i:n]
    out[i] <- sum(dx[i:n] * (axAge - ex[i])^2) / lx[i]
  }
  
  if (m_trunc==0){
    rep(NA,m)
  }else{
    
    if (m == m_trunc){
      out
    }else{
      c(out, rep(NA,m-m_trunc)) 
    }
  }
    
  
  
}

## function to compute standard deviation in ages at deaths
sd_fun <- function(age, mx, dx, lx, ex, ax){
 sqrt(var_fun(age, mx, dx, lx, ex, ax))
  
}




