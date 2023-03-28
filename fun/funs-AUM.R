

## functions for computing the AUM index 
## based on E(T), E(T^2), E(T*Lambda(T))


ETLambdaTCond <- function(age,mx){
  ## length of input and output 
  m <- length(age)
  out <- rep(NA,m)
  ## set 0 death rates at older ages to NA
  mx[mx==0 & age>80] <- NA
  ## check and remove NAs in mx
  age <- age[!is.na(mx)]
  mx <- mx[!is.na(mx)]
  m <- length(age)
  ## upper bound and delta of interval
  age1 <- age + 1
  age1[m] <- age1[m] + 10000
  delta <- age1-age
  ## compute result for each age
  for (i in 1:(m-1)){
    a <- age[i]
    inda <- sum(age <= a)
    # Compute Lambda_T(a) and S_T(a)
    LambdaTa <- Cumhaz(delta,mx,inda-1)+(a-age[inda])*mx[inda]
    STa <- exp(-LambdaTa)
    res <- mx[inda]*(exp(-Cumhaz(delta,mx,inda-1))*exp(age[inda]*mx[inda]))*(Int2(a,age1[inda],mx[inda]) - a*Int1(a,age1[inda],mx[inda]))
    for(j in (inda+1):m){
      term <- (exp(-Cumhaz(delta,mx,j-1))*exp(age[j]*mx[j])) * ( (Cumhaz(delta,mx,j-1)-age[j]*mx[j] - LambdaTa) * Int1(age[j],age1[j],mx[j]) + mx[j]*Int2(age[j],age1[j],mx[j]) )
      if (is.na(term)) term <- 0 
      res <- res + term
    }
    out[i] <- res/STa
  }
  return(out)
}

ETCond <- function(age,mx){
  ## length of input and output 
  m <- length(age)
  out <- rep(NA,m)
  ## set 0 death rates at older ages to NA
  mx[mx==0 & age>80] <- NA
  ## check and remove NAs in mx
  age <- age[!is.na(mx)]
  mx <- mx[!is.na(mx)]
  m <- length(age)
  ## upper bound and delta of interval
  age1 <- age + 1
  age1[m] <- age1[m] + 10000
  delta <- age1-age
  ## compute result for each age
  for (i in 1:(m-1)){
    a <- age[i]
    inda <- sum(age <= a)
    # Compute Lambda_T(a) and S_T(a)
    LambdaTa <- Cumhaz(delta,mx,inda-1)+(a-age[inda])*mx[inda]
    STa <- exp(-LambdaTa)
    res <- ( exp(-Cumhaz(delta,mx,inda-1))*exp(age[inda]*mx[inda]) ) * Int1(age[inda],age1[inda],mx[inda]) 
    for(j in (inda+1):m){
      term <- ( exp(-Cumhaz(delta,mx,j-1))*exp(age[j]*mx[j]) ) * Int1(age[j],age1[j],mx[j]) 
      if (is.na(term)) term <- 0 
      res <- res + term
    }
    out[i] <- res/STa
  }
  return(out)
}


ET2Cond <- function(age,mx){
  ## length of input and output 
  m <- length(age)
  out <- rep(NA,m)
  ## set 0 death rates at older ages to NA
  mx[mx==0 & age>80] <- NA
  ## check and remove NAs in mx
  age <- age[!is.na(mx)]
  mx <- mx[!is.na(mx)]
  m <- length(age)
  ## upper bound and delta of interval
  age1 <- age + 1
  age1[m] <- age1[m] + 10000
  delta <- age1-age
  ## compute result for each age
  for (i in 1:(m-1)){
    a <- age[i]
    inda <- sum(age <= a)
    # Compute Lambda_T(a) and S_T(a)
    LambdaTa <- Cumhaz(delta,mx,inda-1)+(a-age[inda])*mx[inda]
    STa <- exp(-LambdaTa)
    res <- ( exp(-Cumhaz(delta,mx,inda-1))*exp(age[inda]*mx[inda]) ) * Int2(age[inda],age1[inda],mx[inda]) 
    for(j in (inda+1):m){
      term <- ( exp(-Cumhaz(delta,mx,j-1))*exp(age[j]*mx[j]) ) * Int2(age[j],age1[j],mx[j]) 
      if (is.na(term)) term <- 0 
      res <- res + term
    }
    out[i] <- res/STa
  }
  return(out)
}


varTCond <- function(age,mx){
  return(ET2Cond(age,mx) - ETCond(age,mx)^2)
}

CovTLambdaTCond <- function(age,mx){
  return(ETLambdaTCond(age,mx) - ETCond(age,mx)*1)
}

CorrTLambdaTCond <- function(age,mx){
  return(CovTLambdaTCond(age,mx)/sqrt(varTCond(age,mx)))
}


Int1 <- function(a,b,mx){
  if (mx==0 | is.na(mx)){
    res <- 0
  }else{
    res <- exp(-mx*a) * (a+1/mx) - exp(-mx*b) * (b+1/mx)
  }
  return(res)
}

Int2 <- function(a,b,mx){
  if (mx==0 | is.na(mx)){
    res <- 0
  }else{
    res <- -b^2*exp(-mx*b) + a^2*exp(-mx*a) +
      (2/mx) * (exp(-mx*a)*(a+1/mx) - exp(-mx*b)*(b+1/mx))
  }
  return(res) 
} 

Cumhaz <- function(delta,mx,r){
  if(r == 0) {res <- 0}
  if(r > 0) {res <- sum(delta[1:r]*mx[1:r])}
  #print(paste("Cumhaz =",res))
  return(res)
}