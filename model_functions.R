#' initiate all parameters for simulating the data 
#' 
#' @param T the number of time points including the baseline 
#' @param V the number of covariates
#' @param P the number of products (e-cig and c-cig). It can only take 2. 
#' 
#' @example Initiate(nT=5, V=4)
#' 
#' @return a vector of parmeters 
Initiate<-function(nT, V){
  P <- 2
  lambda <- rbind(rep(-2, nT), rep(-3, nT))
  mu <- rbind(rep(1, nT), rep(2, nT))
  alpha <- rbind(sample(c(0.1, 0.2, 0.3), V, replace = TRUE),
                 sample(c(0.1, 0.2, 0.3), V, replace = TRUE))
  eta <- rbind(rep(-0.1, V), 
               rep(-0.2, V))
  
  beta <- rep(2, P)
  gamma <- rep(-2, P)
  delta <- matrix(0, P, P)
  
  simstat <-  list(lambda=lambda,
                   mu=mu,
                   alpha=alpha,
                   eta=eta,
                   beta=beta,
                   gamma=gamma,
                   delta=delta)
  return(simstat)
}

#' Calculate the transitional probabilities 
#' 
#' @param param the values of parameters 
#' @param num the number of observations 
#' @param time the baseline time point in which the transitional probability is evaluated
#' @param Ymat a list storing Y
#' @param Xmat a list storing X
#' 
TransitionProbabilities<-function(param, num=n, time=t, Ymat, Xmat){
  ###initiate values
  for(index in names(param)){
    assign(index, param[[index]])
  }
  
  lambda[,1] <- lambda[,2]
  mu[,1] <- mu[,2]
  
  ######
  PY<-rep(NA,4) #J<- 4
  
  if(Ymat[[num]][time-1,1]==0 & Ymat[[num]][time-1,2]==0){
    PY[1]<-0
    PY[2]<-sum(lambda[2,time], alpha[2,]*Xmat[[num]][time-1,])
    PY[3]<-sum(lambda[1,time], alpha[1,]*Xmat[[num]][time-1,])
    PY[4]<-sum(lambda[,time], delta[2,2], t(alpha)*Xmat[[num]][time-1,])
    
  } else if(Ymat[[num]][time-1,1]==0 & Ymat[[num]][time-1,2]==1){
    PY[1]<-sum(mu[2,time], eta[2,]*Xmat[[num]][time-1,])
    PY[2]<-0
    PY[3]<-sum(lambda[1,time], beta[1], mu[2,time], delta[2,1], (alpha[1,]+eta[2,])*Xmat[[num]][time-1,])
    PY[4]<-sum(lambda[1,time], beta[1], alpha[1,]*Xmat[[num]][time-1,])
    
  } else if(Ymat[[num]][time-1,1]==1 & Ymat[[num]][time-1,2]==0){
    PY[1]<-sum(mu[1,time], eta[1,]*Xmat[[num]][time-1,])
    PY[2]<-sum(lambda[2,time], beta[2], mu[1,time], delta[1,2], (alpha[2,]+eta[1,])*Xmat[[num]][time-1,])
    PY[3]<-0
    PY[4]<-sum(lambda[2,time], beta[2], alpha[2,]*Xmat[[num]][time-1,])
    
  } else if(Ymat[[num]][time-1,1]==1 & Ymat[[num]][time-1,2]==1){
    PY[1]<-sum(mu[1,time], gamma[1], mu[2,time], gamma[2], delta[1,1], t(eta)*Xmat[[num]][time-1,])
    PY[2]<-sum(mu[1,time], gamma[1], eta[1,]*Xmat[[num]][time-1,])
    PY[3]<-sum(mu[2,time], gamma[2], eta[2,]*Xmat[[num]][time-1,])
    PY[4]<-0
  }
  PY=exp(PY)/sum(exp(PY))
  return(PY)
}

#' Simulating N number of data given paramters
#' @param param a list storing all parameters
#' @param N the number of observations 
#' 
#' @return a list including X and Y
#' 
SimulateData<-function(param, N=N){
  #pass the values to parameters
  for(index in names(param)){
    assign(index, param[[index]])
  }
  
  lambda[,1] <- lambda[,2]
  mu[,1] <- mu[,2]
  
  #initiate the X and Y vectors
  nT <- ncol(lambda)
  V <- ncol(alpha)
  P <- 2
  
  Y <- rep(list(matrix(NA, nT, P)), N) #N=2000, T=4, P=2(product)
  X <- rep(list(matrix(NA, nT, V)), N)
  
  covX <- matrix(NA, nrow = V, ncol = V)
  diag(covX) <- 1
  for(i in 1:(ncol(covX)-1)){
    covX[row(covX) == (col(covX) + i)] <- 0.5^i
    covX[row(covX) == (col(covX) - i)] <- 0.5^i
  }
  
  ZERO <- rep(0, V)
  J <- 4
  
  for(n in 1:N){
    for(t in 1:nT){
      X[[n]][t,] <- MASS::mvrnorm(1, ZERO, covX)
      
      if(t>1){
        X[[n]][t,] <- 1/4 * X[[n]][t,] + 3/4 * X[[n]][t-1,]
      } 
      
      #generate the wave 1 
      if(t==1){
        PY <- rep(NA, J)
        
        PY[1]<-0
        PY[2]<-sum(lambda[2,1], alpha[2,]*X[[n]][t,])
        PY[3]<-sum(lambda[1,1], alpha[1,]*X[[n]][t,])
        PY[4]<-sum(lambda[ ,1], delta[1,1], t(alpha)*X[[n]][t,])
        
        PY<-exp(PY)/sum(exp(PY))
      }else{
        PY<-TransitionProbabilities(param=simstat, num=n, time=t, Ymat=Y, Xmat=X)
      }
      
      #generate the outcome
      
      j <- 1
      rand <- runif(1)
      SumPY <- PY[1]
      
      while(SumPY < rand){
        j <- j+1
        SumPY <- SumPY + PY[j]
      }
      
      Y[[n]][t,1] <- (j-1)%/%P
      Y[[n]][t,2] <- (j-1)%%P
    }
  }
  
  return(list(X = X, Y = Y))
}

#' Calulate the likelihood for a given data and a set of parameters
#' 
#' @param param a list of parameters
#' @param data a list of X and Y
#' 
#' @return a numeric value of likelihood

LogLikelihood<-function(param, data){
  P <- 2
  #likelihood after the baseline (t>1)
  X <- data[[1]]
  Y <- data[[2]]
  
  loglike <- 0
  
  namesparam <- sub("[^[:alpha:]]+", "", names(param))
  for(index in unique(namesparam)){
    if(index %in% c("beta", "gamma")){
      tmp <- param[grep(index, names(param))]
      assign(index, tmp)
    } else {
      tmp <- matrix(param[grep(index, names(param))], nrow = 2)
      assign(index, tmp)
    }
  }
  
  #calcualte the likelihood
  N <- length(X)
  for(n in 1:N){
    tempY <- Y[[n]]
    tempX <- X[[n]]
    
    for(t in 2:nT){
      logProb <- c(0,0,0,0)
      
      Y0 <- P*tempY[t-1,1] + tempY[t-1, 2] 
      Y1 <- P*tempY[t  ,1] + tempY[t  , 2] 
      
      
      if(Y0 == 0){
        logProb[2] = sum(lambda[2,t], alpha[2,]*tempX[t-1,])
        logProb[3] = sum(lambda[1,t], alpha[1,]*tempX[t-1,])
        logProb[4] = sum(logProb[2:3], delta[2,2])
        
      } else if(Y0 == 1){
        logProb[1] = sum(mu[2,t], eta[2,]*tempX[t-1,])
        logProb[4] = sum(lambda[1,t], beta[1], alpha[1,]*tempX[t-1,])
        logProb[3] = sum(logProb[c(1,4)], delta[2,1])
        
      } else if(Y0 == 2){
        logProb[1] = sum(mu[1,t], eta[1,]*tempX[t-1,])
        logProb[4] = sum(lambda[2,t], beta[2], alpha[2,]*tempX[t-1,])
        logProb[2] = sum(logProb[c(1,4)], delta[1,2])

      } else if(Y0 == 3){
        logProb[2] = sum(mu[1,t], gamma[1], eta[1,]*tempX[t-1,])
        logProb[3] = sum(mu[2,t], gamma[2], eta[2,]*tempX[t-1,])
        logProb[1] = sum(logProb[2:3], delta[1,1])
      }
      
      loglike = loglike + logProb[Y1+1] - log(sum(exp(logProb))) 
      
    }
  }
  
  return(loglike)
}

#' analyze the data with initial parameters, X, and Y
#' 
#' @param simstat initial values for parameters 
#' @param X a list of N observartions with V covariates 
#' @param Y a list of N observations with P outcomes 
#' 
#' @return 

AnalyzeNew<-function(simstat, X, Y){
  if(length(X)!=length(Y)){
    stop("The lengths of X and Y are not equal.")
  }
  P <- 2
  O <- 2
  J <- 4
  NY <- rep(list(matrix(0, O, O)), P)
  NY0 <- matrix(0, P, O)
  N <- length(X)
  
  #print(sprintf("REPLICATE %1$1.0f", r))
  
  for(p in 1:P){
    for(n in 1:N){
      # numbers of e/cig users at wave 1
      NY0[p, Y[[n]][1, p]+1] <- NY0[p, Y[[n]][1, p]+1] + 1
      for(t in 2:T){
        NY[[p]][Y[[n]][t-1, p]+1, Y[[n]][t, p]+1] <- NY[[p]][Y[[n]][t-1, p]+1, Y[[n]][t, p]+1] + 1   
      }
    }
    # numbers of status trasitions, 0>0, 0>1, 1>0, 1>1, for two products from three waves (i.e., 6000) 
    print(sprintf("p = %1$1.0f %2$1.0f %3$1.0f %4$1.0f %5$1.0f %6$1.0f", p, NY0[p, 2], NY[[p]][1, 1], NY[[p]][1, 2], NY[[p]][2, 1], NY[[p]][2, 2]))
  }
  
  NYY <- matrix(0, J, J)
  for(n in 1:N){
    for(t in 2:T){
      J0 <- P*Y[[n]][t-1, 1] + Y[[n]][t-1, 2]
      J1 <- P*Y[[n]][t, 1] + Y[[n]][t, 2]
      NYY[J0+1,J1+1] <- NYY[J0+1,J1+1]+1
    }
  }
  
  initiation <- matrix(0, P, P)
  cessation <- matrix(0, P, P)
  
  for(J0 in 1:J){
    t <- 2
    n <- 1
    Y[[n]][t-1, 1] <- (J0-1) %/% P
    Y[[n]][t-1, 2] <- (J0-1) %% P
    
    totNYY <- sum(NYY[J0,])
    for(J1 in 1:J){
      print(sprintf("%1$1.0f, %2$1.0f, %3$1.0f, %4$5.3f", J0, J1, NYY[J0,J1], NYY[J0,J1]/totNYY), quote = FALSE)
    }
    
    if(J0 == 1){
      OR <- NYY[J0, 4]*NYY[J0, 1]/(NYY[J0, 3]*NYY[J0, 2]); print(sprintf("OR=%1$5.3f, Delta[11]=%2$5.2f", OR, log(OR)))
      initiation[1, 1] <- NYY[J0, 3]/(NYY[J0, 1]+NYY[J0, 3]);
      initiation[2, 1] <- NYY[J0, 2]/(NYY[J0, 1]+NYY[J0, 2]);
      print(sprintf("initiation[1]=%1$6.4f, initiation[2]=%2$6.4f", initiation[1, 1], initiation[2, 1]))
    } else if(J0 == 2){
      OR <- NYY[J0, 2]*NYY[J0, 3]/(NYY[J0, 1]*NYY[J0, 4]); print(sprintf("OR=%1$5.3f, Delta[10]=%2$5.2f", OR, log(OR)))
      initiation[1, 2] <- NYY[J0, 4]/(NYY[J0, 2]+NYY[J0, 4]);
      cessation[2, 1] <- NYY[J0, 1]/(NYY[J0, 1]+NYY[J0, 2]);
      print(sprintf("initiation[0]=%1$6.4f, cessation[1]=%2$6.4f", initiation[1, 2], cessation[2, 1]));
      print(sprintf("OR=%1$4.2f", initiation[1, 2]*(1-initiation[1, 1])/(initiation[1, 1]*(1-initiation[1, 2]))));
    } else if(J0 == 3){
      OR <- NYY[J0, 2]*NYY[J0, 3]/(NYY[J0, 1]*NYY[J0, 4]); print(sprintf("OR=%1$5.3f, Delta[01]=%2$5.2f", OR, log(OR)))
      cessation[1, 1] <- NYY[J0, 1]/(NYY[J0, 1]+NYY[J0, 3]);
      initiation[2, 2] <- NYY[J0, 4]/(NYY[J0, 3]+NYY[J0, 4]);
      print(sprintf("initiation[0]=%1$6.4f, cessation[1]=%2$6.4f", initiation[2, 2], cessation[1, 1]))
      print(sprintf("OR=%1$4.2f", initiation[2, 2]*(1-initiation[2, 1])/(initiation[2, 1]*(1-initiation[2, 2]))));
    } else if(J0 == 4){
      OR <- NYY[J0, 4]*NYY[J0, 1]/(NYY[J0, 3]*NYY[J0, 2]); print(sprintf("OR=%1$5.3f, Delta[00]=%2$5.2f", OR, log(OR)))
      cessation[1, 2] <- NYY[J0, 2]/(NYY[J0, 2]+NYY[J0, 4]);
      cessation[2, 2] <- NYY[J0, 3]/(NYY[J0, 3]+NYY[J0, 4]);
      print(sprintf("cessation[0]=%1$6.4f, cessation[1]=%2$6.4f", cessation[1, 2], cessation[2, 2]))
      print(sprintf("OR=%1$4.2f", cessation[1, 2]*(1-cessation[1, 1])/(cessation[1, 1]*(1-cessation[1, 2]))));
      print(sprintf("OR=%1$4.2f", cessation[2, 2]*(1-cessation[2, 1])/(cessation[2, 1]*(1-cessation[2, 2]))));
    }
    
  }
  
  simstatvec <- unlist(simstat)
  
  ml <- optim(simstatvec, LogLikelihood, data=list(X, Y), method = "BFGS", hessian=TRUE,
              control = list(fnscale=-1, trace=TRUE, maxit=1000 ))  
  return(ml)
}

