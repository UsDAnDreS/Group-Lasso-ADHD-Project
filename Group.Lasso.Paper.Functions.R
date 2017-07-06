rm(list=ls())

########################
### Needed libraries ###
########################

library(glmnet)
library(matrixcalc)
library(grpreg)
library(boot)

##########################################
## Creates block-diagonal matrix       ###
## with A.list[[i]] as its i-th block  ###
##########################################

block.diag <- function(A.list){
  p <- nrow(A.list[[1]])
  m <- ncol(A.list[[1]])
  if (!is.list(A.list)) return(A.list)
  
  K <- length(A.list)
  
  block_matrix <- matrix(0,K*p,K*m)
  for (i in 1:K){
    block_matrix[1:p + (i-1)*p,1:m + (i-1)*m] <- A.list[[i]]
  }  
  return(block_matrix)
}




#####################
### GENERATES TRANSITION MATRIX FOR THE TIME SERIES
#####################

gen_A = function(  # returns a list of d matrices A_1, ..., A_d
  p    # generate p x p VAR
  ,D=1   # generate VAR(d)
  ,max_eig   # spectral norm of A
  ,edge_density = 0.1 # if different for different lags, provide a vector of length d
  ,nonzero_diag = 1  # ensure all diagonal entries are non-zero; if different for different lags, provide a d-vector of 0 and 1
  ,signs = c(1,-1)  # whether we allow for both + and - elements (c(1,-1)), or just + (1), or just - (-1)
  ,pos.diag = T     # whether all diagonal elements should be positive
  ,stationary = 1  # ensure stationarity
  ,network.family = "random"  # A matrix filled randomly or with specific structure
  ,structure = NULL  # a pre-determined pattern of non-zero elements in generated matrix
){
  
  A = list()
  if (is.null(structure)){
  for (lag in 1:D){
    e_d = ceiling(p^2 * ifelse(length(edge_density) == 1, edge_density, edge_density[lag]))
    diag.ind <- seq(1,p^2,by=p+1)
    temp = sample(setdiff(1:p^2,diag.ind), e_d)
    temp_v = rep(0, p^2); temp_v[temp] = sample(signs,1)
    A[[lag]] = array(temp_v, c(p,p))
    if (nonzero_diag == 1)
      diag(A[[lag]]) = ifelse(pos.diag == T, 1, sample(signs,p,replace=T))
  }
  }
  if(!is.null(structure)){
    for (lag in 1:D){
      A[[lag]]=structure
    }
  }

  ######
  ## Guarantee stationarity of resulting VAR process
  ######
  
  if (stationary){
    A.big = array(0, c(p*D, p*D))
   
    for (i in 1:D)
      A.big[1:p, ((i-1)*p+1):(i*p)] = max_eig*A[[i]]
    
    if (D > 1)
      diag(A.big[(p+1):(p*D), 1:((D-1)*p)]) = 1
    
    temp = max(abs(eigen(A.big)$values))
    count = 0
    while (temp > max_eig){
      count = count+1
      #  print(paste("count:", count, "max_eigen_value:", round(temp, 2), "signal:", round(max(abs(A.big[1:p,])), 2)))
      A.big[1:p,] = A.big[1:p,]*0.95
      temp = max(abs(eigen(A.big)$values))
    }
    
    for (i in 1:D){
      A[[i]] = A.big[1:p, ((i-1)*p+1):(i*p)]
      print(paste("signal reduced to", round(max(abs(A.big[1:p,])), 2)))
    }
  }
  return(list(A=A,
              Signal=round(max(abs(A.big[1:p,])), 2))
  )
}

################
##### Function generating a group of related VAR transition matrices
################


A.setup <- function( # returns K transition matrices
  p,          # number of variables
  D=1,        # order of VAR
  ed,    #  edge density - proportion of non-zero off-diagonal elements in the common matrix
  signs = c(1,-1),  # whether we allow for both + and - elements (c(1,-1)), or just + (1), or just - (-1)  
  pos.diag=pos.diag, # whether diagonal should be positive
  comm,  #  for the case of different matrices - proportion of off-diagonal elements that are not part of common matrix
  max_eig_comm,  #  maximum eigenvalue of transition matrices
  max_eig_ind,  #  maximum eigenvalue of transition matrices
  min_elem, # minimum element value of transition matrices
  K,           # number of entities
  structure = NULL  # pre-determined pattern of non-zeros of common matrix
){
  
  ### Main cases of generated matrices:
  ###    - if comm=0: we generate one single matrix and copy it for each subject;
  ###    - if comm>0: we generate a common matrix and then randomly add non-zero elements to each transition matrix
  ###                 (that way there will be some common off-diagonal zeros, plus some individual effects as well)
  
  min_val <- min_elem
  
    if (comm == 0){
      
      A.true <- list()
      A.ind <- list()
      
      repeat{
        A.obj <- gen_A(p,D=D,ed=ed,signs=signs,max_eig=max_eig_comm,structure=structure)
        Signal <- A.obj$Signal
        A.true <- A.obj$A
        if (Signal > min_val) break;
      } 
      
      A.ind[[1]] <- make.list(matrix(0,p,p),D)
      
      if (K>1){
        for (i in 2:K){
          A.true[[i]] <- A.true
          A.ind[[i]] <- A.ind[[1]]
        }
      }
      A.comm <- A.true
      
      return(list(A.true=A.true,
                  A.comm=A.comm,
                  A.ind=A.ind))
    }
    
    if (comm > 0){ 
      repeat{
        A.Gen <- gen_A(p,D=D,edge_density=ed,signs=signs,max_eig=max_eig_comm,structure=structure)
        if (A.Gen$Signal > min_val) break;
      }
      
      A.true <- list()
      A.ind <- make.list(list(),K)

      for(i in 1:K){
        iter <- 0
        repeat{
          iter <- iter+1
          A11.Gen <- gen_A(p,D=D,edge_density=comm,signs=signs,max_eig=max_eig_ind,nonzero_diag = 0)
          
          A.11.Obj <- list()
          A.big <- array(0, c(p*D, p*D))
          
          for (d in 1:D){
            A.11.Obj[[d]] <- ifelse(A.Gen$A[[d]] != 0,A.Gen$A[[d]],A11.Gen$A[[d]])
            A.big[1:p, ((d-1)*p+1):(d*p)] <- A.11.Obj[[d]]
          }
          
            Signal=round(max(abs(A.big[1:p,])), 2)
            if ((Signal >= min_val) & (Signal <= max_eig)) break;
         }
        A.true[[i]] <- A.11.Obj
        for (d in 1:D) A.ind[[i]][[d]] <- A.11.Obj[[d]] - A.Gen$A[[d]]
      }
    
    return(list(A.true=A.true,
                A.comm=A.Gen$A,
                A.ind=A.ind))
    }
}


##################
### GENERATES DATA WITH PARTICULAR TRANSITION MATRIX AND ERROR MATRIX
###################

require(MASS)
gen_dat = function(  	# returns a p x T matrix of observations {X^1, ..., X^T}
  T = NULL    # number of observed time points
  ,error_sd = NULL # optional, a p-dim vector of the individual sdevs of the p time series
  ,Sigma_error = NULL # input a p x p correlation matrix; otherwise it is set to diag(error_sd) * identity
  ,SNR	# signal-to-noise ratio, used to determine error_sd (= abs(A[i,j])/SNR)
  ,A = NULL	# a list of D adjacency matrices, as returned by gen_A
  ,cut_T = 500*length(A) # time to wait before process reaches stability 
){
  d = length(A)
  p = dim(A[[1]])[1]
  X = array(0, c(p, T+cut_T))
  if(is.null(error_sd))
    error_sd = rep(max(abs(A[[1]]))/SNR, p)
  else if (length(error_sd) == 1)
    error_sd = rep(error_sd, p)
  
  if (is.null(Sigma_error))
    Sigma_error = diag(p)
  Sigma_error = diag(error_sd) %*% Sigma_error %*% diag(error_sd); #print(round(Sigma_error, 4))
  
  X = t(mvrnorm(n = T+cut_T, mu = rep(0, p), Sigma = Sigma_error))
  for (tt in (d+1):(T+cut_T))
    for (lg in 1:d)
      X[,tt] = X[,tt]+A[[lg]] %*% X[,tt-lg]
  
  return(X[,-seq(cut_T)])
}





##############
#### Function calculates all the performance measurements of the estimates:
####             FP,FN,TP,TN,Frobenius difference
##############

Measures.Vec <- function(# returns performance measurements of estimates
  A.est,   # estimates
  A.true   # true matrix
){ 
  if (sum(!A.true) == 0) FP=0
  if (sum(!A.true) != 0)  FP <- sum(!!A.est  & !A.true)/(sum(!A.true))
  if (sum(!!A.true) == 0) FN=0
  if (sum(!!A.true) != 0) FN <- sum(!A.est  & !!A.true)/(sum(!!A.true))
  TN <- 1 - FP
  TP <- 1 - FN
  Frob <- norm(A.est - A.true,type="F")/norm(A.true,type="F")
  return(list(FP=FP,
              FN=FN,
              TP=TP,
              TN=TN,
              Frob=Frob))
}


#######################################
#### Hard-thresholding function #######
#######################################

sparsify <- function(m,a){   
  m1 <- ifelse(abs(m)<a,0,m)
  return(m1)  
}



####################
#### Converts a vectorized matrix (which was stretched into a 1-dimensional vector)
#### back into its original form
####################


ConvertToMatrix.Full <- function(# returns a list of matrices(one matrix per entity)
  vec,        #vectorized matrix
  p           #number of variables
){
  K <- length(vec)/p^2
  A <- list()
  for (i in 1:K){
    A[[i]] <- matrix(vec[(1+(i-1)*(p^2)):(i*p^2)],p,p,byrow=TRUE)
  }
  return(A)
}


##########################################################
## Function connListCalc: calculates a connection list
## for group lasso for p vars per subject 
## It matches corresponding elements of multiple transition matrices into groups
##########################################################

connListCalc <- function(p){  # returns the connection list
  
  connList <- vector("list",2*p^2)
  class(connList) <- "connListObj"
  for (i in 1:(p^2)) connList[[i]] <- as.integer(i+p^2-1)
  for (i in (p^2+1):(2*(p^2))) connList[[i]] <- as.integer(i- p^2-1)
  
  return(connList=connList)
}




#################################################################################
#################################################################################
### IN THAT SECTION I HAVE ALL THE POSSIBLE CRITERIONS                     ######
### FOR PICKING SPARSITY AND FUSION PARAMETERS                             ######
###                                                                        ######
### The parameters for all of them will be:                                ######
###                                                                        ######
###  Est - set of estimates that we calculate the criterion for,           ######
###  X - data matrix,                                                      ######
###  Y - response vector,                                                  ######
###  lambda.path - vector of sparsity parameter values                     ######
###                that estimates from 'Est' correspond to,                ######
###  df.coef - the coefficient for degrees of freedom to be multiplied by  ######
###                                                                        ######
### The functions will return:                                             ######  
###                                                                        ######
###  the estimate that minimizes the criterion,                            ######
###  the corresponding minimum value of the criterion,                     ######
###  set of criterion values for all possible sparsity parameter values,   ######
###  sparsity parameter value corresponding to minimum criterion value,    ######
###  index of that sparsity parameter value in lambda.path,                ######
###  log-likelihood part of the criterion(for the whole lambda path)       ######
###  degrees of freedom of the criterion(for the whole lambda path)        ######
#################################################################################
#################################################################################



################################################
### AIC criterion ##############################
###############################################

AIC <- function(Est,X,Y,lambda.path,df.path,p,K=1,df.coef=2){
  
  AIC <- rep(0,length(lambda.path))
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  n <- nrow(Y)
  t <- n/p + 1
  
  for(i in 1:length(lambda.path)){
    df <- df.path[i]
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- df.coef*df
    AIC[i] <- loglik.part[i] + df.part[i]  
  }
  
  min <- which.min(AIC)
  
  return(list(Est=Est[,min],
              Criter.min=AIC[min],
              Criter=AIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}


################################################
### AICc criterion ##############################
###############################################

AICc <- function(Est,X,Y,p,K=1,lambda.path,df.path){
  
  AIC <- rep(0,length(lambda.path))
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  n <- nrow(Y)
  t <- n/p + 1
  
  for(i in 1:length(lambda.path)){
    df <- df.path[i]
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- 2*df + 2*df*(df+1)/(n-df-1)
    AIC[i] <- loglik.part[i] + df.part[i]  
  }
  
  min <- which.min(AIC)
  
  return(list(Est=Est[,min],
              Criter.min=AIC[min],
              Criter=AIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}


################################################
### BIC criterion ##############################
###############################################

BIC <- function(Est,X,Y,K,p,df.path,lambda.path){
  
  BIC <- rep(0,ncol(Est))
  n <- nrow(Y)
  t <- n/K + 1
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  l <- length(df.path)
  
  for(i in 1:l){ 
    df <- df.path[i]
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- df*log(n)
    BIC[i] <- loglik.part[i] + df*log(n) 
  }
  
  min <- which.min(BIC)
  
  return(list(Est=Est[,min],
              Criter.min=BIC[min],
              Criter=BIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}


##########################################
#### Calculating Matthews coefficient based on TP,FP,TN,FN rates
#### of the estimate
##########################################

Matthews.Coef <- function(TP,FP,TN,FN){
  return(ifelse(TP*TN - FP*FN == 0,0,(TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
}


##################
## MAKES A LIST OF n ELEMENTS
## Elem - type of the elements(matrices,lists, whatever it maybe)
##################

make.list <- function(elem,n){
  res.list <- list()
  for(i in 1:n){
    res.list[[i]] <- elem
  }
  return(res.list)
}


###########################
### Compact Matrix Form Setup function (as in Problem Setup section of the paper)
### Transforms the original K p-dimensional time series of length n.tp into matrices C,X,B
###########################

mat.setup <- function(ts,n.tp,K,p,D=1){
  C.mat <- matrix(0,n.tp-D,K*p)
  for (i in 1:(n.tp-D)){
    for (k in 1:K){
      C.mat[i,((k-1)*p + 1):(k*p)] <- c(ts[((k-1)*p + 1):(k*p),(i+D)])
    }
  }
  
  B.mat <- matrix(0,n.tp-D,K*(D*p))
  for (i in 1:(n.tp-D)){
    for (k in 1:K){
       res <- ts[((k-1)*p + 1):(k*p),i]
       if (D>1){
          for (d in 2:D){
            res <- c(ts[((k-1)*p + 1):(k*p),i+(d-1)],res)
          }
       }
      B.mat[i,((k-1)*(D*p) + 1):(k*(D*p))] <- res
    }
  }
  
  B.list.mat <- list()
  for(k in 1:K){
    B.list.mat[[k]] <- B.mat[,((k-1)*(D*p) + 1):(k*(D*p))]
  }
  
  X.list.mat <- block.diag(B.list.mat)
  
  return(list(C=C.mat,
              B=B.mat,
              X=X.list.mat))
}

########################################
#### OLS estimation for a time series ##
########################################


OLS.tseries <- function(# returns the OLS estimate for transition matrices,
                              #         maximum likelihood estimate of variance
  Data, # time series matrix, rows - variables, columns - time points
  D=1   # order of VAR model for the OLS fit
){   
  p <- nrow(Data)
  t <- ncol(Data)
  
  ### Setting up all the matrices for regression problem
  
  ## Response matrix
  Response <- matrix(c(t(Data[,-c(1:D)])))
  n <- length(Response)
  
  ## Data matrix
  for (d in 1:D){
     if (d==1) X <- matrix(Data[,c((t-d):(D-d+1))],t-D,p,byrow=TRUE)
     if (d>1) X <- cbind(X,matrix(Data[,c((t-d):(D-d+1))],t-D,p,byrow=TRUE))
  }
  X.full <- block.diag(make.list(X,p))
  
  ## Calculating OLS estimate for transition matrices and computing MLE estimate for sigma^2 from residuals
  A.OLS <- solve(t(X.full)%*%X.full)%*%(t(X.full)%*%Response)
  sigma2 <- sum(Response - X.full%*%A.OLS)^2/n
  
  return(list(A.OLS = A.OLS,
              sigma2=sigma2))
}
