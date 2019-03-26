### Executive file for simulation studies 

## Where all the customized functions (for estimation, performance metrics calculation, etc) are contained.
source("MARCH_2019_SUBMISSION_Group_Lasso_Paper_Functions.R")

set.seed(2)

save.results <- TRUE    # Do we save R objects with true/estimated matrices?

### Method you would like to apply.
method.name <- c("Sparse_Lasso",
                 "GroupBridge",
                 "Sparse_Group_Lasso")[3]

### Sparse Group Lasso (SGL) parameters.

alpha <- 0.2    # 0 - sheer group lasso; 1 - sheer sparse lasso.
nlam <- 20      # Grid size for lambda parameter of SGL.
nfolds <- 5     # No. of folds for K-fold CV


p <- 10                                     # number of variables per subject
print(c("p",p))
t.set <- 30                              # number of observed time points per subject
print(c("t.set:",t.set))
K.set <- 5                              # number of subjects per group
print(c("K.set:",K.set))
heter <- c("low","moderate","none")[2]      # level of heterogeneity: low (1%), moderate (2-3%), none (0%)
print(c("heter:",heter))
Thresh2 <- 0.02                             # Threshold to be applied to estimates, in order to get rid of noise
print(c("Thresh2:",Thresh2))
signs <- c(1,-1)                            # Whether effects will be allowed to be positive(1), negative(-1) or both
print(c("Signs:",signs))
SNR <- 1                                    # Signal-to-Noise ratio
print(c("SNR:",SNR))
spectral <- c(0.25,0.4,0.6)[2]              # Maximum eigenvalue allowed in generated VAR transition matrices (either 0.25, 0.4 or 0.6)
print(c("spectral:",spectral))
pos.diag <- T                               # whether all diagonal elements of VAR transition matrices should be positive
print(c("pos.diag=",pos.diag))
n.groups <- 1                               # number of subject groups generated
print(c("n.groups:",n.groups))


rep <- 4                                    # number of simulations to run
print(c("rep",rep))
actual.rep <- 4                            # In case we only need "actual.rep" last replicates
print(c("actual.rep",actual.rep))           # (appropriate when running all "rep" replicates takes too long)


#####
## CRITERIONS AND OTHER PARAMETERS
## Criterions can only be: BIC, AIC and AIC corrected
#####

criter <- "BIC"                                   # selection criterion for final estimate
print(c("criter",criter))

D <- 1                                           # lag order of generated VAR models (this code has been mostly tested for the case of VAR(1) models)
print(c("order:",D))


######
## Generating different settings for different values of p:
## if p <= 20:  5% common off-diagonal edge density, 1% individual density for low heterogeneity, 3% - for moderate.
## if p > 20: 2% common off-diagonal edge density, 1% individual density for low heterogeneity, 2% - for moderate.
######

if (p<=20){
  ed <- 0.05
  print(c("ed:",ed))
  if (heter == "moderate"){
    comm <- 0.03; 
    print(c("comm:",comm))
  }
}
if (p>20){
  ed <- 0.02
  print(c("ed:",ed))
  if (heter == "moderate"){
    comm <- 0.02;
    print(c("comm:",comm))
  }
}
if (heter == "low"){
  comm <- 0.01;
  print(c("comm:",comm))
}

if (heter == "none"){
  comm <- 0;
  print(c("comm:",comm))
}


###############
## Setting the spectral radius for both common and individual component, and the minimum value of non-zero element
###############

if (spectral == 0.6){
  max_eig_comm <- 0.6
  print(c("max_eig_comm:",max_eig_comm))
  max_eig_ind <- 0.6
  print(c("max_eig_ind:",max_eig_ind))
  max_eig <- max(max_eig_comm,max_eig_ind)
  min_elem <- 0.2
  print(c("min_elem:",min_elem))
}

if (spectral == 0.4){
  max_eig_comm <- 0.4
  print(c("max_eig_comm:",max_eig_comm))
  max_eig_ind <- 0.4
  print(c("max_eig_ind:",max_eig_ind))
  max_eig <- max(max_eig_comm,max_eig_ind)
  min_elem <- 0.2
  print(c("min_elem:",min_elem))
}

if (spectral == 0.25){
  max_eig_comm <- 0.25
  print(c("max_eig_comm:",sort(max_eig_comm)))
  max_eig_ind <- 0.25
  print(c("max_eig_ind:",sort(max_eig_ind)))
  max_eig <- max(max_eig_comm,max_eig_ind)
  max_eig.diff <- 0.15
  print(c("max_eig.diff:",max_eig.diff))
  min_elem <- max_eig - max_eig.diff
  print(c("min_elem:",sort(min_elem)))
}


#####
## Parameters for glmnet function: no-intercept model
#####

intercept <- FALSE

## Looping through different numbers K of subjects per group.
for(K in K.set){
  
  ### For VAR of order D we will store the estimates for a single subject in the following (D*p)xp matrix:
  ### A_{(D*p)xp} = [A1_pxp,]
  ###                ...,
  ###               [AD_pxp ]
  
  Final.Est <- make.list(make.list(matrix(0,D*p,p),K),actual.rep)
  A.true.all <- make.list(make.list(matrix(0,D*p,p),K),actual.rep)
  A.comm.all <- make.list(matrix(0,D*p,p),actual.rep)
  A.ind.all <- make.list(make.list(matrix(0,D*p,p),K),actual.rep)
  
  ## Looping through various lengths T of time series.
  for(train in t.set){
    
    t <- train
    ptm <- proc.time()
    
    ######
    ### Initializing vectors of metrics (see Section 3)
    ######
    
    ## Metrics for full estimate A^hat (= A^comm + A^ind)
    FP.Rand.pen <- rep(0,actual.rep)
    FN.Rand.pen <- rep(0,actual.rep)
    TP.Rand.pen <- rep(0,actual.rep)
    TN.Rand.pen <- rep(0,actual.rep)
    Matt.Coef.Rand.pen <- rep(0,actual.rep)
    Frob.Rand.pen <- rep(0,actual.rep)
    
    
    ## Setting up/Creating the folder to save the results into.
    if (save.results){
      if (method.name != "Sparse_Group_Lasso"){
        namedir <- paste(getwd(),"/", method.name, "_Simul_Results_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_Rep=",rep,"_ActualRep=",actual.rep,"_SNR=",SNR,"_criter=",criter,"_heter=",heter,"_signs=",length(signs),"_max_eig_comm=",max_eig_comm,"_max_eig_ind=",max_eig_ind,"_min_elem=",min_elem,sep="")
      }
      if (method.name == "Sparse_Group_Lasso"){
        namedir <- paste(getwd(),"/", method.name, "_Simul_Results_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_Rep=",rep,"_ActualRep=",actual.rep,"_SNR=",SNR,"_criter=CV_alpha=",alpha,"_heter=",heter,"_signs=",length(signs),"_max_eig_comm=",max_eig_comm,"_max_eig_ind=",max_eig_ind,"_min_elem=",min_elem,sep="")
      }
      dir.create(namedir) 
    }
    
    
    ## If we want to generate "rep" replicates, but only do the estimation for last "actual.rep":
    ##    1. We just generate (rep - actual.rep) replicates first.
    if (actual.rep != rep){
      
      for (run in 1:(rep-actual.rep)){
        print(c("run",run))
        
        repeat{
          #####
          ## Generating VAR transition matrices
          #####
          
          A.true <- A.setup(p,
                            D=D,
                            ed=ed,  
                            signs=signs,
                            pos.diag=pos.diag,
                            comm=comm,  
                            max_eig_comm=max_eig_comm, 
                            max_eig_ind=max_eig_ind,
                            min_elem = min_elem,
                            K=K)
          
          A.list <- list()
          
          for (d in 1:D){
            A.full <- list()
            for (k in 1:K) A.full[[k]] <- A.true$A.true[[k]][[d]]
            A.list[[d]] <- block.diag(A.full)
          }
          
          Sigma <- list()
          for(i in 1:K) Sigma[[i]] <- diag(1,p)
          
          Sigma.full <- block.diag(Sigma)
          
          ################################
          ####### DATA GENERATION ########
          ################################
          
          DATA <- gen_dat(T=t,
                          A=A.list,
                          SNR=SNR,
                          Sigma_error=Sigma.full)
          
          print(max(abs(DATA)))
          #  hist(DATA)
          if (max(abs(DATA))<10) break;                  ### making sure generated time series doesn't go overboard with magnitudes of matrix values (happens occasionally)
        }
      }
    }
    
    ##    2. Then, for the remaining "actual.rep" replicates, we generate the data AND conduct the selected estimation method.
    for(run in 1:actual.rep){
      print(c("run",run))
      
      repeat{
        #####
        ## Generating VAR transition matrices
        #####
        
        A.true <- A.setup(p,
                          D=D,
                          ed=ed,  
                          signs=signs,
                          pos.diag=pos.diag,
                          comm=comm,  
                          max_eig_comm=max_eig_comm, 
                          max_eig_ind=max_eig_ind,
                          min_elem = min_elem,
                          K=K)
        A.true.all[[run]] <- A.true$A.true
        A.comm.all[[run]] <- A.true$A.comm
        A.ind.all[[run]] <- A.true$A.ind
        
        A.list <- list()
        
        for (d in 1:D){
          A.full <- list()
          for (k in 1:K) A.full[[k]] <- A.true$A.true[[k]][[d]]
          A.list[[d]] <- block.diag(A.full)
        }
        
        Sigma <- list()
        for(i in 1:K){
          #Sigma[[i]] <- diag(runif(p,0,1))
          Sigma[[i]] <- diag(1,p)
        }
        
        Sigma.full <- block.diag(Sigma)
        
        ################################
        ####### DATA GENERATION ########
        ################################
        
        DATA <- gen_dat(T=t,
                        A=A.list,
                        SNR=SNR,
                        Sigma_error=Sigma.full)
        
        print(max(abs(DATA)))
        #  hist(DATA)
        if (max(abs(DATA))<10) break;                  ### making sure generated time series doesn't go overboard with magnitudes of matrix values (happens occasionally)
      }
      
      
      ### Calculating maximum likelihood estimates of sigma^2 for each subject
      # print("Calculating estimates for sigma_1^2,...,sigma_K^2")
      sigma2 <- rep(0,K)
      sds <- apply(DATA,1,function(x) sd(x))
      
      for (k in 1:K) sigma2[k] <- OLS.tseries(DATA[(k-1)*p + 1:p,],D=D)$sigma2
      
      if (method.name == "Sparse_Lasso"){
        for (k in 1:K){
          # print(k)
          sigma2[k] <- OLS.tseries(DATA[(k-1)*p + 1:p,],D=D)$sigma2
          
          est <- Sparse.tseries(DATA[(k-1)*p + 1:p,],
                                sigma2=sigma2[k],
                                criter=criter)$beta.hat
          
          Final.Est[[run]][[k]] <- matrix(sparsify(est, Thresh2),
                                          byrow=T,
                                          nrow=p)
        }
      }
      
      
      if (method.name != "Sparse_Lasso"){
      # print("Generated")
      
      ### Initializing vectors for final estimates
      
      Group.Est <- make.list(matrix(0,D*p,p),K)
      Sep.Est.Second <-  make.list(matrix(0,D*p,p),K)
      Group.Final <- make.list(matrix(0,D*p,p),K)
      
      A.list <- list()
      
      for (d in 1:D){
        A.full <- list()
        for (k in 1:K) A.full[[k]] <- A.true$A.true[[k]][[d]]
        A.list[[d]] <- block.diag(A.full)
      }
      
      Sigma <- list()
      for(i in 1:K) Sigma[[i]] <- diag(1,p)               # simply setting the error covariance to be identity matrix
      
      Sigma.full <- block.diag(Sigma)
      
      #############################
      ### FULL PROBLEM SETUP   ####
      #############################
      
      M.setup <- mat.setup(DATA,t,K,p,D=D) 
      C.list <- M.setup$C
      X.list <- M.setup$X
      
      ### Vector of group number assignments for group lasso
      group <- c(1:(D*p))
      if (K>1){
        for (j in 2:K) group <- c(group,1:(D*p))
      }
      
      ## Initializing vectors to contain estimates during algorithm iterations
      method.result <- make.list(numeric(D*K*p),p)
      gl.zeros <- list()
      
      for (j in 1:p){
        print(paste("j:",j,sep=""))
        
        ### Initializing response vector and data matrix for standard regression problems for the first stage
        Y <- numeric(K*(t-D))
        for (k in 1:K){
          Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p]
        }
        
        ### Doing group lasso optimization
        D.sigma <- sqrt(diag(c(sapply(sigma2, function(x) return(rep(x,(t-D)))))))
        
        
        ####
        if (method.name == "Sparse_Group_Lasso"){
          seed.stuff <- .Random.seed
          r.cv <- cvSGL(data = list(x=solve(D.sigma) %*% X.list[,order(group)],
                                    y=solve(D.sigma) %*% Y),
                        index=sort(group),
                        type="linear",
                        alpha=alpha,
                        nlam=nlam,
                        nfold=nfolds,
                        standardize = F)
          #verbose=T)
          .Random.seed <- seed.stuff
          
          r.cv$fit$beta
          #est <- r.cv$fit$beta[order(group),]
          est <- matrix(0, 
                        nrow=nrow(r.cv$fit$beta),
                        ncol=ncol(r.cv$fit$beta))
          
          est[order(group),] <- r.cv$fit$beta
          
          cv.err <- r.cv$lldiff
          method.est.out <- est[,which.min(cv.err)]
        }
        if (method.name == "GroupBridge"){
          r <- gBridge(solve(D.sigma) %*% X.list,
                       solve(D.sigma) %*% Y,
                       group=group,
                       family="gaussian")
          
          lambda_G.path <- r$lambda
          est <- r$beta[-1,]
          method.result.df <- r$df
          
          ### Tuning parameter selection
          if (criter == "AIC")
            method.est.out <- AIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,lambda.path=lambda_G.path,df.path=method.result.df)
          if (criter == "BIC")
            method.est.out <- BIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=method.result.df,lambda.path=lambda_G.path)
          if (criter == "AICc")
            method.est.out <- AICc(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=method.result.df,lambda.path=lambda_G.path)
          
          # method.result[[j]] <- sparsify(method.est.out$Est, Thresh2)
          lambda.group <- method.est.out$lambda1
          method.est.out <- method.est.out$Est
        }
        
        #method.result[[j]] <- sparsify(method.est.out$Est, Thresh2)
        method.result[[j]] <- sparsify(method.est.out, Thresh2)
        #lambda.group <- lambda_G.path[which.min(cv.err)]
        
        ### Recording estimates
        for (k in 1:K) Group.Est[[k]][j,] <- method.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)]
      }   
      
      
      for(k in 1:K){
        Group.Final[[k]] <- sparsify(Group.Est[[k]],Thresh2)
        Final.Est[[run]][[k]] <- Group.Final[[k]]
      }
      
    }
      
      
      #########################################################
      ############## CALCULATING PERFORMANCE METRICS   ########
      #########################################################
      
      #######
      #### METRICS FOR FULL ESTIMATE
      #######
      
      true.vec <- NULL
      Rand.pen.vec <- NULL
      for (k in 1:K){
        A.bind <- A.true$A.true[[k]][[1]]
        if (D>1) for (d in 2:D) A.bind <- rbind(A.bind,A.true$A.true[[k]][[d]])
        true.vec <- c(true.vec,A.bind)
        # Rand.pen.vec <- c(Rand.pen.vec,vec(Group.Final[[k]]))
        Rand.pen.vec <- c(Rand.pen.vec,vec(Final.Est[[run]][[k]]))
      }
      
      Measures.Joint <- Measures.Vec(Rand.pen.vec,as.matrix(true.vec))
      FP.Rand.pen[run] <- Measures.Joint$FP
      FN.Rand.pen[run] <- Measures.Joint$FN
      TP.Rand.pen[run] <- Measures.Joint$TP
      TN.Rand.pen[run] <- Measures.Joint$TN
      Matt.Coef.Rand.pen[run] <- Matthews.Coef(TP.Rand.pen[run],
                                               FP.Rand.pen[run],
                                               TN.Rand.pen[run],
                                               FN.Rand.pen[run])
      Frob.Rand.pen[run] <- Measures.Joint$Frob
    }
    
    ###############################################
    #### Saving all the estimates to .rds files ###
    ###############################################
    
    if (save.results){
      saveRDS(Final.Est,paste(namedir,"/Final_Est.rds",sep=""))
      saveRDS(A.true.all,paste(namedir,"/A_true.rds",sep=""))
    }
    
    time.taken <- proc.time() - ptm
    #print("Time taken:")
    #print(time.taken)
    
    ##################################################
    ### Printing all the metrics in a format for latex table:
    ### FP, FN, Matthews for full estimates
    ##################################################
    
    cat("\n")
    cat("\n")
    print(method.name)
    cat("\n")
    cat("\n")
    cat("p=",p," & ",t," & ",K," & ",
        round(mean(FP.Rand.pen),2),"(",round(sd(FP.Rand.pen),2),") & ",
        round(mean(FN.Rand.pen),2),"(",round(sd(FN.Rand.pen),2),") & ",
        round(mean(Matt.Coef.Rand.pen),2),"(",round(sd(Matt.Coef.Rand.pen),2),") & ",
        round(mean(Frob.Rand.pen),2),"(",round(sd(Frob.Rand.pen),2),") ",
        "\\\\",sep="")
    cat("\n")
    cat("\\hline")
    cat('\n') 
    
  }
}



