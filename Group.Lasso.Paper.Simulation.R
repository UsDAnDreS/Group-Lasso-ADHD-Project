# setwd("/home/usdandres/Documents/Study_stuff/George/Group_Lasso_Project/Research18/")
source("Group.Lasso.Paper.Functions.R")

set.seed(2)


p <- 10                                     # number of variables per subject
print(c("p",p))
t.set <- c(30)                             # number of observed time points per subject
print(c("t.set:",t.set))
K.set <- c(20)                              # number of subjects per group
print(c("K.set:",K.set))
heter <- c("low","moderate")[1]             # level of heterogeneity: low (1%) or moderate (2-3%)
print(c("heter:",heter))
Thresh2 <- 0.02                             # Threshold to be applied to estimates get rid of noise
print(c("Thresh2:",Thresh2))
signs <- c(1,-1)                            # Whether effects will be allowed to be positive(1), negative(-1) or both
print(c("Signs:",signs))
SNR <- 2                                    # Signal-to-Noise ratio
print(c("SNR:",SNR))
spectral <- 0.6                             # Maximum eigenvalue allowed in generated VAR transition matrices (either 0.4 or 0.6)
print(c("spectral:",spectral))
pos.diag <- T                               # whether all diagonal elements of VAR transition matrices should be positive
print(c("pos.diag=",pos.diag))
n.groups <- 1                               # number of subject groups generated
print(c("n.groups:",n.groups))
rep <- 1                                   # number of simulations to run
print(c("rep",rep))

#####
## CRITERIONS AND OTHER PARAMETERS
## Criterions can only be: BIC, AIC and AIC corrected
#####

criter <- "BIC"                                   # selection criterion for group lasso (first stage)
print(c("criter",criter))
criter.second <- "BIC"                            # selection criterion for sparse lasso (second stage)
print(c("criter.second:",criter.second))
max.iter <- 20                                    # maximum number of iterations for the two-stage estimation algorithm
print(c("max.iter:",max.iter))
eps <- 0.01                                       # stopping criterion value for two-stage estimation algorithm
print(c("eps:",eps))

D <- 1                                           # lag order of generated VAR models (this code has been mostly tested for the case of VAR(1) models)
print(c("order:",D))


######
## Generating different settings for different values of p:
## if p <= 20:  5% common off-diagonal edge density, 1% individual density of low heterogeneity, 3% - for moderate.
## if p > 20: 2% common off-diagonal edge density, 1% individual density of low heterogeneity, 2% - for moderate.
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

#####
## Parameters for glmnet function: no-intercept standardized model
#####

intercept <- FALSE
standardize <- TRUE


for(K in K.set){
  
  ### For VAR of order D we will store the estimates for a single subject in the following (D*p)xp matrix:
  ### A_{(D*p)xp} = [A1_pxp,]
  ###                ...,
  ###               [AD_pxp ]
  
  Final.Est <- make.list(make.list(matrix(0,D*p,p),K),rep)
  Final.Comm.Est <- make.list(make.list(matrix(0,D*p,p),K),rep)
  Final.Ind.Est <-  make.list(make.list(matrix(0,D*p,p),K),rep)
  A.true.all <- make.list(make.list(matrix(0,D*p,p),K),rep)
  A.comm.all <- make.list(matrix(0,D*p,p),rep)
  A.ind.all <- make.list(make.list(matrix(0,D*p,p),K),rep)
  
  for(train in t.set){

    t <- train
    ptm <- proc.time()
    
    ######
    ### Initializing vectors
    ######
    
    FP.Rand.pen <- rep(0,rep)
    FN.Rand.pen <- rep(0,rep)
    TP.Rand.pen <- rep(0,rep)
    TN.Rand.pen <- rep(0,rep)
    Matt.Coef.Rand.pen <- rep(0,rep)
    
    FP.Rand.pen.comm <- rep(0,rep)
    FN.Rand.pen.comm <- rep(0,rep)
    TP.Rand.pen.comm <- rep(0,rep)
    TN.Rand.pen.comm <- rep(0,rep)
    Matt.Coef.Rand.pen.comm <- rep(0,rep)
    
    FP.Rand.pen.ind <- rep(0,rep)
    FN.Rand.pen.ind <- rep(0,rep)
    TP.Rand.pen.ind <- rep(0,rep)
    TN.Rand.pen.ind <- rep(0,rep)
    Matt.Coef.Rand.pen.ind <- rep(0,rep)
    
    n.iter <- array(0,c(rep,p))
    
    namedir <- paste(getwd(),"/CONSTRAINED_Simul_Results_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_Rep=",rep,"_SNR=",SNR,"_criter=",criter,"_criter-second=",criter.second,"_heter=",heter,"_signs=",length(signs),"_max_eig_comm=",max_eig_comm,"_max_eig_ind=",max_eig_ind,"_min_elem=",min_elem,sep="")
    dir.create(namedir) 
    
    for(run in 1:rep){
      print(c("run",run))
      
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
      
        print("Generated")

      ### Initializing vectors
      
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
        for(i in 1:K){
          #Sigma[[i]] <- diag(runif(p,0,1))
          Sigma[[i]] <- diag(1,p)               # simply setting the error covariance to be identity matrix
        }
        
        Sigma.full <- block.diag(Sigma)
        
        ################################
        ####### DATA GENERATION ########
        ################################
        
        DATA <- gen_dat(T=t,
                        A=A.list,
                        SNR=SNR,
                        Sigma_error=Sigma.full)
        
        ### Calculating maximum likelihood estimates of sigma^2 for each subject
        
           sigma2 <- rep(0,K)
           for (k in 1:K) sigma2[k] <- OLS.tseries(DATA[(k-1)*p + 1:p,],D=D)$sigma2
        
        
        #############################
        ### FULL PROBLEM SETUP   ####
        #############################
        
        M.setup <- mat.setup(DATA,t,K,p,D=D) 

        C.list <- M.setup$C
        B.list <- M.setup$B
        X.list <- M.setup$X
        
        ### vector of group number assignments
        
        group <- c(1:(D*p))
        if (K>1){
          for (j in 2:K){
            group <- c(group,1:(D*p))
          }
        }
        
        
        ## Initializing vectors to contain estimates during algorithm iterations
        
        grouplasso.result_before <- make.list(numeric(D*K*p),p)
        seplasso.result_before <- make.list(numeric(D*K*p),p)
        seplasso.result <- make.list(numeric(D*K*p),p)
        grouplasso.result <- make.list(numeric(D*K*p),p)
        gl.zeros <- list()
        
        for (j in 1:p){
          print(paste("j:",j,sep=""))
          it <- 0
          flag <- 0
          
          while ((it<max.iter) & (flag==0)){ 
            it <- it+1
            #  print(paste("Iter:",it,sep=""))
            
            lasso.freq <- matrix(0,1,D*K*p)  
            grouplasso.freq <- matrix(0,1,D*K*p)
            Y <- numeric(K*(t-D))
            
            Xbeta <- X.list %*% seplasso.result[[j]]
            
            for (k in 1:K){
              Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p] - Xbeta[(k-1)*(t-D) + (1:(t-D))]
            }
            
            ###############################################################
            ## FIRST STAGE: Group lasso estimation of common component ####
            ###############################################################
            
              r <- grpreg(X.list/sqrt(sigma2[j]),
                          Y/sqrt(sigma2[j]),
                          group=group,
                          penalty="grLasso",
                          family="gaussian",
                          intercept=intercept,
                          warn=FALSE)
              
              lambda_G.path <- r$lambda
              est <- r$beta[-1,]
              grouplasso.result.df <- r$df
              
              if (criter == "AIC")
                group.est.out <- AIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,lambda.path=lambda_G.path,df.path=grouplasso.result.df)
              if (criter == "BIC")
                group.est.out <- BIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=grouplasso.result.df,lambda.path=lambda_G.path)
              if (criter == "AICc")
                group.est.out <- AICc(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=grouplasso.result.df,lambda.path=lambda_G.path)
              
              grouplasso.result[[j]] <- sparsify(group.est.out$Est, Thresh2)
              lambda.group <- group.est.out$lambda1
              
              
            ####################################################################
            ## SECOND STAGE: Sparse lasso estimation of individual component ###
            ####################################################################
            
             gl.zeros[[j]] <- (grouplasso.result[[j]] == 0)
            
             Xbeta <- X.list %*% grouplasso.result[[j]]
             for (k in 1:K){
               Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p] - Xbeta[(k-1)*(t-D) + (1:(t-D))]
             }
            
               X.zero <- X.list[,gl.zeros[[j]]]

               r <- glmnet(X.zero/sqrt(sigma2[j]),
                           Y/sqrt(sigma2[j]),
                           family="gaussian",
                           standardize=standardize,
                           intercept=intercept)
            
             
              lambda_SPARS.path <- r$lambda
              est <- r$beta
              sep.df <- r$df

              if (criter.second == "AIC")
                sep.est.out <- AIC(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,df.coef=df.coef,lambda.path=lambda_SPARS.path)
              if (criter.second == "BIC")
                sep.est.out <- BIC(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
              if (criter.second == "AICc")
                sep.est.out <- AICc(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
              
              lambda.sparse <- sep.est.out$lambda1
              seplasso.result[[j]] <- rep(0,D*p*K)
              seplasso.result[[j]][gl.zeros[[j]]] <- sparsify(sep.est.out$Est, Thresh2)

              for (k in 1:K){
                Group.Est[[k]][j,] <- grouplasso.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)]
                Sep.Est.Second[[k]][j,] <- seplasso.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)]
              }
              
            if (it>1){
              #####
              ## STOPPING CRITERION
              #####
                diff_magnitudes <- sum((c(grouplasso.result[[j]],seplasso.result[[j]])-c(grouplasso.result_before[[j]],seplasso.result_before[[j]]))^2)
                # print(diff_magnitudes)
                if (diff_magnitudes < eps){
                  flag <- 1
                  n.iter[run,j] <- it
                }
            }
              grouplasso.result_before <- grouplasso.result  
              seplasso.result_before <- seplasso.result  
          }
          
          if (it == max.iter) n.iter[run,j] <- max.iter
        }
        

        Group.Est.Common <- list()
        Group.Est.Common <- Group.Est
        
        for(k in 1:K){
          Group.Final[[k]] <- Group.Est[[k]] + Sep.Est.Second[[k]]
          Final.Est[[run]][[k]] <- Group.Final[[k]]
          Final.Comm.Est[[run]][[k]] <- Group.Est[[k]] 
          Final.Ind.Est[[run]][[k]] <- Sep.Est.Second[[k]]
          Group.Final[[k]] <- sparsify(Group.Final[[k]],Thresh2)
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
          Rand.pen.vec <- c(Rand.pen.vec,vec(Group.Final[[k]]))
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
      
      #######
      #### METRICS FOR COMMON COMPONENT ESTIMATE
      #######
      
      true.vec.comm <- NULL
      Rand.pen.vec.comm <- NULL
      
      for (k in 1:K){
        A.bind <- A.true$A.comm[[1]]
        if (D>1) for (d in 2:D) A.bind <- rbind(A.bind,A.true$A.comm[[d]])
        true.vec.comm <- c(true.vec.comm,A.bind)
        Rand.pen.vec.comm <- c(Rand.pen.vec.comm,vec(Group.Est[[k]]))
      }
      
      Measures.Joint <- Measures.Vec(Rand.pen.vec.comm,as.matrix(true.vec.comm))
      
      FP.Rand.pen.comm[run] <- Measures.Joint$FP
      FN.Rand.pen.comm[run] <- Measures.Joint$FN
      TP.Rand.pen.comm[run] <- Measures.Joint$TP
      TN.Rand.pen.comm[run] <- Measures.Joint$TN
      Matt.Coef.Rand.pen.comm[run] <- Matthews.Coef(TP.Rand.pen.comm[run],
                                                    FP.Rand.pen.comm[run],
                                                    TN.Rand.pen.comm[run],
                                                    FN.Rand.pen.comm[run])
      
      #######
      #### METRICS FOR INDIVIDUAL COMPONENT ESTIMATE
      #######
      
      true.vec.ind <- NULL
      Rand.pen.vec.ind <- NULL
      
  
          for(k in 1:K){
            A.bind <- A.true$A.ind[[k]][[1]]
            if (D>1) for (d in 2:D) A.bind <- rbind(A.bind,A.true$A.ind[[k]][[d]])
            true.vec.ind <- c(true.vec.ind,A.bind)
            Rand.pen.vec.ind <- c(Rand.pen.vec.ind,vec(Sep.Est.Second[[k]]))
          }
      
      Measures.Joint <- Measures.Vec(Rand.pen.vec.ind,as.matrix(true.vec.ind))
      
      FP.Rand.pen.ind[run] <- Measures.Joint$FP
      FN.Rand.pen.ind[run] <- Measures.Joint$FN
      TP.Rand.pen.ind[run] <- Measures.Joint$TP
      TN.Rand.pen.ind[run] <- Measures.Joint$TN
      Matt.Coef.Rand.pen.ind[run] <- Matthews.Coef(TP.Rand.pen.ind[run],
                                                   FP.Rand.pen.ind[run],
                                                   TN.Rand.pen.ind[run],
                                                   FN.Rand.pen.ind[run])
      
    }
    
      ###############################################
      #### Saving all the estimates to .rds files ###
      ###############################################
      
      saveRDS(Final.Est,paste(namedir,"/Final_Est.rds",sep=""))
      saveRDS(Final.Comm.Est,paste(namedir,"/Common_Est.rds",sep=""))
      saveRDS(Final.Ind.Est,paste(namedir,"/Ind_Est.rds",sep=""))
      saveRDS(A.true.all,paste(namedir,"/A_true.rds",sep=""))
      saveRDS(A.comm.all,paste(namedir,"/A_true_comm.rds",sep=""))
      saveRDS(A.ind.all,paste(namedir,"/A_true_ind.rds",sep=""))
      saveRDS(n.iter,paste(namedir,"/N_Iter.rds",sep=""))
    
      time.taken <- proc.time() - ptm
      #print("Time taken:")
      #print(time.taken)
    
    ##################################################
    ### Printing all the metrics in a format for latex table
    ##################################################
      
    cat("\n")
    cat("\n")
    cat("p=",p," & ",t," & ",K," & ",
        round(mean(FP.Rand.pen),2),"(",round(sd(FP.Rand.pen),2),") & ",
        round(mean(FN.Rand.pen),2),"(",round(sd(FN.Rand.pen),2),") & ",
        round(mean(Matt.Coef.Rand.pen),2),"(",round(sd(Matt.Coef.Rand.pen),2),") & ",
        round(mean(FP.Rand.pen.comm),2),"(",round(sd(FP.Rand.pen.comm),2),") & ",
        round(mean(FN.Rand.pen.comm),2),"(",round(sd(FN.Rand.pen.comm),2),") & ",
        round(mean(Matt.Coef.Rand.pen.comm),2),"(",round(sd(Matt.Coef.Rand.pen.comm),2),") & ",
        round(mean(FP.Rand.pen.ind),2),"(",round(sd(FP.Rand.pen.ind),2),") & ",
        round(mean(FN.Rand.pen.ind),2),"(",round(sd(FN.Rand.pen.ind),2),") & ",
        round(mean(Matt.Coef.Rand.pen.ind),2),"(",round(sd(Matt.Coef.Rand.pen.ind),2),") & ",
        round(mean(n.iter),2),"(",round(sd(n.iter),2),") & ",
        sum(n.iter == max.iter),"\\\\",sep="")
    cat("\n")
    cat("\\hline")
    cat('\n') 
 
  }
}
