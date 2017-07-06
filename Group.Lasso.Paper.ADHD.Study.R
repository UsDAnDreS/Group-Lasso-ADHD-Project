# setwd("/home/usdandres/Documents/Study_stuff/George/Group_Lasso_Project/Research18")
source("Group.Lasso.Paper.Functions.R")

region_names <- c("Left Auditory Cortex",
                  "Right Auditory Cortex",
                  "Stria terminalis",
                  "Left Default Mode Network",
                  "Medial Prefrontal Cortex",
                  "Front Default Mode Network",
                  "Right Default Mode Network",
                  "Occipital Lobe",
                  "Motor Cortex",
                  "Right Dorsolateral Prefrontal Cortex",
                  "Right Polar Frontal Lobe",
                  "Right Parietal Lobe",
                  "Right Inferior Temporal Cortex",
                  "Basal Ganglia",
                  "Left Parietal Lobe",
                  "Left Dorsolateral Prefrontal Cortex",
                  "Left Polar Frontal Lobe",
                  "Left Intraparietal Sulcus",
                  "Right Intraparietal Sulcus",
                  "Left Lateral Occipital Complex",
                  "Primary Visual Cortex",
                  "Right Lateral Occipital Complex",
                  "Dorsal Anterior Cingulate Cortex",
                  "Ventral Anterior Cingulate Cortex",
                  "Right Anterior Insular Cortex",
                  "Left Superior Temporal Sulcus",
                  "Right Superior Temporal Sulcus",
                  "Left Temporoparietal Junction",
                  "Broca Area of Frontal Lobe",
                  "Superior Frontal",
                  "Right Temporoparietal Junction",
                  "Pars Opercularis",
                  "Cerebellum",
                  "Dorsal Posterior Cingulate Cortex",
                  "Left Insular Cortex",
                  "Cingulate Cortex",
                  "Right Insular Cortex",
                  "Left Anterior Intraparietal Sulcus",
                  "Right Anterior Intraparietal Sulcus")


## Uploading ADHD and controls data

initial_path <- getwd()
set.seed(2)
new_path <- paste(initial_path,"/fMRI_data_stuff/ADHD_TwoGroups",sep="")
setwd(new_path)

temp_ADHD = list.files(path=getwd(),pattern=('ADHD=1.*MSDL.*.*detrend=True_standardize=True_TimeSeries.csv'))
temp_Control = list.files(path=getwd(),pattern=('ADHD=0.*MSDL.*.*detrend=True_standardize=True_TimeSeries.csv'))


#####
## Set the number of observed time points
## Then pick all the patients that have at least this many observations
#####

train <- c(120)                                                            # number of observed time points for each subject
print(c("train:",train))

myfiles_ADHD = lapply(temp_ADHD,function(x) read.csv(x,header=F))
Total1 <- length(myfiles_ADHD)
select_ind_1 <- unlist(lapply(myfiles_ADHD,function(x) dim(x)[1]>train))
select_ind_1 <- c(1:Total1)[select_ind_1]

myfiles_Control = lapply(temp_Control,function(x) read.csv(x,header=F))
Total2 <- length(myfiles_Control)
select_ind_2 <- unlist(lapply(myfiles_Control,function(x) dim(x)[1]>train))
select_ind_2 <- c(1:Total2)[select_ind_2]


print(c("Total1:",Total1))                                                # total of ADHD patients with >t time points observed
print(c("Total2:",Total2))                                                # total of control patients with >t time points observed

### Which patient group will be estimated: 1 (ADHD) or 2 (Control)
GN <- 2
print(c("GN:",GN))

#p <- 5
p <- dim(myfiles_ADHD[[select_ind_1[1]]])[2]                            # number of brain regions per patient
print(c("p",p))
Thresh2 <- 0.02                                                         # threshold to be applied to estimates get rid of noise
print(c("Thresh2:",Thresh2))
K <- min(length(select_ind_1),length(select_ind_2))                     # each group has at least K patients
print(c("K:",K))
R <- 30                                                                  # number of bootstrapped samples 
print(c("R:",R))                                                        # (if R == 0 then just make SINGLE RUN for all T time points, no bootstrapping)

bl.len.set <- c(0.3)                                                    # length of generated block bootstraps (as proportion of number of time points)
print(c("bl.len.set",bl.len.set))
offset <- 0                                                             # whether to skip some initial time points of fMRI measurements
print(c("offset:", offset))
d <- 1                                                                  # VAR(d), here the case of d=1 is prioritized
print(c("order:",d))  

#####
## CRITERIONS AND OTHER PARAMETERS
## Criterions can only be: BIC, AIC and AIC corrected
#####

criter.comm <- "BIC"                                                     # selection criterion for group lasso (first stage)
print(c("criter.comm",criter.comm))
criter.ind <- "BIC"                                                      # selection criterion for sparse lasso (second stage)
print(c("criter.ind",criter.ind))
max.iter <- 20                                                          # maximum number of iterations for the two-stage estimation algorithm
print(c("max.iter",max.iter))
eps <- 0.01                                                             # stopping criterion value for two-stage estimation algorithm
print(c("eps:",eps))

#####
## Parameters for glmnet function: no-intercept standardized model
#####

intercept <- FALSE
standardize <- TRUE


###########################################33

### Initializing some vectors
n.iter <- array(0,c(R,p))
t <- train

ptm <- proc.time()


if (R == 0){
  Group.Est <- make.list(matrix(0,p,p),K)
  Sep.Est.Second <-  make.list(matrix(0,p,p),K)
  Group.Final <- make.list(matrix(0,p,p),K)
  Group.Prop <- matrix(0,p,p)
}

if (R > 0){
  Group.Est <- make.list(make.list(matrix(0,p,p),R),K)
  Sep.Est.Second <-  make.list(make.list(matrix(0,p,p),R),K)
  Group.Final <- make.list(make.list(matrix(0,p,p),R),K)
  Group.Prop.Resamp <- make.list(matrix(0,p,p),R)
  Group.Prop <- matrix(0,p,p)
}

for (bl.len in bl.len.set){
  
  print(c("bl.len:",bl.len))
  
  final_path <- paste("ADHD_MSDL_Two_Groups_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_R=",R,"_criter.comm=",criter.comm,"_criter.ind=",criter.ind,"_Thresh=",Thresh2,sep="")
  dir.create(final_path)
  
  ################################
  ####### DATA GENERATION ########
  ################################
  
  selection_1 <- sample(select_ind_1,K)
  selection_2 <- sample(select_ind_2,K)
  
  print(c("GN:",GN))
  
  ##########
  ### Putting together the dataset for the group of patients (ADHD or control)
  ##########
  
  if (GN == 1){
    DATA <- t(myfiles_ADHD[[selection_1[1]]])[1:p, offset + 1:train ]
    for(j in selection_1[-1]){
      DATA <- rbind(DATA,t(myfiles_ADHD[[j]])[1:p, offset + 1:train])
    }
  }
  
  if (GN == 2){
    DATA <- t(myfiles_Control[[selection_2[1]]])[1:p,offset + 1:train]
    for(j in selection_2[-1]){
      DATA <- rbind(DATA,t(myfiles_Control[[j]])[1:p,offset + 1:train])
    }
  }
  
  
  ## Uncomment if you would like to see the time series over each brain region for the first subject
  
  # par(mfrow=c(1,1))
  # for(j in 1:K){ 
  #   if (j == 1){
  #     pdf(paste(final_path,"/GroupType=",GN,"_TimeSeries.pdf",sep=""),width=12)
  #     rbow <- rainbow(p)
  #     for (p1 in c(1:p)){
  #       if (p1 == 1) plot(ts(DATA[(j-1)*p + p1,]),col=rbow[p1],ylab="Region BOLD Signal",ylim=c(-3.3,3.3),main=paste("Subject #",j))
  #       if (p1 > 1) lines(ts(DATA[(j-1)*p + p1,]),col=rbow[p1])
  #     }
  #     dev.off()
  #   }
  # }
  
  ########
  ### Generating R time series block bootstrapped samples (if R>0)
  ########
  
  if (R>0){
    Data.boot <- list()
    for (j in 1:K){
      Data.boot[[j]] <- list()
      tsb <- tsboot(t(DATA[(j-1)*p + 1:p,]),statistic=identity,R=R,l=bl.len*train,sim="geom")
      for (r in 1:R){
        Data.boot[[j]] <- append(Data.boot[[j]],list(matrix(tsb$t[r,],p,dim(DATA)[2],byrow=T)))
      }
    }
  }
  ########
  ### Looping through all bootstrapped samples if R>0
  ### If R == 0 => just a single run for all 'train' data points
  ########
  
  for(i in 1:R){
  
    print(c("Bootstrap:",i))
    
    ######
    ### Compact matrix form setup for the current bootstrap (if R>0)
    ######
    
    if (R>0){
      DATA.boot <- matrix(0,dim(DATA)[1],dim(DATA)[2])
      for (j in 1:K)  DATA.boot[(j-1)*p + 1:p,] <- Data.boot[[j]][[i]]
      Mat.obj.FULL <- mat.setup(DATA.boot,train,K,p,d=1)
    }
    
    ######
    ### Compact matrix form setup for the full set of 'train' data points (if R==0)
    ######
    
    if (R==0){
      Mat.obj.FULL <- mat.setup(DATA,train,K,p,d=1)
    }
    
    C.list <- Mat.obj.FULL$C
    B.list <- Mat.obj.FULL$B
    X.list <- Mat.obj.FULL$X
    
    ### vector of group number assignments
    
    group <- c(1:p)
    for (j in 2:K){
      group <- c(group,1:p)
    }
    
    ## Initializing vectors to contain estimates during algorithm iterations
    
    grouplasso.result_before <- make.list(numeric(K*p),p)
    seplasso.result_before <- make.list(numeric(K*p),p)
    seplasso.result <- make.list(numeric(K*p),p)
    grouplasso.result <- make.list(numeric(K*p),p)
    gl.zeros <- list()
    
    for (j in 1:p){
      print(paste("j:",j,sep=""))
      it <- 0
      flag <- 0
      
      while ((it<max.iter) & (flag==0)){ 
        it <- it+1
        print(paste("Iter:",it,sep=""))
  
        lasso.freq <- matrix(0,1,K*p)  
        grouplasso.freq <- matrix(0,1,K*p)
        Y <- numeric(K*(t-1))
        
        Xbeta <- X.list %*% seplasso.result[[j]]
        
        for (k in 1:K){
          Y[(k-1)*(t-1) + (1:(t-1))] <- C.list[1:(t-1),j + (k-1)*p] - Xbeta[(k-1)*(t-1) + (1:(t-1))]
        }
          
        ###############################################################
        ## FIRST STAGE: Group lasso estimation of common component ####
        ###############################################################
        
          r <- grpreg(X.list,
                      Y,
                      group=group,
                      penalty="grLasso",
                      family="gaussian",
                      intercept=FALSE,
                      warn=FALSE)
          
          lambda_G.path <- r$lambda
          est <- r$beta[-1,]
          grouplasso.result.df <- r$df

          if (criter.comm == "AIC")
            group.est.out <- AIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,lambda.path=lambda_G.path,df.path=grouplasso.result.df)
          if (criter.comm == "BIC")
            group.est.out <- BIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=grouplasso.result.df,lambda.path=lambda_G.path)
          if (criter.comm == "AICc")
            group.est.out <- AICc(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=grouplasso.result.df,lambda.path=lambda_G.path)
          
          grouplasso.result[[j]] <- sparsify(group.est.out$Est, Thresh2)
          lambda.group <- group.est.out$lambda1
        
        
          ###############################################################
          ## SECOND STAGE: Group lasso estimation of common component ###
          ###############################################################
        
          gl.zeros[[j]] <- (grouplasso.result[[j]] == 0)
        
          Xbeta <- X.list %*% grouplasso.result[[j]]
          for (k in 1:K){
            Y[(k-1)*(t-1) + (1:(t-1))] <- C.list[1:(t-1),j + (k-1)*p] - Xbeta[(k-1)*(t-1) + (1:(t-1))]
          }
        
          X.zero <- X.list[,gl.zeros[[j]]]

          r <- glmnet(X.zero,
                      Y,
                      family="gaussian",
                      standardize=standardize,
                      intercept=intercept)
          
          lambda_SPARS.path <- r$lambda
          est <- r$beta
          sep.df <- r$df
          
          if (criter.ind == "AIC")
            sep.est.out <- AIC(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,df.coef=df.coef,lambda.path=lambda_SPARS.path)
          if (criter.ind == "BIC")
            sep.est.out <- BIC(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
          if (criter.ind == "AICc")
            sep.est.out <- AICc(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
          
          lambda.sparse <- sep.est.out$lambda1
          seplasso.result[[j]] <- rep(0,p*K)
          seplasso.result[[j]][gl.zeros[[j]]] <- sparsify(sep.est.out$Est, Thresh2)
          
        if (R>0){
          for (k in 1:K){
            Group.Est[[k]][[i]][j,] <- grouplasso.result[[j]][(k-1)*p+(1:p)]
            Sep.Est.Second[[k]][[i]][j,] <- seplasso.result[[j]][(k-1)*p+(1:p)]
          }
        }
        
        if (R == 0){
          for (k in 1:K){
            Group.Est[[k]][j,] <- grouplasso.result[[j]][(k-1)*p+(1:p)]
            Sep.Est.Second[[k]][j,] <- seplasso.result[[j]][(k-1)*p+(1:p)]
          }
        }

        if (it>1){
          #####
          ## STOPPING CRITERIONS
          #####
            diff_magnitudes <- sum((c(grouplasso.result[[j]],seplasso.result[[j]])-c(grouplasso.result_before[[j]],seplasso.result_before[[j]]))^2)
         #  print(diff_magnitudes)
            if (diff_magnitudes < 0.01){
              flag <- 1
              n.iter[i,j] <- it
            }
        }
        grouplasso.result_before <- grouplasso.result  
        seplasso.result_before <- seplasso.result  
      }
      
      if (it == max.iter) n.iter[i,j] <- max.iter
    }
    
    Group.Est.Common <- list()
    Group.Est.Common <- Group.Est
    
    #######
    ### If we make a single run for all 'train' time points (R==0),
    ### just accumulate estimates for all K subjects
    #######
    
    if (R == 0){
      for(k in 1:K){
        Group.Final[[k]] <- Group.Est[[k]] + Sep.Est.Second[[k]]
        Group.Prop <- Group.Prop + (Group.Est[[k]] != 0)
      }
    }
    
    ######
    ### Making bootstrapped runs - accumulate estimates for all K subjects (k) for each bootstrap (i)
    ######
    
     if (R > 0){
      for(k in 1:K){
        Group.Final[[k]][[i]] <- Group.Est[[k]][[i]] + Sep.Est.Second[[k]][[i]]
        Group.Prop.Resamp[[i]] <- Group.Prop.Resamp[[i]] + (Group.Est[[k]][[i]] != 0)
      }
      Group.Prop.Resamp[[i]] <- Group.Prop.Resamp[[i]]/K
    }
  }
  
  ######
  ### Combine accumulated estimates across bootstrapped samples, 
  ### to calculate proportions of times this effect was selected (Group.Prop)
  ######
  
  if (R>0){ 
    for (i in 1:R) Group.Prop <- Group.Prop + Group.Prop.Resamp[[i]]
    Group.Prop <- Group.Prop/R
  }
  
  ###############################################
  #### Saving all the estimates to .rds files ###
  ###############################################
  
  saveRDS(Group.Prop,paste(final_path,"/Group_Proportions_GN=",GN,".rds",sep=""))
  saveRDS(Group.Est.Common,paste(final_path,"/Common_Estimates_GN=",GN,".rds",sep=""))
  saveRDS(Group.Est,paste(final_path,"/Full_Estimates_GN=",GN,".rds",sep=""))
  saveRDS(Sep.Est.Second,paste(final_path,"/Individ_Estimates_GN=",GN,".rds",sep=""))
  saveRDS(n.iter,paste(final_path,"/NumberIter_GN=",GN,".rds",sep=""))
}

time.taken <- proc.time() - ptm
#print("Time taken:")
#print(time.taken)
