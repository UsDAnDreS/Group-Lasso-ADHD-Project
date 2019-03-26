###########
## Accumulating results across multiple runs with "rep" != "actual.rep".
###########

## Where all the customized functions (for estimation, performance metrics calculation, etc) are contained.
source("MARCH_2019_SUBMISSION_Group_Lasso_Paper_Functions.R")

set.seed(2)

### Method you would like to obtain the results for.
method.name <- c("OurMethod_CONVERGENCE",
                 "Sparse_Lasso",
                 "GroupBridge",
                 "Sparse_Group_Lasso")[4]

# For Sparse Group Lasso.
alpha <- 0.2    # 0 - sheer group lasso; 1 - sheer sparse lasso.

p <- 10                                     # number of variables per subject
print(c("p",p))
t.set <- c(30)                              # number of observed time points per subject
print(c("t.set:",t.set))
K.set <- c(5)                              # number of subjects per group
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
actual.rep <- 2                           # In case we only need "actual.rep" last replicates
print(c("actual.rep",actual.rep))           # (appropriate when running all "rep" replicates takes too long)


#####
## CRITERIONS AND OTHER PARAMETERS
## Criterions can only be: BIC, AIC and AICc
#####

criter <- "BIC"                                   # selection criterion for group lasso (first stage)
print(c("criter",criter))
criter.second <- "BIC"                            # selection criterion for sparse lasso (second stage)
print(c("criter.second:",criter.second))
max.iter <- 30                                    # maximum number of iterations for the two-stage estimation algorithm
print(c("max.iter:",max.iter))
eps <- 0.0001                                       # stopping criterion value for two-stage estimation algorithm
print(c("eps:",eps))

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
## Parameters for glmnet function: no-intercept standardized model
#####

intercept <- FALSE
standardize <- TRUE


K = K.set
train = t.set

Final.Est <- list()
Final.Comm.Est <- list()
Final.Ind.Est <-  list()
A.true.all <- list()
A.comm.all <- list()
A.ind.all <- list()
n.iter <- c()

#####
## In case some are missing, use the [-ind.miss] indexation of rep.set
#####

rep.set <- seq(actual.rep,rep, by=actual.rep)
rep <- length(rep.set)*actual.rep

FP.Rand.pen <- rep(0,rep)
FN.Rand.pen <- rep(0,rep)
TP.Rand.pen <- rep(0,rep)
TN.Rand.pen <- rep(0,rep)
Matt.Coef.Rand.pen <- rep(0,rep)
Frob.Rand.pen <- rep(0,rep)


### If we deal with our two-stage approach - we also need the common and individual component estimation metrics.
if (method.name %in% c("OurMethod_CONVERGENCE")){
  FP.Rand.pen.comm <- rep(0,rep)
  FN.Rand.pen.comm <- rep(0,rep)
  TP.Rand.pen.comm <- rep(0,rep)
  TN.Rand.pen.comm <- rep(0,rep)
  Matt.Coef.Rand.pen.comm <- rep(0,rep)
  Frob.Rand.pen.comm <- rep(0,rep)
  
  FP.Rand.pen.ind <- rep(0,rep)
  FN.Rand.pen.ind <- rep(0,rep)
  TP.Rand.pen.ind <- rep(0,rep)
  TN.Rand.pen.ind <- rep(0,rep)
  Matt.Coef.Rand.pen.ind <- rep(0,rep)
  Frob.Rand.pen.ind <- rep(0,rep)
  
  n.iter <- c()
}



### Scraping R objects (containing transition matrix estimates) from respective folders (depending on the method).
for (r in rep.set){
  if (method.name == "OurMethod_CONVERGENCE"){
    namedir <- paste(method.name, "_Simul_Results_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_Rep=",r,"_ActualRep=",actual.rep,"_eps=",eps,"_SNR=",SNR,"_criter=",criter,"_criter-second=",criter.second,"_heter=",heter,"_signs=",length(signs),"_max_eig_comm=",max_eig_comm,"_max_eig_ind=",max_eig_ind,"_min_elem=",min_elem,sep="")
  }
  if (method.name == "Sparse_Group_Lasso"){
    namedir <- paste(method.name, "_Simul_Results_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_Rep=",r,"_ActualRep=",actual.rep,"_SNR=",SNR,"_criter=CV_alpha=",alpha,"_heter=",heter,"_signs=",length(signs),"_max_eig_comm=",max_eig_comm,"_max_eig_ind=",max_eig_ind,"_min_elem=",min_elem,sep="")
  }
  if (method.name %in% c("Sparse_Lasso", "GroupBridge")){
    namedir <- paste(method.name, "_Simul_Results_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_Rep=",r,"_ActualRep=",actual.rep,"_SNR=",SNR,"_criter=",criter,"_heter=",heter,"_signs=",length(signs),"_max_eig_comm=",max_eig_comm,"_max_eig_ind=",max_eig_ind,"_min_elem=",min_elem,sep="")
  }
  
  Final.Est <- append(Final.Est, readRDS(paste(namedir,"/Final_Est.rds",sep="")))
  A.true.all <- append(A.true.all, readRDS(paste(namedir,"/A_true.rds",sep="")))
  
  if (method.name %in% c("OurMethod_CONVERGENCE")){
    Final.Comm.Est <- append(Final.Comm.Est, readRDS(paste(namedir,"/Common_Est.rds",sep="")))
    Final.Ind.Est <- append(Final.Ind.Est, readRDS(paste(namedir,"/Ind_Est.rds",sep="")))
    
    A.comm.all <- append(A.comm.all, readRDS(paste(namedir,"/A_true_comm.rds",sep="")))
    A.ind.all <- append(A.ind.all, readRDS(paste(namedir,"/A_true_ind.rds",sep="")))
    n.iter <- append(n.iter, readRDS(paste(namedir,"/N_Iter.rds",sep="")))
  }
  
}



#########################################################
############## CALCULATING PERFORMANCE METRICS   ########
#########################################################

for (run in 1:length(Final.Est)){
  
  #######
  #### METRICS FOR FULL ESTIMATE (needed for all methods)
  #######
  
  true.vec <- NULL
  Rand.pen.vec <- NULL
  for (k in 1:K){
    A.bind <- A.true.all[[run]][[k]][[1]]
    if (D>1) for (d in 2:D) A.bind <- rbind(A.bind,A.true.all[[run]][[k]][[d]])
    true.vec <- c(true.vec,A.bind)
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
  
  
  ## Common and individual component metrics - needed only for our two-stage approach.
  
  if (method.name %in% c("OurMethod_CONVERGENCE")){
    #######
    #### METRICS FOR COMMON COMPONENT ESTIMATE
    #######
    
    true.vec.comm <- NULL
    Rand.pen.vec.comm <- NULL
    for (k in 1:K){
      A.bind <- A.comm.all[[run]][[1]]
      if (D>1) for (d in 2:D) A.bind <- rbind(A.bind,A.comm.all[[run]][[d]])
      true.vec.comm <- c(true.vec.comm,A.bind)
      Rand.pen.vec.comm <- c(Rand.pen.vec.comm,vec(Final.Comm.Est[[run]][[k]]))
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
    Frob.Rand.pen.comm[run] <- Measures.Joint$Frob
    
    #######
    #### METRICS FOR INDIVIDUAL COMPONENT ESTIMATE
    #######
    
    true.vec.ind <- NULL
    Rand.pen.vec.ind <- NULL
    for(k in 1:K){
      A.bind <- A.ind.all[[run]][[k]][[1]]
      if (D>1) for (d in 2:D) A.bind <- rbind(A.bind,A.true$A.ind[[run]][[k]][[d]])
      true.vec.ind <- c(true.vec.ind,A.bind)
      Rand.pen.vec.ind <- c(Rand.pen.vec.ind,vec(Final.Ind.Est[[run]][[k]]))
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
    Frob.Rand.pen.ind[run] <- Measures.Joint$Frob
  }
}

####
## Printing out the performance metrics in a LATEX table format.
####

cat("\n")
cat("\n")
print(method.name)
cat("\n")
cat("\n")

## For our two-stage approach - need metrics for full, common and individual matrix estimates.
if (method.name %in% c("OurMethod_CONVERGENCE")){
  cat("p=",p," & ",train," & ",K," & ",
      round(mean(FP.Rand.pen),2),"(",round(sd(FP.Rand.pen),2),") & ",
      round(mean(FN.Rand.pen),2),"(",round(sd(FN.Rand.pen),2),") & ",
      round(mean(Matt.Coef.Rand.pen),2),"(",round(sd(Matt.Coef.Rand.pen),2),") & ",
      round(mean(Frob.Rand.pen),2),"(",round(sd(Frob.Rand.pen),2),") & ",
      round(mean(FP.Rand.pen.comm),2),"(",round(sd(FP.Rand.pen.comm),2),") & ",
      round(mean(FN.Rand.pen.comm),2),"(",round(sd(FN.Rand.pen.comm),2),") & ",
      round(mean(Matt.Coef.Rand.pen.comm),2),"(",round(sd(Matt.Coef.Rand.pen.comm),2),") & ",
      round(mean(Frob.Rand.pen.comm),2),"(",round(sd(Frob.Rand.pen.comm),2),") & ",
      round(mean(FP.Rand.pen.ind),2),"(",round(sd(FP.Rand.pen.ind),2),") & ",
      round(mean(FN.Rand.pen.ind),2),"(",round(sd(FN.Rand.pen.ind),2),") & ",
      round(mean(Matt.Coef.Rand.pen.ind),2),"(",round(sd(Matt.Coef.Rand.pen.ind),2),") & ",
      round(mean(Frob.Rand.pen.ind),2),"(",round(sd(Frob.Rand.pen.ind),2),") & ",
      round(mean(n.iter),2),"(",round(sd(n.iter),2),") & ",
      sum(n.iter == max.iter),"\\\\",sep="")
}

## For other three methods - need metrics for full estimates.
if (!(method.name %in% c("OurMethod_CONVERGENCE"))){
  cat("p=",p," & ",train," & ",K," & ",
      round(mean(FP.Rand.pen),2),"(",round(sd(FP.Rand.pen),2),") & ",
      round(mean(FN.Rand.pen),2),"(",round(sd(FN.Rand.pen),2),") & ",
      round(mean(Matt.Coef.Rand.pen),2),"(",round(sd(Matt.Coef.Rand.pen),2),") & ",
      round(mean(Frob.Rand.pen),2),"(",round(sd(Frob.Rand.pen),2),") ",
      "\\\\",sep="")
}

cat("\n")
cat("\\hline")
cat('\n') 