############
### SECTION 4.
### File working with outputted estimates for ADHD fMRI data study.
#############

## This particular setting will work for T=150, R=100, BIC criterion for both common and individ, Threshold = 0.02
## (output for each setting takes up about 5mb)

## In order to have this code produce the paper figures, you will need to:
##    1) Unpack the 'ADHD200.zip' file, and 
##    2) Have the 'ADHD200' folder in the same directory as this source code.

if (!exists("initial_path")) initial_path <- getwd()
if (exists("initial_path")) setwd(initial_path)

## Where all the customized functions (for estimation, performance metrics calculation, etc) are contained.
source("MARCH_2019_SUBMISSION_Group_Lasso_Paper_Functions.R")

library(plotrix)

### Names of 39 brain regions according to MSDL atlas.
region_names <- c("L Auditory Cortex",
                  "R Auditory Cortex",
                  "Stria terminalis",
                  "L Default Mode Network",
                  "M PreFr Cortex",
                  "Fr Default Mode Network",
                  "R Default Mode Network",
                  "Occipital Lobe",
                  "Motor Cortex",
                  "R Dorsolateral PreFr Cortex",
                  "R Polar Fr Lobe",
                  "R Parietal Lobe",
                  "R Inf Temp Cortex",
                  "Basal Ganglia",
                  "L Parietal Lobe",
                  "L Dorsolateral Prefrontal Cortex",
                  "L Polar Fr Lobe",
                  "L Intraparietal Sulcus",
                  "R Intraparietal Sulcus",
                  "L Lateral Occipital Complex",
                  "Primary Visual Cortex",
                  "R Lateral Occipital Complex",
                  "Dorsal Ant Cingulate Cortex",
                  "Ventral Ant Cingulate Cortex",
                  "R Ant Insular Cortex",
                  "L Sup Temp Sulcus",
                  "R Sup Temp Sulcus",
                  "L Temporoparietal Junction",
                  "Broca Area of Fr Lobe",
                  "Sup Fr",
                  "R Temporoparietal Junction",
                  "Pars Opercularis",
                  "Cerebellum",
                  "Dorsal Post Cingulate Cortex",
                  "L Insular Cortex",
                  "Cingulate Cortex",
                  "R Insular Cortex",
                  "L Ant Intraparietal Sulcus",
                  "R Ant Intraparietal Sulcus")

##  R - Right, 
##  L - Left, 
##  Ant - Anterior, 
##  Post - Posterior,
##  Fr - Frontal, 
##  PreFr - Prefrontal, 
##  M - Medial,
##  Sup - Superior,
##  Inf - Inferior,
##  Temp - Temporal,
##  V - Ventral,
##  IPS - IntraParietal Sulcus,
##  Cing - Cingulate,
##  Ins - Insular,
##  Att - Attention,
##  DMN - Default Mode Network


group_name <- c("ADHD","Control")            ## group names - ADHD and controls
TR <- 2                                      ## repetitions time set to 2s


## Folder where to find both the data and resulting estimates.
initial_path <- getwd()
new_path <- paste(initial_path,"/ADHD200",sep="")
setwd(new_path)

## Importing ADHD and controls data (both the time series and phenotypic characteristics)

temp_ADHD = list.files(path=getwd(),pattern=('ADHD=1_TR=2.*MSDL.*.*detrend=True_standardize=True_TimeSeries.csv'))
temp_Control = list.files(path=getwd(),pattern=('ADHD=0_TR=2.*MSDL.*.*detrend=True_standardize=True_TimeSeries.csv'))

phenotypic_ADHD = list.files(path=getwd(),pattern=('ADHD=1_TR=2.*MSDL.*.*detrend=True_standardize=True_Phenotypic.csv'))
phenotypic_Control = list.files(path=getwd(),pattern=('ADHD=0_TR=2.*MSDL.*.*detrend=True_standardize=True_Phenotypic.csv'))


############
## READING PHENOTYPIC DATA
############

phenotypic <- list()
phen_ADHD = lapply(phenotypic_ADHD,function(x) read.csv(x,header=F))
phen_Control = lapply(phenotypic_Control,function(x) read.csv(x,header=F))


#####
## Set the number of observed time points
## Then pick all the patients that have at least this many observations
#####

train <- c(150)
print(c("train:",train))
t <- train

###############
### TIME SERIES DATA
###############

myfiles_ADHD = lapply(temp_ADHD,function(x) read.csv(x,header=F))
Total1 <- length(myfiles_ADHD)
#unlist(lapply(myfiles_ADHD,function(x) dim(x)[1]))
select_ind_1 <- unlist(lapply(myfiles_ADHD,function(x) dim(x)[1]>=train))
select_ind_1 <- c(1:Total1)[select_ind_1]
phenotypic[[1]] <- phen_ADHD[select_ind_1]

myfiles_Control = lapply(temp_Control,function(x) read.csv(x,header=F))
Total2 <- length(myfiles_Control)
#unlist(lapply(myfiles_Control,function(x) dim(x)[1]))
select_ind_2 <- unlist(lapply(myfiles_Control,function(x) dim(x)[1]>=train))
select_ind_2 <- c(1:Total2)[select_ind_2]
phenotypic[[2]] <- phen_Control[select_ind_2]

###########
##  Setting parameter values
###########

p <- dim(myfiles_ADHD[[select_ind_1[1]]])[2]
print(c("p",p))
Thresh2 <- 0.02                                  ## Hard-threshold value used in the two-stage algorithm estimation, not subject to change for this file
print(c("Thresh2:",Thresh2))
Thresh1 <- 0.02                                  ## Hard-threshold to be applied to the final estimates
print(c("Thresh1:",Thresh1))
Thresh <- 0.75                                   ## Threshold for the group effect consistency: if effects is selected in Thresh2*R bootstrapped samples, it is consistent
print(c("Thresh:",Thresh))
K <- min(length(select_ind_1),length(select_ind_2))  ## K - number of subjects per group
print(c("K:",K))
R <- 100                                          ## number of bootstrapped samples
print(c("R:",R))
n.groups <- 2                                     ## number of experimental group (two in our case: ADHD and controls)
print(c("n.groups:",n.groups))

#####
## SELECTION CRITERIONS
#####

criter.comm <- "BIC"                                 ## for common component selection (first stage)
print(c("criter.comm",criter.comm))
criter.ind <- "BIC"                                  ## for individual component selection (second stage)
print(c("criter.ind",criter.ind))
max.iter <- 30                                                          # maximum number of iterations for the two-stage estimation algorithm
print(c("max.iter",max.iter))
eps <- 0.0001                                                            # stopping criterion value for two-stage estimation algorithm
print(c("eps:",eps))


########

### Initializing some lists to store the estimates of full, common and individual components for each patient.
Group.Prop <- list()
Group.Common.Est <- list()
Group.Prop.Resamp <- list()
Ind.Est <- list()
Ind.Prop.Resamp <- list()
n.iter <- make.list(array(0,c(R,p)),n.groups)

for (GN in 1:n.groups){
  Group.Prop.Resamp[[GN]] <- make.list(matrix(0,p,p),R)
  Group.Prop[[GN]] <- matrix(0,p,p)
  Ind.Prop.Resamp[[GN]] <- make.list(matrix(0,p,p),K)
}



## Reading the estimates saved as .rds files as a result of executing 'MARCH_2019_SUBMISSION_ADHD_Study_Estimation.R' file
## Accumulating them in order to calculate the proportions of times each effect was selected 
## (for both common and individual components)
## Additionally reading the number of iterations needed for convergence of two-stage algorithm

for (GN in 1:n.groups){
  final_path <- paste("ADHD_MSDL_Two_Groups_CONVERGENCE_Thresh=",Thresh2,"_p=",p,"_t=",train,"_K=",K,"_R=",R,"_criter.comm=",criter.comm,"_criter.ind=",criter.ind,"_eps=",eps,"_Thresh=",Thresh2,sep="")
  Group.Common.Est[[GN]] <- readRDS(paste(final_path,"/Common_Estimates_GN=",GN,".rds",sep=""))
  Ind.Est[[GN]] <- readRDS(paste(final_path,"/Individ_Estimates_GN=",GN,".rds",sep=""))
  n.iter[[GN]] <- readRDS(paste(final_path,"/NumberIter_GN=",GN,".rds",sep=""))
  
  for (k in 1:K){
    for (i in 1:R){Ind.Prop.Resamp[[GN]][[k]] <- Ind.Prop.Resamp[[GN]][[k]] + (Ind.Est[[GN]][[k]][[i]] != 0)}
    Ind.Prop.Resamp[[GN]][[k]] <- Ind.Prop.Resamp[[GN]][[k]]/R
  }
  for (i in 1:R){
    for (k in 1:K){Group.Prop.Resamp[[GN]][[i]] <- Group.Prop.Resamp[[GN]][[i]] + (sparsify(Group.Common.Est[[GN]][[k]][[i]],Thresh1) != 0)}
    Group.Prop.Resamp[[GN]][[i]] <- Group.Prop.Resamp[[GN]][[i]]/K
    Group.Prop[[GN]] <- Group.Prop[[GN]] + Group.Prop.Resamp[[GN]][[i]]
  }
  Group.Prop[[GN]] <- Group.Prop[[GN]]/R
}

#####################################################################################################
###  Plotting the proportions matrices for each subject's individual component:
###  in what proportion of bootstrapped samples this individual effect was selected for this subject
###  We don't have any individual effects showing up in more than 20% of subsamples,
###  pointing to homogeneity of the data (all males, mostly teenagers)
######################################################################################################

##########################
### INDIVIDUAL EFFECTS ###
##########################

#pdf("Proportion_Matrices_Individuals.pdf")
for (GN in 1:n.groups){
  for (k in 1:K) color2D.matplot(Ind.Prop.Resamp[[GN]][[k]],
                                 cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                                 show.values=TRUE,
                                 main=paste("Subject ",k," in ",group_name[GN]," group in ",
                                            phenotypic[[GN]][[k]][1,11],
                                            " Sex:",phenotypic[[GN]][[k]][1,15],
                                            " Age:",phenotypic[[GN]][[k]][1,14],sep=""))
}
#dev.off()


#######################################################################
#### NICE GRAPH REPRESENTATION OF TEMPORAL EFFECTS FOR VAR(1) MODEL ###
#######################################################################

library(igraph)

### Assigning brain region names as row and column names for the estimates

for (GN in 1:n.groups){
  rownames(Group.Prop[[GN]]) <- region_names
  colnames(Group.Prop[[GN]]) <- region_names
  for (k in 1:K){
    rownames(Ind.Prop.Resamp[[GN]][[k]]) <- region_names
    colnames(Ind.Prop.Resamp[[GN]][[k]]) <- region_names
  }
}

## Spacings (big empty strings) for brain region names to show up on the graph appropriately
## (might vary depending on the system)

region_names_0 <- paste(region_names,
                        "                                                   ", sep="")

region_names_1 <- paste("                                                                          ", 
                        region_names,
                        sep="")
nodes <- c(region_names_0,region_names_1)


##########################
### GROUP LEVEL EFFECTS ##
##########################

# Thresh <- 0.8

m <- ifelse(Group.Prop[[1]]>Thresh,Group.Prop[[1]],0) * ifelse(Group.Prop[[2]]>Thresh,Group.Prop[[2]],0)

for (GN in 1:n.groups){
  adj.mat <- ifelse(Group.Prop[[GN]]>Thresh,Group.Prop[[GN]],0)
  x <- c(rep(0,p),rep(5.2,p))
  y <- c(seq(p,1,length=p),seq(p,1,length=p))
  from <- NULL
  to <- NULL
  clr <- NULL
  
  for(j in 1:p)
    for(i in 1:p)
      if (adj.mat[i,j] != 0){
        if (i == j) clr <- c(clr,"blue")
        if (i != j){
          if (m[i,j] != 0) clr <- c(clr,"red")
          if (m[i,j] == 0) clr <- c(clr,"green")
        }
        from <- c(from,region_names_0[j])
        to <- c(to,region_names_1[i])
      }
  
  NodeList <- data.frame(nodes, x ,y)
  EdgeList <- data.frame(from, to)
  a <- graph_from_data_frame(vertices = NodeList, d= EdgeList, directed = TRUE)
  
 ## Uncomment if you'd like it saved to a file.
 # pdf(paste(group_name[GN],"_Fancy_graph.pdf",sep=""))
  plot(a,asp=0.2,
       rescale=F,
       xlim=c(-1,6),ylim=c(0,p),
       vertex.size=3,
       vertex.label.cex=0.7,
       vertex.label.degree=pi,
       vertex.label.dist=3.15,
       rescale=FALSE,
       main=paste("Temporal effects for ",group_name[GN]," group    ",sep=""),
       edge.color=clr,
       edge.arrow.size=0.5)
  ## Uncomment if you'd like it saved to a file.
  # dev.off()
}

################################
#### INDIVIDUAL EFFECTS ########
#### (if interested)    ########
################################


# for (GN in 1:n.groups){
#   pdf(paste(group_name[GN],"_Individual_Fancy_graphS.pdf",sep=""))
#   for (k in 1:K){
#     adj.mat <- Ind.Prop.Resamp[[GN]][[k]]
#     x <- c(rep(0,p),rep(5,p))
#     y <- c(seq(p,1,length=p),seq(p,1,length=p))
#     from <- NULL
#     to <- NULL
#     clr <- NULL
#     
#     for(j in 1:p)
#       for(i in 1:p)
#         if (adj.mat[i,j] != 0){
#           from <- c(from,region_names_0[j])
#           to <- c(to,region_names_1[i])
#         }
#     
#     NodeList <- data.frame(nodes, x ,y)
#     if (is.null(from) & is.null(to)) {from=region_names[1]; to=region_names[1]}
#     EdgeList <- data.frame(from, to)
#     a <- graph_from_data_frame(vertices = NodeList, d= EdgeList, directed = TRUE)
#     
#     plot(a,asp=0.2,
#          rescale=F,
#          xlim=c(-1,6),ylim=c(0,p),
#          vertex.size=3,
#          vertex.label.cex=0.7,
#          vertex.label.degree=pi,
#          vertex.label.dist=3.15,
#          rescale=FALSE,
#          main=paste("Subject ",k," in ",group_name[GN]," group in ",
#                     phenotypic[[GN]][[k]][1,11],
#                     " Sex:",phenotypic[[GN]][[k]][1,15],
#                     " Age:",phenotypic[[GN]][[k]][1,14],sep=""),
#          edge.arrow.size=0.5) 
#   }
#   dev.off()
# }


##################################################################################################
### CALCULATING AVERAGE MAGNITUDES OF NON-ZERO EFFECTS:                                   ########
### TOTAL, ONLY DIAGONALS, ONLY OFF-DIAGONALS                                             ########
###  AND                                                                                  ########
### CALCULATING AVERAGE MAGNITUDES OF SELECTED EFFECTS                                    ########
### (ONES THAT HAPPENED TO BE PICKED IN MORE THAN 'Thresh' R BOOTSTRAPPED SAMPLES):      ########
### TOTAL, ONLY DIAGONALS, ONLY OFF-DIAGONALS                                             ########
##################################################################################################

magnitudes <- list()
magnitudes.diag <- list()
magnitudes.offdiag <- list()
magnitudes.selected <- list()
magnitudes.selected.diag <- list()
magnitudes.selected.offdiag <- list()

for (GN in 1:n.groups){ 
  mags <- NULL
  mags.diag <- NULL
  mags.offdiag <- NULL
  mags.selected <- NULL
  mags.selected.diag <- NULL
  mags.selected.offdiag <- NULL
  
  for (k in 1:K)
    for (b in 1:R){
      ### ALL NON-ZEROS
      Z <- Group.Common.Est[[GN]][[k]][[b]]
      mags <- c(mags,abs(Z[Z != 0]))
      mags.diag <- c(mags.diag,abs(diag(Z)[diag(Z) != 0]))
      Z.off <- Z[row(Z) != col(Z)]
      mags.offdiag <- c(mags.offdiag,abs(Z.off[Z.off != 0]))
      ### ONLY SELECTED
      Z <- Group.Common.Est[[GN]][[k]][[b]]
      mags.selected <- c(mags.selected,abs(Z[Group.Prop[[GN]] >= Thresh]))
      mags.selected.diag <- c(mags.selected.diag,abs(diag(Z)[diag(Group.Prop[[GN]]) >= Thresh]))
      Z.off <- Z[row(Z) != col(Z)]
      G.off <- Group.Prop[[GN]][row(Group.Prop[[GN]]) != col(Group.Prop[[GN]])]
      mags.selected.offdiag <- c(mags.selected.offdiag,abs(Z.off[G.off >= Thresh]))
    }
  magnitudes[[GN]] <- mags
  magnitudes.diag[[GN]] <- mags.diag
  magnitudes.offdiag[[GN]] <- mags.offdiag
  magnitudes.selected[[GN]] <- mags.selected
  magnitudes.selected.diag[[GN]] <- mags.selected.diag
  magnitudes.selected.offdiag[[GN]] <- mags.selected.offdiag
}

##################
## PRINTING OUT ##
##################

### ALL NON-ZEROS
for (GN in 1:n.groups){
  print(mean(magnitudes[[GN]]))
  print(sd(magnitudes[[GN]]))
}
for (GN in 1:n.groups){
  print(mean(magnitudes.diag[[GN]]))
  print(sd(magnitudes.diag[[GN]]))
}
for (GN in 1:n.groups){
  print(mean(magnitudes.offdiag[[GN]]))
  print(sd(magnitudes.offdiag[[GN]]))
}

### ONLY SELECTED
for (GN in 1:n.groups){
  print(mean(magnitudes.selected[[GN]]))
  print(sd(magnitudes.selected[[GN]]))
}
for (GN in 1:n.groups){
  print(mean(magnitudes.selected.diag[[GN]]))
  print(sd(magnitudes.selected.diag[[GN]]))
}
for (GN in 1:n.groups){
  print(mean(magnitudes.selected.offdiag[[GN]]))
  print(sd(magnitudes.selected.offdiag[[GN]]))
}

##################################
#### AGE AND GENDER VARIABLES  ###
##################################

ages <- list()
genders <- list()

for (GN in 1:n.groups){
  age <- NULL
  gender <- NULL
  for (k in 1:K){
    age <- c(age,phenotypic[[GN]][[k]][1,14])
    gender <- c(gender,as.character(phenotypic[[GN]][[k]][1,15]))
  }
  ages[[GN]] <- age 
  genders[[GN]] <- gender
}

##################
## PRINTING OUT ##
##################

## AGE

for (GN in 1:n.groups){
  print(mean(ages[[GN]]))
  print(sd(ages[[GN]]))
}

## GENDER
for (GN in 1:n.groups){
  print(table(genders[[GN]]))
}
