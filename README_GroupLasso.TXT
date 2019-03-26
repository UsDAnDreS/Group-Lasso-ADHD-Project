MOST FILES CAN BE FOUND ON PUBLIC ACCESS REPOSITORY:

https://github.com/UsDAnDreS/Group-Lasso-ADHD-Project


####################
##  MARCH_2019_SUBMISSION_Group_Lasso_Paper_Functions.R
####################

Contains all the main functions for data generation, two-stage estimation algorithm, sigma^2 maximum likelihood estimation, performance metrics calculation. Code is very well-documented.


####################
##  MARCH_2019_SUBMISSION_Section_3_TwoStage_Simulation.R
####################

Executive file for the simulation study leading to Table 1 from Section 3 (and Table 1 in supplement as well). Outputs some results in format for Latex tables (FP,FN,Matthews,NFD for full estimate, common and individual components, algorithm convergence metrics), and saves estimates into a folder as .rds files. You can customize such parameters as: number of variables per subject, number of time points, number of subjects per group, etc. Code can be run for VAR(D) models of general lag order D, but has been primarily tested for the VAR(1) models. Whole description can be found in the code comments.


####################
##  MARCH_2019_SUBMISSION_Section_3_Comparative_Simulation_Study.R
####################

Executive file for the comparative simulation study leading to Table 2 from Section 3 (and Table 1 in supplement as well). In particular, for such methods as Sparse Lasso, Group Bridge and Sparse Group Lasso. Outputs some results in format for Latex tables (FP,FN,Matthews,NFD for full estimates), and saves estimates into a folder as .rds files. You can customize such parameters as: number of variables per subject, number of time points, number of subjects per group, etc. Whole description can be found in the code comments.

####################
##  MARCH_2019_SUBMISSION_Section_3_Simulation_Result_Accumulation.R
####################

In case all 50 replicates couldn't be run at once, execution was broken into running a smaller number "actual.rep" at a time (e.g. 10 or 5). This R file compiles the results across multiple such runs. Whole description can be found in the code comments.



####################
## MARCH_2019_SUBMISSION_Section_4_ADHD_Study_Estimation.R
####################

Executive file for the ADHD-200 study from Section 4. Saves resulting estimates (of full transition matrix, of common and individual components) and other characteristics (stability selection proportions, algorithm convergence metrics) into a folder as .rds files. You can customize such parameters as: number of bootstrapped samples to run, which group to estimate, what hard threshold to apply, which tuning parameter selection criterions to use, etc. More details can be found in the code comments.

####################
## MARCH_2019_SUBMISSION_Section_4_ADHD_Output.R
####################

Executive file to obtain graphs and phenotypic summaries corresponding to estimates outputted by 'Group.Lasso.Paper.ADHD.Study.R'.


####################
### ADHD200
####################

Files containing pre-processed brain fMRI time series data (averaged over brain regions according to Multi-Subject Dictionary Learning atlas, MSDL) for ADHD and control patients from ADHD-200 Global Competition (in particular borrowed from Python 'nilearn' package). Additionally it has a folder with a set of ".rds" files
containing the final estimates of joint estimates for this ADHD data, which can be used to output plots via 'Group.Lasso.Paper.ADHD.Output.Plots.R' source code.
