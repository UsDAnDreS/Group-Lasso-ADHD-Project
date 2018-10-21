# Group-Lasso-ADHD-Project

All the files can be found on public access repository:

https://github.com/UsDAnDreS/Group-Lasso-ADHD-Project


####################
##  Group.Lasso.Paper.Functions.R
####################

Contains all the main functions for data generation, two-stage estimation algorithm, sigma^2 maximum likelihood estimation, performance metrics calculation. Code is very well-documented.


####################
##  Group.Lasso.Paper.Simulation.R
####################

Executive file for the simulation study from Section 3. Outputs some results in format for Latex tables (FP,FN,Matthews for full estimate, common and individual components, algorithm convergence metrics), and saves estimates into a folder as .rds files. You can customize such parameters as: number of variables per subject, number of time points, number of subjects per group, etc. Code can be run for VAR(D) models of general lag order D, but has been primarily tested for the VAR(1) models. Whole description can be found in the code comments.

####################
## Group.Lasso.Paper.ADHD.Study.R
####################

Executive file for the ADHD-200 study from Section 4. Saves resulting estimates (of full transition matrix, of common and individual components) and other characteristics (stability selection proportions, algorithm convergence metrics) into a folder as .rds files. You can customize such parameters as: number of bootstrapped samples to run, which group to estimate, what hard threshold to apply, which tuning parameter selection criterions to use, etc. More details can be found in the code comments.

####################
## Group.Lasso.Paper.ADHD.Output.Plots.R
####################

Executive file to obtain graphs and phenotypic summaries corresponding to estimates outputted by 'Group.Lasso.Paper.ADHD.Study.R'.


####################
### ADHD200
####################

Files containing pre-processed brain fMRI time series data (averaged over brain regions according to Multi-Subject Dictionary Learning atlas, MSDL) for ADHD and control patients from ADHD-200 Global Competition (in particular borrowed from Python 'nilearn' package).

