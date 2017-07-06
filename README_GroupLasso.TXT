####################
##  Group.Lasso.Paper.Functions.R
####################

Contains all the main functions for data generation, two-stage estimation algorithm, sigma^2 maximum likelihood estimation, performance metrics calculation. Code is very well-documented.


####################
##  Group.Lasso.Paper.Simulation.R
####################

Executive file for the simulation study from Section 3. You can customize such parameters as: number of variables per subject, number of time points, number of subjects per group, etc. Code can be run for VAR(D) models of general lag order D, but has been primarily tested for the VAR(1) models. Whole description can be found in the code comments.

####################
## Group.Lasso.Paper.ADHD.Study.R
####################

Executive file for the ADHD-200 study from Section 4. You can customize such parameters as: number of bootstrapped samples to run, which group to estimate, what hard threshold to apply, which tuning parameter selection criterions to use, etc. More details can be found in the code comments.


####################
### PA.csv, IL.csv, OH.csv, MI.csv
####################

Files containing Federal Reserve (FRED) seasonally-adjusted monthly time series data on 18 economic indicator variables (all described in the paper) from 1976 up until 2016 for each state. Function 'Load.Data', defined in 'Paper.functions.R' file and executed in 'Real.Data.Four.States.R' file, deals with extracting our time series of interest out of these files.

