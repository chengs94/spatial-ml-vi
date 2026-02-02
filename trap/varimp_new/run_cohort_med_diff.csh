#!/bin/csh
#
#

# Set working directory
cd /home/users/chengsi/Desktop/exposurepred/

# Set local variables
set outfile = 'pred_cohort_med_diff.Rout'.$$
  
# Executable commands
R CMD BATCH /home/users/chengsi/Desktop/exposurepred/pred_cohort.R
