#!/bin/csh
#
#

# Set working directory
cd /home/students/chengsi/Desktop/exposurepred/

# Set local variables
set outfile = 'pred_cohort.Rout'.$$
  
# Executable commands
R CMD BATCH /home/students/chengsi/Desktop/exposurepred/pred_grid.R
