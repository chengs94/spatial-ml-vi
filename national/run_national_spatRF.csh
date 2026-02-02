#!/bin/csh
#
#

# Set working directory
cd /home/users/chengsi/Desktop/exposurepred/national/

# Set local variables
set outfile = 'national_spatRF.Rout'.$$
  
# Executable commands
R-3.6.1 CMD BATCH /home/users/chengsi/Desktop/exposurepred/national/pred_national_spatRF.R
