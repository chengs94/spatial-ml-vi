#!/bin/csh
#
#

# Set working directory
cd /home/students/chengsi/Desktop/exposurepred/

# Set local variables
set outfile = 'trap_cv_2.Rout'.$$
  
# Executable commands
R CMD BATCH /home/students/chengsi/Desktop/exposurepred/model2_fewer.R
