#!/bin/csh
#
#

# Set working directory
cd /home/users/chengsi/Desktop/exposurepred/spatRF/

# Set local variables
set outfile = 'replicate2.Rout'.$$
  
# Executable commands
R-3.6.1 CMD BATCH /home/users/chengsi/Desktop/exposurepred/spatRF/replicate2.R
