#!/bin/csh
#
#

# Set working directory
cd /home/users/chengsi/Desktop/exposurepred/

# Set local variables
set outfile = 'mm.Rout'.$$
  
# Executable commands
R-3.6.1 CMD BATCH /home/users/chengsi/Desktop/exposurepred/model.R
