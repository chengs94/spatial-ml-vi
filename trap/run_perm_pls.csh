#!/bin/csh
#
#

# Set working directory
cd /home/students/chengsi/Desktop/exposurepred/

# Set local variables
set outfile = 'pls_imp.Rout'.$$
  
# Executable commands
R CMD BATCH /home/students/chengsi/Desktop/exposurepred/pls_imp.R
