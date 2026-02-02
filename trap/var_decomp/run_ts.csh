#!/bin/csh
#
#

# Set working directory
cd /home/students/chengsi/Desktop/exposurepred/more

# Set local variables
set outfile = 'ts.Rout'.$$
  
# Executable commands
R CMD BATCH /home/students/chengsi/Desktop/exposurepred/more/two_stage.R
