#!/bin/bash
#
#$ -cwd
#$ -q sheppardlab.q
#$ -S /bin/bash

# Set working directory
cd /home/chengsi/

# Set local variables
set outfile = 'pred_cohort.Rout'.$$
  
# Executable commands
R CMD BATCH /home/chengsi/si/spatRF_code_summary/pred_grid.R
