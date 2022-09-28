#!/bin/bash
#PBS -N a3_nln2                   # Set job name to myjob
#PBS -m ae                       # Email status when job completes
#PBS -M xiaoli@mcw.edu       # Email to this address
#PBS -l nodes=1:ppn=16            # Request 1 node with 8 processors
#PBS -l mem=40gb                  # Request 16gb memory
#PBS -l walltime=72:00:00         # Request 1hr job time
#PBS -j oe                       # Join output with error output
#PBS -V

cd ~/BART/a3/codes/                # Change directory to current working directory
module load gcc/9.2
Rscript bart2_nln.R
