#!/bin/bash
#PBS -N a3_penaft                    # Set job name to myjob
#PBS -m ae                       # Email status when job completes
#PBS -M xiaoli@mcw.edu       # Email to this address
#PBS -l nodes=1:ppn=8            # Request 1 node with 8 processors
#PBS -l mem=16gb                  # Request 16gb memory
#PBS -l walltime=768:00:00         # Request 1hr job time
#PBS -j oe                       # Join output with error output
#PBS -V

cd ~/BART/a3/codes                # Change directory to current working directory

module load gcc/9.2
Rscript real_penaft.R
