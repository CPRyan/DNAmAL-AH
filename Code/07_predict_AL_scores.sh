#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=100g  # Adjust if needed for memory-intensive tasks
#SBATCH -n 4        # Adjust based on parallelism in R
#SBATCH -t 10-0:00
#SBATCH --job-name=predict_al_scores
#SBATCH --mail-type=END,FAIL  # Get email notifications
#SBATCH --mail-user=cpr2139@cumc.columbia.edu
#SBATCH -o Output/Out/predict_AL_scores_%j.out  # SLURM log file goes to Output/Out
#SBATCH --chdir=/ifs/sec/cpc/addhealth/users/cpr/AL_clock  # Set working directory

cd /ifs/sec/cpc/addhealth/users/cpr/AL_clock

# Sets the maximium stack size a process can use - might help fix crashing
ulimit -s unlimited


# Needed to get working on new R
module purge
module load gcc/12.2.0
module load r/4.4.0

export R_LIBS_USER="/ifs/sec/cpc/addhealth/users/cpr/R/x86_64-pc-linux-gnu-library/4.4"
mkdir -p "$R_LIBS_USER"
echo "R path: $(which R)"
R --version
echo "g++: $(which g++)"
g++ --version
echo "libstdc++: $(g++ -print-file-name=libstdc++.so.6)"
strings "$(g++ -print-file-name=libstdc++.so.6)" | grep -F GLIBCXX_3.4.30


# Run the R script from Code/ and store output in Output/Out/
Rscript Code/07_predict_AL_scores.R >> Output/Out/07_predict_AL_scores.Rout