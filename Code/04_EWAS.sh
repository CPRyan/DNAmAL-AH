#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=225g           # Keep at 225 GB, which is already quite high
#SBATCH -n 4                 # 4 CPUs is reasonable, R will use BLAS if enabled
#SBATCH -t 10-0:00           # Up to 10 days
#SBATCH --job-name=EWAS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cpr2139@cumc.columbia.edu
#SBATCH -o Output/Out/EWAS_%j.out
#SBATCH --chdir=/ifs/sec/cpc/addhealth/users/cpr/AL_clock

# Go to working directory
cd /ifs/sec/cpc/addhealth/users/cpr/AL_clock

# Prevent stack overflow (default is often 8192 kB)
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




# Optional: allocate a large tempdir in memory (if tmpfs is large enough)
# This avoids I/O bottlenecks and can reduce RAM disk pressure
export TMPDIR=/tmp/$USER/ewas_$SLURM_JOB_ID
mkdir -p $TMPDIR

# Optional: Print some system info
echo "Starting job on $(hostname) with $SLURM_MEM_PER_NODE MB memory and $SLURM_CPUS_ON_NODE CPUs"
echo "Temporary directory: $TMPDIR"

# Run the R script and append stdout/stderr to Rout
Rscript Code/04_EWAS.R >> Output/Out/04_EWAS.Rout

# Cleanup tempdir
rm -rf $TMPDIR


