#!/bin/bash -l
#
#SBATCH --job-name=estavoir
#SBATCH --output=estavoir_log.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=128
#

module purge
module load gcc/9.3.0
## gsl installed in home directory
export LD_LIBRARY_PATH=$HOME/gsl/lib:$LD_LIBRARY_PATH
#make -f loc.makefile clean
#make -f loc.makefile

start=$SECONDS

./estavoir --threads $SLURM_CPUS_PER_TASK \
--seed 76867 \
--particles 10000 \
--iter 300 \
-p params.txt \
-d data/in/data.tsv \
-G 100 \
-m filter \
--id test \
-x addl-data.txt

duration=$((SECONDS - start))
echo "run time:" $duration
