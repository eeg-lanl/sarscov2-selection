#!/bin/bash -l
#
#SBATCH --job-name=profile-lik
#SBATCH --output=profile-lik_log.%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=64
#
#SBATCH --array=0-50
#

module purge
module load gcc/9.3.0
## gsl installed in home folder
export LD_LIBRARY_PATH=$HOME/gsl/lib:$LD_LIBRARY_PATH
#make -f loc.makefile clean
#make -f loc.makefile

ARGS=($(seq -s ' ' 0.1 0.01 0.6))

start=$SECONDS

./estavoir --threads $SLURM_CPUS_PER_TASK \
    --seed 66867 \
    --particles 10000 \
    --duplicates 10 \
    --iter 300 \
    --fix-par sigma=${ARGS[$SLURM_ARRAY_TASK_ID]} \
    --data-file data/in/data.tsv \
    --param-file params.txt \
    --id test \
    --addl-data-file addl-data.txt

duration=$((SECONDS - start))
echo "run time:" $duration
