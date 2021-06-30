#!/bin/bash -l
#
#SBATCH --job-name=profile-lik
#SBATCH --output=profile-lik_log.%A_%a.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=64
#
#SBATCH --array=0-154
#

## 26 * 5 - 1 = 129

module purge
module load gcc/9.3.0
## gsl installed in home folder
export LD_LIBRARY_PATH=$HOME/gsl/lib:$LD_LIBRARY_PATH
#make -f loc.makefile clean
#make -f loc.makefile


NUM_SIGMA=26
SIGMA_MIN=0.1
SIGMA_DX=0.02
TMAX_MIN=100
TMAX_DX=14

IDX1=$(($SLURM_ARRAY_TASK_ID%$NUM_SIGMA))
IDX2=$(($SLURM_ARRAY_TASK_ID/$NUM_SIGMA))
        
SIGMA=$(awk "BEGIN {print ($SIGMA_MIN + $SIGMA_DX * $IDX1)}")
TMAX=$(awk "BEGIN {print ($TMAX_MIN + $TMAX_DX * $IDX2)}")


start=$SECONDS

./estavoir --threads $SLURM_CPUS_PER_TASK \
    --seed 66867 \
    --particles 10000 \
    --duplicates 10 \
    --iter 300 \
    --fix-par sigma=$SIGMA \
    --data-file data/in/data.tsv \
    --param-file params.txt \
    --end-time $TMAX \
    --id myid_tmax=$TMAX \
    --addl-data-file addl-data.txt

duration=$((SECONDS - start))
echo "run time:" $duration
