#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=00:10:00
#SBATCH --partition=amilan
#SBATCH --job-name=slurm
#SBATCH --output=%x.%j.log  #give slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@"  #log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date  #timestamp

# useful vars
# SLURM_JOB_ID, unset if invoked with: bash scriptname.sbatch
# SLURM_NTASKS, unset if invoked with: bash scriptname.sbatch
cmd="YOUR CMD HERE"
echo $CMD

# run only if invoked with sbatch (useful in testing)
if [ -n "$SLURM_JOB_ID" ]
then
        time eval $cmd
fi
