#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=array-job
#SBATCH --output=%x.%A_%a.log  #like: array-job.ID_1.log

#command line usage: sbatch --array=0-3 jobscript_array.sbatch arg1 arg2 arg3 

#Example: One gzip command for each fastq file.
# 1) copy this script to gzip_array.sbatch, make 'gzip' the command.
# 2) see how many fastq files there are with something like: echo *.fastq | wc -w
# 3) Optionally test if it "sees" the command line args: bash gzip_array.sbatch *.fastq
# 4) Say the answer from (2) is 9. Now run the script with sbatch --array=0-8 gzip_array.sbatch *.fastq

export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
echo "[$0] $SLURM_JOB_NAME $@"  #log the command line
date # timestamp


args=( $@ )  #array of all args
echo "Doing arg index ${SLURM_ARRAY_TASK_ID:=3}"  #testing. use a default value
arg=${args[ $SLURM_ARRAY_TASK_ID ]}

cmd="command $arg"

echo $cmd

if [ -n "$SLURM_JOB_ID" ]  #run with sbatch. use 'bash' to test
then
    time eval $cmd
fi
