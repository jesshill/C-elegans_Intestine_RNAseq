#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=2:00:00 #next time just request 30 min total
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=array-job
#SBATCH --output=%x.%A_%a.log   

# command line:
# sbatch --array=0-23 hisat2_array.sbatch

# use to evaluate jobs
#alias sa='sacct -X --format JobID,JobName,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End'


set -e
module purge

source /curc/sw/anaconda3/latest 
conda activate RNAseq


export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
date # timestamp


args=(../03_output/trimmed_files/Rep_1_CeMbio_intestine_S39_L003_trim_1.fastq.gz           ../03_output/trimmed_files/Rep_2_OP50_intestine_S35_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_CeMbio_to_PA14_intestine_S51_L003_trim_1.fastq.gz   ../03_output/trimmed_files/Rep_2_OP50_to_PA14_intestine_S47_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003_trim_1.fastq.gz  ../03_output/trimmed_files/Rep_2_OP50_to_PA14_whole_worm_S48_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_CeMbio_whole_worm_S40_L003_trim_1.fastq.gz          ../03_output/trimmed_files/Rep_2_OP50_whole_worm_S36_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_OP50_intestine_S33_L003_trim_1.fastq.gz             ../03_output/trimmed_files/Rep_3_CeMbio_intestine_S43_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_OP50_to_PA14_intestine_S45_L003_trim_1.fastq.gz     ../03_output/trimmed_files/Rep_3_CeMbio_to_PA14_intestine_S55_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_OP50_to_PA14_whole_worm_S46_L003_trim_1.fastq.gz    ../03_output/trimmed_files/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_1_OP50_whole_worm_S34_L003_trim_1.fastq.gz            ../03_output/trimmed_files/Rep_3_CeMbio_whole_worm_S44_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_2_CeMbio_intestine_S41_L003_trim_1.fastq.gz           ../03_output/trimmed_files/Rep_3_OP50_intestine_S37_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_2_CeMbio_to_PA14_intestine_S53_L003_trim_1.fastq.gz   ../03_output/trimmed_files/Rep_3_OP50_to_PA14_intestine_S49_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003_trim_1.fastq.gz  ../03_output/trimmed_files/Rep_3_OP50_to_PA14_whole_worm_S50_L003_trim_1.fastq.gz ../03_output/trimmed_files/Rep_2_CeMbio_whole_worm_S42_L003_trim_1.fastq.gz          ../03_output/trimmed_files/Rep_3_OP50_whole_worm_S38_L003_trim_1.fastq.gz)
echo "Doing arg index ${SLURM_ARRAY_TASK_ID:=2}" # testing. use a default value
trim_1=${args[ $SLURM_ARRAY_TASK_ID ]}

 
# script

    echo -e $trim_1
    trim_2=${trim_1/_trim_1/_trim_2}
    root=${trim_1/_trim_1.fastq.gz/}
    root=${root/trimmed_files/tmpdir}
    echo -e $trim_2
    echo -e $root
    
 
cmd="hisat2 \
    --summary-file summary_$SLURM_ARRAY_TASK_ID.txt \
    --no-unal \
    --rna-strandedness RF \
    -p ${SLURM_NTASKS} \
    -x ../ce11_IndexBuild/ce11 \
    -1 ${trim_1} \
    -2 ${trim_2} \
    -S ${root}.sam "

echo $cmd
eval $cmd


