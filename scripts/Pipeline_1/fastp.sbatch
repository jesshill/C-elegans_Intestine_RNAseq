#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=06:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=fastp
#SBATCH --output=%x.%j.log

 
set -e

module purge

source /curc/sw/anaconda3/latest 
conda activate RNAseq


# script

projdir="/projects/jesshill@colostate.edu/RNAseq_pipeline_celegans/"


for R1 in $projdir/01_input/*_R1_001.fastq.gz
do
    echo -e $R1
    R2=${R1/_R1_/_R2_}
    root=${R1/_R1_001.fastq.gz/}
    root=${root/01_input/03_output}
    echo -e $R2
    echo -e $root

    
    fastp -i ${R1} \
    -I ${R2} \
    -o ${root}_trim_1.fastq \
    -O ${root}_trim_2.fastq \
    -h ${root}_report.html \
    -j ${root}_report.json \
    --detect_adapter_for_pe \
    -U \
    --umi_loc=index1 \
    --thread $SLURM_NTASKS \
    -g -x -p


gzip ${root}_trim_1.fastq
gzip ${root}_trim_2.fastq


done


######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASTP VERSION:"
fastp --version
