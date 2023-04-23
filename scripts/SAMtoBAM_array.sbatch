#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=BW_array-job
#SBATCH --output=%x.%A_%a.log   

# command line:
# sbatch --array=0-23 SAMtoBAM_array.sbatch

# use to evaluate jobs
#alias sa='sacct -X --format JobID,JobName,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End'


set -e
module purge

source /curc/sw/anaconda3/latest 
conda activate RNAseq


#### script 1
# Samtools: compress .sam -> .bam

#args=(../03_output/tmpdir/Rep_1_CeMbio_intestine_S39_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_intestine_S41_L003.sam ../03_output/tmpdir/Rep_3_CeMbio_intestine_S43_L003.sam ../03_output/tmpdir/Rep_1_CeMbio_to_PA14_intestine_S51_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_to_PA14_intestine_S53_L003.sam   ../03_output/tmpdir/Rep_3_CeMbio_to_PA14_intestine_S55_L003.sam ../03_output/tmpdir/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003.sam  ../03_output/tmpdir/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003.sam  ../03_output/tmpdir/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003.sam ../03_output/tmpdir/Rep_1_CeMbio_whole_worm_S40_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_whole_worm_S42_L003.sam ../03_output/tmpdir/Rep_3_CeMbio_whole_worm_S44_L003.sam ../03_output/tmpdir/Rep_1_OP50_intestine_S33_L003.sam ../03_output/tmpdir/Rep_2_OP50_intestine_S35_L003.sam ../03_output/tmpdir/Rep_3_OP50_intestine_S37_L003.sam ../03_output/tmpdir/Rep_1_OP50_to_PA14_intestine_S45_L003.sam ../03_output/tmpdir/Rep_2_OP50_to_PA14_intestine_S47_L003.sam ../03_output/tmpdir/Rep_3_OP50_to_PA14_intestine_S49_L003.sam ../03_output/tmpdir/Rep_1_OP50_to_PA14_whole_worm_S46_L003.sam ../03_output/tmpdir/Rep_2_OP50_to_PA14_whole_worm_S48_L003.sam ../03_output/tmpdir/Rep_3_OP50_to_PA14_whole_worm_S50_L003.sam ../03_output/tmpdir/Rep_1_OP50_whole_worm_S34_L003.sam ../03_output/tmpdir/Rep_2_OP50_whole_worm_S36_L003.sam ../03_output/tmpdir/Rep_3_OP50_whole_worm_S38_L003.sam)
#sam=${args[ $SLURM_ARRAY_TASK_ID ]}

#    echo -e $sam

#cmd1="samtools view \
#    -S \
#    -b \
#    $sam > ${sam/.sam/}.bam "

#echo $cmd1
#eval $cmd1





#### script 2
# Samtools: sort .bam -> _sort.bam

#args=(../03_output/samtools/Rep_1_CeMbio_intestine_S39_L003.bam ../03_output/samtools/Rep_1_OP50_to_PA14_whole_worm_S46_L003.bam ../03_output/samtools/Rep_2_OP50_intestine_S35_L003.bam ../03_output/samtools/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003.bam ../03_output/samtools/Rep_1_CeMbio_to_PA14_intestine_S51_L003.bam ../03_output/samtools/Rep_1_OP50_whole_worm_S34_L003.bam ../03_output/samtools/Rep_2_OP50_to_PA14_intestine_S47_L003.bam ../03_output/samtools/Rep_3_CeMbio_whole_worm_S44_L003.bam ../03_output/samtools/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003.bam ../03_output/samtools/Rep_2_CeMbio_intestine_S41_L003.bam ../03_output/samtools/Rep_2_OP50_to_PA14_whole_worm_S48_L003.bam ../03_output/samtools/Rep_3_OP50_intestine_S37_L003.bam ../03_output/samtools/Rep_1_CeMbio_whole_worm_S40_L003.bam ../03_output/samtools/Rep_2_CeMbio_to_PA14_intestine_S53_L003.bam ../03_output/samtools/Rep_2_OP50_whole_worm_S36_L003.bam ../03_output/samtools/Rep_3_OP50_to_PA14_intestine_S49_L003.bam ../03_output/samtools/Rep_1_OP50_intestine_S33_L003.bam ../03_output/samtools/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003.bam ../03_output/samtools/Rep_3_CeMbio_intestine_S43_L003.bam ../03_output/samtools/Rep_3_OP50_to_PA14_whole_worm_S50_L003.bam ../03_output/samtools/Rep_1_OP50_to_PA14_intestine_S45_L003.bam ../03_output/samtools/Rep_2_CeMbio_whole_worm_S42_L003.bam ../03_output/samtools/Rep_3_CeMbio_to_PA14_intestine_S55_L003.bam ../03_output/samtools/Rep_3_OP50_whole_worm_S38_L003.bam)
#bam=${args[ $SLURM_ARRAY_TASK_ID ]}

#    echo -e $bam 

#cmd2="samtools sort \
#    $bam \
#    -o ${bam/.bam/}.sorted.bam "

#echo $cmd2
#eval $cmd2




#### script 3
# Samtools: index _sort.bam -> _sort.bam.bai

#args=(../03_output/samtools/Rep_1_CeMbio_intestine_S39_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_to_PA14_whole_worm_S46_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_intestine_S35_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003.sorted.bam ../03_output/samtools/Rep_1_CeMbio_to_PA14_intestine_S51_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_whole_worm_S34_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_to_PA14_intestine_S47_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_whole_worm_S44_L003.sorted.bam ../03_output/samtools/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_intestine_S41_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_to_PA14_whole_worm_S48_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_intestine_S37_L003.sorted.bam ../03_output/samtools/Rep_1_CeMbio_whole_worm_S40_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_to_PA14_intestine_S53_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_whole_worm_S36_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_to_PA14_intestine_S49_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_intestine_S33_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_intestine_S43_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_to_PA14_whole_worm_S50_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_to_PA14_intestine_S45_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_whole_worm_S42_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_to_PA14_intestine_S55_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_whole_worm_S38_L003.sorted.bam)
#bai=${args[ $SLURM_ARRAY_TASK_ID ]}  

#    echo -e bai 

#cmd3="samtools index \
#    $bai "

#echo $cmd3
#time eval $cmd3



#### script 4
# bamCoverage: Create a .bw file that is normalized. This can be uploaded to IGV or UCSC


args=(../03_output/samtools/Rep_1_CeMbio_intestine_S39_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_to_PA14_whole_worm_S46_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_intestine_S35_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003.sorted.bam ../03_output/samtools/Rep_1_CeMbio_to_PA14_intestine_S51_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_whole_worm_S34_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_to_PA14_intestine_S47_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_whole_worm_S44_L003.sorted.bam ../03_output/samtools/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_intestine_S41_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_to_PA14_whole_worm_S48_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_intestine_S37_L003.sorted.bam ../03_output/samtools/Rep_1_CeMbio_whole_worm_S40_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_to_PA14_intestine_S53_L003.sorted.bam ../03_output/samtools/Rep_2_OP50_whole_worm_S36_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_to_PA14_intestine_S49_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_intestine_S33_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_intestine_S43_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_to_PA14_whole_worm_S50_L003.sorted.bam ../03_output/samtools/Rep_1_OP50_to_PA14_intestine_S45_L003.sorted.bam ../03_output/samtools/Rep_2_CeMbio_whole_worm_S42_L003.sorted.bam ../03_output/samtools/Rep_3_CeMbio_to_PA14_intestine_S55_L003.sorted.bam ../03_output/samtools/Rep_3_OP50_whole_worm_S38_L003.sorted.bam)
bw=${args[ $SLURM_ARRAY_TASK_ID ]}

    echo -e bw

cmd4="bamCoverage \
    -b $bw
    -o ${bw/.sorted.bam/}_sort.bw
    --outFileFormat bigwig
    -p $SLURM_NTASKS
    --normalizeUsing CPM
    --binSize 1
    --filterRNAstrand=forward "

echo $cmd4
time eval $cmd4
