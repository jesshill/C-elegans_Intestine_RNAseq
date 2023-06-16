#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=counts_array-job
#SBATCH --output=%x.%A_%a.log   

# command line:
# sbatch --array=0-23 FeatureCounts_array.sbatch

# use to evaluate jobs
#alias sa='sacct -X --format JobID,JobName,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End'


set -e
module purge

source /curc/sw/anaconda3/latest 
conda activate RNAseq



args=(../03_output/tmpdir/Rep_1_CeMbio_intestine_S39_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_intestine_S41_L003.sam ../03_output/tmpdir/Rep_3_CeMbio_intestine_S43_L003.sam ../03_output/tmpdir/Rep_1_CeMbio_to_PA14_intestine_S51_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_to_PA14_intestine_S53_L003.sam   ../03_output/tmpdir/Rep_3_CeMbio_to_PA14_intestine_S55_L003.sam ../03_output/tmpdir/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003.sam  ../03_output/tmpdir/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003.sam  ../03_output/tmpdir/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003.sam ../03_output/tmpdir/Rep_1_CeMbio_whole_worm_S40_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_whole_worm_S42_L003.sam ../03_output/tmpdir/Rep_3_CeMbio_whole_worm_S44_L003.sam ../03_output/tmpdir/Rep_1_OP50_intestine_S33_L003.sam ../03_output/tmpdir/Rep_2_OP50_intestine_S35_L003.sam ../03_output/tmpdir/Rep_3_OP50_intestine_S37_L003.sam ../03_output/tmpdir/Rep_1_OP50_to_PA14_intestine_S45_L003.sam ../03_output/tmpdir/Rep_2_OP50_to_PA14_intestine_S47_L003.sam ../03_output/tmpdir/Rep_3_OP50_to_PA14_intestine_S49_L003.sam ../03_output/tmpdir/Rep_1_OP50_to_PA14_whole_worm_S46_L003.sam ../03_output/tmpdir/Rep_2_OP50_to_PA14_whole_worm_S48_L003.sam ../03_output/tmpdir/Rep_3_OP50_to_PA14_whole_worm_S50_L003.sam ../03_output/tmpdir/Rep_1_OP50_whole_worm_S34_L003.sam ../03_output/tmpdir/Rep_2_OP50_whole_worm_S36_L003.sam ../03_output/tmpdir/Rep_3_OP50_whole_worm_S38_L003.sam)
sam=${args[ $SLURM_ARRAY_TASK_ID ]}


 
# script

    echo -e $sam
    root=${sam/.sam/}


cmd="featureCounts \
    -p \
    -s 2 \
    -T $SLURM_NTASKS \
    -a ../ce11_IndexBuild/ce11_annotation_ensembl_to_ucsc.gtf \
    -o $root.counts.txt \
    $sam "


echo $cmd
eval $cmd



######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FEATURECOUNTS VERSION:"
featureCounts -v

