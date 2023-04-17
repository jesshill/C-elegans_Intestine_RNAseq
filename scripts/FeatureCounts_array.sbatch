#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=2:00:00
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

args=(../03_output/tmpdir/Rep_1_CeMbio_intestine_S39_L003.sam ../03_output/tmpdir/Rep_2_CeMbio_intestine_S41_L003.sam ../03_output/tmpdir/Rep_3_CeMbio_intestine_S43_L003.sam ../03_outpu$
echo "Doing arg index ${SLURM_ARRAY_TASK_ID:=3}" # testing. use a default value
sam=${args[ $SLURM_ARRAY_TASK_ID ]}
 
# script

    echo -e $sam

cmd="featureCounts \                                                                                                                                                                      
    -p \                                                                                                                                                                                  
    -s 2 \                                                                                                                                                                                
    -T $SLURM_NTASKS \                                                                                                                                                                    
    -a ../ce11_IndexBuild/ce11_annotation_ensembl_to_ucsc.gtf \                                                                                                                           
    -o ../03_output/counts/counts.txt \                                                                                                                              
    $sam "                                                                                                                                                                                
                                                                                                                                                                                          
                                                                                                                                                                                          
echo $cmd
eval $cmd
