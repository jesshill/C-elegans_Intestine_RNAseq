#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=1:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=slurm-BWmerge
#SBATCH --output=%x.%j.log
 
set -e
module purge

source /curc/sw/anaconda3/latest 
conda activate RNAseq
 

# Use bigwigMerge
# Usage: bigWigMerge in1.bw in2.bw .. inN.bw out.bedGraph
 
bigWigMerge ../03_output/04_samtools/Rep_1_CeMbio_intestine_S39_L003_sort.bw           ../03_output/04_samtools/Rep_2_CeMbio_intestine_S41_L003_sort.bw           ../03_output/04_samtools/Rep_3_CeMbio_intestine_S43_L003_sort.bw ../03_output/04_samtools/CeMbio_intestine.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_CeMbio_to_PA14_intestine_S51_L003_sort.bw   ../03_output/04_samtools/Rep_2_CeMbio_to_PA14_intestine_S53_L003_sort.bw   ../03_output/04_samtools/Rep_3_CeMbio_to_PA14_intestine_S55_L003_sort.bw ../03_output/04_samtools/CeMbio_to_PA14_intestine.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_CeMbio_to_PA14_whole_worm_S52_L003_sort.bw  ../03_output/04_samtools/Rep_2_CeMbio_to_PA14_whole_worm_S54_L003_sort.bw  ../03_output/04_samtools/Rep_3_CeMbio_to_PA14_whole_worm_S56_L003_sort.bw ../03_output/04_samtools/CeMbio_to_PA14_whole_worm.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_CeMbio_whole_worm_S40_L003_sort.bw          ../03_output/04_samtools/Rep_2_CeMbio_whole_worm_S42_L003_sort.bw          ../03_output/04_samtools/Rep_3_CeMbio_whole_worm_S44_L003_sort.bw ../03_output/04_samtools/CeMbio_whole_worm.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_OP50_intestine_S33_L003_sort.bw             ../03_output/04_samtools/Rep_2_OP50_intestine_S35_L003_sort.bw             ../03_output/04_samtools/Rep_3_OP50_intestine_S37_L003_sort.bw ../03_output/04_samtools/OP50_intestine.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_OP50_to_PA14_intestine_S45_L003_sort.bw     ../03_output/04_samtools/Rep_2_OP50_to_PA14_intestine_S47_L003_sort.bw     ../03_output/04_samtools/Rep_3_OP50_to_PA14_intestine_S49_L003_sort.bw ../03_output/04_samtools/OP50_to_PA14_intestine.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_OP50_to_PA14_whole_worm_S46_L003_sort.bw    ../03_output/04_samtools/Rep_2_OP50_to_PA14_whole_worm_S48_L003_sort.bw    ../03_output/04_samtools/Rep_3_OP50_to_PA14_whole_worm_S50_L003_sort.bw ../03_output/04_samtools/OP50_to_PA14_whole_worm.bedGraph
bigWigMerge ../03_output/04_samtools/Rep_1_OP50_whole_worm_S34_L003_sort.bw           ../03_output/04_samtools/Rep_2_OP50_whole_worm_S36_L003_sort.bw            ../03_output/04_samtools/Rep_3_OP50_whole_worm_S38_L003_sort.bw ../03_output/04_samtools/OP50_whole_worm.bedGraph

# Use bedgraphtobigwig
# Usage: bedGraphToBigWig in.bedGraph chrom.sizes out.bw
 
bedGraphToBigWig ../03_output/04_samtools/CeMbio_intestine.bedGraph ce11.chrom.sizes ../03_output/04_samtools/CeMbio_intestine.merged.bw
bedGraphToBigWig ../03_output/04_samtools/CeMbio_to_PA14_intestine.bedGraph ce11.chrom.sizes ../03_output/04_samtools/CeMbio_to_PA14_intestine.merged.bw
bedGraphToBigWig ../03_output/04_samtools/CeMbio_to_PA14_whole_worm.bedGraph ce11.chrom.sizes ../03_output/04_samtools/CeMbio_to_PA14_whole_worm.merged.bw
bedGraphToBigWig ../03_output/04_samtools/CeMbio_whole_worm.bedGraph ce11.chrom.sizes ../03_output/04_samtools/CeMbio_whole_worm.merged.bw
bedGraphToBigWig ../03_output/04_samtools/OP50_intestine.bedGraph ce11.chrom.sizes ../03_output/04_samtools/OP50_intestine.merged.bw
bedGraphToBigWig ../03_output/04_samtools/OP50_to_PA14_intestine.bedGraph ce11.chrom.sizes ../03_output/04_samtools/OP50_to_PA14_intestine.merged.bw
bedGraphToBigWig ../03_output/04_samtools/OP50_to_PA14_whole_worm.bedGraph ce11.chrom.sizes ../03_output/04_samtools/OP50_to_PA14_whole_worm.merged.bw
bedGraphToBigWig ../03_output/04_samtools/OP50_whole_worm.bedGraph ce11.chrom.sizes ../03_output/04_samtools/OP50_whole_worm.merged.bw
