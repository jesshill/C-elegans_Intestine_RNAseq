# C. elegans Intestine RNAseq

RNAseq of hand dissected intestine sections from L4 stage N2 worms.

The NEBNext Ultra II directional RNA library prep kit (E7760S) was used in combination with the NEBNext Multiplex Oligos for illumina (UMI adaptors RNA set 1, NEB #E7416)

This kit contains dUTP in the second strand synthesis buffer that allows labeling of the second strand cDNA and subsequent excision with USER Enzyme. Leaving behind only the first strand. 

Paired end sequencing was performed
This kit produces

---

### Scripts for RNAseq pipeline on Alpine

- fastp.sbatch (need to convert into an array job!) 
  > Specify UMI's in the index here

- hisat2_array.sbatch 
  > Specify the RNA strandedness here using 
  > ...

- FeatureCounts_array.sbatch
  > Specify 

--- 

![](pipeline.png)
