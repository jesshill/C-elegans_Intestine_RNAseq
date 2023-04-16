# C. elegans Intestine RNAseq

RNAseq of hand dissected intestine sections

Libraries were prepared with NEBNext Ultra II directional RNA library prep kit, using UMI adapters

---

### Scripts for RNAseq pipeline on Alpine

- fastp.sbatch (need to convert into an array job!)
- hisat2_array.sbatch
- FeatureCounts_array.sbatch

--- 

![](pipeline.png)
