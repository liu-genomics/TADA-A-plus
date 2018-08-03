# Welcome to TADA-A-plus
## Introduction
TADA-A-plus is a extended version of the original TADA-A, which was published in the American Journal of Human Genetics. There are a few significant improvements that are hopefully of general interst to the community studying de novo mutations. These improvements include a mixture binoimal model-based computational framework to estimate the realtive risk of annotations, a efficient pipeline to simulate DNMs given a sample size, a set of base-level mutation rate files and annotations with certain relative risk.

## Relevant annotations.
All of the annotations that are used in TADA-A and TADA-A-plus have a bed file format. There are two major types, one is allele-specific annotations, the other is non-allele-specific annotations. For the former type, each annotation has only one file in the BED format. The later, each has 4 files, each corresponding to as if the mutant allele is A, T, C, or G.

### coding Annotations
#### Nonfunctional synonymous mutations
`/project2/xinhe/TADA-A-plus/other_annotations/coding171029_synonymous_SNV_altA.bed.merge_removing_spidex_lower10pct.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding171029_synonymous_SNV_altC.bed.merge_removing_spidex_lower10pct.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding171029_synonymous_SNV_altG.bed.merge_removing_spidex_lower10pct.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding171029_synonymous_SNV_altT.bed.merge_removing_spidex_lower10pct.bed`
