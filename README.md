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

#### Nonsynonymous mutations
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_nonsynonymous_SNV_altA.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_nonsynonymous_SNV_altC.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_nonsynonymous_SNV_altG.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_nonsynonymous_SNV_altT.bed.merge.bed`

#### LoF mutations
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altA.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altC.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altG.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altT.bed.merge.bed`

#### Mis3 mutations (Polyphen-2 probably damaging mutations)
`/project2/xinhe/TADA-A-plus/other_annotations/coding/Polyphen_HDIV_probably_damaging_altA.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/Polyphen_HDIV_probably_damaging_altC.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/Polyphen_HDIV_probably_damaging_altG.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/Polyphen_HDIV_probably_damaging_altT.bed.merge.bed`

#### CLIP
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT_in_coding_windows.bed`

#### CLIPdb
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT_in_coding_windows.bed`

#### MPC score
`/project2/xinhe/TADA-A-plus/other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altA.sorted.merged_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altC.sorted.merged_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altG.sorted.merged_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altT.sorted.merged_in_coding_windows.bed`

#### Ribosnitch
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG_in_coding_windows.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT_in_coding_windows.bed`

#### RBP-VarDB all RBP binding sites
> All RBP binding sites merged, possibly including sites that are in UTR regions in addtional to those that are in coding regions. 
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge.bed`

#### Ray et al RBP motif hits
`/project2/xinhe/TADA-A-plus/other_annotations/coding/Ray_et_al_RBP_motif_hits/A_Up_Down.best_hit.hg19.bed.merge_in_coding_windows.bed`



#### 
