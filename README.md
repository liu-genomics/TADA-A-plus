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

#### CADD
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altA.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altC.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altG.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altT.bed`

#### GERP
`/project2/xinhe/TADA-A-plus/other_annotations/conservation/Example_windows_coding_gerp_gt2.bed`

#### Spidex
> contains all the based that have a spidex score, no
`/project2/xinhe/TADA-A-plus/data/spidex_public_noncommercial_v1_0.tab_alt_A_lower10pct.bed`
`/project2/xinhe/TADA-A-plus/data/spidex_public_noncommercial_v1_0.tab_alt_C_lower10pct.bed`
`/project2/xinhe/TADA-A-plus/data/spidex_public_noncommercial_v1_0.tab_alt_G_lower10pct.bed`
`/project2/xinhe/TADA-A-plus/data/spidex_public_noncommercial_v1_0.tab_alt_T_lower10pct.bed`

#### RBP_Var_Level2 (Based on Fengbiao Mao's paper)
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altA.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altC.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altG.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altT.bed.merge.bed`

### UTR regions

#### CADD greater than 15
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altA_within_updated_UTRs.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altC_within_updated_UTRs.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altG_within_updated_UTRs.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altT_within_updated_UTRs.bed`

#### CLIP
> Contains all the anotated CLIP regions, not only those that are within UTRs

`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIP/human_combine.merged.bed`

#### CLIPdb
> Contains all the anotated CLIPdb regions, not only those that are within UTRs

`/project2/xinhe/TADA-A-plus/other_annotations/coding/CLIPdb/human_combined.merged.bed`

#### GERP
`/project2/xinhe/TADA-A-plus/other_annotations/conservation/Whole_genome.utr_genename.bed.sorted.merged_gerp_gt2.bed`

#### Ribosnitch
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT.bed`

#### Ray_RBP_motif
> Contains all the anotated Ray RBP motifs, not only those that are within UTRs

`/project2/xinhe/TADA-A-plus/other_annotations/coding/Ray_et_al_RBP_motif_hits/UTR.best_hit.hg19.bed.merge.bed.sorted.merged.bed`

#### RBP-Var level 2
> contain all relevant regions, not only including regoins that are within UTRs

`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altA.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altC.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altG.bed.merge.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altT.bed.merge.bed`

### promoter
#### CADD
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altA.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altC.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altG.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/allele_specific_CADD/whole_genome_SNVs_gt15_altT.bed`

#### GERP
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/data/170820_gerp_gt2.all.sorted.merged.bed`

#### Fetal brain DHS E081 DeepSEA annotation, log odds ratio greater than 0.4
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altA_gt_0.4.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altC_gt_0.4.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altG_gt_0.4.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altT_gt_0.4.bed`

#### Fetal brain DHS E081 DeepSEA annotation, log odds ratio smaller than -0.34
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altA_st_-0.34.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altC_st_-0.34.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altG_st_-0.34.bed`
`/project2/xinhe/TADA-A-plus/other_annotations/deepsea_annotations/180510_merged.E081.allele-specific_altT_st_-0.34.bed`

#### Noonan fetal brain H3K27AC plus Adult brain H3K27AC (This is what I used in the AJHG manuscript)
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus//other_annotations/epigenomic_annotation/Noonan_brain_roadmap_union.bed`

#### Noonan hNSC ATAC-Seq Passage 5, rep1
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/Noonan_hNSC_P5-1_peaks.bed`

#### Noonan hNSC ATAC-Seq Passage 5, rep2
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/Noonan_hNSC_P5-2_peaks.bed`

#### Noonan hNSC ATAC-Seq Passage 15-1
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/Noonan_hNSC_P15-1_peaks.bed`

#### Noonan hNSC ATAC-Seq Passage 50
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/Noonan_hNSC_P15-1_peaks.bed`

#### Noonan hNSC ATAC-Seq Passage 50.
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/Noonan_hNSC_50_peaks.bed`

#### Duan Lab iPSC neuron Day 27
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/GSM1820082_d27_fdr0.05.hot.bed`

#### Duan Lab iPSC neuron Day 33
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/GSM1820083_d33_fdr0.05.hot.bed`

#### Duan Lab iPSC neuron Day 41
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus//other_annotations/epigenomic_annotation/GSM1820085_d411_fdr0.05.hot.bed`

#### Duan Lab iPSC neuron Day 42
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/GSM1820087_d412_fdr0.05.hot.bed`

#### ACC_H3K27AC, PsychENCODE
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/171012_Zhongshan_ACC_H3K27AC_merged_picked.bed`

#### DLPFC_ATAC_seq, PsychENCODE
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/171012_Zhongshan_ACC_H3K27AC_merged_picked.bed`

#### DLPFC_H3K27AC, PsychENCODE
> To get annotations specifically in 2kb promoter regions, need to overlap the annotation with 2kb promoter windows
`/project2/xinhe/TADA-A-plus/other_annotations/epigenomic_annotation/171012_Zhongshan_DLPFC_H3K27AC_merged_picked.bed`



