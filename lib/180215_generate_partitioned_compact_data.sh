# script to print a sbatch file for each R file
# Usage sh script.sh [prifix] [mut_files] [mutrate_scaling_file] [sample_size]
# [prefix] is a direicoty that contains all the partitioned R and sbatch files for one pair of mutation data and annotation data

#!/usr/bin/env bash

MutratePrefix=${5:-../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri}
WindowPrefix=${6:-../data/windows_partition/Example_windows_with_div_score_coding_}
GenePriorFile=${7:-../data/Example_gene_prior.txt}


bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_CLIPdb $2 \
$3 $4 \"../other_annotations/coding/CLIPdb/human_combined.merged.bed\" NA $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_CLIP $2 \
$3 $4 \"../other_annotations/coding/CLIP/human_combine.merged.bed\" NA $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_GERP $2 \
$3 $4 \"../other_annotations/conservation/Example_windows_coding_gerp_gt2.bed\" NA $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_ribosnitch $2 \
$3 $4 NA c\(\"../other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA.bed\",\"../other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC.bed\",\"../other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG.bed\",\"../other_annotations/coding/ribosnitch/hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_CLIPdb_ribosnitch $2 \
$3 $4 NA c\(\"../other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA.bed\",\"../other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC.bed\",\"../other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG.bed\",\"../other_annotations/coding/CLIPdb/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_CLIP_ribosnitch $2 \
$3 $4 NA c\(\"../other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altA.bed\",\"../other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altC.bed\",\"../other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altG.bed\",\"../other_annotations/coding/CLIP/human_combine.merged_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.05.merged.altT.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_nonfunctional_syn $2 \
$3 $4 NA c\(\"../other_annotations/coding/171029_synonymous_SNV_altA.bed.merge_removing_spidex_lower10pct.bed\",\"../other_annotations/coding/171029_synonymous_SNV_altC.bed.merge_removing_spidex_lower10pct.bed\",\"../other_annotations/coding/171029_synonymous_SNV_altG.bed.merge_removing_spidex_lower10pct.bed\",\"../other_annotations/coding/171029_synonymous_SNV_altT.bed.merge_removing_spidex_lower10pct.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_spidex $2 \
$3 $4 NA c\(\"../data/spidex_public_noncommercial_v1_0.tab_alt_A_lower10pct.bed\",\"../data/spidex_public_noncommercial_v1_0.tab_alt_C_lower10pct.bed\",\"../data/spidex_public_noncommercial_v1_0.tab_alt_G_lower10pct.bed\",\"../data/spidex_public_noncommercial_v1_0.tab_alt_T_lower10pct.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_Mis3 $2 \
$3 $4 NA c\(\"../other_annotations/coding/Polyphen_HDIV_probably_damaging_altA.bed.merge.bed\",\"../other_annotations/coding/Polyphen_HDIV_probably_damaging_altC.bed.merge.bed\",\"../other_annotations/coding/Polyphen_HDIV_probably_damaging_altG.bed.merge.bed\",\"../other_annotations/coding/Polyphen_HDIV_probably_damaging_altT.bed.merge.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_LoF $2 \
$3 $4 NA c\(\"../other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altA.bed.merge.bed\",\"../other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altC.bed.merge.bed\",\"../other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altG.bed.merge.bed\",\"../other_annotations/coding/171121_coding_stop_loss_and_gain_SNV_altT.bed.merge.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_CADD $2 \
$3 $4 NA c\(\"../other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altA.bed\",\"../other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altC.bed\",\"../other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altG.bed\",\"../other_annotations/allele_specific_CADD/Example_windows_coding_SNVs_gt15_altT.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_Ray_RBP_motif $2 \
$3 $4 \"../other_annotations/coding/Ray_et_al_RBP_motif_hits/A_Up_Down.best_hit.hg19.bed.merge.bed\" NA $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_RBP_Var_Level2 $2 \
$3 $4 NA c\(\"../other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altA.bed.merge.bed\",\"../other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altC.bed.merge.bed\",\"../other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altG.bed.merge.bed\",\"../other_annotations/coding/RBP-VarDB/RBP.all.bed.merge_overlap_hg19_refGenes_exons.gtf.lg.transc.fa.RNAsnpM3.bed.abspos.p0.1.altT.bed.merge.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile

bash ../lib/180215_print_partitioned_R_script_for_RR_estimate_nonAS_annotations.sh $1_MPC_gt2 $2 \
$3 $4 NA c\(\"../other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altA.sorted.merged.bed\",\"../other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altC.sorted.merged.bed\",\"../other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altG.sorted.merged.bed\",\"../other_annotations/coding/MPC_score/fordist_constraint_official_mpc_values_v2_MPC_gt2_altT.sorted.merged.bed\"\) $MutratePrefix $WindowPrefix $GenePriorFile
# print a Rmd file to calculate relative risk

echo "
\`\`\`{r}
library(data.table)
library(parallel)
source(\"../TADA-A/lib/TADA_Annotation.R\")
\`\`\`

Use coding window file. The previous scaling factors are still applicable. 


\`\`\`{r, cache=TRUE}
prefix_vector <- data.frame(annotation_name = c(\"CADD\",
                                                \"CLIP\",
                                                \"CLIPdb\",
                                                \"CLIPdb_ribosnitch\",
                                                \"CLIP_ribosnitch\",
                                                \"GERP\",
                                                \"LoF\",
                                                \"Mis3\",
                                                \"Nonfunctional_syn\",
                                                \"ribosnitch\",
                                                \"spidex\",
                                                \"Ray_RBP_motif\",
                                                \"RBP_Var_Level2\",
                                                \"MPC_gt2\"), 
                            prefix = c(\"$1_CADD/$1_CADD_\",
                                       \"$1_CLIP/$1_CLIP_\",
                                       \"$1_CLIPdb/$1_CLIPdb_\",
                                       \"$1_CLIPdb_ribosnitch/$1_CLIPdb_ribosnitch_\",
                                       \"$1_CLIP_ribosnitch/$1_CLIP_ribosnitch_\",
                                       \"$1_GERP/$1_GERP_\",
                                       \"$1_LoF/$1_LoF_\",
                                       \"$1_Mis3/$1_Mis3_\",
                                       \"$1_nonfunctional_syn/$1_nonfunctional_syn_\",
                                       \"$1_ribosnitch/$1_ribosnitch_\",
                                       \"$1_spidex/$1_spidex_\",
                                       \"$1_Ray_RBP_motif/$1_Ray_RBP_motif_\",
                                       \"$1_RBP_Var_Level2/$1_RBP_Var_Level2_\",
                                       \"$1_MPC_gt2/$1_MPC_gt2_\"))

report_rr_from_compact_data_in_partitioned_form <- function(prefix_vector){
  report <- data.frame(annotation = NA, log_RR = NA, lower_bound = NA, upper_bound = NA)
  for(i in 1:length(prefix_vector[,1])){
    data_partition <- list()
    test = sprintf(\"%02d\", 0:9)
    for(j in 1:10){
      temp <- readRDS(paste(prefix_vector[i,2], test[j], \".RDS\", sep = \"\"))
      data_partition <- append(data_partition, temp\$base_info)
    }
  temp <- TADA_A_RR_estimate(data = data_partition, selected_annotations = c(1), gene_prior_file = \"$GenePriorFile\", optimization_iteration = 2000, mode = \"single_fast\")
  print(temp)
  report[i,1] <- as.character(prefix_vector[i,1])
  report[i,2:4] <- temp\$rr_report
  }
  return(report)
}

report <- report_rr_from_compact_data_in_partitioned_form(prefix_vector)
\`\`\`

\`\`\`{r}
knitr::kable(report)
\`\`\`

" > $1_RR_estimate.Rmd
