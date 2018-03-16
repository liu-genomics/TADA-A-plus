# script to print a sbatch file for each R file
# Usage sh script.sh [prifix] [mut_files] [mutrate_scaling_file] [sample_size] [nonAS_annotation] [AS_annotation] 
# [prefix] is a direicoty that contains all the partitioned R and sbatch files for one pair of mutation data and annotation data

mkdir -p $1
cd $1

MutratePrefix=${7-../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri}
for i in {00..09} 
do
echo "
### store all data in a compact form for all the genes 
 



library(data.table)
library(parallel)
source(\"../TADA-A/lib/TADA_Annotation.R\")



compact_data_1 <- TADA_A_read_info(mut_files = c(\"$2\"),
                                 window_file = \"../data/windows_partition/Example_windows_with_div_score_coding_$i.with_header.txt\",
                                 mutrate_scaling_files = c(\"$3\"),
                                 sample_sizes = c($4),
                                 gene_prior_file = \"../data/Example_gene_uniform_prior.txt\",
                                 nonAS_noncoding_annotations = c($5),
                                 AS_noncoding_annotations = list($6),
                                 report_proportion = 18665/18665,
                                 chunk_partition_num =1,
                                 node_n = 1,
                                 mutrate_ref_files = c(\"$MutratePrefix.alt_A.mutrate.bw\",
                      \"$MutratePrefix.alt_C.mutrate.bw\",
                      \"$MutratePrefix.alt_G.mutrate.bw\",
                      \"$MutratePrefix.alt_T.mutrate.bw\"),
			MPI = $i
)

saveRDS(compact_data_1,\"$1/$1_$i.RDS\")
quit(\"no\")

" > $1_$i.R
done

for i in *.R
do 
sh ../171116_print_sbatch_24hr_1node_6gb.sh $i $1
done

cd ..
