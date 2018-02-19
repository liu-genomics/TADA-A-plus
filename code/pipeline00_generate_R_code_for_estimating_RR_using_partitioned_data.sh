# MSSNG data, without accounting for div score in mutation rate calibration
bash ../lib/180215_generate_partitioned_compact_data.sh 180216_RR_estimate_DATA_Yuen_2017_ANNO ../data/denovo_db_Yuen_2017_cases_DNM_SNV_cleaned_with_allele_info.txt \
../data/Example_windows_mutrate_scaling_file_for_Yuen_2017_cases_DNM_cleaned.txt 1609

# MSSNG data, with accounting for div score in mutation rate calibration
bash ../lib/180215_generate_partitioned_compact_data.sh 180216_with_div_adjusted_RR_estimate_DATA_Yuen_2017_ANNO ../data/denovo_db_Yuen_2017_cases_DNM_SNV_cleaned_with_allele_info.txt \
../data/Example_windows_mutrate_with_div_score_scaling_file_for_Yuen_NM2017_cases_DNM_cleaned.txt 1609

# SSC 519 trios data, with accouting for div score in mutation rate calibration
bash ../lib/180215_generate_partitioned_compact_data.sh 180216_with_div_adjusted_RR_estimate_Simons_519_new_version_ANNO ../data/Simons_519_new_version_controls_with_allele_info.txt \
../data/Example_windows_mutrate_with_div_score_scaling_file_for_Simons_519_new_cases_DNM.txt 519