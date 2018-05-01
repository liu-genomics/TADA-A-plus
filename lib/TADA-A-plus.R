TADA_A_compare_mutrate_model <- function(mut_file, window_file, sample_size, scale_features, validation_chr = "chr1", mutrate_mode = "regular", genes = "all", validation_level = "window", 
                                         mutrate_ref_files = c("../other_annotations/ERV_mutrate/Example_windows_coding_noncoding_UTRs_ERV.mutrate.standardized.alt_A.mutrate.bw",
                                                               "../other_annotations/ERV_mutrate/Example_windows_coding_noncoding_UTRs_ERV.mutrate.standardized.alt_C.mutrate.bw",
                                                               "../other_annotations/ERV_mutrate/Example_windows_coding_noncoding_UTRs_ERV.mutrate.standardized.alt_G.mutrate.bw",
                                                               "../other_annotations/ERV_mutrate/Example_windows_coding_noncoding_UTRs_ERV.mutrate.standardized.alt_T.mutrate.bw"),
                                         node_n = 1, use_features = TRUE){
  # [mut_file] is the file with DNM infomation in a BED format. 0-based start and 1-based end. e.g., "../data/Example_windows.bed"
  # The first 4 columns must be chr, start, end, and site_index of genomic windows. if [validation_level] is set to be "base", Then need to have two additional columns as ref and alt alleles.
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end", "genename", and "mutrate", followed by features that might affect local background mutation rate.
  # [sample_size] is the number of individuals. 
  # [scale_features] is a vector of the names of features that need to be scaled, this is recommended to apply to continuous features, such as GC content, to make easier the interpretation of the effect sizes of features. 
  # [validation_chr] is the chromosome that is left out during model fitting and will be used to calculate likelihood once parameters are estimated.
  # [mutrate_mode], "regular" means the mutation rate of each window is the total sum of mutation rates over all bases. Under this setting, only in very rare cases a window will have mutation rate
  # equal to 0. The output scaling file would remove these windows, as these windows can't be used in the likelihood estimation. "special" means the mutation rate of each window is the rate of 
  # a specific type of mutations. For example, when analyzing WES data, we would want to use syn mutations to adjust for mutation rates. The mutation rate here should accordingly be the rate of 
  # synnonymous mutations. Under this situation, more windows might have mutation rate as 0, but we still need to get the scaling factor for these windows, as the rate of other mutation types may 
  # not be 0, thus the bases in these windows are informative. If [genes] set to be not "all", then must set [mutrate_mode] to "special".
  # [genes] could be "all", or a list of genes. If set to be "all", all windows will be used to adjust mutation rates. Otherwise, only windows with a gene name that 
  # is included in the list of genes specified by [genes] will be used. 
  # [validation_level] if "window", ll will be calculated at the window level for the validation chromosome. if "base", likelihood will be calculated at the base level. In this situation, [mutrate_ref_files], which is the 
  # [node_n] is the number of nodes used to run a certain chunk of the code, default is 1
  # base level mutation rate ref needs to be provided.
  # the output of this function is the loglikelihood of observing data given the estimated parameters from the training model. THe validation is only performed on the validation
  # chromosome indicated by [validation_chr], the training model didn't use any information from the 
  # [use_features] is TRUE or FALSE. if TRUE, then use all the additional columns of features to account for mutation rates, otherwise, don't use any covariates.
  
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  command <- paste("sed -n '1!p' ", window_file, " | awk {'print $1\"\t\"$2\"\t\"$3\"\t\"$4'} > ", paste(prefix, "_temp_windows.bed", sep = ""), sep = "")
  system(command)
  # get the number of SNVs for each window
  command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", mut_file, " -b ", paste(prefix, "_temp_windows.bed", sep = ""), " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
  system(command)
  # read in the file with the number of SNVs for each window
  coverage <- fread(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- coverage[,1:5]
  colnames(coverage) <- c("chr","start","end","site_index","mut_count")
  # read in the window file that has feature annotations that might affect mutation rates
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  feature_number <- dim(windows)[2] - 6
  # merge [window] and [coverage] to link mutation count with feature annotations
  coverage <- coverage[windows[,-1:-3], on="site_index"]
  if(genes == "all"){
    target_genes <- unique(coverage$genename)
  }else{
    target_genes <- fread(genes, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    target_genes <- target_genes[[1]] 
  }
  if(mutrate_mode == "regular"){ 
    coverage <- coverage[mutrate !=0]
    for(i in 1:length(scale_features)){
      coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
    }
    # write the formula
    f <- paste("mut_count ~ ", paste(colnames(coverage)[8 : (8 + feature_number -1)], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
    #fit mutation rate model using a fixed set of features not including histone modification marks. 
    out.offset <- glm(as.formula(f), family = poisson, data = coverage)
    scaling <- data.table(site_index = coverage$site_index, scaling_factor = out.offset$fitted.values / (2 * coverage$mutrate * sample_size))
  } else if(mutrate_mode == "special"){
    for(i in 1:length(scale_features)){
      coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
    }
    coverage_mutrate_gt0_training <- coverage[mutrate !=0 & chr != validation_chr]
    # only keep windows that are assigned to [target_genes]
    coverage_mutrate_gt0_training <- coverage_mutrate_gt0_training[is.element(genename, target_genes)]
    # consider the case where we don't use any features to adust mutation rates, although features might be already provided in the window file.
    if(use_features){
      f <- paste("mut_count ~ ", paste(colnames(coverage_mutrate_gt0_training)[8 : (8 + feature_number -1)], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
    }else{
      f <- paste("mut_count ~ 1", " + offset(log(2*mutrate*sample_size))", sep = "")
    }
    out.offset <- glm(as.formula(f), family = poisson, data = coverage_mutrate_gt0_training)
    if(use_features){
      scaling_factor <- exp(as.matrix(cbind(1, coverage[,8 : (8 + feature_number -1)])) %*% out.offset$coefficients)
    }else{
      scaling_factor <- exp(out.offset$coefficients)
    }
    coverage <- data.table(coverage, scaling_factor = scaling_factor)
    colnames(coverage)[dim(coverage)[2]] <- "scaling_factor"
    coverage$adjusted_mutrate <- 2 * sample_size * coverage$mutrate * coverage$scaling_factor
    coverage_validation <- coverage[mutrate != 0 & chr == validation_chr & is.element(genename, target_genes)]
    if(validation_level == "window"){
      ll <- sum(mapply(dpois, coverage_validation$mut_count, coverage_validation$adjusted_mutrate, MoreArgs = list(log = TRUE)))
    }else if(validation_level == "base"){
      alt_letters <- c("A", "C", "G", "T")
      window_expansion <- function(table_row){
        start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
        data.frame(table_row[1], start, start+1, paste(table_row["chr"],table_row["genename"],start,sep = "_"), as.numeric(table_row["scaling_factor"]))
      }
      if(node_n != 1){
        coverage_validation_temp <- rbindlist(parApply(cl, coverage_validation, 1, window_expansion))
      }else{
        coverage_validation_temp <- rbindlist(apply(coverage_validation, 1, window_expansion))
      }
      system("echo \"Finished window expansion at the validated chromosome\"")
      system("date")
      colnames(coverage_validation_temp) <- c("chr","start","end","base_ID", "scaling_factor")  
      coverage_validation_temp <- coverage_validation_temp[!duplicated(base_ID)]
      coverage_validation_temp$start <- as.integer(coverage_validation_temp$start)
      coverage_validation_temp$end <- as.integer(coverage_validation_temp$end)
      # write out a bed file to get base-level mutation rates
      fwrite(coverage_validation_temp[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
      # read in allele-specific base-level mutation rate, and record mutations
      # initialize log-likelihood
      ll <- 0
      mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      for(j in 1:length(mutrate_ref_files)){
        command <- paste("../external_tools/bigWigAverageOverBed ", mutrate_ref_files[j], " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
        system(command)
        command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
        system(command)
        coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
        colnames(coverage_noncoding_base_mutrate) <- c("base_ID", "base_mutrate_alt")
        coverage_validation_temp <- coverage_validation_temp[coverage_noncoding_base_mutrate, on = "base_ID"]
        system(paste("echo \"Finished obtaining base-level mutation rate for alt allele ", " now ", ".\"", sep = ""))
        system("date")
        
        mut_allele <- mut[V5 == alt_letters[j]]
        # the subsequent two lines of code are used to prevent a outputing bug when using fwrite. (85000000 to be written as 8.5e7)
        mut_allele$V2 <- as.integer(mut_allele$V2)
        mut_allele$V3 <- as.integer(mut_allele$V3)
        fwrite(mut_allele, paste(prefix, "_temp_mut_allele.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", paste(prefix, "_temp_mut_allele.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""),sep = "")
        system(command)
        base_with_mut <- fread(paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_with_mut <- base_with_mut[,c("V4","V5"), with = FALSE]
        colnames(base_with_mut) <- c("base_ID", "mut_count")
        coverage_validation_temp <- coverage_validation_temp[base_with_mut, on = "base_ID"]
        
        # add to ll
        ll <- ll + sum(mapply(dpois, coverage_validation_temp[base_mutrate_alt != 0]$mut_count, 2 * sample_size * coverage_validation_temp[base_mutrate_alt != 0]$scaling_factor * coverage_validation_temp[base_mutrate_alt != 0]$base_mutrate_alt, list(log = TRUE)))
        coverage_validation_temp <- coverage_validation_temp[,1:5]
      }
    }
  }
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste("echo \"Temp files cleaned and data recording finished!\""))
  system("date")
  return(list(loglikelihood = ll, expected_count = sum(coverage_validation$adjusted_mutrate), observed_count = sum(coverage_validation$mut_count)))
}



# function to calibrate background mutation rate for each DNM study, some features would be shared among studies and others such as GC contents
TADA_A_adjust_mutation_rate_hierarchical <- function(mut_files, window_file, sample_sizes, scale_features, study_specific = "GC_content", scaling_file_names, mutrate_mode = "regular", genes = "all"){
  # [mut_files] is a list of files with DNM infomation in a BED format. 0-based start and 1-based end. e.g., "../data/Example_windows.bed"
  # The first 4 columns must be chr, start, end, and site_index of genomic windows.
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end", "genename", and "mutrate", followed by features that might affect local background mutation rate.
  # [sample_sizes] is a list of sample sizes.
  # [scale_features] is a vector of the names of features that need to be scaled, this is recommended to apply to continuous features, such as GC content, to make easier the interpretation of the effect sizes of features. 
  # [scaling_file_names] is a vector of names of the file that has the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. 
  # [mutrate_mode], "regular" means the mutation rate of each window is the total sum of mutation rates over all bases. Under this setting, only in very rare cases a window will have mutation rate
  # equal to 0. The output scaling file would remove these windows, as these windows can't be used in the likelihood estimation. "special" means the mutation rate of each window is the rate of 
  # a specific type of mutations. For example, when analyzing WES data, we would want to use syn mutations to adjust for mutation rates. The mutation rate here should accordingly be the rate of 
  # synnonymous mutations. Under this situation, more windows might have mutation rate as 0, but we still need to get the scaling factor for these windows, as the rate of other mutation types may 
  # not be 0, thus the bases in these windows are informative. If [genes] set to be not "all", then must set [mutrate_mode] to "special".
  # [genes] could be "all", or a list of genes. If set to be "all", all windows will be used to adjust mutation rates. Otherwise, only windows with a gene name that 
  # is included in the list of genes specified by [genes] will be used. 
  
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  command <- paste("sed -n '1!p' ", window_file, " | awk {'print $1\"\t\"$2\"\t\"$3\"\t\"$4'} > ", paste(prefix, "_temp_windows.bed", sep = ""), sep = "")
  system(command)
  # get the number of SNVs for each window for each mutation dataset
  for(i in 1:length(mut_files)){
    command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", mut_files[i], " -b ", paste(prefix, "_temp_windows.bed", sep = ""), " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
    system(command)
    # read in the file with the number of SNVs for each window
    coverage_tmp <- fread(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if(i == 1){
      coverage <- coverage_tmp[,1:5]
      colnames(coverage) <- c("chr","start","end","site_index", paste("mut_count", i, sep = "_"))
    }
    else{
      coverage <- data.table(coverage, temp = coverage_tmp$V5)
      colnames(coverage)[dim(coverage)[2]] <- paste("mut_count", i, sep = "_")
    }
  }
  # read in the window file that has feature annotations that might affect mutation rates
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  feature_number <- dim(windows)[2] - 6
  # merge [window] and [coverage] to link mutation count with feature annotations
  coverage <- coverage[windows[,-1:-3], on="site_index"]
  if(genes == "all"){
    target_genes <- unique(coverage$genename)
  }else{
    target_genes <- fread(genes, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    target_genes <- target_genes[[1]] 
  }
  if(mutrate_mode == "regular"){ 
    coverage <- coverage[mutrate !=0]
    for(i in 1:length(scale_features)){
      coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
    }
    # write the formula
    f <- paste("mut_count ~ ", paste(colnames(coverage)[8 : (8 + feature_number -1)], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
    #fit mutation rate model using a fixed set of features not including histone modification marks. 
    out.offset <- glm(as.formula(f), family = poisson, data = coverage)
    scaling <- data.table(site_index = coverage$site_index, scaling_factor = out.offset$fitted.values / (2 * coverage$mutrate * sample_size))
  } else if(mutrate_mode == "special"){
    for(i in 1:length(scale_features)){
      coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
    }
    
    for(m in 1:length(mut_files)){
      if(m == 1){
        temp <- data.table(mut_count = coverage[[paste("mut_count", m, sep = "_")]], coverage[,4], coverage[,(5+length(mut_files)):dim(coverage)[2]], sample_size = sample_sizes[m])
        temp <- data.table(temp, study = m)
      }else{
        temp2 <- data.table(mut_count = coverage[[paste("mut_count", m, sep = "_")]], coverage[,4], coverage[,(5+length(mut_files)):dim(coverage)[2]], sample_size = sample_sizes[m])
        temp2 <- data.table(temp2, study = m)
        temp <- rbind(temp, temp2)
        rm(temp2)
      }
    }
    # Factorize the "study" column of temp
    temp$study <- as.factor(temp$study)
    
    temp_gt0 <- temp[mutrate !=0]
    # only keep windows that are assigned to [target_genes]
    temp_gt0 <- temp_gt0[is.element(genename, target_genes)]
    
    
    # currently only assume GC_content is dataset specific, then need to reformat the dataset.
    f <- paste("mut_count ~ ", paste(colnames(temp)[5 : (5 + feature_number -1)], collapse = " + "), " + study +", paste("GC_content", "study", sep = ":"), " + offset(log(2*mutrate*sample_size))", sep = "")
    out.offset <- glm(as.formula(f), family = poisson, data = temp_gt0)
    
    # calculate scaling factors for all windows
    scaling_factor <- exp(model.matrix(as.formula(f), family = poisson, data = temp) %*% out.offset$coefficients)
    scaling <- data.table(site_index = temp$site_index, scaling_factor = scaling_factor)
    colnames(scaling) <- c("site_index", "scaling_factor")
  }
  
  for(n in 1:length(mut_files)){
    fwrite(scaling[((n-1) * dim(coverage)[1] +1) : (n * dim(coverage)[1])], scaling_file_names[n], col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
  # remove intermediate files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  # the return value is the effect sizes of feature annotations on background mutation rates. 
  return(summary(out.offset)$coeff)
}


# The TADA_A_read_info_temp function here don't remove bases without any annotations.
TADA_A_read_info_temp <- function(mut_files = c("../data/Yuen_NM2015_cases_DNM_with_allele_info.txt","../data/Kong_cases_DNM_with_allele_info.txt","../data/Wu_cases_DNM_with_allele_info.txt", "../data/Jiang_cases_DNM_with_allele_info.txt", "../data/Michaelson_cases_DNM_with_allele_info.txt"),
                             window_file = "../data/Example_windows.bed",
                             mutrate_scaling_files = c("../data/Example_windows_mutrate_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../data/Example_windows_mutrate_scaling_file_for_Kong_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Wu_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Jiang_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Michaelson_cases_DNM.txt"),
                             sample_sizes = c(162, 78, 32, 32, 10),
                             gene_prior_file = "../data/Example_gene_prior.txt",
                             nonAS_noncoding_annotations = c("../data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed","../data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed"),
                             AS_noncoding_annotations = NA,
                             report_proportion = 100/18665,
                             chunk_partition_num =1 ,
                             node_n = 6,
                             mutrate_ref_files = c("../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw"),
                             MPI = 1){
  
  # [mut_file] is a vector of files with DNM infomation in a txt format. The first three columns are chromosome, 0-based start and 1-based end, followed by two columns of ref and alt alleles.
  # The code currently only works for SNVs. 
  # The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end" followed by features that might affect local background mutation rate.
  # [sample_sizes] is a vector of the number of individuals in each study. The order of the numbers should match the order of mutation files in [mut_files]
  # [mutrate_scaling_files] is a vector of files that have the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. 
  # each mutrate_scaling_file matches to one pair of window_file and mut_file. The order of files in [mutrate_scaling_file] should match that in [mut_files]
  # [nonAS_noncoding_annotations] a vector of non-allele-specific non-coding annotations. Each element is a name of the file that has one non-coding annotation in BED format. Non-coding annotations that are not overlapped with regions in [window_file] will not be used in model fitting. 
  # [AS_noncoding_annotations] ia NA or a list of vectors of allele specific annotations. i.e., Each type of noncoding annotation is an element in the list. An element is comprised of 4 different bed files, corresponding to noncoding annotatins of 
  # this type based on the alternative allele, A, T, C, or G. e.g., "spidex_lower10pct_alt_A.bed" is a bed file that has genomic intervals representing the union of all bases which, if mutated to an A allele, have a spidex score lower than 10pct of all
  # possible spidex scores. 
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. 
  # [report_proportion] Choose the top X% TADA genes to estimate RR. 
  # [node_n] is the number of nodes used to run a certain chunk of the code, default is 6
  # [mutrate_ref_files] is a vector of mutrate files in the bigwiggle format. These files have base-level mutation rates to a specific allele, A, C, G, T. 
  # [MPI] is the index that will add to temp files, useful when running multipe processes at one time
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, MPI, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  #mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- windows
  # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  for(i in 1:length(mut_files)){
    mutrate_scaling <- fread(mutrate_scaling_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    system(paste("echo \"Finished reading mutrate scaling file ", mutrate_scaling_files[i], ".\"", sep = ""))
    system("date")
    coverage <- coverage[mutrate_scaling, on = "site_index"]
    coverage <- coverage[complete.cases(coverage)] # release memory
    colnames(coverage)[length(colnames(coverage))] <- paste("scaling_factor_study_", i, sep = "")
    rm(mutrate_scaling) # release memory
  }
  
  # get the piror probability of genes.
  gene_prior <- fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # merge gene prior info
  coverage <- coverage[gene_prior, on = "genename"]
  coverage <-coverage[complete.cases(coverage)]
  
  # select genes based on TADA prior probability and [report_proportion]
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(nrow(genes_for_report)*report_proportion)]
    coverage <- coverage[genes_for_report, on = "genename"]
    coverage <-coverage[complete.cases(coverage)]
  }else{
    genes_for_report  <- gene_prior[,1] # choose all the genes in TADA coding table 
  }
  
  # now need to extropolate window mutation file to base level mutation file
  coverage <-coverage[,c("chr","start","end",paste("scaling_factor_study_", seq(1,length(mut_files)), sep = ""),"genename"),with = FALSE]
  coverage$ID <- paste(coverage$genename, coverage$start, sep = "_")
  
  total_rows <- nrow(coverage)
  interval <- floor(total_rows/chunk_partition_num)
  data_bins <- c(rep(seq(1,chunk_partition_num), each = interval),rep(chunk_partition_num, total_rows -interval*chunk_partition_num))
  
  # split into 20 different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long
  coverage <- split(coverage, data_bins)
  
  #funtion to expand windows to bases
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row["chr"],table_row["genename"],start,sep = "_"), table_row["ID"])
  }
  
  options(warn=-1)
  if(node_n != 1){
    cl <- makeCluster(node_n)
  }
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  
  # get nonAS feature number
  if(is.na(nonAS_noncoding_annotations[1])){
    nonAS_feature_number <- 0 
  }else{
    nonAS_feature_number <- length(nonAS_noncoding_annotations)
  }
  # get AS feature number
  if(is.na(AS_noncoding_annotations[1])){
    AS_feature_number <- 0 
  }else{
    AS_feature_number <- length(AS_noncoding_annotations)
  }
  
  # get total feature number
  feature_number = nonAS_feature_number + AS_feature_number
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,4:(4 + feature_number - 1)],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,4:(4 + feature_number - 1)])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  
  # build a list to store data
  data_partition <-list()
  alt_letters <- c("A","C","G","T")
  # here only deals with the noncoding parts that are within 10kb of TSSs of genes. Will deal with 
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    if(node_n != 1){
      coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage[[i]], 1, window_expansion))
    }else{
      coverage_noncoding_for_base_mutrate <- rbindlist(apply(coverage[[i]], 1, window_expansion))
    }
    system(paste("echo \"Finished partitioning base-level coordinates data at Round ", i, ".\"", sep = ""))
    system("date")
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
    
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    # read in allele-specific base-level mutation rate
    for(j in 1:length(mutrate_ref_files)){
      command <- paste("../external_tools/bigWigAverageOverBed ", mutrate_ref_files[j], " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
      system(command)
      command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
      system(command)
      coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
      colnames(coverage_noncoding_base_mutrate) <- c("base_ID",paste("base_mutrate_alt_", alt_letters[j], sep = ""))
      coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[coverage_noncoding_base_mutrate, on = "base_ID"]
      system(paste("echo \"Finished obtaining base-level mutation rate for alt allele ", alt_letters[j], ".\"", sep = ""))
      system("date")
    }
    
    # read in non allele-specific epigenomic annotations
    epi_ID = 1
    if (!is.na(nonAS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in nonAS_noncoding_annotations){
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, sep = "_"))
        coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        system(paste("echo \"Finished reading non-allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
        system("date")
        epi_ID <- epi_ID + 1
      }
    }
    
    # read in allele-specific epigenomic annotations, now the base assignment for such noncoding annotations is based on the 50-bp windows. 
    # Should not be a big problem as distal introns are assigned based on refseq_ID for intron regions, should be consistent with assignment by spidex.
    # though there is a small chance that a gene's distal intron, which is close to the promoter of another gene and within 10 kb from the TSSs of the latter gene, may be mis-assigned to the latter gene.
    # To completely correct for this issue, need to allow epigenetic annotation to take its own gene assignment, which might be necessary in some situations under strict criteria, such as splicing annotaion. 
    
    if (!is.na(AS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in AS_noncoding_annotations){
        for(k in 1:length(epi)){
          command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi[k], " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
          system(command)
          base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
          base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
          colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, alt_letters[k], sep = "_"))
          coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        }
        system(paste("echo \"Finished reading allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
        system("date")
        epi_ID <- epi_ID + 1
      }
    }
    
    
    # now for each study, read in data, and collapse data based on noncoding annotation configuration
    for(m in 1:length(mut_files)){
      mut <- fread(mut_files[m], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      for(letter in alt_letters){
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }
        # a very important step here is removing bases that have adjusted_mutrate_base 0. This happens when the allele of the mutrate we are using is just the reference allele. 
        # By doing this, we make the computation of likelihood valid, also we automatically removed bases with nonAS annotations but with mutant allele the same with ref allele at this step
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage_noncoding_for_base_mutrate_temp[[paste("base_mutrate_alt_", letter, sep = "")]] !=0]
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage[[i]][,c(paste("scaling_factor_study_", m, sep = ""), "genename","ID"), with = FALSE],on = "ID"]
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- coverage_noncoding_for_base_mutrate_temp[, paste("base_mutrate_alt_", letter, sep = ""), with = FALSE] * coverage_noncoding_for_base_mutrate_temp[, paste("scaling_factor_study_", m, sep = ""), with = FALSE] * 2 * sample_sizes[m]
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }
        mut_allele <- mut[V5 == letter]
        # the subsequent two lines of code are used to prevent a outputing bug when using fwrite. (85000000 to be written as 8.5e7)
        mut_allele$V2 <- as.integer(mut_allele$V2)
        mut_allele$V3 <- as.integer(mut_allele$V3)
        fwrite(mut_allele, paste(prefix, "_temp_mut_allele.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", paste(prefix, "_temp_mut_allele.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""),sep = "")
        system(command)
        base_with_mut <- fread(paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_with_mut <- base_with_mut[,c("V4","V5"), with = FALSE]
        colnames(base_with_mut) <- c("base_ID", "mut_count")
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[base_with_mut, on = "base_ID"]
        
        # have to collpase data at this point, otherwise, the I will run out of RAM if process 1000 genes every time. 
        # anno_count <- rep(0, nrow(coverage_noncoding_for_base_mutrate_temp))
        # for(p in 1:feature_number){
        #   anno_count <- anno_count + coverage_noncoding_for_base_mutrate_temp[[(4 + p -1)]]
        # }
        # # remove rows(bases) that don't have any non-coding features, this could save a lot of RAM, so I could use smaller partition number which would greatly accelerate speed when read in data for all genes.
        # coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[anno_count >0,]
        # 
        # first partition by gene for the current chunk
        coverage_noncoding_for_base_mutrate_temp<- split(coverage_noncoding_for_base_mutrate_temp, coverage_noncoding_for_base_mutrate_temp$genename)
        # then partition by feature configuration for each gene in the current chunk
        coverage_noncoding_for_base_mutrate_temp <- sapply(coverage_noncoding_for_base_mutrate_temp, partition_feature, simplify = FALSE)
        # add compact data
        data_partition <- append(data_partition, coverage_noncoding_for_base_mutrate_temp)
        rm(coverage_noncoding_for_base_mutrate_temp) # release memory
        system(paste("echo \"Finished read in mutation data and make them into the compact format for Study ", m, " and allele ", letter, ".\"", sep = ""))
        system("date")
      }
    }
  }
  
  if(node_n != 1){
    stopCluster(cl)
  }
  options(warn = 0)
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste("echo \"Temp files cleaned and data recording finished!\""))
  system("date")
  return(list(base_info = data_partition))
}

