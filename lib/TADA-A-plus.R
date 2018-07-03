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

# functions to read in mutation data that has matched cases and controls. Binomial tests instead of Poisson tests will be performed to estimte the RR of annotations.
TADA_A_read_info_pair <- function(mut_files_1 = c("../data/table.ASCWGS_20180504.WGS1902_hg19_cases_SNV_remove_recurrent_mutations_with_allele_info.txt"),
                                       mut_files_2 = c("../data/table.ASCWGS_20180504.WGS1902_hg19_controls_SNV_remove_recurrent_mutations_with_allele_info.txt"),
                                       window_file = "../data/windows_partition/Example_windows_with_div_score_no_header.bed.temp.partition00.with_header.txt",
                                       gene_prior_file = "../data/Example_gene_prior.txt",
                                       nonAS_noncoding_annotations = c("../data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed","../data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed"),
                                       AS_noncoding_annotations = list(c("../other_annotations/coding/171121_coding_nonsynonymous_SNV_altA.bed.merge.bed", "../other_annotations/coding/171121_coding_nonsynonymous_SNV_altC.bed.merge.bed", "../other_annotations/coding/171121_coding_nonsynonymous_SNV_altG.bed.merge.bed", "../other_annotations/coding/171121_coding_nonsynonymous_SNV_altT.bed.merge.bed")),
                                       mutation_mutrate_adjusting_features_file_1 = NA,
                                       mutation_mutrate_adjusting_features_file_2 = NA,
                                       report_proportion = 100/18665,
                                       MPI = 1){
  
  # [mut_files_1] is a vector of files with cases DNM infomation in a txt format. The first three columns are chromosome, 0-based start and 1-based end, followed by two columns of ref and alt alleles.
  # The code currently only works for SNVs. 
  # The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
  # [mut_files_2] is a vector of files with controls DNM infomation in a txt format. The first three columns are chromosome, 0-based start and 1-based end, followed by two columns of ref and alt alleles.
  # The ordering of files in [mut_files_1] and [mut_files_2] should be consistent.
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end" followed by features that might affect local background mutation rate.
  # [nonAS_noncoding_annotations] a vector of non-allele-specific non-coding annotations. Each element is a name of the file that has one non-coding annotation in BED format. Non-coding annotations that are not overlapped with regions in [window_file] will not be used in model fitting. 
  # [AS_noncoding_annotations] ia NA or a list of vectors of allele specific annotations. i.e., Each type of noncoding annotation is an element in the list. An element is comprised of 4 different bed files, corresponding to noncoding annotatins of 
  # this type based on the alternative allele, A, T, C, or G. e.g., "spidex_lower10pct_alt_A.bed" is a bed file that has genomic intervals representing the union of all bases which, if mutated to an A allele, have a spidex score lower than 10pct of all
  # possible spidex scores. 
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. 
  # [mutation_mutrate_adjusting_features_file_1] is a file each row of which is a mutation bed file, followed by continuous features that might affect mutation rate. The file needs to have
  # features for all the mutations in [mut_files1] and [mut_files2].The variables (e.g. sequencing depth) should come from studies generating the corresponding DNM file in [mut_files1]
  # [mutation_mutrate_adjusting_features_file_2] is a file each row of which is a mutation bed file, followed by continuous features that might affect mutation rate. The file needs to have
  # features for all the mutations in [mut_files1] and [mut_files2].The variables (e.g. sequencing depth) should come from studies generating the corresponding DNM file in [mut_files1]
  # [report_proportion] Choose the top X% TADA genes to estimate RR. 
  # [MPI] is the index that will add to temp files, useful when running multipe processes at one time
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, MPI, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  #mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  fwrite(windows[,c("chr","start","end","genename")], paste(prefix, "_temp_windows_bed.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # only keep bases that have overlap with DNMs, this step will reduce the computation time significantly (If the DNM dataset is relatively smaller compared to the number of bases under consideration). 
  # write a big file that has contains all the DNMs
  command <- paste("echo -n > ", paste(prefix, "_temp_DNM_union.bed", sep = ""), sep = "")
  system(command)
  
  command <- paste("cat", paste(mut_files_1, collapse = " "), paste(mut_files_2, collapse = " "), ">", paste(prefix, "_temp_DNM_union.bed", sep = ""), sep = " ")
  system(command)
  
  command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools intersect -b ", paste(prefix, "_temp_DNM_union.bed", sep = ""), " -a ", paste(prefix, "_temp_windows_bed.bed", sep = ""),
                   " -wa -wb > ", paste(prefix, "_temp_windows_bed_overlap_DNMs.bed", sep = ""), sep = "")
  system(command)
  
  # Read in windows that have DNM overlapped.
  coverage <- fread(paste(prefix, "_temp_windows_bed_overlap_DNMs.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- coverage[,c("V5", "V6", "V7", "V4")]
  colnames(coverage) <- c("chr", "start", "end", "genename")
  coverage <- unique(coverage)
  coverage <- data.table(coverage, ID = paste("base_", seq(1,nrow(coverage)), sep = ""))
  
  
  # # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  # for(i in 1:length(mut_files)){
  #   mutrate_scaling <- fread(mutrate_scaling_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  #   system(paste("echo \"Finished reading mutrate scaling file ", mutrate_scaling_files[i], ".\"", sep = ""))
  #   system("date")
  #   coverage <- coverage[mutrate_scaling, on = "site_index"]
  #   coverage <- coverage[complete.cases(coverage)] # release memory
  #   colnames(coverage)[length(colnames(coverage))] <- paste("scaling_factor_study_", i, sep = "")
  #   rm(mutrate_scaling) # release memory
  # }
  # 
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
  # These information are those necessary to compute log-likelihood in the optimization function, without doing categorization, as we have continuous variables
  # to adjust mutation rates and the number of mutations is not bit. 
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,3:(3 + feature_number - 1)],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,3:(3 + feature_number - 1)])), sum_cases = sum(feature_set$mut_count_1), sum_controls = sum(feature_set$mut_count_2))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  
  # generae a window file to get base level annotations.
  fwrite(coverage[,c(1:3,5)],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # read in non allele-specific epigenomic annotations
  epi_ID = 1
  if (!is.na(nonAS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
    for(epi in nonAS_noncoding_annotations){
      command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
      system(command)
      base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
      colnames(base_in_epi) <- c("ID", paste("Anno",epi_ID, sep = "_"))
      coverage <- coverage[base_in_epi, on = "ID"]
      system(paste("echo \"Finished reading non-allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
      system("date")
      epi_ID <- epi_ID + 1
    }
  }
  
  alt_letters <- c("A","C","G","T")
  
  # read in non-AS annotations
  if (!is.na(AS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
    for(epi in AS_noncoding_annotations){
      for(k in 1:length(epi)){
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi[k], " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("ID", paste("Anno",epi_ID, alt_letters[k], sep = "_"))
        coverage <- coverage[base_in_epi, on = "ID"]
      }
      system(paste("echo \"Finished reading allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
      system("date")
      epi_ID <- epi_ID + 1
    }
  }
  
  
  # build a list to store data
  data_partition <-list()
  
  
  # now for each study, read in data, and collapse data based on noncoding annotation configuration
  for(m in 1:length(mut_files_1)){
    mut_1 <- fread(mut_files_1[m], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mut_2 <- fread(mut_files_2[m], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    for(letter in alt_letters){
      if(!is.na(mutation_mutrate_adjusting_features_file_1) & !is.na(mutation_mutrate_adjusting_features_file_2)){
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_temp <- coverage[, c("genename", "ID", "chr", "start", "end", paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_temp <- coverage[, c("genename", "ID", "chr", "start", "end", paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_temp <- coverage[, c("genename", "ID", "chr", "start", "end", paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }
      }else{
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_temp <- coverage[, c("genename", "ID", paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_temp <- coverage[, c("genename", "ID", paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_temp <- coverage[, c("genename", "ID", paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }
      }
      mut_allele_1 <- mut_1[V5 == letter]
      mut_allele_2 <- mut_2[V5 == letter]
      # the subsequent two lines of code are used to prevent a outputing bug when using fwrite. (85000000 to be written as 8.5e7)
      mut_allele_1$V2 <- as.integer(mut_allele_1$V2)
      mut_allele_1$V3 <- as.integer(mut_allele_1$V3)
      mut_allele_2$V2 <- as.integer(mut_allele_2$V2)
      mut_allele_2$V3 <- as.integer(mut_allele_2$V3)
      fwrite(mut_allele_1, paste(prefix, "_temp_mut_allele_1.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      fwrite(mut_allele_2, paste(prefix, "_temp_mut_allele_2.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      command_1 <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", paste(prefix, "_temp_mut_allele_1.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_mut_allele_1.bed", sep = ""),sep = "")
      command_2 <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", paste(prefix, "_temp_mut_allele_2.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_mut_allele_2.bed", sep = ""),sep = "")
      system(command_1)
      system(command_2)
      base_with_mut_1 <- fread(paste(prefix,"_temp_for_mutrate_overlap_mut_allele_1.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      base_with_mut_1 <- base_with_mut_1[,c("V4","V5"), with = FALSE]
      colnames(base_with_mut_1) <- c("ID", "mut_count_1")
      coverage_temp <- coverage_temp[base_with_mut_1, on = "ID"]
      base_with_mut_2 <- fread(paste(prefix,"_temp_for_mutrate_overlap_mut_allele_2.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      base_with_mut_2 <- base_with_mut_2[,c("V4","V5"), with = FALSE]
      colnames(base_with_mut_2) <- c("ID", "mut_count_2")
      coverage_temp <- coverage_temp[base_with_mut_2, on = "ID"]
      # For binomial test, can't just remove bases without any annotations
      # By contrast, need to remove bases that don't have any mutations either in cases or controls.
      coverage_temp <- coverage_temp[mut_count_1 > 0 | mut_count_2 > 0]  
      
      if(!is.na(mutation_mutrate_adjusting_features_file_1) & !is.na(mutation_mutrate_adjusting_features_file_2)){# if base-level mutation adjusting features are provided
        # read in mutrate adjusting variables
        mutrate_adj_1 <- fread(mutation_mutrate_adjusting_features_file_1[m], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        mutrate_adj_2 <- fread(mutation_mutrate_adjusting_features_file_2[m], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        mutrate_adj <- mutrate_adj_1[mutrate_adj_2, on = c("chr", "start", "end")]
        mutrate_adj_feature_num <- ncol(mutrate_adj_1) - 5
        mutrate_adj <- mutrate_adj[,c(1:3, 6:(5+mutrate_adj_feature_num), ((8+mutrate_adj_feature_num) : (7+2*mutrate_adj_feature_num))), with = FALSE]
        coverage_temp <- mutrate_adj[coverage_temp, on = c("chr", "start", "end")]
        coverage_temp <- unique(coverage_temp)
        coverage_temp <- coverage_temp[,c(-1,-2,-3), with = FALSE]
        # add compact data
        data_partition <- append(data_partition, list(coverage_temp))
      }else{
        coverage_temp<- split(coverage_temp, coverage_temp$genename)
        # then partition by feature configuration for each gene in the current chunk
        coverage_temp <- sapply(coverage_temp, partition_feature, simplify = FALSE)
        # add compact data
        data_partition <- append(data_partition, coverage_temp)
      }
      
      rm(coverage_temp) # release memory
      system(paste("echo \"Finished read in mutation data and make them into the compact format for Study ", m, " and allele ", letter, ".\"", sep = ""))
      system("date")
    }
  }
  
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste("echo \"Temp files cleaned and data recording finished!\""))
  system("date")
  return(list(base_info = data_partition))
}



TADA_A_RR_estimate_binom <-function(data, selected_annotations, gene_prior_file, optimization_iteration = 2000, mode = "regular"){
  #[data] is the [base_info] returned from [TADA_A_reading_info_pair], which contains all the allele specific data across all studies
  #[selected_annotations] is a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
  #[gene_prior_file], #[gene_prior_file], a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  #[mode] is "regular", or "single_fast". "single_fast" is used when estimating RR from only one annotation ([data] only recoreded one annotation) of lots of genes (e.g., all genes), would be 5 times faster.
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  gene_prior$prior <- as.numeric(gene_prior$prior)
  gene_prior$genename = as.character(gene_prior$genename)
  if(mode == "single_fast"){
    further_collapsed_data <- list()
    genes_with_annotation <- unique(names(data))
    
    for(i in 1:length(genes_with_annotation)){
      further_collapsed_data <- append(further_collapsed_data, list(list(`1` = list(feature_vector = 1 , sum_cases= 0, sum_controls = 0))))
    }
    names(further_collapsed_data) <- unique(names(data))
    
    for(i in 1:length(data)){
      Part_of_data = data[i]
      further_collapsed_data[[names(Part_of_data)]][[1]]$sum_cases <- further_collapsed_data[[names(Part_of_data)]][[1]]$sum_cases + Part_of_data[[1]][[1]]$sum_cases
      further_collapsed_data[[names(Part_of_data)]][[1]]$sum_controls <- further_collapsed_data[[names(Part_of_data)]][[1]]$sum_controls + Part_of_data[[1]][[1]]$sum_controls
    }
    data <- further_collapsed_data
  }
  
  # notice the fr function is different from the function that deals with dataset without partition, needs to consider splicing mutation together
  fr <-function(x){ # x is the RRs of selected noncoding annotations. 
    all_rr = x # has an intercept for this case, as mutation rates haven't been adjusted for each dataset. 
    cal_logP_Zg1 <- function(data_partition_element){
      cal_logP_Zg1_level2 <-function(data_partition_element_level2){
        K = exp(all_rr[1] + data_partition_element_level2[[1]][selected_annotations]%*%all_rr[-1])
        data_partition_element_level2[[2]] * (log(K) - log(1+K)) - data_partition_element_level2[[3]] * log(1+K)
      }
      sum(sapply(data_partition_element, cal_logP_Zg1_level2))
    }
    
    cal_logP_Zg0 <- function(data_partition_element){
      cal_logP_Zg0_level2 <-function(data_partition_element_level2){
        K = exp(all_rr[1])
        data_partition_element_level2[[2]] * (log(K) - log(1+K)) - data_partition_element_level2[[3]] * log(1+K)
      }
      sum(sapply(data_partition_element, cal_logP_Zg0_level2))
    }
    
    logP_Zg1 = sapply(data, cal_logP_Zg1)
    logP_Zg0 = sapply(data, cal_logP_Zg0)
    
    logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(logP_Zg1))
    logP_table <- logP_table[gene_prior, on = "genename"]
    logP_table <- logP_table[complete.cases(logP_table)]
    ll_sum1 <- sum(by(logP_table, logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
    ll_sum1
  }
  # Use [optim] to do optimization, for non-splicing_mutations
  feature_number <- length(selected_annotations)
 
  mle <- optim(rep(0.1, (feature_number + 1)), fr ,control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = TRUE)
  
  # get confidence intervals of RR estimates
  fisher_info <- solve(-mle$hessian)
  prop_sigma <- sqrt(c(diag(fisher_info)))
  rr_estimate <- c(mle$par)
  upper<-rr_estimate+1.96*prop_sigma
  lower<-rr_estimate-1.96*prop_sigma
  rr_report <-data.frame(relative_risk = rr_estimate, lower_bound = lower, upper_bournd = upper)
  
  list(mle = mle, rr_report = rr_report)
}




# TADA_A_RR_estimate_binom_v2 is used when base-level mutation rate adjusting featurs were used. 
# base_level mutation mutation-rate adjusting features are different among cases and controls. Otherwise, use [TADA_A_RR_estimate_binom_v3]
TADA_A_RR_estimate_binom_v2 <-function(data, selected_annotations, gene_prior_file, optimization_iteration = 2000, mode = "regular"){
  #[data] is the [base_info] returned from [TADA_A_reading_info_pair], which contains all the allele specific data across all studies. [mutation_mutrate_adjusting_features_file_1] and [mutation_mutrate_adjusting_features_file_2]
  # need to b non-NAs. Otherwise, [TADA_A_RR_estimat_binom] needs to be used. 
  #[selected_annotations] is a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
  #[gene_prior_file], #[gene_prior_file], a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  #[mode] is "regular", or "single_fast". "single_fast" is used when estimating RR from only one annotation ([data] only recoreded one annotation) of lots of genes (e.g., all genes), would be 5 times faster.
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  gene_prior$prior <- as.numeric(gene_prior$prior)
  gene_prior$genename = as.character(gene_prior$genename)
  
  # get the number of mutrate_adjusting_featurs
  mut_adj_feature_num <- which(colnames(data[[1]]) == "genename") - 1
  # get the number of total functional annotations that need to have RR estimated. 
  functional_feature_num <- ncol(data[[1]]) - 2 -2 - mut_adj_feature_num
  data <- rbindlist(data)
  data <- data[gene_prior, on = "genename"]
  data <- data[complete.cases(data)]
  
  data <- data[, c(c(1 : mut_adj_feature_num), mut_adj_feature_num + 2 + selected_annotations, (ncol(data) -2) : ncol(data), (mut_adj_feature_num + 1)) ,with = FALSE] 
  data <- split(data, data$genename)
  # notice the fr function is different from the function that deals with dataset without partition, needs to consider splicing mutation together
  fr <-function(x){ # x is the RRs of selected noncoding annotations. 
    beta0 <- x[1]
    mut_par <- x[2:(mut_adj_feature_num + 1)] * rep( c(1, -1), each = (mut_adj_feature_num/2))
    all_rr <- x[(mut_adj_feature_num + 2): (mut_adj_feature_num + length(selected_annotations) + 1)] # has an intercept for this case, as mutation rates haven't been adjusted for each dataset. 
    
    
    cal_logP <- function(data_partition_element){
      K1 <- exp(beta0 + as.matrix(data_partition_element[, 1:(ncol(data_partition_element) -4), with = FALSE]) %*% c(mut_par, all_rr))
      K1_logP <- data_partition_element$mut_count_1 * (log(K1) - log(K1+1)) - data_partition_element$mut_count_2 * log(1+K1)
      K0 <- exp(beta0 + as.matrix(data_partition_element[, 1:(ncol(data_partition_element) -4), with = FALSE]) %*% c(mut_par, rep(0, length(selected_annotations))))
      K0_logP <- data_partition_element$mut_count_1 * (log(K0) - log(K0+1)) - data_partition_element$mut_count_2 * log(1+K0)
      log(exp(log(data_partition_element[1,]$prior) + sum(K1_logP)) + exp(log((1 - data_partition_element[1,]$prior)) + sum(K0_logP))) 
    }
    
    sum(sapply(data, cal_logP))
  }
  # Use [optim] to do optimization, for non-splicing_mutations
  feature_number <- length(selected_annotations) + mut_adj_feature_num
  
  mle <- optim(rep(0.1, (feature_number + 1)), fr ,control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = TRUE)
  
  # get confidence intervals of RR estimates
  fisher_info <- solve(-mle$hessian)
  prop_sigma <- sqrt(c(diag(fisher_info)))
  rr_estimate <- c(mle$par)
  upper<-rr_estimate+1.96*prop_sigma
  lower<-rr_estimate-1.96*prop_sigma
  rr_report <-data.frame(relative_risk = rr_estimate, lower_bound = lower, upper_bournd = upper)
  
  list(mle = mle, rr_report = rr_report)
}



# TADA_A_RR_estimate_binom_v3 is used when base-level mutation rate adjusting featurs were used and the features for cases and controls are exactly the same. 
TADA_A_RR_estimate_binom_v3 <-function(data, selected_annotations, gene_prior_file, optimization_iteration = 2000, mode = "regular"){
  #[data] is the [base_info] returned from [TADA_A_reading_info_pair], which contains all the allele specific data across all studies. [mutation_mutrate_adjusting_features_file_1] and [mutation_mutrate_adjusting_features_file_2]
  # need to b non-NAs. Otherwise, [TADA_A_RR_estimat_binom] needs to be used. 
  #[selected_annotations] is a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
  #[gene_prior_file], #[gene_prior_file], a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  #[mode] is "regular", or "single_fast". "single_fast" is used when estimating RR from only one annotation ([data] only recoreded one annotation) of lots of genes (e.g., all genes), would be 5 times faster.
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  gene_prior$prior <- as.numeric(gene_prior$prior)
  gene_prior$genename = as.character(gene_prior$genename)
  
  # get the number of mutrate_adjusting_featurs
  mut_adj_feature_num <- (which(colnames(data[[1]]) == "genename") - 1)/2 # there is only one parameter for a pair of covariates that are perfectly correlated. 
  # get the number of total functional annotations that need to have RR estimated. 
  functional_feature_num <- ncol(data[[1]]) - 2 -2 - (which(colnames(data[[1]]) == "genename") - 1)
  data <- rbindlist(data)
  data <- data[gene_prior, on = "genename"]
  data <- data[complete.cases(data)]
  
  data <- data[, c(c(1 : mut_adj_feature_num), 2* mut_adj_feature_num + 2 + selected_annotations, (ncol(data) -2) : ncol(data), (2* mut_adj_feature_num + 1)) ,with = FALSE] 
  data <- split(data, data$genename)
  # notice the fr function is different from the function that deals with dataset without partition, needs to consider splicing mutation together
  fr <-function(x){ # x is the RRs of selected noncoding annotations. 
    beta0 <- x[1]
    mut_par <- x[2:(mut_adj_feature_num + 1)]
    all_rr <- x[(mut_adj_feature_num + 2): (mut_adj_feature_num + length(selected_annotations) + 1)] # has an intercept for this case, as mutation rates haven't been adjusted for each dataset. 
    
    
    cal_logP <- function(data_partition_element){
      K1 <- exp(beta0 + as.matrix(data_partition_element[, 1:(ncol(data_partition_element) -4), with = FALSE]) %*% c(mut_par, all_rr))
      K1_logP <- data_partition_element$mut_count_1 * (log(K1) - log(K1+1)) - data_partition_element$mut_count_2 * log(1+K1)
      K0 <- exp(beta0 + as.matrix(data_partition_element[, 1:(ncol(data_partition_element) -4), with = FALSE]) %*% c(mut_par, rep(0, length(selected_annotations))))
      K0_logP <- data_partition_element$mut_count_1 * (log(K0) - log(K0+1)) - data_partition_element$mut_count_2 * log(1+K0)
      log(exp(log(data_partition_element[1,]$prior) + sum(K1_logP)) + exp(log((1 - data_partition_element[1,]$prior)) + sum(K0_logP))) 
    }
    
    sum(sapply(data, cal_logP))
  }
  # Use [optim] to do optimization, for non-splicing_mutations
  feature_number <- length(selected_annotations) + mut_adj_feature_num
  
  mle <- optim(rep(0.1, (feature_number + 1)), fr ,control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = TRUE)
  
  # get confidence intervals of RR estimates
  fisher_info <- solve(-mle$hessian)
  prop_sigma <- sqrt(c(diag(fisher_info)))
  rr_estimate <- c(mle$par)
  upper<-rr_estimate+1.96*prop_sigma
  lower<-rr_estimate-1.96*prop_sigma
  rr_report <-data.frame(relative_risk = rr_estimate, lower_bound = lower, upper_bournd = upper)
  
  list(mle = mle, rr_report = rr_report)
}


# Function to do simple burden analysis using Fisher exact test. 
simple_burden <- function(file1, file2, annotation){
  # [file1] is the name of cases, in bed format, as long as having the first three columns
  # [file2] is the name of controls in bed format, as long as having the first three columns
  # [annotation] is the name of annotation bed files, need to has four columns, chr, start, end, name. 
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
  command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", file1, " -b ", annotation, " > ", prefix, ".temp.file1.coverageBed", sep = "")
  system(command)
  command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", file2, " -b ", annotation, " > ", prefix, ".temp.file2.coverageBed", sep = "")
  system(command)
  file1_coverage <- fread(paste(prefix, ".temp.file1.coverageBed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  file2_coverage <- fread(paste(prefix, ".temp.file2.coverageBed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  file1_DNM <- read.delim(file1, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  file2_DNM <- read.delim(file2, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  burden <- (sum(file1_coverage$V5)/sum(file2_coverage$V5))/(nrow(file1_DNM)/nrow(file2_DNM))
  p.value <- fisher.test(matrix(c(sum(file1_coverage$V5), sum(file2_coverage$V5), (nrow(file1_DNM) - sum(file1_coverage$V5)), (nrow(file2_DNM) - sum(file2_coverage$V5))), 2, 2), alternative = "greater")$p.value
  system(paste("rm ", prefix, ".temp.file1.coverageBed", sep = ""))
  system(paste("rm ", prefix, ".temp.file2.coverageBed", sep = ""))
  return(list(burden = burden, p.value = p.value, cases_count = sum(file1_coverage$V5), controls_count = sum(file2_coverage$V5)))
}

# Use output from [TADA_A_read_info_pair] when mutation rate adjusting features were used, under this situation, the output will be a list of data.tables, would be easier 
# to calculate burden signal.

TADA_A_RR_estimate_simple_burden<-function(data, gene_list = "all", background = c(117685, 114849)){
  #[data] is the [base_info] returned from [TADA_A_reading_info_pair], which contains all the allele specific data across all studies. [mutation_mutrate_adjusting_features_file_1] and [mutation_mutrate_adjusting_features_file_2]
  # need to b non-NAs. Otherwise, [TADA_A_RR_estimat_binom] needs to be used. 
  #[gene_list],could be "all", or a vector that has a list of genenames from which mutations will be used to calculate burden.
  #[backgroun] is a vector of two numbers, the first is the total number of mutations in cases while the second is the total number of mutations in controls. 
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  #[mode] is "regular", or "single_fast". "single_fast" is used when estimating RR from only one annotation ([data] only recoreded one annotation) of lots of genes (e.g., all genes), would be 5 times faster.
  
  
  data <- rbindlist(data)
  if(gene_list != "all"){
    data <- data[is.element(data$genename, gene_list)]
  }
  # get the number of mutrate_adjusting_featurs
  mut_adj_feature_num <- which(colnames(data) == "genename") - 1
  # get the number of total functional annotations that need to have RR estimated. 
  functional_feature_num <- ncol(data) - 2 -2 - mut_adj_feature_num
  output <- data.table(burden = rep(-1, functional_feature_num), p.value = rep(-1, functional_feature_num), case_count = rep(-1, functional_feature_num), control_count = rep(-1, functional_feature_num))
  for(i in 1:functional_feature_num){
    mut_count_1_with_anno <- sum(data[data[[(mut_adj_feature_num + 2 + i)]] == 1]$mut_count_1)
    mut_count_2_with_anno <- sum(data[data[[(mut_adj_feature_num + 2 + i)]] == 1]$mut_count_2)
    output[i,]$burden <- (mut_count_1_with_anno/mut_count_2_with_anno)/(background[1]/background[2])
    output[i,]$p.value <- fisher.test(matrix(c(mut_count_1_with_anno, mut_count_2_with_anno, (background[1] - mut_count_1_with_anno), (background[2] - mut_count_2_with_anno)),2,2), alternative = "greater")$p.value
    output[i,]$case_count <- mut_count_1_with_anno
    output[i,]$control_count <- mut_count_2_with_anno
  }
  
  list(rr_report = output)
}