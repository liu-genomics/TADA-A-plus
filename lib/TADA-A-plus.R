TADA_A_compare_mutrate_model <- function(mut_file, window_file, sample_size, scale_features, validation_chr = "chr1", mutrate_mode = "regular", genes = "all"){
  # [mut_file] is the file with DNM infomation in a BED format. 0-based start and 1-based end. e.g., "../data/Example_windows.bed"
  # The first 4 columns must be chr, start, end, and site_index of genomic windows.
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
  # the output of this function is the loglikelihood of observing data given the estimated parameters from the training model. THe validation is only performed on the validation
  # chromosome indicated by [validation_chr], the training model didn't use any information from the 
  
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
    
    f <- paste("mut_count ~ ", paste(colnames(coverage_mutrate_gt0_training)[8 : (8 + feature_number -1)], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
    out.offset <- glm(as.formula(f), family = poisson, data = coverage_mutrate_gt0_training)
    scaling_factor <- exp(as.matrix(cbind(1, coverage[,8 : (8 + feature_number -1)])) %*% out.offset$coefficients)
    coverage <- data.table(coverage, scaling_factor = scaling_factor)
    coverage$adjusted_mutrate <- 2 * sample_size * coverage$mutrate * coverage$scaling_factor
    coverage_validation <- coverage[mutrate != 0 & chr == validation_chr & is.element(genename, target_genes)]
    ll <- sum(mapply(dpois, coverage_validation$mut_count, coverage_validation$adjusted_mutrate, MoreArgs = list(log = TRUE)))
  }
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

