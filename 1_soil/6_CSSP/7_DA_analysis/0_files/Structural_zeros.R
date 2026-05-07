Structural_zeros <- function(OTU_table, meta_data, group, ref = NULL,
                             min_reads, min_present_reps){
  if(!(group %in% colnames(meta_data))){
    stop(paste( "Variable", group, "not contained in metadata" ))
  }
  grp_var <- as.data.frame(meta_data)[,group]
  if(!(class(grp_var) %in% c("factor", "character"))){
    stop(paste("Expected character or string as group variable, but recieved", class(grp_var)))
  }
  if(class(grp_var) == "character" & is.null(ref)){
    stop("Reference level must be given when group variable is given as character")
  }
  if(class(grp_var) != "factor"){
    grp_var <- as.factor(grp_var)
    grp_var <- relevel(grp_var, ref = ref)
  }
  grps <- levels(grp_var)
  struc_zero_mat <- matrix(NA, nrow = nrow(OTU_table), ncol = length(grps))
  rownames(struc_zero_mat) <- rownames(OTU_table)
  colnames(struc_zero_mat) <- grps
  for(i in 1:length(grps)){
    idx <- grp_var %in% grps[i]
    struc_zero_mat[,i] <- as.numeric(rowSums(OTU_table[,idx]) == 0)
  }
  
  DA_struc_zero <- list()
  for(i in 2:ncol(struc_zero_mat)){
    struc_zero <- rownames(struc_zero_mat)[rowSums(struc_zero_mat[,c(1, i)]) == 1]
    idx <- grp_var %in% grps[c(1, i)]
    struc_zero_otu_table <- OTU_table[struc_zero,idx]
    rownames(struc_zero_otu_table)[ rowSums(struc_zero_otu_table) > min_reads
                                    & rowSums(struc_zero_otu_table!=0) >= min_present_reps
    ] -> DA_struc_zero[[i-1]]
    if(all(struc_zero_mat[,c(1,i)] == 0)){
      DA_struc_zero[[i-1]] <- character(0)
    }
    
  }
  names(DA_struc_zero) <- grps[-1]
  L <- list(struc_zero_table = struc_zero_mat, struc_zero_DA = DA_struc_zero)
}

Structural_zeros2 <- function(OTU_table, meta_data, group, ref = NULL,
                              min_reads, min_present_reps){
  if(!("library_size" %in% colnames(meta_data))){
    stop("meta_data must contain numeric coulmn named 'library_size'")
  }
  if(!(group %in% colnames(meta_data))){
    stop(paste( "Variable", group, "not contained in metadata" ))
  }
  grp_var <- as.data.frame(meta_data)[,group]
  if(!(class(grp_var) %in% c("factor", "character"))){
    stop(paste("Expected character or string as group variable, but recieved", class(grp_var)))
  }
  if(class(grp_var) == "character" & is.null(ref)){
    stop("Reference level must be given when group variable is given as character")
  }
  if(class(grp_var) != "factor"){
    grp_var <- as.factor(grp_var)
    grp_var <- relevel(grp_var, ref = ref)
  }
  grps <- levels(grp_var)
  struc_zero_mat <- matrix(NA, nrow = nrow(OTU_table), ncol = length(grps))
  rownames(struc_zero_mat) <- rownames(OTU_table)
  colnames(struc_zero_mat) <- grps
  for(i in 1:length(grps)){
    idx <- grp_var %in% grps[i]
    struc_zero_mat[,i] <- as.numeric(rowSums(OTU_table[,idx]) == 0)
  }
  
  RA <- apply(OTU_table, 2, function(x) x/sum(x))
  lib_size <- aggregate(x = meta_data$library_size, by = list(grp_var), mean)
  n_reads <- lib_size$x; names(n_reads) <- as.character(lib_size$Group.1)
  thresh <- 20/min(n_reads)
  
  DA_struc_zero <- list()
  for(i in 2:ncol(struc_zero_mat)){
    struc_zero <- rownames(struc_zero_mat)[rowSums(struc_zero_mat[,c(1, i)]) == 1]
    idx <- grp_var %in% grps[c(1, i)]
    
    if(any(struc_zero_mat[,c(1,i)] != 0)){
      struc_zero_otu_table <- OTU_table[struc_zero,idx]
      struc_zero_otu_table <- RA[struc_zero,idx]
      dt <- data.table(grp = grp_var[idx], t(struc_zero_otu_table))
      mean_RA <- apply(dt[,lapply(.SD, mean), grp][,-1], 2, max)
      reps_present <- apply(dt[,lapply(.SD, function(x) sum(x!=0)), grp][,-1], 2, max)
      
      DA_struc_grp_i <- rownames(struc_zero_otu_table)[ mean_RA > thresh & reps_present>= min_present_reps] 
      DA_struc_zero[[i-1]] <- DA_struc_grp_i
    }
    
    if(all(struc_zero_mat[,c(1,i)] == 0)){
      DA_struc_zero[[i-1]] <- character(0)
    }
    
  }
  names(DA_struc_zero) <- grps[-1]
  L <- list(struc_zero_table = struc_zero_mat, struc_zero_DA = DA_struc_zero)
}

Structural_zeros3 <- function(OTU_table, meta_data, group, ref = NULL,
                              min_reads, min_present_reps){
  if(!("library_size" %in% colnames(meta_data))){
    stop("meta_data must contain numeric coulmn named 'library_size'")
  }
  if(!(group %in% colnames(meta_data))){
    stop(paste( "Variable", group, "not contained in metadata" ))
  }
  grp_var <- as.data.frame(meta_data)[,group]
  if(!(class(grp_var) %in% c("factor", "character"))){
    stop(paste("Expected character or string as group variable, but recieved", class(grp_var)))
  }
  if(class(grp_var) == "character" & is.null(ref)){
    stop("Reference level must be given when group variable is given as character")
  }
  if(class(grp_var) != "factor"){
    grp_var <- as.factor(grp_var)
    grp_var <- relevel(grp_var, ref = ref)
  }
  grps <- levels(grp_var)
  struc_zero_mat <- matrix(NA, nrow = nrow(OTU_table), ncol = length(grps))
  rownames(struc_zero_mat) <- rownames(OTU_table)
  colnames(struc_zero_mat) <- grps
  for(i in 1:length(grps)){
    idx <- grp_var %in% grps[i]
    struc_zero_mat[,i] <- as.numeric(rowSums(OTU_table[,idx]) == 0)
  }
  
  RA <- apply(OTU_table, 2, function(x) x/sum(x))
  lib_size <- aggregate(x = meta_data$library_size, by = list(grp_var), mean)
  n_reads <- lib_size$x; names(n_reads) <- as.character(lib_size$Group.1)
  
  DA_struc_zero <- list()
  for(i in 2:ncol(struc_zero_mat)){
    struc_zero <- rownames(struc_zero_mat)[rowSums(struc_zero_mat[,c(1, i)]) == 1]
    grp_sub <- grps[c(1, i)]
    grp_sub <- grp_sub[order(n_reads[grp_sub])]
    idx <- grp_var %in% grp_sub
    
    if(any(struc_zero_mat[,c(1,i)] != 0)){
      s1 <- n_reads[grp_sub][1]/n_reads[grp_sub][2]
      s <- max(c(s1, 1/s1))
      struc_zero_otu_table <- OTU_table[struc_zero,idx]
      struc_zero_mat_sub <- struc_zero_mat[,grp_sub]
      sz.d <- apply(struc_zero_mat_sub, 1, diff)
      zero_small <- names(sz.d)[sz.d == -1]
      zero_large <- names(sz.d)[sz.d == 1]
      struc_zero_otu_table_small <- struc_zero_otu_table[zero_small,]
      struc_zero_otu_table_large <- struc_zero_otu_table[zero_large,]
      rownames(struc_zero_otu_table_large)[ rowSums(struc_zero_otu_table_large) > min_reads
                                            & rowSums(struc_zero_otu_table_large!=0) >= min_present_reps
                                            ] -> DA_struc_large
      rownames(struc_zero_otu_table_small)[ rowSums(struc_zero_otu_table_small) > min_reads*s
                                            & rowSums(struc_zero_otu_table_small!=0) >= min_present_reps
                                            ] -> DA_struc_small
      DA_struc_zero[[i-1]] <- c(DA_struc_large, DA_struc_small)
    }
    
    if(all(struc_zero_mat[,c(1,i)] == 0)){
      DA_struc_zero[[i-1]] <- character(0)
    }
    
  }
  names(DA_struc_zero) <- grps[-1]
  L <- list(struc_zero_table = struc_zero_mat, struc_zero_DA = DA_struc_zero)
}

Structural_zeros_simplified <- function(OTU_table, group){
  grp_var <- as.data.frame(meta_data)[,group]
  if(!(class(group) %in% c("factor", "character"))){
    stop(paste("Expected character or string as group variable, but recieved", class(group)))
  }
  grps <- levels(group)
  struc_zero_mat <- matrix(NA, nrow = nrow(OTU_table), ncol = 2)
  rownames(struc_zero_mat) <- rownames(OTU_table)
  colnames(struc_zero_mat) <- grps
  for(i in 1:2){
    idx <- group %in% grps[i]
    struc_zero_mat[,i] <- as.numeric(rowSums(OTU_table[,idx]) == 0)
  }
  
  DA_struc_zero <- rownames(struc_zero_mat)[apply(struc_zero_mat, 1, diff) != 0]
  
  L <- list(struc_zero_table = struc_zero_mat, struc_zeros = DA_struc_zero)
  
  return(L)
}
