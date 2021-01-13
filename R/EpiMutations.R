#' @export
EpiMutations<-function(case_samples, control_panel, fd, method = "manova", chr = NULL, start = NULL, end = NULL, epi_params = epi_params(), bump_cutoff =  0.1, min_cpg = 3, verbose = TRUE)
{
  
  # Identify type of input and extract required data:
  #	* betas
  #	* sample's classification
  #	* feature annotation
  if(is.null(case_samples))
  {
    stop("The argument 'case_samples' must be introduced")
    
  }
  if(is.null(control_panel))
  {
    stop("The argument 'case_samples' must be introduced")
    
  }
  if(!class(case_samples)[1] %in% c("matrix", "data.frame")){
    stop("'case_samples' must be a matrix or data frame")
  }
  if(!class(control_panel)[1] %in% c("matrix", "data.frame")){
    stop("'control_panel' must be a matrix or data frame")
  }
  fd <- .fd_cols(fd)
  
  if(!is.null(start) & !is.null(end)){
    if(is.null(chr)){
      stop("Argument 'chr' must be inroduced with 'start' and 'end' parameters")
    }
    if(start > end){
      stop("'start' cannot be higher than 'end'")
    }
    
  }
  if(!is.null(start) & is.null(end) | is.null(start) & !is.null(end)){
    stop("'start' and 'end' arguments must be introduced together")
    
  }
  avail <- c("manova", "mlm", "isoforest", "mahdistmcd", "barbosa", "qn")
  method <- charmatch(method, avail)
  method <- avail[method]
  if(is.na(method)) stop("Invalid method was selected'")
  if(method == "barbosa" & ncol(case_samples) <= 1){
    warning("More than 1 case sample must be introduced to use barbosa")
  }
  if(verbose) message("Selected epimutation detection method '", method, "'")
  #if(method %in% c()) {
  #	stop("Method not implemented yet")
  #}
  
  if(!is.null(chr)){
    if(!is.null(start) & !is.null(end)){
      fd <- fd[fd$seqnames %in% chr & fd$start >= start & fd$end <= end,]
    }else{
      fd <- fd[fd$seqnames %in% chr,]
    }
    case_samples <- case_sample[rownames(fd),]
    control_panel <- control_panel[rownames(fd),]
    betas <- betas[rownames(fd),]
    #pd <- pd[which(rownames(pd) %in% colnames(betas)),]
  }
  
  
  # Identify cases and controls
  cas_sam <- colnames(case_samples)
  ctr_sam <- colnames(control_panel)
  
  # Differentiate between methods that required region detection that the ones
  # that finds outliers to identify regions
  if(method %in% c("manova", "mlm", "mahdistmcd", "isoforest")) {
    if(verbose) message(paste0("Selected method '", method, "' required of 'bumphunter'"))
    # Prepare model to be evaluated
    rst <- do.call(rbind, lapply(cas_sam, function(case){
      samples_names <- c(ctr_sam, case)
      status <- samples_names == case
      status <- as.data.frame(status, row.names = samples_names)
      betas <- cbind(control_panel, case_samples[,case,drop=FALSE])
      model <- stats::model.matrix(~status, status)
      # Run bumphunter for region partitioning
      bumps <- bumphunter::bumphunter(object = betas, design = model,
                                      pos = fd$start, chr = fd$seqnames, cutoff = bump_cutoff)$table
      
      suppressWarnings(
        if(!is.na(bumps) | nrow(bumps) != 0){
          
          bumps <- bumps[bumps$L >= min_cpg, ]
          bumps$sz <- bumps$end - bumps$start
          # bumps <- bumps[bumps$sz < length(ctr_sam), ] # <--------------- TODO
          if(verbose) message(paste0(nrow(bumps), " candidate regions were found for case sample '", case, "'"))
          
          # Identify outliers according to selected method
          bump_out  <- do.call(rbind, lapply(seq_len(nrow(bumps)), function(ii) {
            bump <- bumps[ii, ]
            beta_bump <- betas_from_bump(bump, fd, betas)
            if(method == "mahdistmcd") {
              dst <- epi_mahdistmcd(beta_bump, epi_params$mahdistmcd$nsamp)
              threshold <- sqrt(qchisq(p = 0.975, df = ncol(beta_bump)))
              outliers <- which(dst$statistic >= threshold)
              outliers <- dst$ID[outliers]
              return(res_mahdistmcd(case, bump, beta_bump, outliers))
            } else if(method == "mlm") {
              sts <- epi_mlm(beta_bump, model)
              return(res_mlm(bump, beta_bump, sts, case))
            } else if(method == "manova") {
              sts <- epi_manova(beta_bump, model, case)
              return(res_manova(bump, beta_bump, sts, case))
            } else if(method == "isoforest") {
              sts <- epi_isoforest(beta_bump, case)
              return(res_isoforest(bump, beta_bump, sts, case))
            }
            
          }))
        })
    }))
    if(nrow(rst) == 0){
      rst <- NA
    }
  }else if(method == "barbosa") {
    # Compute reference statistics
    if(verbose) message("Calculating statistics from reference distribution required by Barbosa et. al. 2019")
    bctr_min <- apply(control_panel, 1, min, na.rm = TRUE)
    bctr_max <- apply(control_panel, 1, max, na.rm = TRUE)
    bctr_mean <- apply(control_panel, 1, mean, na.rm = TRUE)
    bctr_prc <- suppressWarnings(apply(control_panel, 1, quantile, probs = c(0.01, 0.99), na.rm = TRUE))
    bctr_pmin <- bctr_prc[1, ]
    bctr_pmax <- bctr_prc[2, ]
    rm(bctr_prc)
    #case <- betas[ , cas_sam[1], drop=FALSE]
    
    # Run region detection
    rst <- do.call(rbind, lapply(cas_sam, function(case) {
      x <- epi_barbosa(case_samples[ , case, drop=FALSE], fd, bctr_min, bctr_max, bctr_mean, 
                       bctr_pmin, bctr_pmax, window_sz, min_cpg, epi_params$barbosa$offset_mean, epi_params$barbosa$offset_abs)
      if(nrow(x) != 0){
        x$sample <- case 
      }else{
        x <- NA
      }
      x
    }))
    # rst$epi_id <- sapply(seq_len(nrow(rst)), function(ii) paste0("epi_", method, "_", ii))
    # colnames(rst) <- c("chromosome", "start", "end", "sz", "cpg_n", "cpg_ids", 
    # 	   "outlier_score", "outlier_significance", "outlier_direction", 
    # 	   "sample", "epi_id")
    # rownames(rst) <- seq_len(nrow(rst))
    # 
    # return(rst[ , c(11, 10, 1:9)])
  } else { # if(method == "qn") {
    nbetas <- qn_norm(cbind(control_panel, case_samples), qn = TRUE)
    regions <- qn_bump(nbetas[,cas_sam], fd, window = epi_params$qn$window_sz, cutoff = bump_cutoff)
    rst <- do.call(rbind, lapply(cas_sam, function(case) {
      x <- qn_outlier(case, regions, nbetas, fd, min_cpg, epi_params$qn$qn_th)
      if(nrow(x) != 0){
        x$sample <- case 
      }else{
        x <- NA
      }
      x
    }))
    rst <- rst[rst$outlier_direction != "", ]
  }
  suppressWarnings(
    if(is.na(rst)){
      return(message("No outliers found"))
      
    }else{
      rst$epi_id <- sapply(seq_len(nrow(rst)), function(ii) paste0("epi_", method, "_", ii))
      colnames(rst) <- c("chromosome", "start", "end", "sz", "cpg_n", "cpg_ids", 
                         "outlier_score", "outlier_significance", "outlier_direction", 
                         "sample", "epi_id")
      rownames(rst) <- seq_len(nrow(rst))
      rst <- rst[ , c(11, 10, 1:9)]
      if(method == "manova"){
        pvalue_cutoff <- epi_params$manova$pvalue_cutoff
      }else if(method == "mlm"){
        pvalue_cutoff <- epi_params$mlm$pvalue_cutoff
        
      }
      rst <- filter_results(rst, method, pvalue_cutoff, epi_params$isoforest$outlier_score_cutoff)
      rst <- tibble::as_tibble(rst)
      return(rst)
    })
}

# Helper functions

.fd_cols <- function(fd){
  if(class(fd) != "data.frame"){
    stop("'fd' must be a data frame")
  }
  seqnames_field <- c("seqnames", "seqname",
                      "chromosome", "chrom",
                      "chr", "chromosome_name",
                      "seqid")
  start_field <- "start"
  end_field <- c("end", "stop")
  strand_field <- c("strand")
  
  
  if(seqnames_field[1] %in% colnames(fd) | seqnames_field[2] %in% colnames(fd) | seqnames_field[3] %in% colnames(fd) | seqnames_field[4] %in% colnames(fd) | seqnames_field[5] %in% colnames(fd) | seqnames_field[6] %in% colnames(fd) | seqnames_field[7] %in% colnames(fd)){
    seqnames_pos <- which(colnames(fd) %in% seqnames_field)
    seqnames <- fd[,seqnames_pos]
  }else{
    stop("In feature dataset chromosome column name must be specified as: 'seqnames', 'seqname', 'chromosome', 'chrom', 'chr', 'chromosome_name' or 'seqid'")
  }
  if(start_field %in% colnames(fd)){
    start_pos <- which(colnames(fd) %in% start_field)
    start <- fd[,start_pos]
  }else{
    stop("In feature dataset start column name must be specified as: 'start'")
  }
  if(end_field[1] %in% colnames(fd) | end_field[2] %in% colnames(fd)){
    end_pos <- which(colnames(fd) %in% end_field)
    end <- fd[,end_pos]
  }else{
    stop("In feature dataset end column name  must be specified as: 'end' or 'stop'")
  }
  if( strand_field %in% colnames(fd)){
    strand_pos <- which(colnames(fd) %in% strand_field)
    strand <- fd[,strand_pos]
  }else{
    stop("In feature dataset strand column name  must be specified as: 'strand'")
  }
  fd_rownames <- rownames(fd)
  fd <- data.frame("seqnames" = seqnames, "start" = start, "end" = end, "strand" = strand)
  rownames(fd) <- fd_rownames
  return(fd)
}
