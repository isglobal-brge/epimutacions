#' @export
EpiMutations<-function(case_sample, control_panel = NULL, method = "manova", chr = NULL, min_cpg = 3, verbose = TRUE)
{
  window_sz = 10
  offset_mean = 0.15
  offset_abs = 0.1
  
  bump_cutoff <-  0.1
  nsamp <- "deterministic"
  qn_th <- 3
  
  # Identify type of input and extract required data:
  #	* betas
  #	* sample's classification
  #	* feature annotation
  
  if(is.null(case_sample))
  {
    stop("The argument 'case_sample' must be introduced")
    
  }else{
    
    if(ncol(case_sample) != 1){
      stop("The argument 'case_sample' must have 1 sample")
    }
    if(class(case_sample) == "GenomicRatioSet") {
      if(verbose) message("Input of type 'GeomicRatioSet'")
      if (is.null(control_panel)){
        #Combine control panel and case sample
        set <- minfi::combineArrays(grs_control_panel, case_sample,
                                    outType = c("IlluminaHumanMethylation450k",
                                                "IlluminaHumanMethylationEPIC",
                                                "IlluminaHumanMethylation27k"),
                                    verbose = TRUE)
      }else{
        if(class(control_panel)!= "GenomicRatioSet"){
          stop("The type of the arguments 'control_panel' and 'case_sample' must be the same,'GeomicRatioSet'")
        }
        set <- minfi::combineArrays(control_panel, case_sample,
                                    outType = c("IlluminaHumanMethylation450k",
                                                "IlluminaHumanMethylationEPIC",
                                                "IlluminaHumanMethylation27k"),
                                    verbose = TRUE)
      }
      betas <- minfi::getBeta(set)
      pd <- as.data.frame(SummarizedExperiment::colData(set))
      fd <- as.data.frame(SummarizedExperiment::rowRanges(set))
      rownames(fd) <- rownames(set) 
      }else if(class(case_sample) == "ExpressionSet") {
        if(verbose) message("Input of type 'ExpressionSet")
        if (is.null(control_panel)){
          #Combine control panel and case sample
          set <- a4Base::combineTwoExpressionSet(es_control_panel, case_sample)
      }else{
        if(class(control_panel)!= "ExpressionSet"){
          stop("The type of the arguments 'control_panel' and 'case_sample' must be the same,'ExpressionSet'")
        }
        set <- a4Base::combineTwoExpressionSet(control_panel, case_sample)
      }
      betas <- Biobase::exprs(set)
      pd <- Biobase::pData(set)
      fd <- Biobase::fData(set)
    } else {
      stop("Input data 'case_sample' must be a 'GenomicRatioSet' or an 'ExpressionSet'")
    }
  }
  # Identify the method to be used
  if(!("start" %in%  colnames(fd)) & !("seqnames" %in%  colnames(fd))){
    stop("In feature data variable 'seqnames' and 'start' must be introduced specifying the chromosome and start position")
  }
  avail <- c("manova", "mlm", "isoforest", "mahdistmcd", "barbosa", "qn")
  method <- charmatch(method, avail)
  method <- avail[method]
  if(is.na(method)) stop("Invalid method was selected'")
  if(verbose) message("Selected epimutation detection method '", method, "'")
  #if(method %in% c()) {
  #	stop("Method not implemented yet")
  #}

  # Identify cases and controls
  cas_sam <- colnames(case_sample)
  pd$status <- ifelse(rownames(pd) == cas_sam, 1, 0)
  ctr_sam <- rownames(pd)[pd$status == 1]
  
  # Differentiate between methods that required region detection that the ones
  # that finds outliers to identify regions
  if(method %in% c("manova", "mlm", "mahdistmcd", "isoforest")) {
    if(verbose) message(paste0("Selected method '", method, "' required of 'bumphunter'"))
      # Prepare model to be evaluated
      model <- stats::model.matrix(~status, pd)
      # Select chromosome
      if(!is.null(chr)){
       fd <- fd[fd$seqnames == chr,]
       betas <- betas[rownames(fd),]
      }
      
      
      # Run bumphunter for region partitioning
        bumps <- bumphunter::bumphunter(object = betas, design = model,
                                        pos = fd$start, chr = fd$seqnames, cutoff = bump_cutoff)$table
        if(!is.na(bumps)){
      
      bumps <- bumps[bumps$L >= min_cpg, ]
      bumps$sz <- bumps$end - bumps$start
      #bumps <- bumps[bumps$sz < length(ctr_sam), ] # <--------------- TODO
      if(verbose) message(paste0(nrow(bumps), " candidate regions were found for case sample '", cas_sam, "'"))
      
      # Identify outliers according to selected method
      rst <- do.call(rbind, lapply(seq_len(nrow(bumps)), function(ii) {
        bump <- bumps[ii, ]
        beta_bump <- betas_from_bump(bump, fd, betas)
        
        if(method == "mahdistmcd") {
          dst <- epi_mahdistmcd(beta_bump, nsamp)
          threshold <- sqrt(qchisq(p = 0.975, df = ncol(beta_bump)))
          outliers <- which(dst$statistic >= threshold)
          outliers <- dst$ID[outliers]
          return(res_mahdistmcd(cas_sam, bump, beta_bump, outliers))
        } else if(method == "mlm") {
          sts <- epi_mlm(beta_bump, model)
          return(res_mlm(bump, beta_bump, sts, cas_sam))
        } else if(method == "manova") {
          sts <- epi_manova(beta_bump, model, cas_sam)
          return(res_manova(bump, beta_bump, sts, cas_sam))
        } else if(method == "isoforest") {
          sts <- epi_isoforest(beta_bump, cas_sam)
          return(res_isoforest(bump, beta_bump, sts, cas_sam))
        }
        
      }))
        }else{
          rst <- NA
          warning("No bumps found!")
        }
  }else if(method == "barbosa") {
    # Compute reference statistics
    if(verbose) message("Calculating statistics from reference distribution required by Barbosa et. al. 2019")
    bctr_min <- apply(betas[ , ctr_sam], 1, min, na.rm = TRUE)
    bctr_max <- apply(betas[ , ctr_sam], 1, max, na.rm = TRUE)
    bctr_mean <- apply(betas[ , ctr_sam], 1, mean, na.rm = TRUE)
    bctr_prc <- suppressWarnings(apply(betas[ , ctr_sam], 1, quantile, probs = c(0.01, 0.99), na.rm = TRUE))
    bctr_pmin <- bctr_prc[1, ]
    bctr_pmax <- bctr_prc[2, ]
    rm(bctr_prc)
    #case <- betas[ , cas_sam[1], drop=FALSE]
    
    # Run region detection
      x <- epi_barbosa(betas[ , cas_sam, drop=FALSE], fd, bctr_min, bctr_max, bctr_mean, 
                       bctr_pmin, bctr_pmax, window_sz, min_cpg, offset_mean, offset_abs)
      x$sample <- cas_sam
      x
    
    # rst$epi_id <- sapply(seq_len(nrow(rst)), function(ii) paste0("epi_", method, "_", ii))
    # colnames(rst) <- c("chromosome", "start", "end", "sz", "cpg_n", "cpg_ids", 
    # 	   "outlier_score", "outlier_significance", "outlier_direction", 
    # 	   "sample", "epi_id")
    # rownames(rst) <- seq_len(nrow(rst))
    # 
    # return(rst[ , c(11, 10, 1:9)])
  } else { # if(method == "qn") {
    nbetas <- qn_norm(betas, qn = TRUE)
    regions <- qn_bump(nbetas[ , cas_sam], fd, window = window_sz, cutoff = bump_cutoff)
    x <- qn_outlier(cas_sam, regions, nbetas, fd, min_cpg, qn_th)
    x$sample <- cas_sam
    x
    rst <- rst[rst$outlier_direction != "", ]
  }
  if(is.na(rst)){
    return(rst)
  }else{
    rst$epi_id <- sapply(seq_len(nrow(rst)), function(ii) paste0("epi_", method, "_", ii))
    colnames(rst) <- c("chromosome", "start", "end", "sz", "cpg_n", "cpg_ids", 
                       "outlier_score", "outlier_significance", "outlier_direction", 
                       "sample", "epi_id")
    rownames(rst) <- seq_len(nrow(rst))
    #rst <- rst[ , c(11, 10, 1:9)]
    rst <- rst[apply(rst, 1, function(x) !all(is.na(x))),]
    return(rst)
  }
}
  