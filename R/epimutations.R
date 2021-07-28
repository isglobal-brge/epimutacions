#' @title Epimutations analysis based on outlier detection methods
#' @description The function identifies  Differentially Methylated Regions
#' in a case sample by comparing it against a control panel. 
#'  @param case_samples a GenomicRatioSet object containing the case samples.
#' See the constructor function \link[minfi]{GenomicRatioSet}, \link[minfi]{makeGenomicRatioSetFromMatrix}. 
#' @param control_panel a GenomicRatioSet object containing the control panel (control panel).
#' @param method a character string naming the outlier detection method to be used. 
#' This can be set as: \code{"manova"}, \code{"mlm"}, \code{"isoforest"}, \code{"mahdistmcd"}, 
#' \code{"quantile"} and \code{"beta"}. 
#' The default is \code{"manova"}. 
#' For more information see \strong{Details}. 
#' @param chr a character string containing the sequence names to be analysed. The default value is \code{NULL}. 
#' @param start an integer specifying the start position. The default value is \code{NULL}.
#' @param end an integer specifying the end position. The default value is \code{NULL}.
#' @param epi_params the parameters for each method. See the function \link[epimutations]{epi_parameters}.  
#' @param maxGap the maximum location gap used in \link[bumphunter]{bumphunter} method. 
#' @param bump_cutoff a numeric value of the estimate of the genomic profile above the 
#' cutoff or below the negative of the cutoff will be used as candidate regions. 
#' @param min_cpg an integer specifying the minimum CpGs number in a DMR.  
#' @param verbose logical. If TRUE additional details about the procedure will provide to the user. 
#' The default is TRUE. 
#' @details The function compares a case sample against a control panel to identify epimutations in the given 
#' sample. First, the DMRs are identified using the \link[bumphunter]{bumphunter} approach. 
#' After that, CpGs in those DMRs are tested in order to detect regions
#' with CpGs being outliers.  For that, different outlier detection methods can be selected:  
#'  * Multivariate Analysis of Variance (\code{"manova"}). \link[stats]{manova}
#'  * Multivariate Linear Model (\code{"mlm"})
#'  * Isolation Forest (\code{"isoforest"}) \link[isotree]{isolation.forest}
#'  * Robust Mahalanobis Distance (\code{"mahdistmcd"}) \link[robustbase]{covMcd}
#'  * Quantile distribution (\code{"quantile"})
#'  * Beta (\code{"beta"})
#'  
#' We defined candidate epimutation regions (found in candRegsGR) based on the 450K 
#' array design. As CpGs are not equally distributed along the genome, only CpGs closer
#' to other CpGs can form an epimutation. More information can be found in candRegsGR documentation.
#' 
#' @return The function returns an object of class tibble containing the outliers regions.  
#' The results are composed by the following columns: 
#' * \code{epi_id}: systematic name for each epimutation identified. It provides the name 
#' of the used anomaly detection method. 
#' * \code{sample}: the name of the sample containing the epimutation. 
#' * \code{chromosome}, \code{start} and \code{end}: indicate the location of the epimutation.
#' * \code{sz}: the window's size of the event.
#' * \code{cpg_n}: the number of CpGs in the epimutation.
#' * \code{cpg_ids}: the names of CpGs in the epimutation.
#' * \code{outlier_score}: 
#'    * For method \code{manova} it provides the approximation to F-test and the Pillai score, separated by \code{/}.
#'    * For method \code{mlm} it provides the approximation to F-test and the R2 of the model, separated by \code{/}.
#'    * For method \code{isoforest} it provides the magnitude of the outlier score.
#'    * For method \code{beta} it provides the mean outlier p-value.
#'    * For methods \code{quantile} and \code{mahdistmcd} it is filled with NA.
#' * \code{outlier_direction}: indicates the direction of the outlier with \code{"hypomethylation"} and \code{"hypermethylation"}
#'    * For \code{manova}, \code{mlm}, \code{isoforest}, and \code{mahdistmcd} it is computed from the values obtained from bumphunter.
#'    * For \code{quantile} it is computed from the location of the sample in the reference distribution (left vs. right outlier).
#'    * For method \code{beta} it return a NA.
#' * \code{pvalue}: 
#'    * For methods \code{manova}, \code{mlm}, and \code{isoforest} it provides the p-value obtained from the model.
#'    * For method \code{quantile}, \code{mahdistmcd} and \code{beta} is filled with NA.    
#' * \code{adj_pvalue}: for methods with p-value (\code{manova} and \code{mlm} adjusted p-value with Benjamini-Hochberg based on the total number of regions detected by Bumphunter.
#' * \code{epi_region_id}: Name of the epimutation region as defined in \code{candRegsGR}.
#' * \code{CRE}: cREs (cis-Regulatory Elements) as defined by ENCODE overlapping the epimutation region. Different cREs are separated by ;.
#' * \code{CRE_type}: Type of cREs (cis-Regulatory Elements) as defined by ENCODE. Different type are separeted by , and different cREs are separated by ;.
#' @examples 
#' \dontrun{
#' library(epimutacions)
#' 
#' # The data for this example is available in epimutacionsData (ExperimentHub) package
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' query(eh, c("epimutacionsData"))
#' methy <- eh[["EH6692"]]
#' 
#' #Find epimutations in GSM2562701 sample of methy dataset
#' 
#' case_samples <- methy[,"GSM2562701"]
#' control_panel <- methy[,-51]
#' epimutations(case_samples, control_panel, method = "manova")
#' }
#' @export
epimutations <- function(case_samples, control_panel,
                         method = "manova", 
                         chr = NULL, start = NULL, end = NULL, 
                         epi_params = epi_parameters(), 
                         maxGap = 1000, bump_cutoff =  0.1, 
                         min_cpg = 3, verbose = TRUE)
{
  
  # 1. Inputs check and data extraction
  ## Inputs check
  if(is.null(case_samples)){
    stop("The argument 'case_samples' must be introduced")
    
  }
  if(is.null(control_panel)){
    stop("The argument 'case_samples' must be introduced")
  }

  if(class(case_samples) != "GenomicRatioSet"){
    stop("'case_samples' must be of class 'GenomicRatioSet'") 
  }
  if(class(control_panel) != "GenomicRatioSet"){
    stop("'control_panel' must be of class 'GenomicRatioSet'. 
         To create a 'GenomicRatioSet' object use 'makeGenomicRatioSetFromMatrix'
         function from minfi package")
  }
  if(minfi::annotation(case_samples)[1] != minfi::annotation(control_panel)[1] & minfi::annotation(case_samples)[2] != minfi::annotation(control_panel)[2]){
    stop("The annotation of 'case_samples' and 'control_panel' must be the same")
  }
  
  if(!is.null(start) & !is.null(end)){
    if(is.null(chr)){
      stop("Argument 'chr' must be inroduced with 'start' and 'end' parameters")
    }
    if(length(start) != length(end) & length(chr) != length(start)){
      stop("'start' and 'end' length must be same")
    }
    for(i in 1:length(start)){
      if(start[i] > end[i]){
        stop("'start' cannot be higher than 'end'")
      }
      
    }
    
  }
  if(!is.null(start) & is.null(end) | is.null(start) & !is.null(end)){
    stop("'start' and 'end' arguments must be introduced together")
    
  }
  
  avail <- c("manova", "mlm", "isoforest", "mahdistmcd", "quantile", "beta")
  method <- charmatch(method, avail)
  method <- avail[method]
  if(is.na(method)) stop("Invalid method was selected'")
  
  if(verbose) message("Selected epimutation detection method '", method, "'")
  
  
  ## Extract required data:
      #	* feature annotation
      #	* betas
      #	* sample's classification
  
  ### Feature annotation
  fd <- as.data.frame(GenomicRanges::granges(case_samples))
  rownames(fd) <- rownames(case_samples)
  ### Betas
  betas_case <- minfi::getBeta(case_samples)
  betas_case <- betas_case[rownames(fd),,drop =FALSE]
  betas_control <- minfi::getBeta(control_panel)
  betas_control <- betas_control[rownames(fd),]
  
  ### Select CpGs in the specified in the arguments 'chr', 'start' and 'end'
  
  if(!is.null(chr)){
    if(!is.null(start) & !is.null(end)){
      fd_split <- NULL
      for(i in 1:length(chr)){
        fd_split <- rbind(fd_split, fd[fd$seqnames %in% chr[i] & fd$start >= start[i] & fd$end <= end[i],])
      }
      fd <- fd_split
      rm(fd_split)
      
    }else{
      fd <- fd[fd$seqnames %in% chr,]
    }
    if(nrow(fd) == 0){
      stop("No CpG was found in the specified region")
    }
    
    betas_case <- betas_case[rownames(fd),,drop=FALSE]
    betas_control <- betas_control[rownames(fd),]
  }
  
  ### Identify case and control samples
  cas_sam <- colnames(betas_case)
  ctr_sam <- colnames(betas_control)
  
  
  # 2. Epimutations definition (using different methods)
  ## Methods that need bumphunter ("manova", "mlm", "mahdistmcd" and "isoforest")
  if(method %in% c("manova", "mlm", "mahdistmcd", "isoforest")) {
    if(verbose) message(paste0("Selected method '", method, "' required of 'bumphunter'"))
    # Prepare model to be evaluated
    rst <- do.call(rbind, lapply(cas_sam, function(case){
      samples_names <- c(ctr_sam, case)
      status <- samples_names == case
      status <- as.data.frame(status, row.names = samples_names)
      betas <- cbind(betas_control, betas_case[,case,drop=FALSE])
      model <- stats::model.matrix(~status, status)
      
      # Run bumphunter for region partitioning
      bumps <- bumphunter::bumphunter(object = betas, design = model,
                                      pos = fd$start, chr = fd$seqnames, 
                                      maxGap = maxGap,
                                      cutoff = bump_cutoff)$table
      suppressWarnings(
        if(!is.na(bumps)){
          
          ## Homogeneize output of bumphunter to epimutacions naming
          bumps$chromosome <- bumps$chr
          bumps$cpg_n <- bumps$L
          bumps$sz <- bumps$end - bumps$start
          
          bumps <- bumps[bumps$L >= min_cpg, ]
          # bumps <- bumps[bumps$sz < length(ctr_sam), ] # <--------------- TODO
          if(verbose) message(paste0(nrow(bumps), " candidate regions were found for case sample '", case, "'"))
          if(nrow(bumps) != 0){
          # Identify outliers according to selected method
          bumps  <- do.call(rbind, lapply(seq_len(nrow(bumps)), function(ii){
            bump <- bumps[ii, ]
            beta_bump <- betas_from_bump(bump, fd, betas)
            if(method == "mahdistmcd") {
              dst <- try(epi_mahdistmcd(beta_bump, epi_params$mahdistmcd$nsamp), silent = TRUE)
              if(class(dst) != "try-error"){
                threshold <- sqrt(stats::qchisq(p = 0.999975, df = ncol(beta_bump)))
                outliers <- which(dst$statistic >= threshold)
                outliers <- dst$ID[outliers] 
                x <- res_mahdistmcd(case, bump, beta_bump, outliers)
              }
            } else if(method == "mlm") {
              sts <- try(epi_mlm(beta_bump, model), silent = TRUE)
              #sts <- epi_mlm(beta_bump, model)
              x <- res_mlm(bump, beta_bump, sts, case)
            } else if(method == "manova") {
              sts <- epi_manova(beta_bump, model, case)
              x <- res_manova(bump, beta_bump, sts, case)
            } else if(method == "isoforest") {
              sts <- epi_isoforest(beta_bump, case, epi_params$isoforest$ntrees)
              x <- res_isoforest(bump, beta_bump, sts, case, epi_params$isoforest$outlier_score_cutoff)
            } 
          }))
          ## Filter using the adjusted p-value calculated from regions identified in each sample ("manova and "mlm")
          if(method == "manova" & !is.null(bumps)){
            bumps <- filter_manova(bumps, epi_params$manova$pvalue_cutoff)
          }
          if(method == "mlm" & !is.null(bumps)){
            bumps <- filter_mlm(bumps, epi_params$mlm$pvalue_cutoff)
          }
       
          }
        })
      #Add a row filled by NAs for the samples with any epimutations
      if(is.null(nrow(bumps)) || nrow(bumps) == 0){
        bumps <- data.frame(chromosome = 0,
                            start = 0,
                            end = 0,
                            sz = NA,
                            cpg_n = NA,
                            cpg_ids = NA,
                            outlier_score = NA,
                            outlier_direction = NA,
                            pvalue = NA,
                            adj_pvalue = NA,
                            sample = case)
      }
      bumps
    }))
    ## Methods that do not need bumphunter ("quantile" and "beta")
  }else if(method == "quantile") {
      # Compute reference statistics
    if(verbose) message("Calculating statistics from 'quantile' method")
    #bctr_min <- suppressWarnings(apply(betas_control, 1, min, na.rm = TRUE))
    #bctr_max <- suppressWarnings(apply(betas_control, 1, max, na.rm = TRUE))
    #bctr_mean <- suppressWarnings(apply(betas_control, 1, mean, na.rm = TRUE))
    #bctr_prc <- suppressWarnings(apply(betas_control, 1, quantile, probs = c(0.999975, 0.000025), na.rm = TRUE))
    if(verbose) message("Using quantiles ", epi_params$quantile$qinf, " and ", epi_params$quantile$qsup)
    bctr_prc <- suppressWarnings(apply(betas_control, 1, quantile, probs = c(epi_params$quantile$qinf, epi_params$quantile$qsup), na.rm = TRUE))
    bctr_pmin <- bctr_prc[1, ]
    bctr_pmax <- bctr_prc[2, ]
    rm(bctr_prc)
    #case <- betas[ , cas_sam[1], drop=FALSE]
      # Run region detection
    rst <- do.call(rbind, lapply(cas_sam, function(case) {
      #x <- epi_quantile(betas_case[ , case, drop = FALSE], fd, bctr_min, bctr_max, bctr_mean, 
      #                 bctr_pmin, bctr_pmax, window_sz, min_cpg, epi_params$quantile$offset_mean, epi_params$quantile$offset_abs)
      x <- epi_quantile(betas_case[ , case, drop = FALSE], fd, bctr_pmin, bctr_pmax, epi_params$quantile$window_sz, min_cpg, epi_params$quantile$offset_abs)
      if(is.null(x) || nrow(x) == 0){
        x <- data.frame(
          chromosome = 0,
          start = 0,
          end = 0,
          sz = NA,
          cpg_n = NA,
          cpg_ids = NA,
          outlier_score = NA,
          outlier_direction = NA,
          pvalue = NA,
          adj_pvalue = NA,
          sample = case
        )
      } else {
        x$sample <- case 
      }
      x
    }))
  } else if(method == "beta") {
    
    # Get Beta distribution params
    message("Computing beta distribution parameters")
    beta_params <- getBetaParams(t(betas_control))
    beta_mean <- rowMeans(betas_control, na.rm = TRUE)
    
    message("Defining Regions")
    rst <- do.call(rbind, lapply(cas_sam, function(case) {
      x <- epi_beta(beta_params, beta_mean, 
                    betas_case[ , case, drop = FALSE], 
                    GenomicRanges::makeGRangesFromDataFrame(fd[rownames(betas_control),]),
                    epi_params$beta$pvalue_cutoff, 
                    epi_params$beta$diff_threshold, min_cpg, maxGap)
      if(nrow(x) == 0){
        x <- data.frame(
          chromosome = 0,
          start = 0,
          end = 0,
          sz = NA,
          cpg_n = NA,
          cpg_ids = NA,
          outlier_score = NA,
          outlier_direction = NA,
          pvalue = NA,
          adj_pvalue = NA,
          sample = case
        )
      } else {
        x$sample <- case 
      }
      x
    }))
  }   
      # 3. Prepare the output and addition of CREs
      ## Prepare the output
      rst$epi_id <- sapply(seq_len(nrow(rst)), function(ii) paste0("epi_", method, "_", ii))
      rownames(rst) <- seq_len(nrow(rst))
      rst <- rst[ , c(12, 11, 1:10)]
      
      ## Add CREs and epi_region_id
      rst$CRE_type <- rst$CRE <- rst$epi_region_id <- NA
      rst_c <- rst
      rst_c <- tryCatch({
        rstGR <- GenomicRanges::makeGRangesFromDataFrame(rst)
        ensembldb::seqlevelsStyle(rstGR) <- "UCSC" ## Ensure chromosomes have the same format
        over <- GenomicRanges::findOverlaps(rstGR, candRegsGR)
        
        rst$epi_region_id[S4Vectors::from(over)] <- names(candRegsGR[S4Vectors::to(over)])
        rst$CRE[S4Vectors::from(over)] <- candRegsGR[S4Vectors::to(over)]$CRE
        rst$CRE_type[S4Vectors::from(over)] <- candRegsGR[S4Vectors::to(over)]$CRE_type
        
        rst
      }, error = function(e) { rst })
      
      ## Convert rst into a tibble class
      rst <- tibble::as_tibble(rst_c)
      return(rst)
    }
  #})

