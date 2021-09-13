#' @title Generates a GRanges object 
#' @description  This function makes a GRanges object from a \code{GenomicRatioSet}.
#' @param methy a \code{GenomicRatioSet} object containing the control and
#' case samples used in \link[epimutacions]{epimutations} or 
#' \link[epimutacions]{epimutations_one_leave_out} function.  
#' @param cpg_ids a character string specifying the name of the
#' CpGs in the DMR of interest. 
#' @return  The function returns a GRanges object containing the beta values and the
#' genomic ranges of the CpGs of interest. 
#' @importFrom methods is
create_GRanges_class <- function(methy, cpg_ids){
  
  #Identify type of input and extract required data:
  # * The object class  must be 'GenomicRatioSet'
  # * beta values and genomic ranges
  if(is(methy,"GenomicRatioSet")) {
    betas <- minfi::getBeta(methy)
    fd <- as.data.frame(SummarizedExperiment::rowRanges(methy))
    rownames(fd) <- rownames(betas)
  } else {
    stop("Input data 'methy' must be a 'GenomicRatioSet'")
  }
  
  #CpG names
  cpg_ids <- unlist(strsplit(cpg_ids, ","))
  #Genomic ranges
  fd <- fd[cpg_ids,]
  fd <- cols_names(fd, cpg_ids_col = FALSE) #common cols names (epi_plot)
  #Beta values matrix
  betas <- betas[cpg_ids,]
  rownames(fd) <- rownames(betas)
  #Generate the GenomicRanges class object
  gr <- GenomicRanges::makeGRangesFromDataFrame(fd)
  S4Vectors::values(gr) <- betas
  return(gr)
}
#' @title Sets common column names in a data frame
#' @description  Sets common column names in a given data frame 
#' containing the CpGs genomic ranges or a DMR (result of
#' \link[epimutacions]{epimutations} or 
#' \link[epimutacions]{epimutations_one_leave_out} function).
#' @param x a data frame containing the genomic ranges or 
#' a DMR (a row of the results of \link[epimutacions]{epimutations} or 
#' \link[epimutacions]{epimutations_one_leave_out} function). 
#' @param cpg_ids_col a boolean, if TRUE the input data frame contains
#' the CpGs names column. 
#' @return The function returns a data frame containing the column names
#' to carry out the analysis without any error. 
#' 

cols_names <- function(x, cpg_ids_col = FALSE){
  
  #Accepted names for 'dmr' column names
  seqnames_field <- c("seqnames", "seqname",
                      "chromosome", "chrom",
                      "chr", "chromosome_name",
                      "seqid")
  start_field <- "start"
  end_field <- c("end", "stop")
  strand_field <- c("strand")
  cpg_field <- c("cpg_ids", "cpgnames", "cpg")
  sample_field <- "sample"
  
  #Identify the column name for each argument
  if(seqnames_field[1] %in% colnames(x) | seqnames_field[2] %in% colnames(x) | seqnames_field[3] %in% colnames(x) | seqnames_field[4] %in% colnames(x) | seqnames_field[5] %in% colnames(x) | seqnames_field[6] %in% colnames(x) | seqnames_field[7] %in% colnames(x)){
    seqnames_pos <- which(colnames(x) %in% seqnames_field)
    seqnames <- x[,seqnames_pos]
  }else{
    stop("Chromosome column name must be specified as: 'seqnames', 'seqname', 'chromosome', 'chrom', 'chr', 'chromosome_name' or 'seqid'")
  }
  if(start_field %in% colnames(x)){
    start_pos <- which(colnames(x) %in% start_field)
    start <- x[,start_pos]
  }else{
    stop("Start column name  must be specified as: 'start'")
  }
  if(end_field[1] %in% colnames(x) | end_field[2] %in% colnames(x)){
    end_pos <- which(colnames(x) %in% end_field)
    end <- x[,end_pos]
  }else{
    stop("End column name  must be specified as: 'end' or 'stop'")
  }
  
  if(cpg_ids_col == TRUE){
    
    if(sample_field %in% colnames(x)){
      sample_pos <- which(colnames(x) %in% sample_field)
      sample <- x[,sample_pos]
    }else{
      stop("Sample column name  must be specified as: 'sample'")
    }
    
    if(cpg_field[1] %in% colnames(x) | cpg_field[2] %in% colnames(x) |cpg_field[3] %in% colnames(x)){
      cpg_pos <- which(colnames(x) %in% cpg_field)
      cpg_ids <- x[,cpg_pos]
    }else{
      stop("CpGs column name  must be specified as: 'cpg_ids', 'cpgnames' or 'cpg'")
    }
    df <- data.frame("sample" = sample, "seqnames" = seqnames, "start" = start, "end" = end, "cpg_ids" = cpg_ids)
  }else{
    if( strand_field %in% colnames(x)){
      strand_pos <- which(colnames(x) %in% strand_field)
      strand <- x[,strand_pos]
    }else{
      stop("In feature dataset strand column name  must be specified as: 'strand'")
    }
    df <- data.frame("seqnames" = seqnames, "start" = start, "end" = end, "strand" = strand)
  }
  #establish common column names
  x_rownames <- rownames(x)
  rownames(df) <- x_rownames
  return(df)
}

#' @title Computes beta values, standard deviation and mean 
#' to plot the epimutation
#' @description  Computes the beta values, population mean and 
#'  1, 1.5, and 2 standard deviations from the mean of the distribution
#'  necessary to plot the epimutations.
#' @param gr  a GRanges object obtained from \link[epimutacions]{create_GRanges_class}
#' function. 
#' @return The function returns a list containing the melted beta values, 
#' the population mean and 1, 1.5, and 2 standard deviations 
#' from the mean of the distribution.
#' 

betas_sd_mean <- function(gr){
  
  df <- as.data.frame(gr)
  colnames(df) <- c("seqnames", "start", "end", "width", "strand",colnames(GenomicRanges::elementMetadata(gr)))
  #Compute: 
  # * beta values
  # * Population mean 
  # * 1, 1.5, and 2 standard deviations from the mean of the distribution
  betas <- as.data.frame(S4Vectors::values(gr))
  mean <- rowMeans(betas)
  sd <- apply(betas,1,sd)
  sd_1_lower <- abs(mean - sd)
  sd_1_upper <- mean + sd
  sd_1.5_lower <- mean - 1.5*sd
  sd_1.5_lower <- ifelse(sd_1.5_lower > 0, sd_1.5_lower, 0)
  sd_1.5_upper <- mean +  1.5*sd
  sd_2_lower <- mean - 2*sd
  sd_2_lower <- ifelse(sd_2_lower > 0, sd_2_lower, 0)
  sd_2_upper <- mean +  2*sd
  
  sd <- cbind(df[,c("seqnames", "start","end","width","strand")], sd_1_lower, sd_1_upper, sd_1.5_lower,sd_1.5_upper,sd_2_lower,sd_2_upper)
  mean <- cbind(df[,c("seqnames", "start","end","width","strand")], mean)
  
  
  #Melt beta values, mean and sd object (necessary for the ggplot)
  
  beta_values <- reshape2::melt(df, id = c("seqnames", "start", "end", "width", "strand"))
  mean <- reshape2::melt(mean, id = c("seqnames", "start", "end", "width", "strand","mean"))
  sd <- reshape2::melt(sd, id = c("seqnames", "start", "end", "width", "strand","sd_1_lower", "sd_1_upper", "sd_1.5_lower","sd_1.5_upper","sd_2_lower","sd_2_upper"))
  
  #Create the output list
  output <- list("beta_values" = beta_values, "mean" = mean, "sd" = sd)
  return(output)
}

#' @title UCSC gene annotations
#' @description  UCSC gene annotations for a given genome assembly. 
#' @param genome  genome asambly. Can be set as: 
#' \code{'hg38'}, \code{'hg19'} and \code{'hg18'}. 
#' @return The function returns gene annotations for the specified genome assembly.
#' 


UCSC_annotation <- function(genome = "hg19"){
  if(genome == "hg19" & requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")){
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if (genome == "hg38" & requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if(genome == "hg18" & requireNamespace("TxDb.Hsapiens.UCSC.hg18.knownGene")){
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene
  } else{
    warning("Genes are not shown since TxDb database is not installed in you computer")
    txdb <- NULL
  }
  
  if(!is.null(txdb)){
    all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
    all_genes$symbol <- AnnotationDbi::mapIds(Homo.sapiens::Homo.sapiens, 
                                              keys = all_genes$gene_id,
                                              keytype = "ENTREZID",
                                              column = "SYMBOL")
  }else{
    all_genes <- NULL
  }
  
  return(all_genes)
}

#' @title UCSC annotation
#' @description  UCSC annotations for CpG Islands, H3K27Ac and H3K4Me3
#'  for a given genome assembly and genomic coordinates. 
#' @param genome  genome asambly. Can be set as: 
#' \code{'hg38'}, \code{'hg19'} and \code{'hg18'}. 
#' @param chr 
#' @param chr a character string containing the sequence names to be analysed.
#' @param from,to scalar, specifying the range of genomic coordinates.
#' Note that \code{from} cannot be larger than \code{to}. 
#' @return \code{UCSC_regulation} returns
#' a list containing CpG Islands, H3K27Ac and H3K4Me3 tacks. 
#' 
UCSC_regulation <- function(genome, chr, from, to){
  
  #cpgIslands
  cpgIslands <- suppressWarnings(Gviz::UcscTrack(genome = genome, chromosome = chr,
                                                 track = "CpG Island", from = from,
                                                 to = to, trackType = "AnnotationTrack",
                                                 start = "chromStart", end = "chromEnd",
                                                 id = "name", shape = "box",
                                                 fill = "#FA9114", name = "CpG",
                                                 background.title = "#9D5D10", 
                                                 rotation.title = 0))
  
  #H3K27Ac, H3K4Me3 and H3K27Me3
  mySession <-  rtracklayer::browserSession("UCSC")
  rtracklayer::genome(mySession) <- "hg19"
  granges <- GenomicRanges::GRanges(chr, 
                                    IRanges::IRanges(from, to))
  
  #H3K27Ac
  H3K27Ac <- suppressWarnings(rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, 
                                                    track = "Layered H3K27Ac",
                                                    range = granges)))
  H3K27Ac$seqnames <- chr
  value <- H3K27Ac$value
  H3K27Ac <- GenomicRanges::makeGRangesFromDataFrame(H3K27Ac)
  S4Vectors::values(H3K27Ac) <- value 
  H3K27Ac <- Gviz::DataTrack(H3K27Ac, 
                             type = "hist", window = "auto",
                             col.histogram = "darkblue",
                             fill.histogram = "darkblue", data = "X", 
                             name = "H3K27Ac", chr = chr,
                             background.title = "#C0E4B0")
  
  
  #H3K4Me3
  H3K4Me3 <- suppressWarnings(rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, 
                                                                                track = "Layered H3K4Me3",
                                                                                range = granges)))
  H3K4Me3$seqnames <- chr
  value <- H3K4Me3$value
  H3K4Me3 <- GenomicRanges::makeGRangesFromDataFrame(H3K4Me3)
  S4Vectors::values(H3K4Me3) <- value 
  H3K4Me3 <- Gviz::DataTrack(H3K4Me3, 
                             type = "hist", window = "auto",
                             col.histogram = "darkred",
                             fill.histogram = "darkred", data = "X", 
                             name = "H3K4Me3", chr = chr,
                             background.title = "#C0E4B0")

  #H3K27Me3
  
  if(genome == "hg19"){
    ah <- AnnotationHub::AnnotationHub()
    H3K27Me3 <- AnnotationHub::query(ah , c("UCSC", "H3K27me3", "hg19"))
    H3K27Me3 <- Gviz::DataTrack(H3K27Me3[["AH23260"]], 
                                type = "hist", window = "auto",
                                col.histogram = "darkgreen",
                                fill.histogram = "darkgreen", data = "score", 
                                name = "H3K27Me3", chr = chr, start = from, 
                                end = to,
                                background.title = "#C0E4B0")
  }

  
  return(list("cpgIslands" = cpgIslands, "H3K27Ac" = H3K27Ac, "H3K4Me3" = H3K4Me3, "H3K27Me3" = H3K27Me3))
}