dmr_cols <- function(dmr){
  seqnames_field <- c("seqnames", "seqname",
                      "chromosome", "chrom",
                      "chr", "chromosome_name",
                      "seqid")
  start_field <- "start"
  end_field <- c("end", "stop")
  cpg_field <- c("cpg_ids", "cpgnames", "cpg")
  sample_field <- "sample"
  
  if(sample_field %in% colnames(dmr)){
    sample_pos <- which(colnames(dmr) %in% sample_field)
    sample <- dmr[,sample_pos]
  }else{
    stop("Chromosome column name must be specified as: 'seqnames', 'seqname', 'chromosome', 'chrom', 'chr', 'chromosome_name' or 'seqid'")
  }
  if(seqnames_field[1] %in% colnames(dmr) | seqnames_field[2] %in% colnames(dmr) | seqnames_field[3] %in% colnames(dmr) | seqnames_field[4] %in% colnames(dmr) | seqnames_field[5] %in% colnames(dmr) | seqnames_field[6] %in% colnames(dmr) | seqnames_field[7] %in% colnames(dmr)){
    seqnames_pos <- which(colnames(dmr) %in% seqnames_field)
    seqnames <- dmr[,seqnames_pos]
  }else{
    stop("Sample column name  must be specified as: 'sample'")
  }
  if(start_field %in% colnames(dmr)){
    start_pos <- which(colnames(dmr) %in% start_field)
    start <- dmr[,start_pos]
  }else{
    stop("Start column name  must be specified as: 'start'")
  }
  if(end_field[1] %in% colnames(dmr) | end_field[2] %in% colnames(dmr)){
    end_pos <- which(colnames(dmr) %in% end_field)
    end <- dmr[,end_pos]
  }else{
    stop("End column name  must be specified as: 'end' or 'stop'")
  }
  if(cpg_field[1] %in% colnames(dmr) | cpg_field[2] %in% colnames(dmr) |cpg_field[3] %in% colnames(dmr)){
    cpg_pos <- which(colnames(dmr) %in% cpg_field)
    cpg_ids <- dmr[,cpg_pos]
  }else{
    stop("CpGs column name  must be specified as: 'cpg_ids', 'cpgnames', 'cpg'")
  }
  
  dmr <- data.frame("sample" = sample,"seqnames" = seqnames, "start" = start, "end" = end, "cpg_ids" = cpg_ids)
  return(dmr)
}