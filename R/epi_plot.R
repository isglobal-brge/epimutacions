epi_plot <- function(dmr, methy, genome = "hg19", from = NULL, to = NULL){
  #Comprobations
  ## NULL elements
  if(is.null(dmr)){
    stop("The argument 'dmr' must be introduced")
  }
  if(is.null(methy)){
    stop("The argument 'beta' must be introduced")
  }
  if(is.null(genome)){
    stop("The argument 'genome' must be introduced")
  }
  ##Dimensions 
  if(nrow(dmr) > 1){
    warning("more than one DMR introduced (nrow > 1) only the first element will be used")
    dmr <- dmr[1,]
  }
  ##Genome posible options
  if(genome != "hg19" & genome != "hg18"){
    stop("Argument 'genome' must be 'hg19' or 'hg18'")
  }
  ##From and to arguments
  if(is.null(from) & !is.null(to) | !is.null(from) & is.null(to)){
    stop("Arguments 'from' and 'to' must be provided together")	
  }
  if(!is.null(from) & !is.null(to)){
    if(from > to){
      stop("The value of argument 'from' must be smaller than  'to'")	
    }
  }
  
  #Select dmr colnames and columns
  dmr <- .dmr_cols(dmr)
  
  #Generate GenomicRanges object with the beta values in the DMR 
  gr <- create_GRanges_class(methy, dmr[,"cpg_ids"])
  df <- as.data.frame(gr)
  #Generate all the necessary objects for the plot
  betas <- as.data.frame(S4Vectors::values(gr))
  mean <- rowMeans(betas)
  sd <- apply(betas,1,sd)
  sd_1_lower <- mean - sd
  sd_1_upper <- mean + sd
  sd_1.5_lower <- mean - 1.5*sd
  sd_1.5_upper <- mean +  1.5*sd
  sd_2_lower <- mean - 2*sd
  sd_2_upper <- mean +  2*sd
  
  sd <- cbind(df[,c("seqnames", "start","end","width","strand")], sd_1_lower, sd_1_upper, sd_1.5_lower,sd_1.5_upper,sd_2_lower,sd_2_upper)
  mean <- cbind(df[,c("seqnames", "start","end","width","strand")], mean)
  
  
  #Melt
  beta_values<- reshape2::melt(df, id = c("seqnames", "start", "end", "width", "strand"))
  mean <- reshape2::melt(mean, id = c("seqnames", "start", "end", "width", "strand","mean"))
  sd <- reshape2::melt(sd, id = c("seqnames", "start", "end", "width", "strand","sd_1_lower", "sd_1_upper", "sd_1.5_lower","sd_1.5_upper","sd_2_lower","sd_2_upper"))
  ###
  beta_values$status <- ifelse(beta_values$variable == dmr$sample, dmr$sample, "control")
  beta_values$color <- ifelse(beta_values$status == dmr$sample,"red","black")
  beta_values$lines <- ifelse(beta_values$status == "control","longdash","solid")
  names <- beta_values[beta_values$variable == dmr$sample,]
  names$id <- names(gr)
  
  #Plot epimutations
  
  plot_betas <- ggplot2::ggplot() + 
    ggplot2::geom_line(data = beta_values, ggplot2::aes(x = start, y = value, group = variable, color = status), linetype = beta_values$lines) +
    ggplot2::geom_point(data = beta_values, ggplot2::aes(x = start, y = value, group = variable, color = status), color = beta_values$color) 
  
  plot_sd <- plot_betas +
    ggplot2::geom_ribbon(data = sd, ggplot2::aes(x = start, ymin = sd_2_lower, ymax = sd_2_upper), fill = "gray39", alpha = 0.4) +
    ggplot2::geom_ribbon(data = sd, ggplot2::aes(x = start, ymin = sd_1.5_lower, ymax = sd_1.5_upper), fill = "gray40", alpha = 0.4) +
    ggplot2::geom_ribbon(data = sd, ggplot2::aes(x = start, ymin = sd_1_lower, ymax = sd_1_upper), fill = "gray98", alpha = 0.4)
  
  plot_mean <-  plot_sd +
    ggplot2::geom_line(data = mean, ggplot2::aes(x = start, y = mean, color = "mean")) +
    ggplot2::geom_point(data = mean, ggplot2::aes(x = start, y = mean), show.legend = TRUE)
  
  plot_cpg_names <- plot_mean +
    ggrepel::geom_text_repel() + 
    ggplot2::annotate(geom="text", x=names$start, y=names$value + 0.05, label=names$id, color="black")
  
  plot <- plot_cpg_names + 
    #ylim(0,1) +  
    ggplot2::scale_colour_manual(name="Status",values=c("black","red","darkblue"))
  
  #Plot gene annotations
  
  genes <- UCSC_annotation(genome)
  ideo_track <- Gviz::IdeogramTrack(genome = genome, chromosome = dmr$seqnames)
  genome_track <- Gviz::GenomeAxisTrack()
  gene_track <- Gviz::GeneRegionTrack(genes, 
                                      chromosome = dmr$seqnames,
                                      name = "Genes",
                                      transcriptAnnotation = "symbol",
                                      background.title = "#7EA577")				
  
  #Plot window
  
  dev.new(width = 1080, height = 1350, unit = "px")
  p1 <- plot
  if(is.null(from) & is.null(to)){
    p2 <- grid::grid.grabExpr(Gviz::plotTracks(list(ideo_track,genome_track,gene_track), 
                                               from = dmr$start - 100000, 
                                               to = dmr$end + 100000, 
                                               add = TRUE))
  }else{
    p2 <- grid::grid.grabExpr(Gviz::plotTracks(list(ideo_track,
                                                    genome_track,
                                                    gene_track), 
                                               from = from, 
                                               to = to, add = TRUE))
  }
  
  gridExtra::grid.arrange(grobs = list(p1,p2), row=2)
}

# Helper functions

.dmr_cols <- function(dmr){
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
    stop("Sample column name  must be specified as: 'sample'")
  }
  if(seqnames_field[1] %in% colnames(dmr) | seqnames_field[2] %in% colnames(dmr) | seqnames_field[3] %in% colnames(dmr) | seqnames_field[4] %in% colnames(dmr) | seqnames_field[5] %in% colnames(dmr) | seqnames_field[6] %in% colnames(dmr) | seqnames_field[7] %in% colnames(dmr)){
    seqnames_pos <- which(colnames(dmr) %in% seqnames_field)
    seqnames <- dmr[,seqnames_pos]
  }else{
    stop("Chromosome column name must be specified as: 'seqnames', 'seqname', 'chromosome', 'chrom', 'chr', 'chromosome_name' or 'seqid'")
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
    stop("CpGs column name  must be specified as: 'cpg_ids', 'cpgnames' or 'cpg'")
  }
  
  dmr_rownames <- rownames(dmr)
  dmr <- data.frame("sample" = sample,"seqnames" = seqnames, "start" = start, "end" = end, "cpg_ids" = cpg_ids)
  rownames(dmr) <- dmr_rownames
  return(dmr)
}