#' @title  Plot a given epimutation and locate it along the genome
#' @description This function plots a given epimutation
#' and UCSC annotations for the specified genomic region.  
#' @param dmr epimutation obtained as a result of \link[epimutacions]{epimutations}
#' function. 
#' @param methy a GenomicRatioSet object containing the information
#' of control and case samples used for the analysis in the \link[epimutacions]{epimutations}
#' function. See the constructor function \link[minfi]{GenomicRatioSet},
#' \link[minfi]{makeGenomicRatioSetFromMatrix}.
#' @param genome a character string specifying the genome of reference. 
#' It can be set as \code{"hg38"},\code{"hg19"} and \code{"hg18"}. 
#' The default is \code{"hg19"}. 
#' @param regulation a boolean. If TRUE UCSC annotations for CpG Islands, H3K27Ac and H3K4Me3
#' are plotted. The default is FALSE. The running process when \code{regulation} is TRUE 
#' can take several minutes.
#' @param from,to scalar, specifying the range of genomic coordinates 
#' for the plot of gene annotation region.
#' If \code{NULL} the plotting ranges are derived from the individual track. 
#' Note that \code{from} cannot be larger than \code{to}. 
#' @details 
#' The tracks are plotted vertically. Each track is separated by different background
#' colour and a section title. The colours and titles are preset and cannot be set by 
#' the user. 
#' 
#' Note that if you want to see the UCSC annotations maybe you need to take a bigger
#' genomic region. 
#' 
#' @return The function returns a plot divided in two parts: 
#' * ggplot graph including the individual with 
#' the epimutation in red, the control samples in dashed black lines and
#' population mean in blue. Grey shaded regions indicate 1, 1.5 and 2 standard
#' deviations from the mean of the distribution. 
#' * UCSC gene annotations for the specified genomic region. 
#' 
#' @examples 
#' \dontrun{
#' library(epimutacions)
#' data(methy)
#' 
#' #Find epimutations in GSM2562701 sample of methy dataset
#' 
#' case_sample <- methy[,"GSM2562701"]
#' control_panel <- methy[,-51]
#' 
#' results <- epimutations(case_sample, control_panel, method = "manova")
#' 
#' plot_epimutations(as.data.frame(results[1,]), methy)
#' }
#' 
#' @export
plot_epimutations <- function(dmr, methy, genome = "hg19", regulation = FALSE, from = NULL, to = NULL){
  
  # Identify type of input and extract required data:
  #	* Input parameters: class, not null (if requered), nrow/ncols... 
  #	* sample's classification
  #	* feature annotation
  
  ## NULL arguments
  if(is.null(dmr)){
    stop("The argument 'dmr' must be introduced")
  }
  if(is.null(methy)){
    stop("The argument 'beta' must be introduced")
  }
  if(is.null(genome)){
    stop("The argument 'genome' must be introduced")
  }
  ##Unique DMR  
  if(nrow(dmr) > 1){
    warning("more than one DMR introduced (nrow > 1) only the first element will be used")
    dmr <- dmr[1,]
  }
  ##Genome assembly
  if(genome != "hg38" & genome != "hg19" & genome != "hg18"){
    stop("Argument 'genome' must be 'hg38', 'hg19' or 'hg18'")
  }
  ##Epimutation start('from') and end('to') possitions
   ## * 'from' and 'to' introduced together
  if(is.null(from) & !is.null(to) | !is.null(from) & is.null(to)){
    stop("Arguments 'from' and 'to' must be provided together")	
  }
   ## * 'from' is smaller than 'to' 
  if(!is.null(from) & !is.null(to)){
    if(from > to){
      stop("The value of argument 'from' must be smaller than 'to'")	
    }
  }
  
  # DMR column names must be always the same (set the common column names)
  dmr  <- cols_names(dmr, cpg_ids_col = TRUE)  #epi_plot
  
  # Set 'from' and 'to' arguments value
  if(is.null(from) & is.null(to)){
    from <- dmr$start - 1000
    to <- dmr$end + 1000
  }
  
  #Generate GenomicRanges object to contain in the same object:
  ## * Genomic ranges of each CpG in the DMR
  ## * Beta values 
  gr <- create_GRanges_class(methy, dmr[,"cpg_ids"]) #epi_plot
  betas_sd_mean <- betas_sd_mean(gr) #epi_plot
  
  #Generate variables in 'beta_values' data frame containing:
  # * status: case sample name/'control'
  # * color: 'red' for case sample and 'black' for control sample
  # * lines: 'longdash' for controls and 'solid' for case and population mean
  
  betas_sd_mean$beta_values$status <- ifelse(betas_sd_mean$beta_values$variable == dmr$sample, dmr$sample, "control")
  betas_sd_mean$beta_values$color <- ifelse(betas_sd_mean$beta_values$status == dmr$sample,"red","black")
  betas_sd_mean$beta_values$lines <- ifelse(betas_sd_mean$beta_values$status == "control","longdash","solid")
  
  #Generate a variable with the CpGs names
  names <- betas_sd_mean$beta_values[betas_sd_mean$beta_values$variable == dmr$sample,]
  names$id <- names(gr)
  
  #Plot epimutations
  
  plot_betas <- ggplot2::ggplot() + 
    ggplot2::geom_line(data = betas_sd_mean$beta_values, ggplot2::aes(x = start, y = value, group = variable, color = status), linetype = betas_sd_mean$beta_values$lines) +
    ggplot2::geom_point(data = betas_sd_mean$beta_values, ggplot2::aes(x = start, y = value, group = variable, color = status), color = betas_sd_mean$beta_values$color)
  
  plot_sd <- plot_betas +
    ggplot2::geom_ribbon(data = betas_sd_mean$sd, ggplot2::aes(x = start, ymin = sd_2_lower, ymax = sd_2_upper), fill = "gray39", alpha = 0.4) +
    ggplot2::geom_ribbon(data = betas_sd_mean$sd, ggplot2::aes(x = start, ymin = sd_1.5_lower, ymax = sd_1.5_upper), fill = "gray40", alpha = 0.4) +
    ggplot2::geom_ribbon(data = betas_sd_mean$sd, ggplot2::aes(x = start, ymin = sd_1_lower, ymax = sd_1_upper), fill = "gray98", alpha = 0.4)
  
  plot_mean <-  plot_sd +
    ggplot2::geom_line(data = betas_sd_mean$mean, ggplot2::aes(x = start, y = mean, color = "mean")) +
    ggplot2::geom_point(data = betas_sd_mean$mean, ggplot2::aes(x = start, y = mean), show.legend = TRUE)
  
  plot_cpg_names <- plot_mean +
    ggrepel::geom_text_repel() + 
    ggplot2::annotate(geom="text", x=names$start, y=names$value + 0.05, label=names$id, color="black")
  
  plot <- plot_cpg_names + 
    #ggplot2::lims(y = c(0,1)) +  
    ggplot2::scale_colour_manual(name="Status", values=c("black","red","darkblue")) +
    ggplot2::theme_bw() + 
    ggplot2::labs(x = "Coordinates") + 
    ggplot2::labs(y = "DNA methylation level")
  
  #Plot gene annotations
  
  genes <- UCSC_annotation(genome) #epi_plot
  ideo_track <- Gviz::IdeogramTrack(genome = genome, chromosome = dmr$seqnames)
  genome_track <- Gviz::GenomeAxisTrack()
  gene_track <- Gviz::GeneRegionTrack(genes, 
                                      chromosome = dmr$seqnames,
                                      name = "Genes",
                                      transcriptAnnotation = "symbol",
                                      background.title = "#8F913A",
                                      rotation.title = 0)		
  
  tracks_Highlight <- Gviz::HighlightTrack(trackList = list(genome_track, gene_track),
                                           start = dmr$start, end = dmr$end,
                                           chromosome = dmr$seqnames,
                                           col = "#7EA577", fill = "#C6D7C3",
                                           alpha = 0.4,
                                           inBackground = FALSE)
  if(regulation == TRUE){
    
      sz <- to - from
      if(sz > 200000){
        stop("The region is too large (> 200kb) to download the annotations from 'UCSC'")
      }
      annotation <- UCSC_regulation(genome, dmr$seqnames, from, to)

    if(genome ==  "hg19"){
      tracks_Highlight <- Gviz::HighlightTrack(trackList = list(genome_track, gene_track, 
                                                                annotation$cpgIslands,
                                                                annotation$H3K4Me3,
                                                                annotation$H3K27Ac, 
                                                                annotation$H3K27Me3),
                                                                start = dmr$start, end = dmr$end,
                                                                chromosome = dmr$seqnames,
                                                                col = "#7EA577", 
                                                                fill = "#C6D7C3",
                                                                alpha = 0.4,
                                                                inBackground = FALSE)
      
    }else{
      tracks_Highlight <- Gviz::HighlightTrack(trackList = list(genome_track, gene_track, 
                                                                annotation$cpgIslands,
                                                                annotation$H3K4Me3,
                                                                annotation$H3K27Ac),
                                               start = dmr$start, end = dmr$end,
                                               chromosome = dmr$seqnames,
                                               col = "#7EA577", fill = "#C6D7C3",
                                               alpha = 0.4,
                                               inBackground = FALSE)
    }
    
  }
  
  #Plot window
  
  dev.new(width = 1080, height = 1350, unit = "px")
  p1 <- plot
  p2 <- grid::grid.grabExpr(Gviz::plotTracks(list(ideo_track, tracks_Highlight), 
                                               from = from, 
                                               to = to, add = TRUE))
  gridExtra::grid.arrange(grobs = list(p1,p2), row=2)
}



