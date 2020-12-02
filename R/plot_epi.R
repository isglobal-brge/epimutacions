#' @title  Creates DMRs and outlier CpGs plot
#' @description This function plots DMRs, outlier CpGs, 
#' methylation in beta values and UCSC annotations for a specified epi-mutation
#' @param methy a \code{GenomicRatioSet} or \code{ExpressionSet}
#' @param epi_res \link{\code{epimutacions} function result
#' @param sam A character specifying the sample to plot
#' @param chr A character specifying the chromosome of the epi-mutation to be plotted  
#' @param genome The genome of reference. 
#' It can be set as \code{"hg18"} and \code{"hg19"}. The default is \code{"hg19"}. 
#' @param from,to Scalar, specifying the range of genomic coordinates of the plot.
#' If \code{NULL} the plotting ranges are derived from the individual track. 
#' Note that from cannot be larger than to. 
#' 
#' @details 
#' The tracks are plotted vertically. Each track is separated by different background
#' colour and a section title. 
#' 
#' Note that if you want to see the UCSC annotations maybe you need to take a bigger
#' genomic region. However, if there are many DMRs can be difficult to see the 
#' methylation values plot for each CpGs, in this case, you need to take a smaller genomic region. 
#' The mentioned can be adjusted using \code{from} and \code{to} parameters. 
#' 
#' @return A genomic graphic including different tracks specified by the user. 
#' 
#' @examples  
#' 
#' \dontrun{
#' data(methy)
#' # Find epi-mutations in a specific sample
#' 
#' epi_manova <- epimutacions(methy, method = "manova")
#' 
#' # Plot the identified epi-mutations
#'
#' plot_epi(methy, epi_res = epi_manova, sam = "GSM2562699", chr = "chr7", genome = "hg19") 
#' }
#' 
#' @export
plot_epi <- function(methy, epi_res, sam, chr, genome = "hg19", from = NULL, to = NULL) {
	if(length(chr) != 1) {
		stop("Argument 'chr' must be of length one (aka. 'chr10')")
	}
	if(length(sam) != 1 | !sam %in% epi_res$sample) {
		stop("Argument 'sam' must be of length one and present in 'epi_res'")
	}
	if(genome != "hg19" & genome != "hg18"){
		stop("Argument 'genome' must be 'hg19' or 'hg18'")
	}
	if(is.null(from) & !is.null(to) | !is.null(from) & is.null(to)){
		stop("Arguments 'from' and 'to' must be provided together")	
	}
	if(!is.null(from) & !is.null(to)){
		if(from > to){
		stop("The value of argument 'from' must be smaller than  'to'")	
		}
	}
	
	epi_res <- epi_res[epi_res$sample == sam & epi_res$chromosome == chr, ]
	
	ideo_track <- Gviz::IdeogramTrack(genome = genome, chromosome = chr)
	genome_track <- Gviz::GenomeAxisTrack()
	
	dmrs_track <- Gviz::AnnotationTrack(epi_res, name = "DMRs")
	Gviz::displayPars(dmrs_track) <- list(background.title = "#AAAADD", background.panel = "#EEEEFF", col = "#AAAADD", fill = "#AAAADD")
	
	gr <- create_GRanges_class(methy, epi_res, sam, chr)
	
	cpgs_track <- Gviz::AnnotationTrack(gr, name = "CpGs")
	Gviz::displayPars(cpgs_track) <- list(background.title = "#DDAAAA", background.panel = "#FFEEEE", col = "#DDAAAA", fill = "#DDAAAA")
	

	betas_track <- Gviz::DataTrack(gr, name = "Beta values", start  = min(start(gr)), end = max(end(gr)))
	groups <- ifelse(names(GenomicRanges::mcols(gr)) == sam, "case", "control")
	Gviz::displayPars(betas_track) <- list(background.panel = "#F5F8F9", col = c("#E98963","#6CAF95"), fill = c("#E98963","#6CAF95"), groups = groups, type = c("a", "p"), legend = FALSE)
	
	genes <- UCSC_annotation(genome)
	
	if(is.null(genes)){
		if(is.null(from) & is.null(to)){
			Gviz::plotTracks(list(ideo_track, genome_track, dmrs_track, cpgs_track, betas_track), from = min(start(gr)), to = max(end(gr)) + 2)
		}else{
			Gviz::plotTracks(list(ideo_track, genome_track, dmrs_track, cpgs_track, betas_track), from = from, to = to + 2)	
		}
		
	}else{
		
		gene_track <- Gviz::GeneRegionTrack(genes, chromosome = chr, start = min(start(gr)),  end = max(end(gr)), name = "Genes")
		Gviz::displayPars(gene_track) <- list(background.title = "#7EA577", background.panel = "#EEFFEE", transcriptAnnotation = "symbol")
		
		
		if(is.null(from) & is.null(to)){
			Gviz::plotTracks(list(ideo_track, genome_track, dmrs_track, cpgs_track, betas_track, gene_track), from = min(start(gr)), to = max(end(gr))+2)
		}else{
			Gviz::plotTracks(list(ideo_track, genome_track, dmrs_track, cpgs_track, betas_track, gene_track), from = from, to = to+2)	
		}
	}

	
	
	
	#if(method %in% c("barbosa", "mahdistmcd")) {
	#	Gviz::plotTracks(list(ideoTrack, gatrack, atrack, qtrack))	
	#} else if(method %in% c("manova", "mlm")) {
	#	qtrack = Gviz::DataTrack(epi_res, data = -log10(epi_res$outlier_significance),
	#							 name = "Significance [-log10(pval)]", genome = genome)
	#	Gviz::plotTracks(list(ideoTrack, gatrack, atrack, qtrack))
	#} else if(method == "isoforest") {
	#	qtrack = Gviz::DataTrack(epi_res, data = epi_res$outlier_score,
	#							 name = "Score", genome = genome)
	#	Gviz::plotTracks(list(ideoTrack, gatrack, atrack, qtrack))
	#}
	
}
