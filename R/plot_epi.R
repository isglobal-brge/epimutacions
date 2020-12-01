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
			Gviz::plotTracks(list(ideo_track, genome_track, dmrs_track, cpgs_track, betas_track, gene_track), from = min(start(gr)), to = max(end(gr)))
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
