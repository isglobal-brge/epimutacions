UCSC_annotation <- function(genome = "hg19"){
	if(genome == "hg19" & requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")){
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
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
