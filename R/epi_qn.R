qn_norm <- function(betas, fd, qn = TRUE) {
	# if(median_na) {
	# 	filled <- t(apply(betas, 1, function(row) {
	# 		m <- median(row, na.rm = TRUE)
	# 		row[is.na(row)] <- m
	# 		row
	# 	}))
	# } else {
	# 	filled <- betas
	# }
	
	if(qn) { # identify the center as median and the scale as Qn score
		center <- apply(betas, 1, median, na.rm = TRUE)
		vari   <- apply(betas, 1, function(row) {
			robustbase::Qn(row[!is.na(row)])
		})
	} else { # identify the center as mean and the scañe as sd
		center <- apply(betas, 1, mean, na.rm = TRUE)
		vari   <- apply(betas, 1, sd, na.rm = TRUE)
	}
	
	scaled <- t(scale(t(betas), center = center, scale = vari))
	
	return(scaled)
}

qn_bump <- function(nbetas, fd, window_sz = 1000, cutoff = 0.05) {
	chr <- as.character(fd$seqnames)
	pos <- as.numeric(fd$start)
	cl <- bumphunter::clusterMaker(chr, pos, maxGap = window_sz)
	reg <- lapply(as.data.frame(nbetas), function(sam) {
		suppressMessages(bumphunter::regionFinder(sam, chr, pos, cl, cutoff = cutoff))
	})
	names(reg) <- colnames(nbetas)
	return(reg)
}

qn_outlier <- function(case, reg, nbetas, fd, min_cpg = 3, th = 3) {
	reg <- reg[[case]]
	reg <- reg[reg$L >= min_cpg, ]
	
	
	do.call(rbind, lapply(seq_len(nrow(reg)), function(ii) {
		sel <- reg[ii, ]
		
		cpgs <- rownames(fd[ fd$seqnames == sel$chr & fd$start >= sel$start & fd$end <= sel$end, ])
		
		# Check outilers into the top side (right side)
		cpg_top <- cpgs[nbetas[cpgs, case] > th]
		cpg_top_cnt <- sum(nbetas[cpgs, case] > th)
		cpg_top_mean <- mean(nbetas[cpgs, case][nbetas[cpgs, case] > th])
		cpg_top_mean <- ifelse(is.nan(cpg_top_mean), NA, cpg_top_mean)
		
		# Check outliers into the bottom side (left side)
		cpg_btm <- cpgs[nbetas[cpgs, case] < (th * -1)]
		cpg_btm_cnt <- sum(nbetas[cpgs, case] < (th * -1))
		cpg_btm_mean <- mean(nbetas[cpgs, case][nbetas[cpgs, case] < (th * -1)])
		cpg_btm_mean <- ifelse(is.nan(cpg_btm_mean), NA, cpg_btm_mean)
		
		outlier_direction <- ""
		outlier_score <- NA
		if(length(cpg_top) == 0 & length(cpg_btm) == 0) {
			cpg_labels <- ""
		} else if(length(cpg_top) != 0 & length(cpg_btm) != 0) {
			stop(paste0("For patient '", patient, "' found region (", sel$chr, "-", sel$start, ":", sel$end, ") with CpGs as right and left outliers."))
		} else if(length(cpg_top) != 0) {
			cpg_labels <- paste(cpg_top, collapse = ",", sep = "")
			outlier_score <- cpg_top_mean
			outlier_direction <- "hypermethylation"
		} else {
			cpg_labels <- paste(cpg_btm, collapse = ",", sep = "")
			outlier_score <- cpg_btm_mean
			outlier_direction <- "hypomethylation"
		}
		
		data.frame(chromosome = sel$chr, start = sel$start, end = sel$end,
			length = sel$end - sel$start, N_CpGs = length(cpgs), 
			CpG_ids = cpg_labels,
			outlier_score = outlier_score,
			outlier_significance = NA,
			outlier_direction = outlier_direction
		)
	}))
}


