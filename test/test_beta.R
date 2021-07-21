library(epimutacions)
library(SummarizedExperiment)
library(minfi)
library(BiocParallel)

data("methy")

methy.fil <- methy[rowMeans(is.na(getBeta(methy))) < 0.05, ]


case_samples <- methy.fil[,methy.fil$status == "case"]
control_samples <- methy.fil[,methy.fil$status == "control"]

epi_mvo <- epimutations(case_samples, control_samples, method = "manova")
mvo.fil <- subset(epi_mvo, pvalue < 0.05/40408)
epi_bar <- epimutations(case_samples[1:7000, ], control_samples[1:7000, ], method = "barbosa")

epi_bdist <- epimutations(case_samples, control_samples, method = "beta")
par <- epi_parameters()
par$beta$diff_threshold <- 0.3
epi_bdist <- epimutations(case_samples, control_samples, method = "beta", epi_params = par)

par$manova$pvalue_cutoff <- 1e-50
epi_mvo2 <- epimutations(case_samples, control_samples, method = "manova", epi_params = par)

bd <- makeGRangesFromDataFrame(subset(epi_bdist, sample == "GSM2562701"))
man <- makeGRangesFromDataFrame(subset(mvo.fil, sample == "GSM2562701"))
findOverlaps(man, bd)

epi_bar <- epimutations(case_samples[1:70, ], control_samples[1:70, ], method = "barbosa")


## No funciona
epimutations_one_leave_out(methy[, 1:10], BPPARAM = SnowParam(2))
microbenchmark::microbenchmark(epimutations_one_leave_out(methy[, 1:10], BPPARAM = MulticoreParam(20)),
                               epimutations_one_leave_out(methy[, 1:10]), times = 2)

