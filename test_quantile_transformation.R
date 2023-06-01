### Test approximation to normalize quantile

library(tidyverse)
library(preprocessCore)
library(minfi)


## Generate data ####
n <- 1000
meanA <- 0.1
sdA <- 5e-4
meanB <- 0.15
sdB <- 5e-4

simulateData <- function(meanA, meanB, sd, n, n_outs, diff_outs){
  
  termA <- (meanA * (1 - meanA)) / sd
  distA <- rbeta(n,  meanA * (termA - 1), (1 - meanA) * (termA - 1))
  
  termB <- (meanB * (1 - meanB)) / sd
  distB <- rbeta(n,  meanB * (termB - 1), (1 - meanB) * (termB - 1))
  ## Añadir outliers
  distB <- c(distB, runif(n_outs, meanB + diff_outs - 0.1, meanB + diff_outs + 0.1))
  
  list(A = distA, B = distB)  
  
}


vals <- getBeta(methy["cg22469274", ])
quantiles <-  quantile(getBeta(brge_methy["cg22469274", ]), probs =  seq(0.05, 0.95, 0.05))

normalizeQuantile <- function(vals, quantiles){
  
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center]
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(valsf, ncol = 1), quantiles)
  moddf <- data.frame(Q = as.vector(valsQ), f = as.vector(valsf))

  diff <- as.vector(valsf) - valsQ
  mod <- lm(Q ~ f,  data = moddf)
  bpars <- getBetaParams(matrix(valsf, ncol = 1))
  ps <- pbeta(vals, bpars[1], bpars[2])
  quants <- pmax(quants, 1 - quants)
  vals_out <- (1 - quants)/0.5*(predict(mod, data.frame(f = as.vector(vals)))) + (quants - 0.5)/0.5*(vals - mean(diff))
  
  vals_out[vals_out < 1e-3] <- 1e-3
  vals_out[vals_out > 1-1e-3] <- 1-1e-3
  vals_out

  }



normalizeBeta <- function(vals, beta_params){
  
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center]
  bpars <- getBetaParams(matrix(vals, ncol = 1))
  
  ps <- pbeta(vals, bpars[1], bpars[2])
  vals_pit <- qbeta(ps, beta_params[1], beta_params[2])
  
}

#### New distribution higher ####
vals_high <- simulateData(0.1, 0.2, 5e-4, 1000, 10, 0.1)

qA <- minfi::logit2(quantile(vals_high$A, seq(0.05, 0.95, 0.05)))
distB_beta <- normalizeBeta(vals_high$B, getBetaParams(matrix(vals_high$A, ncol = 1)))

data.frame(Distribution = c(vals_high$A, vals_high$B, distB_beta), 
           Group = rep(c("A", "B", "B normalize"), c(lengths(vals_high), length(distB_beta)))) %>%
  ggplot(aes(x = Group, y = Distribution)) +
  geom_boxplot() +
  theme_bw()


data.frame(Distribution = c(vals_high$A, vals_high$B, distB_beta), 
           Group = rep(c("A", "B", "B normalize"), c(lengths(vals_high), length(distB_beta)))) %>%
  ggplot(aes(color = Group, x = Distribution)) +
  geom_density() +
  theme_bw()


#### New distribution higher ####
vals_low <- simulateData(0.2, 0.1, 5e-4, 1000, 10, 0.1)

qA <- minfi::logit2(quantile(vals_low$A, seq(0.05, 0.95, 0.05)))
distB_norm <- normalizeBeta(vals_low$B, getBetaParams(matrix(vals_low$A, ncol = 1)))

data.frame(Distribution = c(vals_low$A, vals_low$B, distB_norm), 
           Group = rep(c("A", "B", "B normalize"), c(lengths(vals_low), length(distB_norm)))) %>%
  ggplot(aes(x = Group, y = Distribution)) +
  geom_boxplot() +
  theme_bw()


data.frame(Distribution = c(vals_low$A, vals_low$B, distB_norm), 
           Group = rep(c("A", "B", "B normalize"), c(lengths(vals_low), length(distB_norm)))) %>%
  ggplot(aes(color = Group, x = Distribution)) +
  geom_density() +
  theme_bw()



## Real datasets ####
library(ExperimentHub)
library(meffil)
library(limma)
library(brgedata)
library(epimutacions)
eh <- ExperimentHub()


### Load data
reference_panel <- eh[["EH6691"]]
methy <- eh[["EH6690"]]
data("brge_methy", package = "brgedata")
brge_methy


## Process ref_methy
ref_methy <- preprocessFunnorm(reference_panel)
ref_methy$dataset <- "Reference"
assay(ref_methy, "M") <- getM(ref_methy)

## Combine datasets
methy$dataset <- "methy"
brge_methy$dataset <- "BRGE"

com_cpgs <- intersect(rownames(methy), rownames(brge_methy))
comb <- makeGenomicRatioSetFromMatrix(cbind(getBeta(methy[com_cpgs, ]), getBeta(brge_methy[com_cpgs, ])))
comb$Dataset <- rep(c("methy", "BRGE"), c(ncol(methy), ncol(brge_methy)))

## Explore global methylation patterns ####
pc_comb <- meffil.methylation.pcs(getBeta(comb), probe.range = 40000, full.obj = TRUE)

pc_comb.vars <- pc_comb$sdev^2/sum(pc_comb$sdev^2)
pc_comb_df <- data.frame(pc_comb$x[, 1:10], dataset = comb$Dataset)

ggplot(pc_comb_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point() +
  scale_x_continuous(name = paste0("PC1 (", round(pc_comb.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_comb.vars[2]*100, 1), "%)")) +
  theme_bw()


## Get CpGs with differences in methylation
mod <- model.matrix(~ Dataset, colData(comb))
comb_fit <- lmFit(getBeta(comb), mod) %>% eBayes()
comb_res <- topTable(comb_fit, coef = 2, n = Inf)

top_cpgs <- rownames(comb_res)[1:9]
getBeta(comb[top_cpgs, ]) %>%
  t() %>%
  data.frame() %>%
  mutate(dataset = comb$Dataset,
         sample = colnames(comb)) %>%
  gather(CpG, Methylation, 1:9) %>%
  ggplot(aes(x = dataset, y = Methylation)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  facet_wrap(~ CpG)

ref_quantiles <- rowQuantiles(getBeta(brge_methy), probs =  c(0.005, 0.995, seq(0.05, 0.95, 0.05)))

beta_params <- getBetaParams(t(getBeta(brge_methy[com_cpgs, ])))
mod_methy <- integrateReferenceQuantile(getBeta(methy[com_cpgs, ]), beta_params)



plot(getBeta(methy[top_cpgs[1], ]), mod_methy[1, ])


getBeta(comb[com_cpgs, ]) %>%
  cbind(mod_methy) %>%
  t() %>%
  data.frame() %>%
  mutate(dataset = c(comb$Dataset, rep("methy_mod", ncol(methy)))) %>%
  gather(CpG, Methylation, 1:9) %>%
  ggplot(aes(x = dataset, y = Methylation)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  facet_wrap(~ CpG)

# Compute modified CpGs ####
mod_methy <- integrateReferenceQuantile(getBeta(methy[com_cpgs, ]), ref_quantiles[com_cpgs, ])


comb_mod <- makeGenomicRatioSetFromMatrix(cbind(mod_methy, getBeta(brge_methy[com_cpgs, ])))
comb_mod$Dataset <- rep(c("methy", "BRGE"), c(ncol(mod_methy), ncol(brge_methy)))


getBeta(comb[strsplit(epi_ori$cpg_ids[2], ",")[[1]], ]) %>%
  cbind(mod_methy[strsplit(epi_ori$cpg_ids[2], ",")[[1]],]) %>%
  t() %>%
  data.frame() %>%
  mutate(dataset = c(comb$Dataset, rep("methy_mod", ncol(methy)))) %>%
  gather(CpG, Methylation, 1:6) %>%
  ggplot(aes(x = dataset, y = Methylation)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  facet_wrap(~ CpG)

## Test ####
normalizeQuantile(getBeta(methy["cg22469274", ]), 
                  quantile(getBeta(brge_methy["cg22469274", ]), probs =  seq(0.05, 0.95, 0.05)))




## Explore global methylation patterns ####
pc_comb_mod <- meffil.methylation.pcs(getBeta(comb_mod), probe.range = 40000, full.obj = TRUE)

pc_comb_mod.vars <- pc_comb_mod$sdev^2/sum(pc_comb_mod$sdev^2)
pc_comb_mod_df <- data.frame(pc_comb_mod$x[, 1:10], dataset = comb_mod$Dataset)

ggplot(pc_comb_mod_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point() +
  scale_x_continuous(name = paste0("PC1 (", round(pc_comb_mod.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_comb_mod.vars[2]*100, 1), "%)")) +
  theme_bw()

## Compare epimutations
mod_methy_gset <- makeGenomicRatioSetFromMatrix(mod_methy, pData = colData(methy))

epi_ori <- epimutations( methy[com_cpgs,methy$status == "case"], methy[com_cpgs,methy$status == "control"],
                         method = "quantile")
epi_mod <- epimutations( mod_methy_gset[com_cpgs,mod_methy_gset$status == "case"], 
                         mod_methy_gset[com_cpgs,mod_methy_gset$status == "control"],
                         method = "quantile")

epi_integrate <- epimutations( methy[, methy$status == "case"], methy[, methy$status == "control"],
                               method = "quantile_reference", quantile_reference = ref_quantiles)
  epi_integrate <- subset(epi_integrate, cpg_n >= 3)

library(cowplot)
plot_grid(plot_epimutations(epi_ori, methy) + ggtitle("original"), 
  plot_epimutations(epi_ori, mod_methy_gset) + ggtitle("modified"), ncol = 2)

lapply(c(1:3, 5), function(i){
  print(i)
  plot_grid(plot_epimutations(epi_ori[i, ], methy) + ggtitle("original"), 
            plot_epimutations(epi_ori[i, ], mod_methy_gset) + ggtitle("modified"), ncol = 2)
  
})

lapply(c(1, 3, 4), function(i){
  print(i)
  plot_grid(plot_epimutations(epi_mod[i, ], methy) + ggtitle("original"), 
            plot_epimutations(epi_mod[i, ], mod_methy_gset) + ggtitle("modified"), ncol = 2)
  
})

lapply(c(1:2, 4), function(i){
  print(i)
  plot_grid(plot_epimutations(epi_integrate[i, ], methy) + ggtitle("original"), 
            plot_epimutations(epi_integrate[i, ], mod_methy_gset) + ggtitle("modified"), ncol = 2)
  
})


## Test in simulations ####
library(ramr)
load("simulated_GRS.Rdata")
ref_quantiles <- rowQuantiles(getBeta(sim.refGRS), probs =  c(0.005, 0.995, 0.50), na.rm = TRUE)
ref_betas <- epimutacions:::getBetaParams(t(getBeta(sim.refGRS)))

epi_pca <- epimutations(sim.refGRS[, 1:10], sim.newGRS[, 1:10], 
                            method = "quantile", pca_correction = TRUE)


epi_refbeta <- epimutations(sim.newGRS, sim.newGRS, 
                                 quantile_reference = list(quantiles = ref_quantiles, beta_params = ref_betas), method = "quantile_reference")

ramr.data <- rowRanges(sim.newGRS)
mcols(ramr.data) <- getBeta(sim.newGRS)

ref.data <- rowRanges(sim.refGRS)
mcols(ref.data) <- getBeta(sim.refGRS)

plotAMR(ref.data, sim.refGRS$sample, amrs.unique[1])

sel_cpgs <- names(subsetByOverlaps(rowRanges(sim.newGRS), amrs.unique[1]))

epi_refbeta_mini <- epimutations(sim.newGRS[sel_cpgs, ], sim.newGRS[sel_cpgs, ], 
                            quantile_reference = list(quantiles = ref_quantiles, beta_params = ref_betas), method = "quantile_reference")

