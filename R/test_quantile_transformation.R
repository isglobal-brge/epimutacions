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




normalizeQuantile <- function(vals, quantiles){
  
  ### Select values from quantiles 5-95%
  valsm <- minfi::logit2(vals)
  vals_f <- valsm[valsm > quantile(valsm, 0.05) & valsm < quantile(valsm, 0.95)]
  
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(vals_f, ncol = 1), quantiles)
  mod <- lm(valsQ ~ vals_f)
  valsm_out <- predict(mod, data.frame(vals_f = valsm))
  vals_out <- minfi::ilogit2(valsm_out)
  }

#### New distribution higher ####
vals_high <- simulateData(0.1, 0.2, 5e-4, 1000, 10, 0.1)

qA <- minfi::logit2(quantile(vals_high$A, seq(0.05, 0.95, 0.05)))
distB_norm <- normalizeQuantile(vals_high$B, qA)

data.frame(Distribution = c(vals_high$A, vals_high$B, distB_norm), 
           Group = rep(c("A", "B", "B normalize"), c(lengths(vals_high), length(distB_norm)))) %>%
  ggplot(aes(x = Group, y = Distribution)) +
  geom_boxplot() +
  theme_bw()


data.frame(Distribution = c(vals_high$A, vals_high$B, distB_norm), 
           Group = rep(c("A", "B", "B normalize"), c(lengths(vals_high), length(distB_norm)))) %>%
  ggplot(aes(color = Group, x = Distribution)) +
  geom_density() +
  theme_bw()


#### New distribution higher ####
vals_low <- simulateData(0.2, 0.1, 5e-4, 1000, 10, 0.1)

qA <- minfi::logit2(quantile(vals_low$A, seq(0.05, 0.95, 0.05)))
distB_norm <- normalizeQuantile(vals_low$B, qA)

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
