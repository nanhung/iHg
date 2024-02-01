#
rm(list = ls())

#
load("outputs/iHgMice_mcmc.Rdata")
load("outputs/iHgRat_mcmc.Rdata")
load("outputs/iHgHuman")
no_sample <- 20
sample_iters <- sample(seq_len(dim(mice_mcmc)[1]), no_sample)
sample_mice_mcmc <- mice_mcmc[sample_iters, , ]
nd2 <- dim(sample_mice_mcmc)[3]
dim(sample_mice_mcmc) <- c(4 * no_sample, nd2)
dim(sample_mice_mcmc)
