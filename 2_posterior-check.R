# Load packages -------------------------------------------------------------------------
library(data.table) # fread
library(purrr) # map
library(ggplot2)
library(coda)

out <- c("outputs/iHgRat_6734.out",
         "outputs/iHgRat_4880.out",
         "outputs/iHgRat_5916.out")

data <- out |> map(fread) |> map(as.data.frame)

n_chains <- length(data)
sample_number <- dim(data[[1]])[1]
dim <- c(sample_number, n_chains, dim(data[[1]])[2])
n_iter <- dim(data[[1]])[1]
n_param <- dim(data[[1]])[2]
x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
for (i in 1:n_chains) {
  x[, i, ] <- as.matrix(data[[i]][1:n_iter, ])
}
dimnames(x)[[3]] <- names(data[[1]])
pars_name <- dimnames(x)[[3]]
dim(x)

# Posterior check
j <- 10002:20001
mnt <- monitor(x[j, , pars_name[-1]], digit = 6, print = T)

# Plot
x <- read.delim("outputs/iHgRat_check_5916.out")
pdf(file = "outputs/validation.pdf", width = 9, height = 9)
ggplot(data = x) + 
  geom_point(aes(x=Data, y=Prediction)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(slope = 1)
dev.off()

# Save final population parameters 
save_directory <- "outputs"
file_name <- paste("MCMCfinal_", format(Sys.time(), "%Y-%m-%d"), ".RData", sep = "")
str <- which(pars_name == "lnPLC(1)")
end <- which(pars_name == "V_lnKbrnC(1)")
pars <- pars_name[str:end]
mcmc_out <- x[j, , pars]
l <- dim(mcmc_out)[1] * dim(mcmc_out)[2]
dim(mcmc_out) <- c(l, length(pars))
colnames(mcmc_out) <- pars
save(mcmc_out, file = file.path(save_directory, file_name))

