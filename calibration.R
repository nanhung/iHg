library(data.table)
library(purrr)
library(rstan)
out <- c("outputs/iHgRat_3365.out",
         "outputs/iHgRat_4880.out",
         "outputs/iHgRat_5916.out",
         "outputs/iHgRat_6734.out")
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
dim(x)

x[seq(50001, 100001, 10), ,-1] |> monitor()

iHgRat_mcmc <- x[seq(50001, 100001, 10), ,]
save(iHgRat_mcmc, file = "iHgMice_mcmc.RData")

load("outputs/iHgMice_mcmc.Rdata")
no_sample <- 10
sample_iters <- sample(1:dim(iHgMice_mcmc)[1], no_sample)
sample_iHgMice_mcmc <- iHgMice_mcmc[sample_iters, ,]
nd2 <- dim(sample_iHgMice_mcmc)[3]
dim(sample_iHgMice_mcmc) <- c(4*no_sample, nd2)

for(iter in seq(dim(sample_iHgMice_mcmc)[1])){
  head(sample_iHgMice_mcmc, iter) |> head() |>
    write.table(file="MCMC.check.dat", row.names=F, sep="\t")  
  model <- "iHgMiceBW.model"
  input <- "iHgMice.MCMC.check.in"
  RMCSim::mcsim(model = model, input = input, dir = "modeling")
  out <- read.delim("MCMC.check.out")
  out$iter <- iter
  if (iter==1) X <- out
    else X <- rbind(X, out)
}
X |> tail()

  
MCMCfinal <- do.call(rbind, list(X[[3]], X[[2]], X[[1]]))
MCMCfinal |> tail(1) |> write.table(file="MCMC.check.dat", row.names=F, sep="\t")  
model <- "gPYR_analytic_ss.model"
#if(j %in% c(1:4)){
  input <- paste0("gPYR_analytic_5met_check_", cohort[j], "_Total.mcmc.in")
#} else {
#  input <- paste0("gPYR_analytic_3met_check_", cohort[j], "Total.mcmc.in")
#}  
RMCSim::mcsim(model = model, input = input, dir = "MCSim")

