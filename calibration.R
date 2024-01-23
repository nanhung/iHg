library(data.table)
library(purrr)
library(rstan)
out <- c("outputs/iHgMice_3365.out",
         "outputs/iHgMice_4880.out",
         "outputs/iHgMice_5916.out",
         "outputs/iHgMice_6734.out")
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
