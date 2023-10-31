# Load packages -------------------------------------------------------------------------
library(data.table) # fread
library(purrr) # map
library(ggplot2)
library(coda)
library(rstan)

MC_dens<-function(x, param, dprior = NULL, dprior.arg, n = 1000, lty = 2, 
                  ref_posterior = NULL){
  par(mar=c(2,2,3,1))
  df <- as.data.frame(x[,,param])
  dens <- apply(df, 2, density)
  if (is.character(param)==TRUE){
    plot(NA, 
         xlim=range(sapply(dens, "[", "x")), 
         ylim=range(sapply(dens, "[", "y")), 
         xlab="", ylab="", main = param)
  } else {
    plot(NA, 
         xlim=range(sapply(dens, "[", "x")), 
         ylim=range(sapply(dens, "[", "y")), 
         xlab="", ylab="", main = dimnames(x)[[3]][param])
  }
  mapply(lines, dens, col=1:length(dens))
  
  if (is.null(dprior) == FALSE){
    mapply(lines, dens, col=1:length(dens))
    plot.range <- seq(par("usr")[1], par("usr")[2], length.out = n)
    prior.dens <- do.call(dprior, c(list(x=plot.range), dprior.arg))
    lines(plot.range, prior.dens, lty=lty)
  } 
  if (is.null(ref_posterior) == FALSE){
    x1<-ref_posterior[1]
    x2<-ref_posterior[2]
    y1<-par("usr")[3]
    y2<-par("usr")[4]
    polygon(x = c(x1, x1, x2, x2), 
            y = c(y1, y2, y2, y1), 
            col = rgb(1, 0, 0, 0.1), lty = 2, border = "grey")
  }
}

out <- c("outputs/iHgRat_3365.out",
         "outputs/iHgRat_6734.out",
         "outputs/iHgRat_4880.out",
         "outputs/iHgRat_5916.out")

out <- c("outputs/iHgMice_3365.out",
         "outputs/iHgMice_6734.out",
         "outputs/iHgMice_4880.out",
         "outputs/iHgMice_5916.out")

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
str <- 1 #floor(dim(x)[1]/2) # + 30001
end <- dim(x)[1]  
iters <- str:end 

mnt <- monitor(x[, , pars_name[-1]], digit = 4, print = F)
print(mnt)

#
select_pars <- c("M_lnPLC(1)", "M_lnPKC(1)", "M_lnPBrnC(1)", "M_lnPBrnC(1)",
                 "M_lnPRestC(1)", "M_lnKabsC(1)", "M_lnKunabsC(1)", "M_lnKbileC(1)",
                 "M_lnKurineC(1)", "M_lnKbrnC(1)")
 
par(mfrow = c(2,5))
for (i in 1:length(select_pars)){
  MC_dens(x[iters,,], param = select_pars[i])
}


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

