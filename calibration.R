#
library(data.table)
library(purrr)
library(rstan)
library(ggplot2)
library(dplyr)


#
names(pdfFonts())
my_font <- "sans"
my_theme <- theme(
  plot.title       = element_text(family = my_font, size = 14, 
                                  face = "bold", hjust = 0.5, 
                                  margin = ggplot2::margin(0, 0, 10, 0)),
  axis.title.y     = element_text(face = "bold", size = 10, family = my_font, 
                                  margin = ggplot2::margin(0, 10, 0, 0)),
  axis.text.y      = element_text(color = "black"), 
  axis.ticks.y     = element_line(color = "black"),
  axis.text.x      = element_text(color = "black"),
  axis.ticks.x     = element_line(color = "black"),
  axis.title.x     = element_text(face = "bold", size = 10, 
                                  family = my_font, 
                                  margin = ggplot2::margin(10, 0, 0, 0)),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border     = element_blank(),
  panel.background = element_blank(),
  plot.caption     = element_text(hjust = 0, family = my_font, size = 8)
  )

#
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

x[seq(50001, 100001, 10), , -1] |> monitor(warmup = 0)

mice_mcmc <- x[seq(50001, 100001, 10), , ]
save(mice_mcmc, file = "outputs/iHgMice_mcmc.RData")

#
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

x[seq(50001, 100001, 10), , -1] |> monitor(warmup = 0)

rat_mcmc <- x[seq(50001, 100001, 10), , ]
save(rat_mcmc, file = "outputs/iHgRat_mcmc.RData")

#
out <- c("outputs/iHgHuman_3365.out",
  "outputs/iHgHuman_4880.out",
  "outputs/iHgHuman_5916.out",
  "outputs/iHgHuman_6734.out")
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

x[, , -1] |> monitor(warmup = 0)

human_mcmc <- x
save(human_mcmc, file = "outputs/iHgHuman_mcmc.RData")

#
rm(list = ls())

#
load("outputs/iHgMice_mcmc.Rdata")
no_sample <- 10
sample_iters <- sample(seq_len(dim(mice_mcmc)[1]), no_sample)
sample_mice_mcmc <- mice_mcmc[sample_iters, , ]
nd2 <- dim(sample_mice_mcmc)[3]
dim(sample_mice_mcmc) <- c(4 * no_sample, nd2)
dim(sample_mice_mcmc)

model <- "iHgMiceBW.model"
if (!file.exists("mcsim.iHgMiceBW.model.exe")) {
   RMCSim::makemcsim(model, dir = "modeling")
}
for (iter in seq(dim(sample_mice_mcmc)[1])){
  head(sample_mice_mcmc, iter) |> tail(1) |>
    write.table(file = "MCMC.check.dat", row.names = FALSE, sep = "\t")
  input <- "iHgMice.MCMC.check.in"
  RMCSim::mcsim(model = model, input = input, dir = "modeling")
  out <- read.delim("MCMC.check.out")
  out$iter <- iter
  if (iter == 1) xx <- out
  else xx <- rbind(xx, out)
}
xx |> tail()

xx$Data[xx$Data == -1] <- NA

xx |> filter(Simulation == 1 & Time > 0) |> 
  ggplot() +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(Output_Var ~ Simulation, scales = "free") + my_theme

