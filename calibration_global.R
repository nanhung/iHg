#
library(scales)
library(dplyr)
library(ggplot2)
library(viridis)

#
load("outputs/iHgMice_mcmc.Rdata")
no_sample <- 20
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
mice_xx <- xx

#
load("outputs/iHgRat_mcmc.Rdata")
no_sample <- 20
sample_iters <- sample(seq_len(dim(rat_mcmc)[1]), no_sample)
sample_rat_mcmc <- rat_mcmc[sample_iters, , ]
nd2 <- dim(sample_rat_mcmc)[3]
dim(sample_rat_mcmc) <- c(4 * no_sample, nd2)
dim(sample_rat_mcmc)
model <- "iHgRatBW.model"
if (!file.exists("mcsim.iHgRatBW.model.exe")) {
  RMCSim::makemcsim(model, dir = "modeling")
}
for (iter in seq(dim(sample_rat_mcmc)[1])){
  head(sample_rat_mcmc, iter) |> tail(1) |>
    write.table(file = "MCMC.check.dat", row.names = FALSE, sep = "\t")
  input <- "iHgRat.MCMC.check.in"
  RMCSim::mcsim(model = model, input = input, dir = "modeling")
  out <- read.delim("MCMC.check.out")
  out$iter <- iter
  if (iter == 1) xx <- out
  else xx <- rbind(xx, out)
}
xx$Output_Var |> unique()
rat_xx <- xx

#
load("outputs/iHgHuman_mcmc.RData")
no_sample <- 20
sample_iters <- sample(seq_len(dim(human_mcmc)[1]), no_sample)
sample_human_mcmc <- human_mcmc[sample_iters, , ]
nd2 <- dim(sample_human_mcmc)[3]
dim(sample_human_mcmc) <- c(4 * no_sample, nd2)
dim(sample_human_mcmc)
model <- "iHgHumanBW.model"
if (!file.exists("mcsim.iHgHumanBW.model.exe")) {
  RMCSim::makemcsim(model, dir = "modeling")
}
for (iter in seq(dim(sample_human_mcmc)[1])){
  head(sample_human_mcmc, iter) |>
    tail(1) |>
    write.table(file = "MCMC.check.dat", row.names = FALSE, sep = "\t")
  input <- "iHgHuman.MCMC.check.in"
  RMCSim::mcsim(model = model, input = input, dir = "modeling")
  out <- read.delim("MCMC.check.out")
  out$iter <- iter
  if (iter == 1) xx <- out
  else xx <- rbind(xx, out)
}
xx$Output_Var |> unique()
human_xx <- xx

#
mice_xx$species <- "Mice"
rat_xx$species <- "Rat"
human_xx$species <- "Human"

x <- do.call(rbind, list(mice_xx, rat_xx, human_xx))
x <- x |>
  mutate(conc = ifelse(Output_Var == "Aurine", "Urine",
                                 ifelse(Output_Var == "CKU", "Kidney",
                                        ifelse(Output_Var == "CBrnU", "Brain",
                                               ifelse(Output_Var == "CBldU", "Blood", ifelse(Output_Var == "CLU", "Liver", "Feces"))))))
x$Data[x$Data == -1] <- NA
x <- x[complete.cases(x), ]
x <- x |> mutate(ratio = Data/Prediction) |>
  mutate(acceptance = ifelse(ratio > 3.3, "Outside 3X difference", 
                             ifelse(ratio < 0.333, "Outside 3X difference", 
                                    "Between 3X difference")))
rng <- range(x$Data)

#
# define plotting element
set_theme <- theme(
  legend.position  = "top",
  axis.text.y      = element_text(color = "black"),
  axis.ticks.y     = element_line(color = "black"),
  axis.text.x      = element_text(color = "black"),
  axis.ticks.x     = element_line(color = "black"),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  axis.title       = element_blank(),
  legend.title     = element_blank(),  
  panel.background = element_blank()
)

pdf(file = "calibration_global.pdf", height = 5, width = 13)
x |> ggplot()+ 
  facet_wrap(~species) + 
  scale_x_log10(
    #lim = rng,
    breaks = trans_breaks("log10", function(x) 10^x, n = 5),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(
    #lim = rng,
    breaks = trans_breaks("log10", function(x) 10^x, n = 5),
    labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = -0.52, lty = 2) +
  geom_abline(slope = 1, intercept = 0.52, lty = 2) +
  geom_point(aes(x = Data, y = Prediction, shape = conc, color = acceptance), 
             alpha =0.6) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, end = 0.8) + 
  set_theme
dev.off()
