# packag
library(data.table)
library(purrr)
library(rstan)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

# posterior check
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
x[seq(1, 10001, 2), , -1] |> monitor(warmup = 0)

# Save to RData
human_mcmc <- x[seq(1, 10001, 2), , ]
save(human_mcmc, file = "outputs/iHgHuman_mcmc.RData")

# tidy
rm(list = ls())

# data manipulate (randon sample 20 iters from 4 chains)
load("outputs/iHgHuman_mcmc.Rdata")
no_sample <- 20
sample_iters <- sample(seq_len(dim(human_mcmc)[1]), no_sample)
sample_human_mcmc <- human_mcmc[sample_iters, , ]
nd2 <- dim(sample_human_mcmc)[3]
dim(sample_human_mcmc) <- c(4 * no_sample, nd2)
dim(sample_human_mcmc)

# posterior predictive simulation
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

# ouput manipulate
xx <- xx |>
  mutate(conc = ifelse(Output_Var == "Aurine", "Urine",
                       ifelse(Output_Var == "CKU", "Kidney",
                              ifelse(Output_Var == "CBrnU", "Brain",
                                     ifelse(Output_Var == "CBldU", "Blood",
                                            ifelse(Output_Var == "CLU", "Liver",
                                                   "Feces"))))))
xx <- xx |>
  mutate(label = ifelse(Simulation == 1, "IV: 0.025 ug Hg/kg",
                        ifelse(Simulation == 2, "Oral: 0.09375 ug Hg/kg", "")))
xx$Data[xx$Data == -1] <- NA
adj_level <- xx$label |> unique()
xx$label <- factor(xx$label, level = adj_level)
xx |> tail()

# define plotting element
set_theme <- theme(
  axis.text.y      = element_text(color = "black"),
  axis.ticks.y     = element_line(color = "black"),
  axis.text.x      = element_text(color = "black"),
  axis.ticks.x     = element_line(color = "black"),
  axis.line.x      = element_line(color = "black"),
  axis.line.y      = element_line(color = "black"),
  legend.key       = element_blank(),
  axis.title       = element_blank(),
  panel.background = element_blank()
)
p1 <- xx |>
  filter(Simulation == 1 & Time > 0) |>
  ggplot() +
  scale_y_log10(lim = c(10^-2, 10^0),
                breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ label, scales = "free") +
  theme_bw() +
  set_theme
p2 <- xx |>
  filter(Simulation == 2 & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ label, scales = "free") +
  theme_bw() +
  set_theme

# add the title and axis label
title <- ggdraw() +
  draw_label(
    "Human",
    fontface = "bold",
    x = 0,
    size = 18,
    hjust = 0
  ) +
  theme(plot.margin = margin(0, 0, 0, 1))
xlab <- ggdraw() +
  draw_label(
    "Time (hr)",
    fontface = "bold", size = 14, hjust = 0,
  ) + theme(plot.margin = margin(0, 0, 0, 1))
ylab <- ggdraw() +
  draw_label(
    "Amount (ug) / Concenthumanion (ug/mL)",
    fontface = "bold", size = 14, vjust = 0, angle = 90
  ) + theme(plot.margin = margin(0, 0, 0, 1))

# plot
pdf(file = "calibration_human.pdf", height = 6, width = 18)
plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(p1, p2, nrow = 1, labels = c("A", "B"), rel_widths = c(0.5, 0.5)),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)
  ),
  nrow = 1, rel_widths = c(0.02, 1)
)
dev.off()
