# packages
library(data.table)
library(purrr)
library(rstan)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

# posterior check
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

# Save to RData
mice_mcmc <- x[seq(50001, 100001, 10), , ]
save(mice_mcmc, file = "outputs/iHgMice_mcmc.RData")

# tidy
rm(list = ls())

# data manipulate (randon sample 20 iters from 4 chains)
load("outputs/iHgMice_mcmc.Rdata")
no_sample <- 20
sample_iters <- sample(seq_len(dim(mice_mcmc)[1]), no_sample)
sample_mice_mcmc <- mice_mcmc[sample_iters, , ]
nd2 <- dim(sample_mice_mcmc)[3]
dim(sample_mice_mcmc) <- c(4 * no_sample, nd2)
dim(sample_mice_mcmc)

# posterior predictive simulation
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

# ouput manipulate
xx <- xx |> 
  mutate(conc = ifelse(Output_Var == "AUrine", "Urine",
                                 ifelse(Output_Var == "CKU", "Kidney",
                                        ifelse(Output_Var == "CBrnU", "Brain",
                                               ifelse(Output_Var == "CBldU", "Blood", "Liver")))))
xx <- xx |> 
  mutate(label = ifelse(Simulation == 1, "IV: 400 ug Hg/kg",
                                  ifelse(Simulation == 2, "Oral water: 1,000 ug Hg/kg",
                                         ifelse(Simulation == 3, "Oral gavage: 3,000 ug Hg/kg/d",
                                                ifelse(Simulation == 4, "Oral gavage: 15,000 ug Hg/kg/d",
                                                       ifelse(Simulation == 5, "Oral gavage: 75,000 ug Hg/kg/d",
                                                              ifelse(Simulation == 6, "Oral gavage: 925 ug Hg/kg/d",
                                                                     ifelse(Simulation == 7, "Oral gavage: 3,695 ug Hg/kg/d",
                                                                            "Oral gavage: 14,775 ug/kg/d"))))))))
xx$Data[xx$Data == -1] <- NA
adj_level <- xx$label |> unique()
xx$label <- factor(xx$label,
                   level = adj_level)
xx |> tail()
#
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
p1 <- xx |> filter(Simulation == 1 & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ label, scales = "free") +
  theme_bw() +
  set_theme
p2 <- xx |> filter(Simulation == 2 & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ label, scales = "free") +
  theme_bw() +
  set_theme
p3 <- xx |> filter(Simulation %in% c(3:5) & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ label, scales = "free") +
  theme_bw() +
  set_theme
p4 <- xx |> filter(Simulation %in% c(6:8) & Time > 0) |>
  ggplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_line(aes(x = Time, y = Prediction, group = iter), color = "grey") +
  geom_point(aes(x = Time, y = Data)) +
  facet_grid(conc ~ label, scales = "free") +
  theme_bw() +
  set_theme

# add the title and axis label
title <- ggdraw() +
  draw_label(
    "Mice",
    fontface = "bold",
    x = 0,
    size = 18,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 1)
  )
xlab <- ggdraw() +
  draw_label(
    "Time (hr)",
    fontface = "bold", size = 14, hjust = 0,
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )
ylab <- ggdraw() +
  draw_label(
    "Amount (ug) / Concentration (ug/mL)",
    fontface = "bold", size = 14, vjust = 0, angle = 90
  ) + theme(
    plot.margin = margin(0, 0, 0, 1)
  )

# plot
pdf(height = 11, width = 18)
plot_grid(
  ylab,
  plot_grid(
    title,
    plot_grid(
      plot_grid(p1, p2, nrow = 2, labels = c("A", "B"),
                rel_heights = c(3 / 4, 1 / 4)),
      plot_grid(
        p3, p4, nrow = 2,
        labels = c("C", "D")
      ),
      nrow = 1, rel_widths = c(0.33, 0.66)
    ),
    xlab, nrow = 3, rel_heights = c(0.05, 1, 0.05)),
  nrow = 1, rel_widths = c(0.02, 1)
)
dev.off()