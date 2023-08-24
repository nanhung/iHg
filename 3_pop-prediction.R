# Load packages -------------------------------------------------------------------------
library(RMCSim)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)

# Data
Tissues <- c(rep("Kidney", 3), rep("Liver", 3),
            rep("Kidney", 3), rep("Liver", 3),
            rep("Kidney", 3), rep("Liver", 3), rep("Brain", 3))
Time <- rep(c(1428, 2868,	4308), 7)
Conc <- c(22.95, 33.67, 47.8, 0.23, 0.14, 0.14,
          100.87, 97.19, 86.08, 0.25, 0.38, 0.44,
          101.87, 122.19, 92.78, 1.00, 1.62, 1.86, 0.06, 0.03, 0.01)
Experiment <- c(rep("Exp 1 - 312 ug/kg", 6),
               rep("Exp 2 - 1250 ug/kg", 6),
               rep("Exp 3 - 5000 ug/kg", 9))

dat <- data.frame(Tissues, Time, Conc, Experiment)

# Load MCMCfinal
MCMC_file_to_be_load <- list.files(paste0("./outputs/"), pattern='MCMCfinal')
load(paste0("outputs/",MCMC_file_to_be_load))
parms <- colnames(mcmc_out)

# Setpoint
i <- sample(dim(mcmc_out)[1], 100) # random select 100 draws
tmp.x <- mcmc_out[i, parms] 
tmp.x |> write.table(file="poppred.dat", row.names=T, sep="\t")
model <- "iHgRatBW.model"
input <- "iHgRat_MCMC_setpts.in" 
mcsim(model = model, input = input, dir = "MCSim")

# Pop prediction
pred_pop <- fread("poppred.out") |> as.data.frame()
dim(pred_pop)
str_1.1 <- which(names(pred_pop) == "CKU_1.1")

Conc.med <- apply(pred_pop[str_1.1:dim(pred_pop)[2]], 2, quantile, probs=0.5) |> as.numeric()
Conc.upper <- apply(pred_pop[str_1.1:dim(pred_pop)[2]], 2, quantile, probs=0.975) |> as.numeric()
Conc.lower <- apply(pred_pop[str_1.1:dim(pred_pop)[2]], 2, quantile, probs=0.025) |> as.numeric()
Time <- rep(seq(0, 4320, 1), 9)
Tissues <- rep(rep(c("Kidney", "Liver", "Brain"), each=4321), 3)
Experiment <- c(rep(c("Exp 1 - 312 ug/kg", 
                      "Exp 2 - 1250 ug/kg", 
                      "Exp 3 - 5000 ug/kg"), each=12963))
pred <- data.frame(Experiment, Tissues, Time, Conc.med, Conc.lower, Conc.upper)

# plot
#pdf(file = "outputs/pop-prediction.pdf", width = 9, height = 9)
png(file = "outputs/pop-prediction.png", width = 3000, height = 1800, res = 300)
pred |> 
  ggplot() + geom_line(aes(x=Time, y=Conc.med)) + 
  geom_line(aes(x=Time, y= Conc.upper), lty="dotted", lwd=0.5) +
  geom_line(aes(x=Time, y= Conc.lower), lty="dotted", lwd=0.5) +
  geom_point(data = dat, aes(x=Time, y= Conc)) +
  facet_grid(Experiment~Tissues) +
  scale_y_log10() +
  theme_bw() + 
  labs(x="Time, hr", y="Concentration") +
  theme(
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    axis.text = element_text(size = 12, colour = "black", face = "bold"))
dev.off()  

# Housekeeping
system("rm *.out *.dat *.exe")
