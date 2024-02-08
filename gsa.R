# Tidy
rm(list=ls())

# Packages
pacman::p_load(pksensi, ggplot2)

# Parameter range
load("outputs/iHgHuman_mcmc.RData")
out_vars <- dimnames(human_mcmc)[[3]] 
nd1 <- dim(human_mcmc)[1] * 4
nd2 <- dim(human_mcmc)[3]
dim(human_mcmc) <- c(nd1,nd2)
colnames(human_mcmc) <- out_vars
str <- which(colnames(human_mcmc) == "M_lnPLC(1)")
end <- which(colnames(human_mcmc) == "M_lnKbrnC(1)")
pop_min <- human_mcmc[,str:end] |> apply(2, min)
pop_max <- human_mcmc[,str:end] |> apply(2, max)
pop_max - pop_min

# Model Compilation
mName <- "iHgHumanBW"
file.copy("modeling/iHgHumanBW.model", "iHgHumanBW.model") # issue to compile file in subfolder
compile_model(mName, application = "mcsim")

# Parameter space generating
params <- c("lnPLC", "lnPKC", "lnPBrnC", "lnPRestC", 
  "lnKabsC", "lnKunabsC", "lnKbileC", "lnKurineC", "lnKbrnC")
q <- rep("qunif", 9)
q.arg <- list(list(pop_min[1], pop_max[1]), 
  list(pop_min[2], pop_max[2]), 
  list(pop_min[3], pop_max[3]), 
  list(pop_min[4], pop_max[4]), 
  list(pop_min[5], pop_max[5]), 
  list(pop_min[6], pop_max[6]), 
  list(pop_min[7], pop_max[7]), 
  list(pop_min[8], pop_max[8]), 
  list(pop_min[9], pop_max[9]) 
)
set.seed(1234)
x <- rfast99(params = params, n = 10000, q = q, q.arg = q.arg, replicate = 1)

# Single dose simulation
conditions <- c("BW0 = 64;", "BWgrowth = 0;", "Growthrate = 0;", 
  "sex = 1;", "TChng = 0.5;", "PDose = PerDose(0.09375, 24,  0, 0.05);",
  "IVDose = PerDose(0.0 , 24,  0, 0.003);",
  " expowk =  PerDose(1.0, 168,  0, 0.05);", 
  "expodur = PerDose(1.0, 1850, 0, 0.05);", "Drink = 0.0;")
vars <- c("AUCCL", "AUCCK", "AUCCBrn")
times <- c(12)
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
check(out)

# Plot
set_theme <- theme(
  legend.position  = "none",
  axis.text.x      = element_blank(),
  axis.ticks.x     = element_blank(),
  axis.title       = element_blank(),
)
pdf()
heat_check(out, order = c("first order", "total order"), show.all = T, text = T) + 
  ggtitle("Human") + theme_bw() + set_theme
dev.off()
