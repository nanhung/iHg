#
library(pksensi)

# 
load("outputs/iHgHuman_mcmc.RData")
out_vars <- dimnames(human_mcmc)[[3]] 
nd1 <- dim(human_mcmc)[1] * 4
nd2 <- dim(human_mcmc)[3]
dim(human_mcmc) <- c(nd1,nd2)
colnames(human_mcmc)
str <- which(colnames(human_mcmc) == "M_lnPLC(1)")
end <- which(colnames(human_mcmc) == "M_lnKbrnC(1)")
pop_min <- human_mcmc[,str:end] |> apply(2, min)
pop_max <- human_mcmc[,str:end] |> apply(2, max)

#
mName <- "iHgHumanBW"
file.copy("modeling/iHgHumanBW.model", "iHgHumanBW.model") # issue to compile file in subfolder
compile_model(mName, application = "mcsim")

#
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
x <- rfast99(params = params, n = 512, q = q, q.arg = q.arg, replicate = 10)







conditions <- c("mgkg_flag = 1",
                "OralExp_APAP = NDoses(2, 1, 0, 0, 0.001)",
                "OralDose_APAP_mgkg = 20.0")
vars <- c("lnCPL_APAP_mcgL", "lnCPL_AG_mcgL", "lnCPL_AS_mcgL")
times <- c(0.1, 0.5, 1, 2, 3, 4, 6, 8, 12)
out <- solve_mcsim(x, mName = mName,
                   params = params, 
                   time = times, 
                   vars = vars,
                   condition = conditions, 
                   rtol = 1e-7, atol = 1e-9)
heat_check(out, order = "total order", show.all = T, times = c(8))
