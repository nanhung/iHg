load("validation_rat.RData")
load("validation_mice.RData")
load("validation_human.RData")

x |> ggplot() + 
  geom_point(aes(x=ratio, y=label)) +
  scale_x_log10() +
  facet_grid(species~organs, scales = "free_y") +
  geom_vline(xintercept = 1, lty=2) +
  labs(x="Prediction / Observation", y="Study") +
  theme(strip.background = element_blank())
        