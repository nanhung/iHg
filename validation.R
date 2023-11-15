load("validation_rat.RData")
load("validation_mice.RData")
load("validation_human.RData")

rate_text <- do.call(rbind, list(validation_rat, validation_mice, validation_human)) |>
  select(species, organs, accept_rate) |> distinct(species, organs, accept_rate) 

do.call(rbind, list(validation_rat, validation_mice, validation_human)) |> 
  ggplot() + 
  geom_vline(xintercept = 0.33, lty=3, col=2) +
  geom_vline(xintercept = 3, lty=3, col=2) +
  geom_vline(xintercept = 1, lty=2, col=2) +
  geom_point(aes(x=ratio, y=label), col="grey", alpha=0.4) +
  scale_x_log10() +
  facet_grid(species~organs, scales = "free_y") +
  geom_label(data = rate_text, mapping = aes(x = 0.01, y = 0.9, label = round(accept_rate, 2))) +
  labs(x="Prediction / Observation", y="Study") +
  theme(strip.background = element_blank())
        