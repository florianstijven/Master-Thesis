set.seed(1)
library(tidyverse)
library(latex2exp)
load(file = "results_ovarian_100k.RData")

#make artificial data set to plot the four plots in a square grid
sens_data_marg_helper = sens_data_marg %>%
  group_by(ordering) %>%
  mutate(assumptions = "No Assumptions") %>%
  ungroup()
sens_data_marg_helper = sens_data_marg_helper %>%
  bind_rows(
    sens_data_marg %>% filter(tau_s0s1 > 0, tau_s0t0 > 0, tau_s0t1 > 0,
                              tau_s1t0 > 0, tau_s1t1 > 0, tau_t0t1 > 0) %>%
      group_by(ordering) %>%
      mutate(assumptions = "Monotonicity") %>%
      ungroup()
  )
sens_data_marg_helper = sens_data_marg_helper %>%
  bind_rows(
    sens_data_marg %>% filter(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0))) %>%
      group_by(ordering) %>%
      mutate(assumptions = "Weaker CA") %>%
      ungroup()
  )
sens_data_marg_helper = sens_data_marg_helper %>%
  bind_rows(
    sens_data_marg %>% filter(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0)),
                              tau_s0s1 > 0, tau_s0t0 > 0, tau_s0t1 > 0,
                              tau_s1t0 > 0, tau_s1t1 > 0, tau_t0t1 > 0) %>%
      group_by(ordering) %>%
      ungroup() %>%
      mutate(assumptions = "Monotonicity + Weaker CA")
  )

sens_data_marg_helper = left_join(x = sens_data_marg_helper,
                                  y = sens_data_marg_helper %>%
                                    group_by(ordering, assumptions) %>%
                                    summarise(n = n()) %>%
                                    ungroup(),
                                  by = c("ordering", "assumptions")
)
sens_data_marg_helper = sens_data_marg_helper %>% mutate(icc = 1 - exp(-2*minfo))

sens_data_marg_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions, ordering, n) %>%
  summarise(Range = paste0("[", round(min(sp_rho), 3), ", ",
                           round(max(sp_rho), 3), "]"),
            perc = paste0("[", round(quantile(sp_rho, 0.01), 3), ", ",
                          round(quantile(sp_rho, 0.99), 3), "]"),
            median = median(sp_rho)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(ordering, assumptions) %>%
  knitr::kable(format = "latex")

sens_data_marg %>%
  group_by(ordering) %>%
  summarise(mean_0 = mean(tau_s0t0), mean_1 = mean(tau_s1t1))

png(filename = "Figures/ovarian_results_exploration_ca_100k_rho.png", width = 6.6, height = 4.95,
    units = "cm", pointsize = 10, res = 500)
sens_data_marg %>%
  mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0)),
                      "Yes", "No")) %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = ca, y = sp_rho)) +
  scale_x_continuous(limits = c(-1, 1), 
                     name = TeX("$\\frac{\\tau_{S_0T_1} + \\tau_{S_1T_0}}{2}$")) +
  scale_y_continuous(name = TeX("$\\rho_s")) +
  geom_vline(xintercept = 0.816) +
  geom_point(alpha = 0.3, size = 0.03) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2), legend.background = element_blank())
dev.off()

# png(filename = "Figures/ovarian_results_exploration_wa.png", width = 13.2, height = 9.9, 
#     units = "cm", pointsize = 10, res = 500)
sens_data_marg %>%
  mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
  mutate(Monotonicity = ifelse(pmin(tau_s0s1, tau_t0t1) > 0,
                               "Yes", "No")) %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = wa, y = sp_rho, color = Monotonicity)) +
  scale_x_continuous(limits = c(-1, 1), 
                     name = TeX("$\\frac{\\tau_{S_0S_1} + \\tau_{T_0T_1}}{2}$")) +
  scale_y_continuous(name = TeX("$\\rho_s")) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2), legend.background = element_blank())
# dev.off()
