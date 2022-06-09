set.seed(1)
library(tidyverse)
library(latex2exp)
load(file = "ovarian_sens_results.RData")

model_comparison %>% 
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 2))

#make artificial data set to plot the four plots in a square grid
sens_data_marg_helper = sens_data_marg %>%
  group_by(ordering) %>%
  slice_sample(n = 1000) %>%
  mutate(assumptions = "No Assumptions") %>%
  ungroup()
sens_data_marg_helper = sens_data_marg_helper %>%
  bind_rows(
    sens_data_marg %>% filter(tau_s0s1 > 0, tau_s0t0 > 0, tau_s0t1 > 0,
                              tau_s1t0 > 0, tau_s1t1 > 0, tau_t0t1 > 0) %>%
      group_by(ordering) %>%
      slice_sample(n = 1000) %>%
      mutate(assumptions = "Monotonicity") %>%
      ungroup()
  )
sens_data_marg_helper = sens_data_marg_helper %>%
  bind_rows(
    sens_data_marg %>% filter(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0))) %>%
      group_by(ordering) %>%
      slice_sample(n = 1000) %>%
      mutate(assumptions = "Weaker CA") %>%
      ungroup()
  )
sens_data_marg_helper = sens_data_marg_helper %>%
  bind_rows(
    sens_data_marg %>% filter(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0)),
                              tau_s0s1 > 0, tau_s0t0 > 0, tau_s0t1 > 0,
                              tau_s1t0 > 0, tau_s1t1 > 0, tau_t0t1 > 0) %>%
      group_by(ordering) %>%
      slice_sample(n = 1000) %>%
      ungroup() %>%
      mutate(assumptions = "Monotonicity + Weaker CA")
  )

#R_H
png(filename = "Figures/ovarian_results_sens_ord.png", width = 12, height = 9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = 1 - exp(-2*minfo))) +
  scale_x_continuous(limits = c(0.5, 1), name = TeX("$R_H^2$")) +
  geom_histogram(fill = "gray", color = "black") +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()

png(filename = "Figures/ovarian_results_sens_no_ord.png", width = 12, height = 9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "No Ordering") %>%
  ggplot(aes(x = 1 - exp(-2*minfo))) +
  scale_x_continuous(limits = c(0.5, 1), name = TeX("$R_H^2$")) +
  geom_histogram(fill = "gray", color = "black") +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()

#Spearman's rho
png(filename = "Figures/ovarian_results_sens_ord_rho.png", width = 12, height = 9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = sp_rho)) +
  scale_x_continuous(limits = c(0.5, 1), name = TeX("$\rho_S$")) +
  geom_histogram(fill = "gray", color = "black") +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()

png(filename = "Figures/ovarian_results_sens_no_ord_rho.png", width = 12, height = 9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "No Ordering") %>%
  ggplot(aes(x = sp_rho)) +
  scale_x_continuous(limits = c(0.5, 1), name = TeX("$\rho_s$")) +
  geom_histogram(fill = "gray", color = "black") +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()


sens_data_marg_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions, ordering) %>%
  summarise(Range = paste0("[", round(min(1 - exp(-2*minfo)), 3), ", ",
                           round(max(1 - exp(-2*minfo)), 3), "]"),
            SD = sd(1 - exp(-2*minfo)),
            perc = paste0("[", round(quantile(1 - exp(-2*minfo), 0.01), 3), ", ",
                          round(quantile(1 - exp(-2*minfo), 0.99), 3), "]"),
            median = median(1 - exp(-2*minfo))
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(ordering, assumptions) %>%
  knitr::kable(format = "latex")

sens_data_marg_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions, ordering) %>%
  summarise(Range = paste0("[", round(min(sp_rho), 3), ", ",
                           round(max(sp_rho), 3), "]"),
            SD = sd(sp_rho),
            perc = paste0("[", round(quantile(sp_rho, 0.01), 3), ", ",
                          round(quantile(sp_rho, 0.99), 3), "]"),
            median = median(sp_rho)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(ordering, assumptions) %>%
  knitr::kable(format = "latex")

sens_data_marg_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions, ordering) %>%
  summarise(Range = paste0("[", round(min(kendall), 3), ", ",
                           round(max(kendall), 3), "]"),
            SD = sd(kendall),
            perc = paste0("[", round(quantile(kendall, 0.01), 3), ", ",
                          round(quantile(kendall, 0.99), 3), "]"),
            median = median(kendall)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(ordering, assumptions) %>%
  knitr::kable(format = "latex")

# sens_data_marg_helper %>% 
#   pivot_longer(names_to = "extra", values_to = "value", cols = 4:9) %>%
#   filter(!(extra %in% c("tau_s0t0", "tau_s1t1"))) %>%
#   filter(assumptions == "Monotonicity + Weaker CA") %>%
#   ggplot(aes(x = value, y = sp_rho)) +
#   geom_point() +
#   facet_grid(ordering~extra)

# sens_data_marg_helper %>% 
#   pivot_longer(names_to = "extra", values_to = "value", cols = 4:9) %>%
#   filter(!(extra %in% c("tau_s0t0", "tau_s1t1"))) %>%
#   ggplot(aes(x = value, y = sp_rho, color = assumptions)) +
#   geom_point(alpha = 0.1, size = 0.5) +
#   facet_grid(ordering~extra)

#this plot may be less useful
#I think these plots could be used to illustrate another way to look
#at additional assumptions. Instead of prespecifying all the assumptions
#we could also just look at the relation between ICA on the one hand,
#and the unidentifiable parameters on the other hand

# sens_data_marg_helper %>% 
#   mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
#   ggplot(aes(x = ca, y = 1 - exp(-2*minfo))) +
#   scale_x_continuous(limits = c(-1, 1), 
#                      name = TeX("$\\frac{\\tau_{S_0T_1} + \\tau_{S_1T_0}}{2}$")) +
#   scale_y_continuous(name = TeX("$\\R_H^2")) +
#   geom_point(alpha = 0.3, size = 1) +
#   theme_bw() +
#   facet_wrap(facets = "assumptions")
# 
# sens_data_marg_helper %>% 
#   mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
#   ggplot(aes(x = wa, y = 1 - exp(-2*minfo))) +
#   scale_x_continuous(limits = c(-1, 1), 
#                      name = TeX("$\\frac{\\tau_{S_0T_0} + \\tau_{S_1T_1}}{2}$")) +
#   scale_y_continuous(name = TeX("$\\R_H^2")) +
#   geom_point(alpha = 0.3, size = 1) +
#   theme_bw() +
#   facet_wrap(facets = "assumptions")

sens_data_marg %>%
  group_by(ordering) %>%
  summarise(mean_0 = mean(tau_s0t0), mean_1 = mean(tau_s1t1))

png(filename = "Figures/ovarian_results_exploration_ca.png", width = 12, height = 9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg %>%
  mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0)),
                      "Yes", "No")) %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = ca, y = 1 - exp(-2*minfo), color = WCA)) +
  scale_x_continuous(limits = c(-1, 1), 
                     name = TeX("$\\frac{\\tau_{S_0T_1} + \\tau_{S_1T_0}}{2}$")) +
  scale_y_continuous(name = TeX("$\\R_H^2")) +
  geom_vline(xintercept = 0.75) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2), legend.background = element_blank())
dev.off()

png(filename = "Figures/ovarian_results_exploration_wa.png", width = 12, height = 9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg %>%
  mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
  mutate(Monotonicity = ifelse(pmin(tau_s0s1, tau_t0t1) > 0,
                      "Yes", "No")) %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = wa, y = 1 - exp(-2*minfo), color = Monotonicity)) +
  scale_x_continuous(limits = c(-1, 1), 
                     name = TeX("$\\frac{\\tau_{S_0S_1} + \\tau_{T_0T_1}}{2}$")) +
  scale_y_continuous(name = TeX("$\\R_H^2")) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2), legend.background = element_blank())
dev.off()

# 
# sens_data_marg %>%
#   mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
#   filter(ordering == "No Ordering") %>%
#   ggplot(aes(x = ca, y = sp_rho)) +
#   scale_x_continuous(limits = c(-1, 1), 
#                      name = TeX("$\\frac{\\tau_{S_0T_1} + \\tau_{S_1T_0}}{2}$")) +
#   scale_y_continuous(name = TeX("$\\rho_s")) +
#   geom_vline(xintercept = 0.80) +
#   geom_point(alpha = 0.3) +
#   theme_bw()

