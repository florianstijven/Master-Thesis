set.seed(1)
#load packages
library(tidyverse)
library(latex2exp)
#load data obtained from sensitivity analysis
load(file = "Code/Ovarian Cancer Case Study/ovarian_sens_results.RData")

#model fitting results
model_comparison %>% 
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 2))

#recode ICC is a variable
sens_data_marg = sens_data_marg %>% mutate(icc = 1 - exp(-2*minfo))

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

#R_H
png(filename = "Figures/ovarian_results_sens_ord.png", width = 13.2, height = 9.9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "Ordering") %>%
  mutate(assumptions = paste0(assumptions, " (n = ", n, ")")) %>%
  ggplot(aes(x = icc)) +
  coord_cartesian(xlim = c(0.5, 1)) +
  scale_x_continuous(name = TeX("$R_h^2$")) +
  geom_histogram(fill = "gray", color = "black", 
                 binwidth = 0.025, boundary = 1) +
  scale_y_continuous(limits = c(0,4000)) +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()



png(filename = "Figures/ovarian_results_sens_no_ord.png", width = 13.2, height = 9.9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "No Ordering") %>%
  mutate(assumptions = paste0(assumptions, " (n = ", n, ")")) %>%
  ggplot(aes(x = 1 - exp(-2*minfo))) +
  coord_cartesian(xlim = c(0.5, 1)) +
  scale_x_continuous(name = TeX("$R_h^2$")) +
  geom_histogram(fill = "gray", color = "black", 
                 binwidth = 0.025, boundary = 1) +
  scale_y_continuous(limits = c(0,4000)) +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()

#Spearman's rho
png(filename = "Figures/ovarian_results_sens_ord_rho.png", width = 13.2, height = 9.9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "Ordering") %>%
  mutate(assumptions = paste0(assumptions, " (n = ", n, ")")) %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  geom_histogram(fill = "gray", color = "black", 
                 binwidth = 0.05, boundary = 1) +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()

png(filename = "Figures/ovarian_results_sens_no_ord_rho.png", width = 13.2, height = 9.9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg_helper %>%
  filter(ordering == "No Ordering") %>%
  ggplot(aes(x = sp_rho)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(name = TeX("$\\rho_s$")) +
  geom_histogram(fill = "gray", color = "black", 
                 binwidth = 0.05, boundary = 1) +
  theme_bw() +
  facet_wrap(facets = ~assumptions)
dev.off()

#TABLES
sens_data_marg_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions, ordering, n) %>%
  summarise(Range = paste0("[", round(min(icc), 3), ", ",
                           round(max(icc), 3), "]"),
            perc = paste0("[", round(quantile(icc, 0.01), 3), ", ",
                          round(quantile(icc, 0.99), 3), "]"),
            median = median(icc)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(ordering, assumptions) %>%
  knitr::kable(format = "latex")

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

sens_data_marg_helper %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  group_by(assumptions, ordering, n) %>%
  summarise(Range = paste0("[", round(min(kendall), 3), ", ",
                           round(max(kendall), 3), "]"),
            perc = paste0("[", round(quantile(kendall, 0.01), 3), ", ",
                          round(quantile(kendall, 0.99), 3), "]"),
            median = median(kendall)
  ) %>%
  mutate_if(.predicate = is.numeric, .funs = ~round(x = .x, digits = 3)) %>%
  arrange(ordering, assumptions) %>%
  knitr::kable(format = "latex")


#extra explorative graphs
sens_data_marg %>%
  group_by(ordering) %>%
  summarise(mean_0 = mean(tau_s0t0), mean_1 = mean(tau_s1t1))

png(filename = "Figures/ovarian_results_exploration_ca.png", width = 13.2, height = 9.9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg %>%
  mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
  mutate(WCA = ifelse(pmin(abs(tau_s0t0), abs(tau_s1t1)) > pmax(abs(tau_s0t1), abs(tau_s1t0)),
                      "Yes", "No")) %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = ca, y = 1 - exp(-2*minfo), color = WCA)) +
  scale_x_continuous(limits = c(-1, 1), 
                     name = TeX("$\\frac{\\tau_{S_0T_1} + \\tau_{S_1T_0}}{2}$")) +
  scale_y_continuous(name = TeX("$\\R_h^2")) +
  geom_vline(xintercept = 0.816) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2), legend.background = element_blank())
dev.off()

# sens_data_marg %>%
#   mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
#   filter(ordering == "Ordering") %>%
#   group_by(interval = (ca %/% 0.1)*0.1 + 0.05) %>%
#   summarise(min = min(icc)) %>%
#   ungroup() %>%
#   ggplot(aes(x = interval, y = min)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(breaks = seq(-1, 1, 0.10))
# 
# sens_data_marg %>%
#   filter(ordering == "Ordering") %>%
#   group_by(interval = (pmax(tau_s0t1, tau_s1t0) %/% 0.1)*0.1 + 0.05) %>%
#   summarise(min = min(icc)) %>%
#   ungroup() %>%
#   ggplot(aes(x = interval, y = min)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(breaks = seq(-1, 1, 0.10))

png(filename = "Figures/ovarian_results_exploration_wa.png", width = 13.2, height = 9.9, 
    units = "cm", pointsize = 10, res = 500)
sens_data_marg %>%
  mutate(ca = 0.5*(tau_s0t1 + tau_s1t0), wa = 0.5*(tau_s0s1 + tau_t0t1)) %>%
  mutate(Monotonicity = ifelse(pmin(tau_s0s1, tau_t0t1) > 0,
                      "Yes", "No")) %>%
  filter(ordering == "Ordering") %>%
  ggplot(aes(x = wa, y = 1 - exp(-2*minfo), color = Monotonicity)) +
  scale_x_continuous(limits = c(-1, 1), 
                     name = TeX("$\\frac{\\tau_{S_0S_1} + \\tau_{T_0T_1}}{2}$")) +
  scale_y_continuous(name = TeX("$\\R_h^2")) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2), legend.background = element_blank())
dev.off()


