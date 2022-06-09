library(ggplot2)
library(latex2exp)
library(ggpubr)
p = ggplot() +  
  geom_polygon(mapping = aes(x = c(0, 0, 1), 
                             y = c(1, 0, 1)), fill = "gray") +
  geom_abline(slope = 1) +
  annotate(geom = "text", x = 0.35, y = 0.7, label = TeX(string = "$S_0^* < T_0$")) +
  scale_x_continuous(breaks = NULL, name = TeX("Time to non-terminal event, $S_0^*$"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = NULL, name = TeX("Time to terminal event, $T_0$"),
                     expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm")))) +
  ggtitle(TeX("$f_U (s_0^*, t_0)$"))
p
x = seq(0, 1, 0.001)
d = dweibull(x = x, shape = 1.4, scale = 0.27)
q = ggplot() +
  scale_x_continuous(breaks = NULL, expand = c(0,0), name = NULL) +
  scale_y_continuous(breaks = NULL, expand = c(0,0), name = TeX("$T_0 = \\infty$")) +
  geom_line(mapping = aes(x = x, y = d)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.line.x = element_blank(),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm")))) +
  coord_flip() +
  ggtitle(TeX("f_{\\infty}(t_0)"))
q

ggarrange(p, q, ncol = 2, widths = c(2, 1))
