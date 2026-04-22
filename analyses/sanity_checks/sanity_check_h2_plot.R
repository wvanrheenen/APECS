library(tidyverse)
library(wesanderson)


df = read.table("./analyses/simPed_h2sim.results.out", header=T, sep="\t")

# h2 by K
pdf("simPed_h2sim_by_prevalence.pdf", w=3, h=3)
filter(df, lambda==2, k==2, K > 0.001, K < 0.2) %>%
  mutate(Kf = as.factor(K),
         h2f = as.factor(h2)) %>%
  ggplot(., aes(x=h2f, y=h2sim, fill=Kf)) + 
    geom_boxplot(lwd=0.4, fatten=1, outlier.size = 0.5) +
    theme_bw() +
    theme(legend.position="top",
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(color = "black")) +
    xlab(expression(italic("h"^2))) +
    scale_y_continuous("estimated heritability", 
                       breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       limits = c(-0.05,1)) +
  scale_fill_manual(values = wes_palette(n=5, name="Darjeeling1"), 
                      name = "K")
dev.off()

# h2 by k
pdf("simPed_h2sim_by_N_generations.pdf", w=3, h=3)
filter(df, K == 0.1, lambda == 2 ) %>%
  mutate(kf = as.factor(k),
         h2f = as.factor(h2)) %>%
  ggplot(., aes(x=h2f, y=h2sim, fill=kf)) + 
  geom_boxplot(lwd=0.4, fatten=1, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position="top",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  xlab(expression(italic("h"^2))) +
  scale_y_continuous("estimated heritability", 
                     breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.05,1)) +
  scale_fill_manual(values = wes_palette(n=4, name="Darjeeling1"),  
                    name = "generations")
dev.off()
# h2 by lambda
pdf("simPed_h2sim_by_N_offspring.pdf", w=3, h=3)
filter(df, K == 0.1, k == 2 ) %>%
  mutate(lf = as.factor(lambda),
         h2f = as.factor(h2)) %>%
  ggplot(., aes(x=h2f, y=h2sim, fill=lf)) + 
  geom_boxplot(lwd=0.4, fatten=1, outlier.size = 0.5) +
  theme_bw() +
  theme(legend.position="top",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  xlab(expression(italic("h"^2))) +
  scale_y_continuous("estimated heritability", 
                     breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.05,1)) +
  scale_fill_manual(values = wes_palette(n=4, name="Darjeeling1"),  
                    name = "Mean offspring")
dev.off()