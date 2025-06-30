# ATB
# supp figure 2
# conjugation rates

#load packages
library("tidyverse")
library("readxl")
library("ggtext")
library("cowplot")
library("ggpubr")
library("rstatix")
library("scales")

# part A - lactose utilization
lactose <- data.frame(partner = c("F128+", "F128+"), interaction_type = c("mutualism", "facilitation"),
                      cfus = c(0.85, 0.85))

partA <- lactose %>%
  ggplot(aes(x = interaction_type, y = cfus)) +
  geom_point(size = 2, stroke = 1.5)+
  ylab("Percent lactose positive *S. enterica*<br>at the end of co-culture growth in liquid") +
  theme_bw(base_size = 16)+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  ylim(0, 10)+
  annotate("text", x = 1.5, y = 1.5, label = "LOD")+
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown(), strip.text = element_markdown())

# part B - conjugation on plates
dilution_rate = 10^-5
plated_volume_mL = 0.1
data_from_LC <- data.frame(media = c("lacmet", "lacmet", "lacmet", "lacmet", "lacmet",
                                     "glucmet", "glucmet", "glucmet", "glucmet", "glucmet"),
                           citrate_num = c(178, 183, 207, 226, 95, 412, 428, 239, 220, 202),
                           citrate_tet_num = c(27, 32, 38, 37, 12, 26, 25, 31, 19, 17))

stats <- data_from_LC %>%
  mutate(media = ifelse(media == "lacmet", "facilitation","competition")) %>%
  mutate(total_rate = citrate_num / (plated_volume_mL * dilution_rate),
         tet_rate = citrate_tet_num / (plated_volume_mL * dilution_rate),
         conj_rate = tet_rate / total_rate) %>%
  t_test(conj_rate ~ media) %>%
  add_significance("p") %>%
  add_xy_position(x = "media") %>%
  mutate(label = paste0(scientific(p,digits = 2), " (", p.signif, ")"))

partB <- data_from_LC %>%
  mutate(media = ifelse(media == "lacmet", "facilitation","competition")) %>%
  mutate(total_rate = citrate_num / (plated_volume_mL * dilution_rate),
         tet_rate = citrate_tet_num / (plated_volume_mL * dilution_rate),
         conj_rate = tet_rate / total_rate) %>%
  ggplot(aes(x = media, y = conj_rate * 100)) +
  stat_summary(aes(fill = media), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5)+
  ylab("Percent plasmid positive *S. enterica*<br>at the end of co-culture growth on agar") +
  stat_pvalue_manual(stats, label = "label", tip.length = 0.0, size = 3) +
  theme_bw(base_size = 16)+
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown(), strip.text = element_markdown())

# part C - liquid conjugation
liquid <- data.frame(partner = c("F128+", "F128+"), interaction_type = c("mutualism", "facilitation"),
                      cfus = c(0.85, 0.85))

partC <- liquid %>%
  ggplot(aes(x = interaction_type, y = cfus)) +
  geom_point(size = 2, stroke = 1.5)+
  ylab("Percent plasmid positive *S. enterica*<br>at the end of co-culture growth in liquid") +
  theme_bw(base_size = 16)+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  ylim(0,10)+
  annotate("text", x = 1.5, y = 1.5, label = "LOD")+
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown(), strip.text = element_markdown())


# final figure
suppfigure2 <- plot_grid(partA, partB, partC, ncol = 3, label_size = 26, labels = c("A", "B", "C"))

png(here::here("figures", "final-figs", "imgs", "supp-figure-2.png"), res = 300, width = 3000, height = 1500)
suppfigure2
dev.off()
