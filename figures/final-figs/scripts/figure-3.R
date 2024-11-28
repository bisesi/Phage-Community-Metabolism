# ATB
# figure 3
# partner dynamics

# load packages
library("tidyverse")
library("readxl")
library("patchwork")
library("ggtext")
library("ggpubr")
library("cowplot")
library("chromatographR")
library("broom")
library("ggVennDiagram")
library("ggfortify")
library("ggpubfigs")

# part A - community composition
date <- "6November2023"
no_infection <- FALSE

source(here::here("data-generation", "experimental", "clean-mutualism-tecan-data.R"))

partA <- all_data_corrected %>%
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  mutate(M_corrected_OD = ifelse(partner %in% c("M0104", "S0240 M0104"), abs(OD - E_corrected_OD - S_corrected_OD - 0.05), 0),
         M_corrected_OD = ifelse(strain == "M0104", OD - 0.08, M_corrected_OD)) %>%
  pivot_longer(cols = c("E_corrected_OD", "S_corrected_OD", "M_corrected_OD"), names_to = "species") %>%
  mutate(species = case_when(species == "E_corrected_OD" ~ "*E. coli*",
                             species == "S_corrected_OD" ~ "*S. enterica*",
                             species == "M_corrected_OD" ~ "*M. extorquens*")) %>%
  filter(partner != "JS002" | is.na(partner)) %>%
  filter(!strain %in% c("S0240", "M0104")) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(is.na(partner) ~ "mono",
                             partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("mono", "*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(value = case_when(partner == "*S. enterica*" & species == "*M. extorquens*" ~ 0,
                           partner == "*M. extorquens*" & species == "*S. enterica*" ~ 0,
                           partner == "mono" & species == "*S. enterica*" ~ 0,
                           TRUE ~ value),
         value = na_if(value, 0)) %>%
  filter(well != "G8") %>% 
  filter(hour < 125) %>% filter(partner != "mono") %>%
  ggplot(aes(x = hour, y = value, color = species, alpha = well)) +
  geom_smooth(method = "loess", span = 0.25, se = FALSE) +
  facet_grid(partner~strain) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("*E. coli*" = "#004488", "*S. enterica*" = "#DDAA33", "*M. extorquens*" = "#BB5566")) +
  ylab("inferred species-specific OD600") +
  xlab("hour") +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_markdown())

# part B - composition
# S
# M
# SM
date <- "6November2023"
source(here::here("data-generation", "experimental", "clean-mutualism-tecan-data.R"))
nov6 <- all_data_corrected %>% 
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  mutate(M_corrected_OD = ifelse(partner %in% c("M0104", "S0240 M0104"), abs(OD - E_corrected_OD - S_corrected_OD - 0.05), 0),
         M_corrected_OD = ifelse(strain == "M0104", OD - 0.08, M_corrected_OD)) %>%
  pivot_longer(cols = c("E_corrected_OD", "S_corrected_OD", "M_corrected_OD"), names_to = "species") %>%
  mutate(species = case_when(species == "E_corrected_OD" ~ "*E. coli*",
                             species == "S_corrected_OD" ~ "*S. enterica*",
                             species == "M_corrected_OD" ~ "*M. extorquens*")) %>%
  mutate(date = date)

date <- "20November2023"
source(here::here("data-generation", "experimental", "clean-mutualism-tecan-data.R"))
nov20 <- all_data_corrected %>%
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  mutate(M_corrected_OD = ifelse(partner %in% c("M0104", "S0240 M0104"), abs(OD - E_corrected_OD - S_corrected_OD - 0.05), 0),
         M_corrected_OD = ifelse(strain == "M0104", OD - 0.08, M_corrected_OD)) %>%
  pivot_longer(cols = c("E_corrected_OD", "S_corrected_OD", "M_corrected_OD"), names_to = "species") %>%
  mutate(species = case_when(species == "E_corrected_OD" ~ "*E. coli*",
                             species == "S_corrected_OD" ~ "*S. enterica*",
                             species == "M_corrected_OD" ~ "*M. extorquens*")) %>%
  mutate(date = date)

partB <- nov6 %>% filter(well != "G8") %>% rbind(., nov20) %>% filter(!strain %in% c("S0240", "M0104")) %>% filter(is.na(partner) | partner != "JS002") %>%
  group_by(well, date) %>% filter(OD == max(OD)) %>% slice_head(n = 3) %>% ungroup() %>%
  mutate(value = case_when(is.na(partner) & species == "*S. enterica*" ~ 0,
                           partner == "S0240" & species == "*M. extorquens*" ~ 0,
                           partner == "M0104" & species == "*S. enterica*" ~ 0,
                           TRUE ~ value)) %>%
  select(well, OD, strain, partner, species, value, date) %>%
  pivot_wider(names_from = species, values_from = value) %>%
  mutate(percent_partner = (`*S. enterica*` + `*M. extorquens*`) / (`*S. enterica*` + `*M. extorquens*` + `*E. coli*`)) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(is.na(partner) ~ "mono",
                             partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("mono", "*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(comparison = "composition") %>% rename(value = percent_partner) %>%
  filter(partner != "mono") %>% ggplot(aes(x = strain, y = value, color = strain)) +
  geom_point(size = 2, stroke = 1.5) +
  stat_summary(aes(fill = strain), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  facet_wrap(~partner, scales = "free_y", ncol = 3) +
  ylab("percent partner") +
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "F128+<br>M13+" = "#0F2080")) +
  theme_bw(base_size = 16) + theme(axis.text = element_markdown(), axis.title.x = element_blank(),
                                   legend.position = "none", strip.background = element_blank(),
                                   strip.text = element_markdown())

# part C - function
partC <- nov6 %>% filter(well != "G8") %>% rbind(., nov20) %>% filter(!strain %in% c("S0240", "M0104")) %>% filter(is.na(partner) | partner != "JS002") %>%
  group_by(well, date) %>% filter(OD == max(OD)) %>% slice_head(n = 3) %>% ungroup() %>%
  mutate(value = case_when(is.na(partner) & species == "*S. enterica*" ~ 0,
                           partner == "S0240" & species == "*M. extorquens*" ~ 0,
                           partner == "M0104" & species == "*S. enterica*" ~ 0,
                           TRUE ~ value)) %>%
  select(well, OD, strain, partner, species, value, date) %>%
  pivot_wider(names_from = species, values_from = value) %>%
  mutate(total_OD = `*S. enterica*` + `*M. extorquens*` + `*E. coli*`) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(is.na(partner) ~ "mono",
                             partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("mono", "*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(comparison = "function") %>% rename(value = total_OD) %>% filter(partner != "mono") %>% 
  ggplot(aes(x = strain, y = value, color = strain)) +
  geom_point(size = 2, stroke = 1.5) +
  stat_summary(aes(fill = strain), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  facet_wrap(~partner, scales = "free_y", ncol = 3) +
  ylab("max OD600") +
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "F128+<br>M13+" = "#0F2080")) +
  theme_bw(base_size = 16) + theme(axis.text = element_markdown(), axis.title.x = element_blank(),
                                   legend.position = "none", strip.background = element_blank(),
                                   strip.text = element_markdown())

# final figure
left <- plot_grid(partB, partC, ncol = 1, label_size = 26, labels = c("B", "C"))
figure3 <- plot_grid(partA, left, ncol = 2, label_size = 26, labels = c("A", ""), rel_widths = c(0.8,1))

png(here::here("figures", "final-figs", "imgs", "figure-3.png"), res = 300, width = 4000, height = 2000)
figure3
dev.off()






