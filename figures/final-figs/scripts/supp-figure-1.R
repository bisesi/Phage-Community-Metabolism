# ATB
# supp figure 1
# COMETS modeling

#load packages
library("tidyverse")
library("readxl")
library("ggtext")
library("cowplot")
library("ggpubr")
library("rstatix")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

# part A - ratio of amino acids
setwd(here::here("fba-data", "figure-1", "cobra"))
reactions <- read_csv("core_reactions.csv") %>% select(-c("...1"))
amino_acids <- c("ala__L_c", "arg__L_c", "asn__L_c", "asp__L_c", "cys__L_c", "gln__L_c", "glu__L_c", "gly_c", "his__L_c", "ile__L_c",
                                       "leu__L_c", "lys__L_c", "met__L_c", "pro__L_c", "ser__L_c", "thr__L_c", "trp__L_c", "tyr__L_c", "val__L_c", "phe__L_c")

bolded <- c("glycine", "arginine", "glutamine", "glutamate", "asparagine", "aspartate", "histidine", "cysteine")

partA <- reactions %>% filter(Chemical %in% amino_acids) %>%
  inner_join(., reactions %>% filter(Chemical %in% amino_acids) %>% group_by(reaction) %>%
               summarize(total_aa = sum(Coefficient)), by = "reaction") %>% 
  mutate(percent_aa = Coefficient / total_aa) %>%
  mutate(reaction = case_when(reaction == "host" ~ "*E. coli*",
                              reaction == "plasmid" ~ "F128",
                              reaction == "virus" ~ "M13")) %>%
  mutate(Chemical = case_when(Chemical == "ala__L_c" ~ "alanine",
                              Chemical == "arg__L_c" ~ "arginine",
                              Chemical == "asn__L_c" ~ "asparagine",
                              Chemical == "asp__L_c" ~ "aspartate",
                              Chemical == "cys__L_c" ~ "cysteine",
                              Chemical == "gln__L_c" ~ "glutamine",
                              Chemical == "glu__L_c" ~ "glutamate",
                              Chemical == "gly_c" ~ "glycine",
                              Chemical == "his__L_c" ~ "histidine",
                              Chemical == "ile__L_c" ~ "isoleucine",
                              Chemical == "leu__L_c" ~ "leucine",
                              Chemical == "lys__L_c" ~ "lysine",
                              Chemical == "met__L_c" ~ "methionine",
                              Chemical == "pro__L_c" ~ "proline",
                              Chemical == "ser__L_c" ~ "serine",
                              Chemical == "thr__L_c" ~ "threonine",
                              Chemical == "trp__L_c" ~ "tryptophan",
                              Chemical == "tyr__L_c" ~ "tyrosine",
                              Chemical == "val__L_c" ~ "valine",
                              Chemical == "phe__L_c" ~ "phenylalanine")) %>%
  mutate(Chemical = ifelse(Chemical %in% bolded, paste0("***", Chemical, "***"), Chemical)) %>%
  ggplot(aes(x = fct_reorder(Chemical, percent_aa), y = percent_aa, fill = reaction)) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 16) +
  coord_flip() + ylab("percent of total aa stoichiometry") + 
  theme(axis.title.y = element_blank(), axis.text = element_markdown(), legend.position = "none", 
        legend.text = element_markdown(), legend.title = element_blank()) +
  scale_fill_manual(values = c("*E. coli*" = "black", "F128" = "#117733", "M13" = "#88CCEE"))

# part B - time series
setwd(here::here("fba-data", "supp-figure-1"))
es <- read_csv("ES_biomass_range.csv") %>% select(-c("...1")) %>% rename(E0 = iJO1366, S0 = STM_v1_0) %>%
  pivot_longer(cols = c(E0, S0)) %>% mutate(interaction = "*S. enterica*")
em <- read_csv("EM_biomass_range.csv") %>% select(-c("...1")) %>% rename(E0 = iJO1366, M0 = Mrubrum_AM1_hprKO) %>%
  pivot_longer(cols = c(E0, M0)) %>% mutate(interaction = "*M. extorquens*")
esm <- read_csv("ESM_biomass_range.csv") %>% select(-c("...1")) %>% rename(E0 = iJO1366, S0 = STM_v1_0, M0 = Mrubrum_AM1_hprKO) %>%
  pivot_longer(cols = c(E0, S0, M0)) %>% mutate(interaction = "*S. enterica*<br>+<br>*M. extorquens*")

setwd(here::here("fba-data", "figure-1", "comets"))
e <- read_csv("E_biomass_range.csv") %>% select(-c("...1")) %>% rename(E0 = iJO1366)

bounds <- e %>% mutate(optimized = ifelse(optimized == "virus" & lower_bound == 0, "host", optimized)) %>% pull(lower_bound) %>% unique()
plasmid_bound <- bounds[10]
virus_bound <- bounds[18]

e <- e %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  pivot_longer(cols = E0) %>% mutate(interaction = "monoculture")

partB <- rbind(es,em,esm) %>%
  mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  rbind(., e) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(name = case_when(name == "E0" ~ "*E. coli*",
                          name == "S0" ~ "*S. enterica*",
                          name == "M0" ~ "*M. extorquens*")) %>%
  mutate(interaction = factor(interaction, levels = c("monoculture", "*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  ggplot(aes(x = cycle, y = log10(value), color = name)) +
  ylab("log10(simulated biomass)") +
  scale_color_manual(values = c("*E. coli*" = "#004488", "*S. enterica*" = "#DDAA33", "*M. extorquens*" = "#BB5566")) +
  geom_line(linewidth = 2) + facet_grid(interaction~optimized, scales = "free") + theme_bw(base_size = 16) +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_markdown())

# part C - percent partner
partC <- rbind(es,em,esm) %>%
  mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(name = case_when(name == "E0" ~ "*E. coli*",
                          name == "S0" ~ "*S. enterica*",
                          name == "M0" ~ "*M. extorquens*")) %>% pivot_wider(names_from = name, values_from = value) %>% 
  filter(cycle == max(cycle)) %>%
  mutate(across(c(`*S. enterica*`, `*M. extorquens*`), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(percent_partner = (`*S. enterica*` + `*M. extorquens*`) / (`*S. enterica*` + `*M. extorquens*` + `*E. coli*`)) %>%
  mutate(interaction = factor(interaction, levels = c("*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  ggplot(aes(x = optimized, y = percent_partner, color = optimized)) +
  geom_point(size = 2, stroke = 1.5) + facet_wrap(~interaction) +
  scale_color_manual(values = c("uninf" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16)+
  ylab("percent partner") +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(), strip.text = element_markdown())

# part D - total biomass
partD <- rbind(es,em,esm) %>%
  mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(name = case_when(name == "E0" ~ "*E. coli*",
                          name == "S0" ~ "*S. enterica*",
                          name == "M0" ~ "*M. extorquens*")) %>% pivot_wider(names_from = name, values_from = value) %>% 
  filter(cycle == max(cycle)) %>%
  mutate(across(c(`*S. enterica*`, `*M. extorquens*`), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(total_yield = (`*S. enterica*` + `*M. extorquens*` + `*E. coli*`)) %>%
  mutate(interaction = factor(interaction, levels = c("*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  ggplot(aes(x = optimized, y = log10(total_yield), color = optimized)) +
  geom_point(size = 2, stroke = 1.5) + facet_wrap(~interaction) +
  scale_color_manual(values = c("uninf" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16)+
  ylab("log10(max biomass)") +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(), strip.text = element_markdown())

# part E - what are the partners using
setwd(here::here("fba-data", "supp-figure-1"))
es_flux <- read_csv("ES_fluxes_range.csv") %>% select(., c(lower_bound, optimized, cycle, species, contains("EX_"))) %>%
  pivot_longer(cols = -c("lower_bound", "optimized", "cycle", "species"), values_to = "flux") %>%
  mutate(interaction = "*S. enterica*")
em_flux <- read_csv("EM_fluxes_range.csv") %>% select(., c(lower_bound, optimized, cycle, species, contains("EX_"))) %>%
  pivot_longer(cols = -c("lower_bound", "optimized", "cycle", "species"), values_to = "flux") %>%
  mutate(interaction = "*M. extorquens*")
esm_flux <- read_csv("ES_fluxes_range.csv") %>% select(., c(lower_bound, optimized, cycle, species, contains("EX_"))) %>%
  pivot_longer(cols = -c("lower_bound", "optimized", "cycle", "species"), values_to = "flux") %>%
  mutate(interaction = "*S. enterica*<br>+<br>*M. extorquens*")

compounds_SM <- c("EX_ac_e", "EX_cit_e", "EX_for_e", "EX_succ_e", "EX_hxa_e") 

compounds_E <- c("EX_hxa_e", "EX_ac_e", "EX_for_e", "EX_glyclt_e", "EX_succ_e", "EX_cit_e") 

shared <- c("acetate", "hexanoate", "formate")

uptake <- rbind(es_flux, esm_flux, em_flux) %>%
  filter(species != "iJO1366") %>%
  mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  filter(flux != 0) %>%
  group_by(species, optimized, name, interaction) %>%
  filter(flux == min(flux)) %>%
  filter(name %in% compounds_SM) %>%
  slice_head(n = 1) %>% ungroup() %>% filter(flux < 0) %>%
  dplyr::select(name, interaction, flux, species, optimized) %>%
  ungroup() %>% pivot_wider(names_from = optimized, values_from = flux) %>%
  mutate(across(c(host, plasmid, `virus with plasmid`), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(log_plasmid = log2(abs(plasmid) /abs(host)),
         log_virusplasmid = log2(abs(`virus with plasmid`) / abs(host))) %>%
  dplyr::select(-c(host, plasmid, `virus with plasmid`)) %>%
  pivot_longer(cols = c("log_plasmid", "log_virusplasmid"), names_to = "type") %>%
  mutate(type = ifelse(type == "log_plasmid", "F128+", "F128+<br>M13+")) %>%
  mutate(name = case_when(name == "EX_ac_e" ~ "acetate",
                          name == "EX_succ_e" ~ "succinate",
                          name == "EX_for_e" ~ "formate",
                          name == "EX_cit_e" ~ "citrate",
                          name == "EX_hxa_e" ~ "hexanoate")) %>%
  filter(interaction != "*S. enterica*<br>+<br>*M. extorquens*") %>%
  mutate(interaction = factor(interaction, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  mutate(name = ifelse(name %in% shared, paste0("***", name, "***"), name))

partE <- uptake %>%
  ggplot(aes(x = fct_reorder(name, value), y = value, fill = type)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_wrap(~interaction, scales = "free") + coord_flip() +
  scale_fill_manual(values = c("F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16)+
  ylab("log2(fold-change in partner max uptake)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(axis.text = element_markdown(), axis.title.y = element_blank(),
        legend.position = "none", strip.background = element_blank(), strip.text = element_markdown())
  
# part F - was is E spitting out
secretion <- rbind(es_flux, esm_flux, em_flux) %>%
  filter(species == "iJO1366") %>%
  mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  filter(flux != 0) %>%
  group_by(species, optimized, name, interaction) %>%
  filter(flux == max(flux)) %>%
  filter(name %in% compounds_E) %>%
  slice_head(n = 1) %>% ungroup() %>% filter(flux > 0) %>%
  dplyr::select(name, interaction, flux, species, optimized) %>%
  ungroup() %>% pivot_wider(names_from = optimized, values_from = flux) %>%
  mutate(across(c(host, plasmid, `virus with plasmid`), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(log_plasmid = log2(plasmid / host),
         log_virusplasmid = log2(`virus with plasmid` / host)) %>%
  dplyr::select(-c(host, plasmid, `virus with plasmid`)) %>%
  pivot_longer(cols = c("log_plasmid", "log_virusplasmid"), names_to = "type") %>%
  mutate(type = ifelse(type == "log_plasmid", "F128+", "F128+<br>M13+")) %>%
  filter(interaction != "*S. enterica*<br>+<br>*M. extorquens*") %>%
  mutate(interaction = factor(interaction, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  mutate(name = case_when(name == "EX_hxa_e" ~ "hexanoate",
                          name == "EX_ac_e" ~ "acetate",
                          name == "EX_for_e" ~ "formate",
                          name == "EX_glyclt_e" ~ "glycolate",
                          name == "EX_succ_e" ~ "succinate",
                          name == "EX_cit_e" ~ "citrate")) %>%
  mutate(name = ifelse(name %in% shared, paste0("***", name, "***"), name))

partF <- secretion %>%
  ggplot(aes(x = fct_reorder(name, value), y = value, fill = type)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_wrap(~interaction, scales = "free") + coord_flip() +
  scale_fill_manual(values = c("F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16)+
  ylab("log2(fold-change in *E. coli* max secretion)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(axis.text = element_markdown(), axis.title.y = element_blank(), axis.title.x = element_markdown(),
        legend.position = "none", strip.background = element_blank(), strip.text = element_markdown())

# final figure
right <- plot_grid(partC, partD, ncol = 2, label_size = 26, labels = c("C", "D"))
left <- plot_grid(partA, partB, ncol = 2, label_size = 26, labels = c("A", "B"))
top <- plot_grid(left, right, ncol = 1, rel_heights = c(1, 0.6))
bottom <- plot_grid(partF, partE, ncol = 2, labels = c("E", "F"), label_size = 26)
suppfigure1 <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 0.5))

png(here::here("figures", "final-figs", "imgs", "supp-figure-1.png"), res = 300, width = 4500, height = 4500)
suppfigure1
dev.off()

