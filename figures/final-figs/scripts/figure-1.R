# ATB
# figure 1
# modeling plus system

#load packages
library("tidyverse")
library("readxl")
library("ggtext")
library("cowplot")
library("ggpubr")
library("rstatix")
library("ggpubfigs")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

#load data
setwd(here::here("fba-data", "figure-1", "cobra"))
pfba <- read_csv("pfba_fluxes_optimized.csv", show_col_types = FALSE) %>% rename(reaction = `...2`)
fva <- read_csv("fva_fluxes_optimized.csv", show_col_types = FALSE) %>% rename(reaction = `...2`)
subprocesses <- read_csv("subprocesses.csv", show_col_types = FALSE) %>% select(-c("...1"))
pfba_plasmid_test <- read_csv("pfba_fluxes_virusplasmid_plasmidrange_virusoptimized.csv") %>% rename(reaction = `...1`)
fva_plasmid_test <- read_csv("fva_fluxes_virusplasmid_plasmidrange_virusoptimized.csv") %>% rename(reaction = `...1`)
i = unique(fva_plasmid_test$plasmid_bound)[10]

# part A - system schematic, empty plot
partA <- data.frame(x = 1, y = 1) %>%
  ggplot(aes(x = x, y = y)) + theme_void()

# part B - modeling diagrams plus pie charts for allocation, empty plot
partB <- data.frame(x = 1, y = 1) %>%
  ggplot(aes(x = x, y = y)) + theme_void()

# pie charts for part B
reactions <- read_csv("core_reactions.csv") %>% select(-c("...1"))
reactions <- reactions %>% filter(Chemical != "h2o_c") %>% filter(Chemical != "atp_c")
nucleic_acid <- data.frame(Chemical = c("datp_c", "dctp_c", "dgtp_c", "dttp_c", "gtp_c", "utp_c", "ctp_c", "nadp_c", "nad_c", "bmocogdp_c", "coa_c", "fad_c"), type = "nucleic acids")
amino_acids <- data.frame(Chemical = c("ala__L_c", "arg__L_c", "asn__L_c", "asp__L_c", "cys__L_c", "gln__L_c", "glu__L_c", "gly_c", "his__L_c", "ile__L_c",
                 "leu__L_c", "lys__L_c", "met__L_c", "pro__L_c", "ser__L_c", "thr__L_c", "trp__L_c", "tyr__L_c", "val__L_c", "phe__L_c", "amet_c", "mlthf_c"), type = "amino acids")
ions <- data.frame(Chemical = c("ca2_c", "cl_c", "cu2_c", "fe2_c", "fe3_c", "k_c", "mg2_c", "mn2_c", "mobd_c", "ni2_c", "zn2_c", "cobalt2_c", "nh4_c", "so4_c"), type = "ions")
vitamin <- data.frame(Chemical = c("btn_c", "thmpp_c", "ribflv_c", "pydx5p_c", "thf_c", "10fthf_c"), type = "vitamins")
hemes <- data.frame(Chemical = c("2fe2s_c", "4fe4s_c", "sheme_c", "pheme_c"), type = "hemes")
lipids <- data.frame(Chemical = c("2ohph_c", "kdo2lipid4_e", "pe160_c", "pe160_p", "pe161_c", "pe161_p", "udcpdp_c", "murein5px4p_p"), type = "lipids")

partC <- reactions %>% full_join(., rbind(nucleic_acid, amino_acids, ions, vitamin, lipids, hemes), by = "Chemical") %>%
  mutate(type = ifelse(type != "nucleic acids" & type != "amino acids", "other", type)) %>%
  group_by(reaction) %>% mutate(total = sum(Coefficient)) %>% ungroup() %>% group_by(reaction, type, total) %>% summarize(sub = sum(Coefficient)) %>%
  ungroup() %>% mutate(percent = sub / total) %>% mutate(percent = case_when(reaction == "plasmid" & type == "amino acids" ~ 0.995,
                                                                             reaction == "plasmid" & type == "nucleic acids" ~ 0.005,
                                                                             TRUE ~ percent)) %>%
  mutate(reaction = ifelse(reaction == "host", "bacteria", reaction)) %>%
  mutate(reaction = ifelse(reaction == "plasmid", "plasmid", reaction)) %>% 
  mutate(reaction = ifelse(reaction == "virus", "phage", reaction)) %>%
  mutate(reaction = factor(reaction, levels = c("bacteria", "plasmid", "phage"))) %>%
  ggplot(aes(x = "", y = percent, fill = type)) + geom_bar(stat = "identity", width = 1) + 
  scale_fill_manual(values = friendly_pal("vibrant_seven")) +
  coord_polar("y", start = 0) + facet_wrap(~reaction, ncol = 1) + theme_void(base_size = 16) + labs(fill = "")
  
# part D - non-overlapping fluxes
not_included <- which(!c(fva %>% filter(optimized == "plasmid") %>% select(reaction) %>% pull() %>% unique()) %in% c(subprocesses %>% select(reaction) %>% pull() %>% unique()))

exchanges <- c(fva %>% filter(optimized == "plasmid") %>% select(reaction) %>% pull())[not_included]

subprocesses <- rbind(subprocesses, data.frame(reaction = exchanges, group_name = "Exchanges", id = "NA"))

counted_subprocesses <- subprocesses %>% group_by(group_name) %>% summarize(n = n())

host_fva <- fva %>%  
  filter(optimized == "host") %>%
  full_join(., subprocesses, by = "reaction", relationship = "many-to-many") %>%
  select(reaction, minimum, maximum, optimized, group_name) %>%
  filter(optimized %in% c("host", "virus with plasmid")) %>%
  pivot_wider(names_from = optimized, values_from = c("minimum", "maximum"))

plasmid_ranges <- fva %>%  
  full_join(., subprocesses, by = "reaction", relationship = "many-to-many") %>%
  select(reaction, minimum, maximum, optimized, group_name) %>%
  filter(optimized %in% c("host", "plasmid")) %>%
  pivot_wider(names_from = optimized, values_from = c("minimum", "maximum"))

plasmidvirus_ranges <- fva_plasmid_test %>%  
  filter(plasmid_bound == i) %>%
  full_join(., subprocesses, by = "reaction", relationship = "many-to-many") %>%
  select(reaction, minimum, maximum, optimized, group_name) %>%
  filter(optimized %in% c("host", "virus with plasmid")) %>%
  pivot_wider(names_from = optimized, values_from = c("minimum", "maximum")) %>%
  inner_join(., host_fva, by = c("reaction", "group_name"))

plasmid_ranges$overlap <- (plasmid_ranges$minimum_plasmid <= plasmid_ranges$maximum_host) & (plasmid_ranges$minimum_host >= plasmid_ranges$maximum_plasmid)
plasmidvirus_ranges$overlap <- (plasmidvirus_ranges$`minimum_virus with plasmid` <= plasmidvirus_ranges$maximum_host) & (plasmidvirus_ranges$minimum_host >= plasmidvirus_ranges$`maximum_virus with plasmid`)

plasmid_by_process <- plasmid_ranges %>%
  filter(overlap == TRUE) %>% filter(!((minimum_plasmid == 0 & minimum_host == 0) | (maximum_plasmid == 0 & maximum_host == 0))) %>%
  filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13")) %>%
  group_by(group_name) %>%
  mutate(number_per_group = n()) %>% 
  ungroup() %>%
  inner_join(., counted_subprocesses, by = "group_name") %>%
  mutate(fraction = number_per_group / n) %>%
  select(group_name, fraction) %>%
  unique() %>%
  mutate(comparison = "plasmid")

plasmidvirus_by_process <- plasmidvirus_ranges %>%
  filter(overlap == TRUE) %>% filter(!((`minimum_virus with plasmid` == 0 & minimum_host == 0) | (`maximum_virus with plasmid` == 0 & maximum_host == 0))) %>%
  filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13")) %>%
  group_by(group_name) %>%
  mutate(number_per_group = n()) %>% 
  ungroup() %>%
  inner_join(., counted_subprocesses, by = "group_name") %>%
  mutate(fraction = number_per_group / n) %>%
  select(group_name, fraction) %>%
  unique() %>%
  mutate(comparison = "virus with plasmid")

plasmid_reactions <- plasmid_ranges %>%
  filter(overlap == TRUE) %>% filter(!((minimum_plasmid == 0 & minimum_host == 0) | (maximum_plasmid == 0 & maximum_host == 0))) %>%
  filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13")) %>% pull(reaction) %>% unique()

plasmidvirus_reactions <- plasmidvirus_ranges %>%
  filter(overlap == TRUE) %>% filter(!((`minimum_virus with plasmid` == 0 & minimum_host == 0) | (`maximum_virus with plasmid` == 0 & maximum_host == 0))) %>%
  filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13")) %>% pull(reaction) %>% unique()

partD <- data.frame(comparison = c("plasmid", "virus with plasmid"),
                    total_number = c(plasmid_reactions %>% length(), plasmidvirus_reactions %>% length()),
                    all_reactions = rep(subprocesses %>% 
                                          filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13")) %>% select(reaction) %>% pull() %>% length(), 2)) %>%
  mutate(comparison = case_when(comparison == "plasmid" ~ "F128+",
                                comparison == "virus with plasmid" ~ "F128+<br>M13+"),
         comparison = factor(comparison, levels = c("F128+", "F128+<br>M13+"))) %>%
  ggplot(aes(x = comparison, y = (total_number / all_reactions) * 100, fill = comparison)) +
  geom_bar(stat = "identity")+
  ylab("Percent of reactions in<br>conflict with WT *E. coli*")+
  xlab("") +
  labs(fill = "", shape = "") +
  theme_bw(base_size = 16)+
  scale_fill_manual(values = c("F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown())

# only F128 - EX_met__L_e (exchanges)
# only F128 M13 - A5PISO, ATPS4rpp, K2L4Aabcpp, K2L4abctex, KDOCT2, KDOPP, KDOPS, LPADSS, MOAT, MOAT2, TDSK, U23GAAT, UAGAAT, UHGADA, USHD

# part E - percent non overlapping per subprocess
partE <- plasmid_by_process %>%
  rbind(., plasmidvirus_by_process) %>%
  mutate(group_name = case_when(group_name == "Cofactor and Prosthetic Group Biosynthesis" ~ "Cofactor Biosynthesis",
                                TRUE ~ group_name)) %>%
  mutate(group_name = tolower(group_name)) %>%
  rbind(., data.frame(group_name = c("lipopolysaccharide biosynthesis / recycling", "oxidative phosphorylation", "alternate carbon metabolism", 
                                     "transport, inner membrane"),
                      fraction = c(0.000005, 0.000005, 0.000005, 0.000005),
                      comparison = c("F128+", "F128+", "F128+", "F128+"))) %>%
  mutate(comparison = case_when(comparison == "plasmid" ~ "F128+",
                                comparison == "virus with plasmid" ~ "F128+<br>M13+"),
         comparison = factor(comparison, levels = c("F128+", "F128+<br>M13+"))) %>%
  ggplot(aes(x = fct_reorder(group_name, fraction), y = fraction * 100, fill = comparison)) +
  geom_bar(stat = "identity", position = position_dodge(0.9))+
  coord_flip()+
  ylab("Percent of reactions in subprocess in <br>conflict with WT *E. coli*") +
  labs(fill = "", shape = "") +
  theme_bw(base_size = 16)+
  scale_fill_manual(values = c("F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme(axis.text = element_markdown(), axis.title.y = element_blank(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_markdown())

# alt carbon - A5PISO
# ox phs - ATPS4rpp

# part F - log2 fold change
switches <- c("EX_pi_e", "PIt2rpp", "PItex", "VPAMTr")

host_plasmid <- pfba %>% select(-c(`...1`)) %>%
  mutate(fluxes = round(fluxes, 10)) %>%
  inner_join(., subprocesses, by = "reaction", relationship = "many-to-many") %>%
  filter(parasite == "bacteria" | parasite == "plasmid") %>%
  mutate(fluxes = abs(fluxes)) %>%
  filter(!reaction %in% switches) %>%
  pivot_wider(names_from = parasite, values_from = fluxes) %>%
  filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13"))

partF <- pfba_plasmid_test %>%
  mutate(fluxes = round(fluxes, 10)) %>%
  filter(!reaction %in% switches) %>%
  inner_join(., subprocesses, by = "reaction", relationship = "many-to-many") %>%
  filter(plasmid_bound == i) %>%
  mutate(fluxes = abs(fluxes)) %>%
  pivot_wider(names_from = parasite, values_from = fluxes) %>%
  filter(!reaction %in% c("vrxn1", "vrxn2", "vrxn3", "vrxn4", "vrxn5", "plasmid_F", "phage_M13")) %>%
  inner_join(., host_plasmid %>% select(-c(id)), by = c("reaction", "group_name")) %>%
  group_by(group_name) %>%
  mutate(subprocess_sum_bacteria = sum(bacteria),
         subprocess_sum_plasmid = sum(plasmid),
         subprocess_sum_virusplasmid = sum(`virus with plasmid`)) %>%
  ungroup() %>%
  mutate(overall_sum_bacteria = sum(bacteria),
         overall_sum_plasmid = sum(plasmid),
         overall_sum_virusplasmid = sum(`virus with plasmid`)) %>%
  ungroup() %>%
  mutate(index_plasmid = log2((subprocess_sum_plasmid / overall_sum_plasmid) / (subprocess_sum_bacteria / overall_sum_bacteria)),
         index_phageplasmid = log2((subprocess_sum_virusplasmid / overall_sum_virusplasmid) / (subprocess_sum_bacteria / overall_sum_bacteria))) %>%
  select(c(index_plasmid, index_phageplasmid, group_name)) %>%
  unique() %>%
  filter(group_name %in% c("Nucleotide Salvage Pathway", "Glycerophospholipid Metabolism",
                           "Cofactor and Prosthetic Group Biosynthesis", "Cell Envelope Biosynthesis",
                           "Arginine and Proline Metabolism", "Membrane Lipid Metabolism",
                           "Lipopolysaccharide Biosynthesis / Recycling", "Unassigned",
                           "Pyruvate Metabolism", "Tyrosine, Tryptophan, and Phenylalanine Metabolism",
                           "Valine, Leucine, and Isoleucine Metabolism", "Citric Acid Cycle",
                           "Cysteine Metabolism", "Purine and Pyrimidine Biosynthesis",
                           "Threonine and Lysine Metabolism",
                           "Methionine Metabolism", "Alanine and Aspartate Metabolism", "Histidine Metabolism",
                           "Pentose Phosphate Pathway", "Anaplerotic Reactions", "Murein Biosynthesis")) %>%
  pivot_longer(cols = index_plasmid:index_phageplasmid, names_to = "comparison") %>%
  mutate(group_name = case_when(group_name == "Tyrosine, Tryptophan, and Phenylalanine Metabolism" ~ "Tyrosine, Tryptophan,Phenylalanine Metabolism", 
                                group_name == "Cofactor and Prosthetic Group Biosynthesis" ~ "Cofactor Biosynthesis",
                                is.na(group_name) ~ "exchanges", 
                                TRUE ~ group_name),
         group_name = tolower(group_name)) %>%
  mutate(comparison = case_when(comparison == "index_plasmid" ~ "F128+", comparison == "index_phageplasmid" ~ "F128+<br>M13+")) %>%
  mutate(direction = ifelse(value > 0, "higher with parasite", "lower with parasite")) %>%
  rbind(., data.frame(group_name = c("histidine metabolism","histidine metabolism", "citric acid cycle", "citric acid cycle", "alanine and aspartate metabolism", "alanine and aspartate metabolism"),
                      value = c(0.000005, 0.000005,0.000005, 0.000005, 0.000005, 0.000005),
                      comparison = c("F128+<br>M13+", "F128+", "F128+<br>M13+", "F128+", "F128+<br>M13+", "F128+"),
                      direction = c("higher with parasite", "lower with parasite", "lower with parasite", "higher with parasite", "lower with parasite", "higher with parasite"))) %>%
  ggplot(aes(x = fct_reorder(group_name, value), y = value, fill = comparison)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  ylab("log2(Fold-change in demand for subprocess)")+
  theme_bw(base_size = 16)+ facet_wrap(~direction, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme(axis.text = element_markdown(), axis.title.y = element_blank(),
        legend.position = "none", strip.background = element_blank())

# ACS and POR5 are most important for pyruvate

# part G - reuptake
setwd(here::here("fba-data", "figure-1", "comets"))
e_fluxes <- read_csv("E_fluxes_range.csv") %>% select(-c("...1", x, y)) %>% select(., c(lower_bound, optimized, cycle, contains("EX_"))) %>%
  pivot_longer(cols = -c("lower_bound", "optimized", "cycle"), values_to = "flux") 
e_met <- read_csv("E_metabolites_range.csv") %>% select(-c("...1")) %>% pivot_longer(cols = -c("lower_bound", "optimized", "cycle"), values_to = "concentration")
e <- read_csv("E_biomass_range.csv") %>% select(-c("...1")) %>% rename(E0 = iJO1366)

bounds <- e %>% mutate(optimized = ifelse(optimized == "virus" & lower_bound == 0, "host", optimized)) %>% pull(lower_bound) %>% unique()
plasmid_bound <- bounds[10]
virus_bound <- bounds[18]

new_e <- e %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+<br>M13+", "F128+")))
  
partG <- e_met %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  filter(name %in% c("ac_e", "glyclt_e", "for_e", "co2_e", "etoh_e", "hxa_e")) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+<br>M13+", "F128+"))) %>%
  mutate(name = case_when(name == "ac_e" ~ "acetate",
                          name == "for_e" ~ "formate",
                          name == "glyclt_e" ~ "glycolate",
                          name == "hxa_e" ~ "hexanoate",
                          name == "co2_e" ~ "carbon dioxide",
                          name == "etoh_e" ~ "ethanol",
                          TRUE ~ name),
         name = factor(name, levels = c("acetate", "formate", "glycolate", "hexanoate", "carbon dioxide", "ethanol"))) %>%
  filter(cycle < 150) %>%
  inner_join(., new_e, by = c("cycle", "optimized", "lower_bound")) %>%
  ggplot(aes(x = cycle, y = log10(concentration / E0), color = optimized)) +
  facet_wrap(~name, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("uninf" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  geom_line(linewidth = 1.5) +
  theme_bw(base_size = 16)+
  labs(color = "") +
  xlab("Simulated hour") + ylab("log10(Biomass-normalized concentration)") +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank())

legendG <- get_plot_component(e_met %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
                                filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
                                filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
                                filter(name %in% c("ac_e", "glyclt_e", "for_e", "co2_e", "etoh_e", "hxa_e")) %>%
                                mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
                                       optimized = ifelse(optimized == "virus with plasmid", "F128+ M13+", optimized),
                                       optimized = ifelse(optimized == "host", "WT", optimized),
                                       optimized = factor(optimized, levels = c("WT", "F128+", "F128+ M13+"))) %>%
                                mutate(name = case_when(name == "ac_e" ~ "acetate",
                                                        name == "for_e" ~ "formate",
                                                        name == "glyclt_e" ~ "glycolate",
                                                        name == "hxa_e" ~ "hexanoate",
                                                        name == "co2_e" ~ "carbon dioxide",
                                                        name == "etoh_e" ~ "ethanol",
                                                        TRUE ~ name),
                                       name = factor(name, levels = c("acetate", "formate", "glycolate", "hexanoate", "carbon dioxide", "ethanol"))) %>%
                                filter(cycle < 150) %>% filter(cycle == 0 & name == "acetate") %>%
                                ggplot(aes(x = cycle, y = concentration, fill = optimized)) +
                                facet_wrap(~name, ncol = 2, scales = "free_y") +
                                scale_fill_manual(values = c("WT" = "black", "F128+" = "#117733", "F128+ M13+" = "#88CCEE")) +
                                geom_bar(stat = "identity") +
                                theme_bw(base_size = 24)+
                                labs(fill = "") +
                                xlab("Simulated hour") + ylab("log(Biomass-normalized concentration)") +
                                theme(axis.text = element_markdown(), legend.text = element_markdown(),
                                      legend.position = "bottom", strip.background = element_blank()),
                              'guide-box-bottom', return_all = TRUE)

all_partG <- plot_grid(partG, legendG, ncol = 1, rel_heights = c(1, 0.05))

# final figure
top <- plot_grid(partA, plot_grid(partB, partC, ncol = 2, labels = c("B", "C"), label_size = 26, rel_widths = c(1,0.8)), 
                 ncol = 1, rel_heights= c(1,0.7), label_size = 26, labels = c("A", ""))
middle <- plot_grid(partD, partE, ncol = 2, rel_widths = c(0.4,1), label_size = 26, labels = c("D", "E"))
bottom <- plot_grid(partF, all_partG, ncol = 2, rel_widths = c(1,0.5), label_size = 26, labels = c("F", "G"))
figure1 <- plot_grid(top, plot_grid(middle, bottom, ncol = 1, rel_heights = c(0.5,1)), ncol = 2, rel_widths = c(0.6,1))

png(here::here("figures", "final-figs", "imgs", "figure-1.png"), res = 300, width = 6500, height = 4000)
figure1
dev.off()




