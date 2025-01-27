# ATB
# figure 4
# partner metabolomics

# load packages
library("tidyverse")
library("readxl")
library("ggtext")
library("cowplot")
library("ggpubr")
library("rstatix")
library("scales")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

# part A - metabolomics volcano plots
midlog <- read_csv(here::here("experimental-data", "metabolomics", "spent-media-metabolomics-cfus.csv")) %>%
  janitor::clean_names()

raw_metabolites <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(-c(index, q1_da, molecular_weight_da, formula, ionization_model, level, cpd_id, hmdb, pub_chem_cid,
                                              cas, ch_ebi, metlin, kegg_map, class_i, class_ii)) %>%
  mutate(e1_midlog = e1_midlog / (midlog %>% filter(status == "e1_midlog") %>% dplyr::select(cfu) %>% pull()),
         e2_midlog = e2_midlog / (midlog %>% filter(status == "e2_midlog") %>% dplyr::select(cfu) %>% pull()),
         e3_midlog = e3_midlog / (midlog %>% filter(status == "e3_midlog") %>% dplyr::select(cfu) %>% pull()),
         p1_midlog = p1_midlog / (midlog %>% filter(status == "p1_midlog") %>% dplyr::select(cfu) %>% pull()),
         p2_midlog = p2_midlog / (midlog %>% filter(status == "p2_midlog") %>% dplyr::select(cfu) %>% pull()),
         p3_midlog = p3_midlog / (midlog %>% filter(status == "p3_midlog") %>% dplyr::select(cfu) %>% pull()),
         pm1_midlog = pm1_midlog / (midlog %>% filter(status == "pm1_midlog") %>% dplyr::select(cfu) %>% pull()),
         pm2_midlog = pm2_midlog / (midlog %>% filter(status == "pm2_midlog") %>% dplyr::select(cfu) %>% pull()),
         pm3_midlog = pm3_midlog / (midlog %>% filter(status == "pm3_midlog") %>% dplyr::select(cfu) %>% pull()),
         s_e1_depleted = s_e1_depleted / (midlog %>% filter(status == "s_e1_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_e2_depleted = s_e2_depleted / (midlog %>% filter(status == "s_e2_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_e3_depleted = s_e3_depleted / (midlog %>% filter(status == "s_e3_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_p1_depleted = s_p1_depleted / (midlog %>% filter(status == "s_p1_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_p2_depleted = s_p2_depleted / (midlog %>% filter(status == "s_p2_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_p3_depleted = s_p3_depleted / (midlog %>% filter(status == "s_p3_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_pm1_depleted = s_pm1_depleted / (midlog %>% filter(status == "s_pm1_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_pm2_depleted = s_pm2_depleted / (midlog %>% filter(status == "s_pm2_depleted") %>% dplyr::select(cfu) %>% pull()),
         s_pm3_depleted = s_pm3_depleted / (midlog %>% filter(status == "s_pm3_depleted") %>% dplyr::select(cfu) %>% pull()),
         m_pm1_depleted = m_pm1_depleted / (midlog %>% filter(status == "m_pm1_depleted") %>% dplyr::select(cfu) %>% pull()),
         m_pm2_depleted = m_pm2_depleted / (midlog %>% filter(status == "m_pm2_depleted") %>% dplyr::select(cfu) %>% pull()),
         m_pm3_depleted = m_pm3_depleted / (midlog %>% filter(status == "m_pm3_depleted") %>% dplyr::select(cfu) %>% pull()))

class_names <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(c(compounds, class_i, class_ii))

metabolites <- read_excel(here::here("experimental-data", "metabolomics", "metware-initial-report.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(-c(formula, level, cas, index, fdr))

aggregated <- raw_metabolites %>% 
  dplyr::select(-c(qc01, qc02, qc03)) %>%
  pivot_longer(-c(compounds), names_to = "name", values_to = "intensity") %>%
  mutate(replicate = case_when(grepl("3", name) == TRUE ~ 3,
                               grepl("2", name) == TRUE ~ 2,
                               grepl("1", name) == TRUE ~ 1)) %>%
  mutate(name = case_when(name %in% c("e1_midlog", "e2_midlog", "e3_midlog") ~ "uninfected",
                          name %in% c("p1_midlog", "p2_midlog", "p3_midlog") ~ "plasmid",
                          name %in% c("pm1_midlog", "pm2_midlog", "pm3_midlog") ~ "infected",
                          name %in% c("s_e1_depleted", "s_e2_depleted", "s_e3_depleted") ~ "s uninfected",
                          name %in% c("s_p1_depleted", "s_p2_depleted", "s_p3_depleted") ~ "s plasmid",
                          name %in% c("s_pm1_depleted", "s_pm2_depleted", "s_pm3_depleted") ~ "s infected",
                          name %in% c("m_pm1_depleted", "m_pm2_depleted", "m_pm3_depleted") ~ "m infected")) %>%
  pivot_wider(names_from = replicate, names_glue = "replicate_{replicate}", values_from = intensity) %>%
  dplyr::mutate(replicate_1 = log(replicate_1),
                replicate_2 = log(replicate_2),
                replicate_3 = log(replicate_3)) %>%
  dplyr::mutate(replicate_1 = ifelse(is.infinite(replicate_1), 0, replicate_1),
                replicate_2 = ifelse(is.infinite(replicate_2), 0, replicate_2),
                replicate_3 = ifelse(is.infinite(replicate_3), 0, replicate_3))

p_vals_all_comparisons <- data.frame()
peaks <- aggregated %>% dplyr::select(compounds) %>% pull() %>% unique()
for (i in peaks){
  avg_replicates <- gather(aggregated %>% filter(compounds == i) %>% dplyr::select(-compounds), key, val, -name) %>%
    do(data.frame(tidy(pairwise.t.test(.$val, .$name, p.adjust.method="BH", 
                                       paired=FALSE, pool.sd=FALSE)))) %>%
    dplyr::mutate(compounds = as.character(i))
  p_vals_all_comparisons <- rbind(p_vals_all_comparisons, avg_replicates)
}

p_vals_all_comparisons$p.adj <- p_vals_all_comparisons$p.value

no_transformation <- raw_metabolites %>% 
  dplyr::select(-c(qc01, qc02, qc03)) %>%
  pivot_longer(-c(compounds), names_to = "name", values_to = "intensity") %>%
  dplyr::mutate(replicate = case_when(grepl("3", name) == TRUE ~ 3,
                                      grepl("2", name) == TRUE ~ 2,
                                      grepl("1", name) == TRUE ~ 1)) %>%
  dplyr::mutate(name = case_when(name %in% c("e1_midlog", "e2_midlog", "e3_midlog") ~ "uninfected",
                                 name %in% c("p1_midlog", "p2_midlog", "p3_midlog") ~ "plasmid",
                                 name %in% c("pm1_midlog", "pm2_midlog", "pm3_midlog") ~ "infected",
                                 name %in% c("s_e1_depleted", "s_e2_depleted", "s_e3_depleted") ~ "s uninfected",
                                 name %in% c("s_p1_depleted", "s_p2_depleted", "s_p3_depleted") ~ "s plasmid",
                                 name %in% c("s_pm1_depleted", "s_pm2_depleted", "s_pm3_depleted") ~ "s infected",
                                 name %in% c("m_pm1_depleted", "m_pm2_depleted", "m_pm3_depleted") ~ "m infected")) %>%
  pivot_wider(names_from = replicate, names_glue = "replicate_{replicate}", values_from = intensity) 

prep_for_comparisons <- p_vals_all_comparisons %>% 
  filter((group1 == "m infected" & group2 == "infected") | (group2 == "m infected" & group1 == "infected") | (group1 == "s infected" & group2 == "infected") |
           (group2 == "s infected" & group1 == "infected") | (group1 == "s uninfected" & group2 == "uninfected") | (group2 == "s uninfected" & group1 == "uninfected") |
           (group1 == "s plasmid" & group2 == "plasmid") | (group2 == "s plasmid" & group1 == "plasmid")) %>%
  dplyr:: mutate(denominator = ifelse(group1 %in% c("uninfected", "infected", "plasmid"), group1, group2)) %>%
  dplyr::mutate(numerator = ifelse(group1 == denominator, group2, group1)) %>%
  dplyr::select(-c(group1, group2)) %>% mutate(comparison = paste0(denominator, "_", numerator)) %>%
  pivot_longer(denominator:numerator, names_to = "type", values_to = "name")

compounds_of_interest <- prep_for_comparisons %>% dplyr::select(compounds) %>% pull() %>% unique()
partner_frame <- data.frame()

plasmid <- prep_for_comparisons %>% filter(comparison == "plasmid_s plasmid")
sinfect <- prep_for_comparisons %>% filter(comparison == "infected_s infected")
uninfected <- prep_for_comparisons %>% filter(comparison == "uninfected_s uninfected")
minfect <- prep_for_comparisons %>% filter(comparison == "infected_m infected")
for (i in compounds_of_interest){
  slice_intensities <- no_transformation %>% filter(compounds == i)
  sliced_pvals <- plasmid %>% filter(compounds == i)
  length_test <- plasmid %>% filter(compounds == i) %>% dim()
  if (length_test[1] == 0) {next}
  joined <- inner_join(slice_intensities, sliced_pvals, by = c("compounds", "name"))
  to_plasmid <- joined %>% dplyr::select(-c(type)) %>% 
    mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
    pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(`s plasmid` / plasmid)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
  partner_frame <- rbind(partner_frame, to_plasmid)
}

for (i in compounds_of_interest){
  slice_intensities <- no_transformation %>% filter(compounds == i)
  sliced_pvals <- uninfected %>% filter(compounds == i)
  length_test <- uninfected %>% filter(compounds == i) %>% dim()
  if (length_test[1] == 0) {next}
  joined <- inner_join(slice_intensities, sliced_pvals, by = c("compounds", "name"))
  to_uninfected <- joined %>% dplyr::select(-c(type)) %>% 
    mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
    pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(`s uninfected` / uninfected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
  partner_frame <- rbind(partner_frame, to_uninfected)
}

for (i in compounds_of_interest){
  slice_intensities <- no_transformation %>% filter(compounds == i)
  sliced_pvals <- sinfect %>% filter(compounds == i)
  length_test <- sinfect  %>% filter(compounds == i) %>% dim()
  if (length_test[1] == 0) {next}
  joined <- inner_join(slice_intensities, sliced_pvals, by = c("compounds", "name"))
  to_sinfect <- joined %>% dplyr::select(-c(type)) %>% 
    mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
    pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(`s infected` / infected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
  partner_frame <- rbind(partner_frame, to_sinfect)
}

for (i in compounds_of_interest){
  slice_intensities <- no_transformation %>% filter(compounds == i)
  sliced_pvals <- minfect %>% filter(compounds == i)
  length_test <- minfect  %>% filter(compounds == i) %>% dim()
  if (length_test[1] == 0) {next}
  joined <- inner_join(slice_intensities, sliced_pvals, by = c("compounds", "name"))
  to_minfect <- joined %>% dplyr::select(-c(type)) %>% 
    mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
    pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(`m infected` / infected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
  partner_frame <- rbind(partner_frame, to_minfect)
}

partA <- partner_frame %>%
  dplyr::mutate(partner = ifelse(comparison == "infected_m infected", "*M. extorquens*", "*S. enterica*")) %>%
  dplyr::mutate(parasite = case_when(comparison == "infected_m infected" | comparison == "infected_s infected" ~ "F128+<br>M13+",
                                     comparison == "uninfected_s uninfected" | comparison == "infected_s infected" ~ "uninf",
                                     comparison == "plasmid_s plasmid" ~ "F128+"),
                partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  dplyr::mutate(colors = ifelse(-log10(as.numeric(p.adj)) <= -log10(0.05), "non sig", "sig")) %>%
  dplyr::mutate(shape = ifelse(compounds == "Lactic acid", "shape", "no shape")) %>%
  dplyr::mutate(parasite = factor(parasite, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(log2fc = ifelse(log2fc < 0 & is.infinite(log2fc), -18, log2fc)) %>%
  mutate(log2fc = ifelse(log2fc > 0 & is.infinite(log2fc), 18, log2fc)) %>%
  ggplot(aes(x = as.numeric(log2fc), y = -log10(as.numeric(p.adj)), color = colors, shape = shape)) +
  geom_point(size = 2) +
  ylab("-log10(adjusted p-value)") +
  xlab("log2(fold-change)") +
  xlim(-20, 20)+
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  geom_vline(xintercept = -17, color = "black", linetype = "dashed")+
  geom_vline(xintercept = 17, color = "black", linetype = "dashed")+
  facet_grid(parasite~partner)+
  scale_color_manual(values = c("sig" = "black", "non sig" = "darkgrey"))+
  theme_bw(base_size = 16)+
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_markdown())

# part B - heat map of depleted compounds
partB <- partner_frame %>% 
  mutate(logp = -log10(as.numeric(p.adj)),
         log2fc = as.numeric(log2fc)) %>%
  filter(logp >= 1.3) %>%
  inner_join(., class_names, by = "compounds") %>% filter(log2fc < 0) %>%
  dplyr::mutate(partner = ifelse(comparison == "infected_m infected", "*M. extorquens*", "*S. enterica*")) %>%
  dplyr::mutate(parasite = case_when(comparison == "infected_m infected" | comparison == "infected_s infected" ~ "F128+<br>M13+",
                                     comparison == "uninfected_s uninfected" | comparison == "infected_s infected" ~ "uninf",
                                     comparison == "plasmid_s plasmid" ~ "F128+"),
                partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  dplyr::mutate(parasite = factor(parasite, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  group_by(comparison, class_i, parasite, partner) %>% summarize(n = n()) %>% ungroup() %>% group_by(comparison) %>% 
  mutate(total = sum(n)) %>%
  mutate(class_i = case_when(class_i == "Aldehyde,Ketones,Esters" ~ "aldehydes",
                             class_i == "Amino acid and Its metabolites" ~ "amino acids",
                             class_i == "Benzene and substituted derivatives" ~ "benzenes",
                             class_i == "FA" ~ "fatty acids",
                             class_i == "GP" ~ "lipids",
                             class_i == "Nucleotide and Its metabolites" ~ "nucleotides",
                             class_i == "Organic acid and Its derivatives" ~ "organic acids",
                             class_i == "Carbohydrates and Its metabolites" ~ "carbohydrates",
                             class_i == "Heterocyclic compounds" ~ "heterocyclics",
                             class_i == "Others" ~ "other")) %>%
  filter(class_i %in% c("aldehydes", "amino acids", "nucleotides", "organic acids")) %>%
  ggplot(aes(x = fct_reorder(class_i, n), y = (n / total) * 100)) +
  geom_bar(stat = "identity") + facet_grid(parasite~partner, scale = "free") + coord_flip() + ylab("percent of significant compounds") +
  theme_bw(base_size = 16) + theme(strip.background = element_blank(), axis.text.x = element_markdown(), axis.title.y = element_blank(),
                                   strip.text = element_markdown())

# part C - lactate depletion
pvals <- partner_frame %>% filter(compounds == "Lactic acid") %>%
  mutate(partner = ifelse(grepl("_s", comparison) == TRUE, "*S. enterica*", "none")) %>%
  mutate(partner = ifelse(grepl("_m", comparison) == TRUE, "*M. extorquens*", partner)) %>%
  mutate(baseline = case_when(comparison %in% c("infected_s infected", "infected_m infected") ~ "phage",
                              comparison == "plasmid_s plasmid" ~ "plasmid",
                              comparison == "uninfected_s uninfected" ~ "uninfected"))

partC <- raw_metabolites %>% filter(compounds == "Lactic acid") %>% 
  pivot_longer(cols = c(e1_midlog:m_pm3_depleted)) %>% 
  dplyr::select(-c(qc01, qc02, qc03)) %>% 
  mutate(baseline = case_when(name %in% c("e1_midlog", "e2_midlog", "e3_midlog") ~ "uninfected", 
                   name %in% c("p1_midlog", "p2_midlog", "p3_midlog") ~ "plasmid", 
                   name %in% c("pm1_midlog", "pm2_midlog", "pm3_midlog") ~ "phage", 
                   name %in% c("s_e1_depleted", "s_e2_depleted", "s_e3_depleted") ~ "uninfected", 
                  name %in% c("s_p1_depleted", "s_p2_depleted", "s_p3_depleted") ~ "plasmid", 
                  name %in% c("s_pm1_depleted", "s_pm2_depleted", "s_pm3_depleted") ~ "phage", 
                  name %in% c("m_pm1_depleted", "m_pm2_depleted", "m_pm3_depleted") ~ "phage")) %>%
  mutate(partner = ifelse(grepl("s_", name) == TRUE, "*S. enterica*", "none")) %>%
  mutate(partner = ifelse(grepl("m_", name) == TRUE, "*M. extorquens*", partner)) %>%
  group_by(baseline, partner) %>% mutate(rep = 1:n()) %>% ungroup() %>%
  dplyr::select(-c(name)) %>%
  pivot_wider(names_from = partner, values_from = value) %>%
  mutate(log2fc_s = log2(`*S. enterica*` / none)) %>%
  mutate(log2fc_m = log2(`*M. extorquens*` / none)) %>% dplyr::select(baseline, rep, log2fc_s, log2fc_m) %>%
  pivot_longer(cols = c(log2fc_s:log2fc_m)) %>% rename(log2fc = value) %>% 
  mutate(partner = ifelse(name == "log2fc_s", "*S. enterica*", "*M. extorquens*")) %>% dplyr::select(-c(name)) %>%
  inner_join(., pvals %>% select(p.adj, partner, baseline), by = c("partner", "baseline")) %>%
  mutate(parasite = case_when(baseline == "phage" ~ "F128+<br>M13+",
                              baseline == "plasmid" ~ "F128+",
                              baseline == "uninfected" ~ "uninf")) %>%
  dplyr::mutate(parasite = factor(parasite, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  dplyr::mutate(partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  mutate(label = case_when(p.adj > 0.05 ~ "ns",
                           p.adj < 0.05 & p.adj > 0.01 ~ "*",
                           p.adj < 0.01 ~ "**")) %>%
  ggplot(aes(x = parasite, y = log2fc, color = partner)) +
  geom_point(size = 2, stroke = 1.5) +
  ylab("log2(fold-change) in lactate") +
  geom_text(data = . %>% 
        group_by(baseline, partner, parasite) %>% 
        slice(1), aes(x = parasite, label = label, vjust = -1, color = "black", size = 5)) +
  facet_wrap(~partner, scales = "free_x") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  stat_summary(aes(fill = partner), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  scale_color_manual(values = c("*S. enterica*" = "#DDAA33", "*M. extorquens*" = "#BB5566")) +
  theme_bw(base_size = 16)+
  ylim(-12,2)+
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown(), text = element_text(color = "black"))

# part D - lactate supp
date <- "3May2024"
no_infection = TRUE

source(here::here("data-generation", "experimental", "clean-mutualism-tecan-data.R"))

stats <- all_data_corrected %>% inner_join(., plate_layout %>% select(well, supplement), by = "well") %>%
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  mutate(M_corrected_OD = ifelse(partner %in% c("M0104", "S0240 + M0104"), abs(OD - E_corrected_OD - S_corrected_OD - 0.05), 0),
         M_corrected_OD = ifelse(strain == "M0104", OD - 0.08, M_corrected_OD)) %>%
  filter(!is.na(partner)) %>%
  filter(partner != "JS002") %>%
  group_by(well) %>% filter(OD == max(OD)) %>% slice_head() %>% ungroup() %>%
  mutate(percent_M = case_when(partner == "M0104" ~ M_corrected_OD / (E_corrected_OD + M_corrected_OD),
                               partner == "S0240" ~ NA,
                               partner == "S0240 + M0104" ~ M_corrected_OD / (S_corrected_OD + E_corrected_OD + M_corrected_OD)),
         percent_S = case_when(partner == "S0240" ~ S_corrected_OD / (E_corrected_OD + S_corrected_OD),
                               partner == "M0104" ~ NA,
                               partner == "S0240 + M0104" ~ S_corrected_OD / (S_corrected_OD + E_corrected_OD + M_corrected_OD))) %>%
  select(partner, strain, supplement, percent_M, percent_S) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 + M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(supplement = ifelse(is.na(supplement), "std", "lact+"),
         supplement = factor(supplement, levels = c("std", "lact+"))) %>%
  pivot_longer(percent_M:percent_S, names_to = "species") %>%
  group_by(species, partner) %>%
  drop_na() %>%
  t_test(value ~ supplement) %>%
  add_significance("p") %>%
  add_xy_position(x = "supplement") %>%
  mutate(label = paste0(scientific(p,digits = 2), " (", p.signif, ")"))

partD <- all_data_corrected %>% inner_join(., plate_layout %>% select(well, supplement), by = "well") %>%
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  mutate(M_corrected_OD = ifelse(partner %in% c("M0104", "S0240 + M0104"), abs(OD - E_corrected_OD - S_corrected_OD - 0.05), 0),
         M_corrected_OD = ifelse(strain == "M0104", OD - 0.08, M_corrected_OD)) %>%
  filter(!is.na(partner)) %>%
  filter(partner != "JS002") %>%
  group_by(well) %>% filter(OD == max(OD)) %>% slice_head() %>% ungroup() %>%
  mutate(percent_M = case_when(partner == "M0104" ~ M_corrected_OD / (E_corrected_OD + M_corrected_OD),
                               partner == "S0240" ~ NA,
                               partner == "S0240 + M0104" ~ M_corrected_OD / (S_corrected_OD + E_corrected_OD + M_corrected_OD)),
         percent_S = case_when(partner == "S0240" ~ S_corrected_OD / (E_corrected_OD + S_corrected_OD),
                               partner == "M0104" ~ NA,
                               partner == "S0240 + M0104" ~ S_corrected_OD / (S_corrected_OD + E_corrected_OD + M_corrected_OD))) %>%
  select(partner, strain, supplement, percent_M, percent_S) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 + M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(supplement = ifelse(is.na(supplement), "std", "lact+"),
         supplement = factor(supplement, levels = c("std", "lact+"))) %>%
  pivot_longer(percent_M:percent_S, names_to = "species") %>%
  ggplot(aes(x = supplement, y = value, color = species)) + 
  stat_summary(aes(fill = species), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5) +
  facet_wrap(~partner, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 16) +
  ylab("percent partner") +
  stat_pvalue_manual(stats, label = "label", tip.length = 0.0, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = c("percent_S" = "#DDAA33", "percent_M" = "#BB5566")) +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown())

# final figure
top <- plot_grid(partA, partB, ncol = 2, rel_widths = c(0.85,1), label_size = 26, labels = c("A", "B"))
bottom <- plot_grid(partC, partD, ncol = 2, rel_widths = c(0.7, 1), label_size = 26, labels = c("C", "D"))
figure4 <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1,0.85))

png(here::here("figures", "final-figs", "imgs", "figure-4.png"), res = 300, width = 3200, height = 2400)
figure4
dev.off()





