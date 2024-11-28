# ATB
# figure 4
# partner metabolomics

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
  ggplot(aes(x = fct_reorder(class_i, n), y = (n / total) * 100)) +
  geom_bar(stat = "identity") + facet_grid(parasite~partner, scale = "free") + coord_flip() + ylab("percent of significant compounds") +
  theme_bw(base_size = 16) + theme(strip.background = element_blank(), axis.text.x = element_markdown(), axis.title.y = element_blank(),
                                   strip.text = element_markdown())

# part C - lactate depletion
no_data <- data.frame(compounds = c("Lactic acid"), log2fc = c(0.005, 0.005), partner = c("M0"), parasite = c("uninf", "F128+"))

partC <- partner_frame %>%
  dplyr::mutate(partner = ifelse(comparison == "infected_m infected", "*M. extorquens*", "*S. enterica*")) %>%
  dplyr::mutate(parasite = case_when(comparison == "infected_m infected" | comparison == "infected_s infected" ~ "F128+<br>M13+",
                                     comparison == "uninfected_s uninfected" | comparison == "infected_s infected" ~ "uninf",
                                     comparison == "plasmid_s plasmid" ~ "F128+"),
                partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  dplyr::mutate(colors = ifelse(-log10(as.numeric(p.adj)) <= -log10(0.05), "non sig", "sig")) %>%
  dplyr::mutate(shape = ifelse(compounds == "Lactic acid", "shape", "no shape")) %>%
  dplyr::mutate(parasite = factor(parasite, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  select(compounds, log2fc, partner, parasite) %>%
  filter(compounds == "Lactic acid") %>%
  rbind(., no_data) %>%
  ggplot(aes(x = parasite, y = as.numeric(log2fc), fill = partner)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  ylab("lactate depletion") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("*S. enterica*" = "#DDAA33", "*M. extorquens*" = "#BB5566")) +
  theme_bw(base_size = 16)+
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown())

# part D - galK/pss bioassay
date <- "18June2024"

OD_path <- generate_paths("od", date, "csv")
JFP_path <- generate_paths("jfp", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
RFP_path <- generate_paths("rfp", date, "csv")
cfus <- generate_paths("cfu", date, "xlsx")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, OD_path), "OD") %>% filter(!is.na(cycle))
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, CFP_path), "CFP") %>% filter(!is.na(cycle))
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, YFP_path), "YFP") %>% filter(!is.na(cycle))
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, RFP_path), "RFP") %>% filter(!is.na(cycle))
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, "plate_layout.csv")) %>%
  janitor::clean_names()
CFU <- load_pfu_data(here::here("experimental-data", "tecan-data","auxotroph-bioassays", date, cfus)) %>%
  filter(!is.na(cfu))

all_data_plate1 <- OD %>%
  inner_join(., CFP %>% select(cycle, well, CFP), by = c("well", "cycle")) %>%
  inner_join(., RFP %>% select(cycle, well, RFP), by = c("well", "cycle")) %>%
  inner_join(., YFP %>% select(cycle, well, YFP), by = c("well", "cycle")) %>%
  inner_join(., plate_layout %>% filter(plate_number == 1), by = c("well")) %>%
  rename(partner = media) %>%
  filter(strain != "none") %>%
  mutate(partner = ifelse(partner == "none", NA, partner))

partD <- CFU %>% select(well, cfu, plate, plate_number) %>%
  inner_join(., plate_layout, by = c("well", "plate_number")) %>%
  rename(partner = media) %>%
  filter(plate_number == 1) %>%
  pivot_wider(names_from = plate, values_from = cfu) %>%
  mutate(total_biomass = E + all_S,
         total_S_ratio = all_S / (E + all_S),
         percent_of_S_as_KO = KO / all_S) %>%
  pivot_longer(cols = total_biomass:percent_of_S_as_KO, names_to = "comparison") %>%
  filter(comparison == "percent_of_S_as_KO") %>%
  filter(!is.na(partner)) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(partner == "pps WT" ~ "∆*pps*",
                             partner == "galK WT" ~ "∆*galK*"),
         partner = factor(partner, levels = c("∆*pps*", "∆*galK*"))) %>%
  ggplot(aes(x = strain, y = value, color = strain)) + 
  stat_summary(aes(fill = strain), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5) +
  facet_wrap(~partner, ncol = 2, scales = "free_y") +
  theme_bw(base_size = 16) +
  ylab("mutant relative fitness") +
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "F128+<br>M13+" = "#0F2080")) +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown())

# part E - lactate supp
date <- "3May2024"
no_infection = TRUE

source(here::here("data-generation", "experimental", "clean-mutualism-tecan-data.R"))

partE <- all_data_corrected %>% inner_join(., plate_layout %>% select(well, supplement), by = "well") %>%
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
  ylab("% partner") +
  scale_color_manual(values = c("percent_S" = "#DDAA33", "percent_M" = "#BB5566")) +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown())

# final figure
top <- plot_grid(partA, partB, ncol = 2, rel_widths = c(0.7,1), label_size = 26, labels = c("A", "B"))
bottom <- plot_grid(partC, partD, partE, ncol = 3, rel_widths = c(0.5,0.85,1), label_size = 26, labels = c("C", "D", "E"))
figure4 <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1,0.7))

png(here::here("figures", "final-figs", "imgs", "figure-4.png"), res = 300, width = 4000, height = 2500)
figure4
dev.off()





