# ATB
# figure 2
# E metabolomics

# load packages
library("tidyverse")
library("readxl")
library("ggtext")
library("cowplot")
library("ggpubr")
library("rstatix")
library("ggVennDiagram")
library("pmp")
library("scales")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

# load metabolomics data
midlog <- read_csv(here::here("experimental-data", "metabolomics", "spent-media-metabolomics-cfus.csv")) %>%
  janitor::clean_names()

raw_metabolites <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(-c(index, q1_da, molecular_weight_da, formula, ionization_model, level, cpd_id, hmdb, pub_chem_cid,
                                              cas, ch_ebi, metlin, kegg_map, class_i, class_ii)) %>%
  dplyr::mutate(e1_midlog = e1_midlog / (midlog %>% filter(status == "e1_midlog") %>% dplyr::select(cfu) %>% pull()),
                e2_midlog = e2_midlog / (midlog %>% filter(status == "e2_midlog") %>% dplyr::select(cfu) %>% pull()),
                e3_midlog = e3_midlog / (midlog %>% filter(status == "e3_midlog") %>% dplyr::select(cfu) %>% pull()),
                p1_midlog = p1_midlog / (midlog %>% filter(status == "p1_midlog") %>% dplyr::select(cfu) %>% pull()),
                p2_midlog = p2_midlog / (midlog %>% filter(status == "p2_midlog") %>% dplyr::select(cfu) %>% pull()),
                p3_midlog = p3_midlog / (midlog %>% filter(status == "p3_midlog") %>% dplyr::select(cfu) %>% pull()),
                pm1_midlog = pm1_midlog / (midlog %>% filter(status == "pm1_midlog") %>% dplyr::select(cfu) %>% pull()),
                pm2_midlog = pm2_midlog / (midlog %>% filter(status == "pm2_midlog") %>% dplyr::select(cfu) %>% pull()),
                pm3_midlog = pm3_midlog / (midlog %>% filter(status == "pm3_midlog") %>% dplyr::select(cfu) %>% pull()))

class_names <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(c(compounds, class_i, class_ii))

# part A - representative growth curves of monoculture
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
  filter(is.na(partner) & !strain %in% c("S0240", "M0104", "JS002")) %>% filter(species == "E_corrected_OD") %>%
  mutate(species = case_when(species == "E_corrected_OD" ~ "*E. coli*",
                             species == "S_corrected_OD" ~ "*S. enterica*",
                             species == "M_corrected_OD" ~ "*M. extorquens*")) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  filter(hour < 90) %>%
  ggplot(aes(x = hour, y = value, color = species, alpha = well)) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE) +
  facet_wrap(~strain) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("*E. coli*" = "#004488", "*S. enterica*" = "#DDAA33", "*M. extorquens*" = "#BB5566")) +
  ylab("monoculture OD600") +
  xlab("hour") +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_markdown())

# part B - growth rate and productivity thanks to infection
date <- "6November2023"
source(here::here("data-generation", "experimental", "clean-mutualism-tecan-data.R"))
nov6 <- all_data_corrected %>% 
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  mutate(M_corrected_OD = ifelse(partner %in% c("M0104", "S0240 M0104"), abs(OD - E_corrected_OD - S_corrected_OD - 0.05), 0),
         M_corrected_OD = ifelse(strain == "M0104", OD - 0.08, M_corrected_OD)) %>%
  pivot_longer(cols = c("E_corrected_OD", "S_corrected_OD", "M_corrected_OD"), names_to = "species") %>%
  filter(is.na(partner) & !strain %in% c("S0240", "M0104", "JS002")) %>% filter(species == "E_corrected_OD") %>%
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
  filter(is.na(partner) & !strain %in% c("S0240", "M0104", "JS002")) %>% filter(species == "E_corrected_OD") %>%
  mutate(species = case_when(species == "E_corrected_OD" ~ "*E. coli*",
                             species == "S_corrected_OD" ~ "*S. enterica*",
                             species == "M_corrected_OD" ~ "*M. extorquens*")) %>%
  mutate(date = date)

rates <- rbind(nov6, nov20) %>% 
  group_by(well, date, strain, partner) %>%
  summarize(model_fit = fit_baranyi(hour, log(value), 5),
            fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>% ungroup() %>% filter(fit_variable %in% c("growth_rate"))

stats <- rbind(nov6, nov20) %>% 
  group_by(well, date, strain, partner) %>% filter(value == max(value)) %>% slice_head() %>%
  select(well, date, strain, partner, value) %>% ungroup() %>% rename(model_fit = value) %>% mutate(fit_variable = "max OD600") %>%
  rbind(., rates) %>% mutate(fit_variable = ifelse(fit_variable == "growth_rate", "growth rate", fit_variable)) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  group_by(fit_variable) %>%
  tukey_hsd(model_fit ~ strain) %>%
  add_xy_position(x = "strain") %>%
  mutate(y.position = case_when((fit_variable == "growth rate" & group1 == "uninf" & group2 == "F128+<br>M13+") ~ 0.77,
                                (fit_variable == "growth rate" & group1 == "uninf" & group2 == "F128+") ~ 0.72,
                                 (fit_variable == "growth rate" & group1 == "F128+" & group2 == "F128+<br>M13+") ~ 0.82,
                                 (fit_variable == "max OD600" & group1 == "uninf" & group2 == "F128+") ~ 0.36,
                                 (fit_variable == "max OD600" & group1 == "uninf" & group2 == "F128+<br>M13+") ~ 0.38,
                                 (fit_variable == "max OD600" & group1 == "F128+" & group2 == "F128+<br>M13+") ~ 0.40,
                                 TRUE ~ y.position)) %>%
  mutate(label = paste0(scientific(p.adj,digits = 2), " (", p.adj.signif, ")"))

partB <- rbind(nov6, nov20) %>% 
  group_by(well, date, strain, partner) %>% filter(value == max(value)) %>% slice_head() %>%
  select(well, date, strain, partner, value) %>% ungroup() %>% rename(model_fit = value) %>% mutate(fit_variable = "max OD600") %>%
  rbind(., rates) %>% mutate(fit_variable = ifelse(fit_variable == "growth_rate", "growth rate", fit_variable)) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  ggplot(aes(x = strain, y = model_fit, color = strain)) +
  geom_point(size = 2, stroke = 1.5) +
  stat_pvalue_manual(stats, label = "label", tip.length = 0.0, size = 3) +
  stat_summary(aes(fill = strain), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  facet_wrap(~fit_variable, scales = "free") +
  ylab("monoculture estimate") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = c("uninf" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16) + theme(axis.text = element_markdown(), axis.title.x = element_blank(),
                                   legend.position = "none", strip.background = element_blank())

# part C - PCA
effm <- raw_metabolites %>% select(compounds, c(e1_midlog:pm3_midlog, qc01:qc03)) 

for_pca <- as.matrix(effm[,-c(1)])

qc <- which(grepl("qc", colnames(for_pca)) == TRUE)
indices <- which(grepl("qc", colnames(for_pca)) != TRUE)
name_vector <- c(rep("E", 3), rep("P", 3), rep("PM", 3), rep("qc", 3))

for_pca <- for_pca[rowSums(is.na(for_pca)) != ncol(for_pca), ]

subset_data <- pqn_normalisation(for_pca, classes=name_vector[c(indices,qc)], qc_label="qc")
subset_data <- mv_imputation(subset_data, method="KNN", k=5, rowmax=0.5,
                             colmax=0.5, maxp=NULL, check_df=FALSE)
subset_data <- glog_transformation(subset_data, classes=name_vector[c(indices,qc)], qc_label="qc")
clean_QCs <- subset_data[, -grep("qc", colnames(subset_data))]
clean_QCs <- prcomp(t(clean_QCs), center=TRUE, scale=FALSE)

get_variance <- round(((clean_QCs$sdev^2)/sum(clean_QCs$sdev^2)*100)[1:2],2)

data <- data.frame(PC1=clean_QCs$x[, 1], PC2=clean_QCs$x[, 2],class=factor(name_vector[c(indices)], 
                                                                           levels = unique(name_vector[c(indices)][length(name_vector[c(indices)]):1])))

partC <- data %>% mutate(class = case_when(class == "E" ~ "uninf",
                                           class == "P" ~ "F128+",
                                           class == "PM" ~ "F128+ M13+")) %>%
  mutate(class = factor(class, levels = c("uninf", "F128+", "F128+ M13+"))) %>%
                           ggplot(aes(x = PC1, y = PC2, color = class)) + geom_point(size = 4) +
  xlab(paste0("PC1 (", get_variance[1], "%)")) +
  ylab(paste0("PC2 (", get_variance[2], "%)")) +
  scale_color_manual(values = c("uninf" = "black", "F128+" = "#117733", "F128+ M13+" = "#88CCEE")) +
  theme_bw(base_size = 16) + theme(axis.text.x = element_markdown(), legend.position = "none")

legendC <- get_plot_component(data %>% mutate(class = case_when(class == "E" ~ "uninf",
                                                                class == "P" ~ "F128+",
                                                                class == "PM" ~ "F128+ M13+")) %>%
                                mutate(class = factor(class, levels = c("uninf", "F128+", "F128+ M13+"))) %>%
                                ggplot(aes(x = PC1, y = PC2, color = class)) + geom_point(size = 4) +
                                xlab(paste0("PC1 (", get_variance[1], "%)")) +
                                ylab(paste0("PC2 (", get_variance[2], "%)")) +
                                scale_color_manual(values = c("uninf" = "black", "F128+" = "#117733", "F128+ M13+" = "#88CCEE")) +
                                theme_bw(base_size = 16) + theme(axis.text.x = element_markdown(), legend.position = "bottom", legend.title = element_blank()),
                                'guide-box-bottom', return_all = TRUE)

# part D - venn diagram of detected compounds
e <- raw_metabolites %>% select(compounds, e1_midlog, e2_midlog, e3_midlog) %>% 
  pivot_longer(cols = c(e1_midlog:e3_midlog)) %>% filter(value > 0) %>% pull(compounds) %>% unique()

f <- raw_metabolites %>% select(compounds, p1_midlog, p2_midlog, p3_midlog) %>% 
  pivot_longer(cols = c(p1_midlog:p3_midlog)) %>% filter(value > 0) %>% pull(compounds) %>% unique()

fm <- raw_metabolites %>% select(compounds, pm1_midlog, pm2_midlog, pm3_midlog) %>% 
  pivot_longer(cols = c(pm1_midlog:pm3_midlog)) %>% filter(value > 0) %>% pull(compounds) %>% unique()

metabolites <- list(uninf = e, F128 = f, F128_M13 = fm)

partD <- ggVennDiagram(metabolites, label_alpha=0, category.names = c("uninf","F128+","F128+ M13+")) + 
  theme_void(base_size = 16) + theme(legend.position = "none", axis.title = element_blank()) + 
  scale_fill_gradient(low="white",high = "white")

# part E - E, F, FM volcano plots
aggregated <- raw_metabolites %>% 
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
    do(data.frame(tidy(pairwise.t.test(.$val, .$name, p.adjust.method="none", 
                                       paired=FALSE, pool.sd=FALSE)))) %>%
    dplyr::mutate(compounds = as.character(i))
  p_vals_all_comparisons <- rbind(p_vals_all_comparisons, avg_replicates)
}

p_vals_all_comparisons$p.adj <- p.adjust(p_vals_all_comparisons$p.value, method = "BH")

prep_for_comparisons <- p_vals_all_comparisons %>% filter(group1 %in% c("infected", "plasmid", "uninfected") & group2 %in% c("infected", "plasmid", "uninfected")) %>%
  mutate(denominator = ifelse(group1 == "uninfected", group1, group2)) %>% filter(denominator == "uninfected") %>% dplyr::mutate(numerator = group2) %>%
  dplyr::select(-c(group1, group2)) %>% dplyr::mutate(comparison = paste0(denominator, "_", numerator)) %>%
  pivot_longer(denominator:numerator, names_to = "type", values_to = "name")

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


compounds_of_interest <- prep_for_comparisons %>% dplyr::select(compounds) %>% pull() %>% unique()
metabolism_ecoli_frame <- data.frame()
for (i in compounds_of_interest){
  slice_intensities <- no_transformation %>% filter(compounds == i)
  sliced_pvals <- prep_for_comparisons %>% filter(compounds == i)
  joined <- inner_join(slice_intensities, sliced_pvals, by = c("compounds", "name"))
  if (("uninfected_infected" %in% (joined %>% dplyr::select(comparison) %>% pull() %>% unique())) == TRUE & ("uninfected_plasmid" %in% (joined %>% dplyr::select(comparison) %>% pull() %>% unique())) == FALSE){
    to_infected <- joined %>% filter(comparison == "uninfected_infected") %>% dplyr::select(-c(type)) %>% 
      dplyr::mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
      pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(infected / uninfected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
    new_frame <- to_infected
  } 
  if (("uninfected_plasmid" %in% (joined %>% dplyr::select(comparison) %>% pull() %>% unique())) == TRUE & ("uninfected_infected" %in% (joined %>% dplyr::select(comparison) %>% pull() %>% unique())) == FALSE){
    to_plasmid <- joined %>% filter(comparison == "uninfected_plasmid") %>% dplyr::select(-c(type)) %>% 
      dplyr::mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
      pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(plasmid / uninfected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
    new_frame <- to_plasmid
  }
  if (("uninfected_plasmid" %in% (joined %>% dplyr::select(comparison) %>% pull() %>% unique())) == TRUE & ("uninfected_infected" %in% (joined %>% dplyr::select(comparison) %>% pull() %>% unique())) == TRUE){
    to_plasmid <- joined %>% filter(comparison == "uninfected_plasmid") %>% dplyr::select(-c(type)) %>% 
      dplyr::mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
      pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(plasmid / uninfected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
    to_infected <- joined %>% filter(comparison == "uninfected_infected") %>% dplyr::select(-c(type)) %>% 
      dplyr::mutate(average = (replicate_1 + replicate_2 + replicate_3) / 3) %>% dplyr::select(compounds, name, p.adj, average, comparison) %>% 
      pivot_wider(names_from = "name", values_from = "average") %>% dplyr::mutate(log2fc = log2(infected / uninfected)) %>% dplyr::select(compounds, p.adj, log2fc, comparison)
    new_frame <- rbind(to_infected, to_plasmid)
  }
  metabolism_ecoli_frame <- rbind(metabolism_ecoli_frame, new_frame)
}

partE <- metabolism_ecoli_frame %>%
  mutate(comparison = ifelse(comparison == "uninfected_plasmid", "plasmid", "virus with plasmid")) %>%
  mutate(comparison = ifelse(comparison == "plasmid", "F128+", "F128+<br>M13+"),
         comparison = factor(comparison, levels = c("F128+", "F128+<br>M13+"))) %>%
  mutate(colors = ifelse(-log10(as.numeric(p.adj)) <= -log10(0.05), "non sig", "sig")) %>%
  mutate(shape = ifelse(compounds == "Lactic acid", "shape", "no shape")) %>%
  mutate(log2fc = ifelse(log2fc < 0 & is.infinite(log2fc), -18, log2fc)) %>%
  mutate(log2fc = ifelse(log2fc > 0 & is.infinite(log2fc), 18, log2fc)) %>%
  ggplot(aes(x = as.numeric(log2fc), y = -log10(as.numeric(p.adj)), color = colors, shape = shape)) +
  geom_point(size = 2) +
  ylab("-log10(adjusted p-value)") +
  xlab("log2(fold-change)") +
  xlim(-20,20)+
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  geom_vline(xintercept = -17, color = "black", linetype = "dashed")+
  geom_vline(xintercept = 17, color = "black", linetype = "dashed")+
  facet_wrap(~comparison, ncol = 2)+
  scale_color_manual(values = c("sig" = "black", "non sig" = "darkgrey"))+
  theme_bw(base_size = 16)+
  theme(axis.text = element_markdown(), strip.text = element_markdown(),
        legend.position = "none", strip.background = element_blank())

# part F - number of compound significantly enriched by class
partF <- metabolism_ecoli_frame %>% 
  mutate(logp = -log10(as.numeric(p.adj)),
         log2fc = as.numeric(log2fc)) %>%
  filter(logp >= 1.3) %>%
  inner_join(., class_names, by = "compounds") %>% filter(log2fc > 0) %>%
  mutate(comparison = ifelse(comparison == "uninfected_plasmid", "plasmid", "virus with plasmid")) %>%
  mutate(comparison = ifelse(comparison == "plasmid", "F128+", "F128+<br>M13+"),
         comparison = factor(comparison, levels = c("F128+", "F128+<br>M13+"))) %>%
  mutate(compounds = tolower(compounds)) %>%
  group_by(comparison, class_i) %>% summarize(n = n()) %>%
  ungroup() %>% group_by(comparison) %>% mutate(total = sum(n)) %>%
  mutate(class_i = case_when(class_i == "Aldehyde,Ketones,Esters" ~ "aldehydes",
                                                    class_i == "Amino acid and Its metabolites" ~ "amino acids",
                                                    class_i == "Benzene and substituted derivatives" ~ "benzenes",
                                                    class_i == "FA" ~ "fatty acids",
                                                    class_i == "GP" ~ "lipids",
                                                    class_i == "Nucleotide and Its metabolites" ~ "nucleotides",
                                                    class_i == "Organic acid and Its derivatives" ~ "organic acids",
                                                    class_i == "Carbohydrates and Its metabolites" ~ "carbohydrates",
                                                    class_i == "Heterocyclic compounds" ~ "heterocyclics")) %>%
  filter(class_i %in% c("aldehydes", "amino acids", "nucleotides", "organic acids")) %>%
  ggplot(aes(x = fct_reorder(class_i, n), y = (n / total) * 100)) +
  geom_bar(stat = "identity") + facet_wrap(~comparison) + coord_flip() + ylab("percent of significant compounds") +
  theme_bw(base_size = 16) + theme(strip.background = element_blank(), axis.text.x = element_markdown(), axis.title.y = element_blank(),
                                   strip.text = element_markdown())
  
# part G - HPLC
hplc <- read_csv(here::here("experimental-data", "hplc", "24April2024", "aggregated_signals_24April2024.csv")) %>% dplyr::select(-c("...1")) %>%
  pivot_longer(cols = c(-time), names_to = "sampleID") %>%
  dplyr::mutate(class = case_when(grepl("std", sampleID) == TRUE ~ "standard",
                                  grepl("blank", sampleID) == TRUE ~ "blank",
                                  sampleID %in% c("e1", "e2", "e3", "f1", "f2", "f4", "fm1", "fm3", "fm4") ~ "sample",
                                  TRUE ~ sampleID))

partG <- hplc %>% filter(class == "sample") %>% 
  dplyr::mutate(strain = case_when(sampleID %in% c("e1", "e2", "e3") ~ "E0224",
                                   sampleID %in% c("f1", "f2", "f4") ~ "E0224 F+",
                                   sampleID %in% c("fm1", "fm3", "fm4") ~ "E0224 F+ M13+")) %>%
  dplyr::mutate(strain = case_when(strain == "E0224" ~ "uninf",
                                   strain == "E0224 F+" ~ "F128+",
                                   strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
                strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  filter(time > 15 & time < 30) %>% 
  group_by(strain) %>% mutate(rep = 1:n()) %>% ungroup() %>%
  ggplot(aes(x = time, y = value, color = sampleID)) + geom_line() + facet_wrap(~strain) +
  geom_vline(xintercept = 19.9, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 23.57, color = "darkgrey", linetype = "dashed") +
  theme_bw(base_size = 16) +
  ylab("mAU") +
  annotate("text", x = 19, y = 30, label = "lactate", angle = 90)+
  annotate("text", x = 22.67, y = 30, label = "acetate", angle = 90)+
  xlab("retention time (minutes)") +
  theme(axis.text = element_markdown(), strip.text = element_markdown(),
        legend.position = "none", strip.background = element_blank())

# final figure
top <- plot_grid(partA, partB, labels = c("A", "B"), label_size = 26, ncol = 2)
middle <- plot_grid(plot_grid(partC, legendC, ncol = 1, rel_heights = c(1,0.1)), partD, partE, ncol = 3, labels = c("C", "D", "E"), label_size = 26, rel_widths = c(0.7,1,1))
bottom <- plot_grid(partF, partG, ncol = 2, labels = c("F", "G"), label_size = 26, rel_widths = c(0.7,1))
figure2 <- plot_grid(top, middle, bottom, ncol = 1, rel_heights = c(0.9,0.9,1))

png(here::here("figures", "final-figs", "imgs", "figure-2.png"), res = 300, width = 4250, height = 3250)
figure2
dev.off()

