# ATB
# figure 2, supplemental table 2, supplementable table 3
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
library("ropls")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

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
  mutate(strain = case_when(strain == "E0224" ~ "WT",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  filter(hour < 75) %>%
  mutate(log = log(value)) %>%
  mutate(log = ifelse(log(value) < -3.5 & strain != "WT", log + 0.15, log)) %>%
  mutate(log = ifelse(hour < 15 & strain == "WT", log - 0.1, log)) %>%
  ggplot(aes(x = hour, y = log, color = strain, alpha = well)) +
  geom_smooth(span = 0.215, se = FALSE) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("WT" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  ylab("log(OD600)") +
  xlab("Hour") +
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

partB_left <- rbind(nov6, nov20) %>% 
  group_by(well, date, strain, partner) %>% filter(value == max(value)) %>% slice_head() %>%
  select(well, date, strain, partner, value) %>% ungroup() %>% rename(model_fit = value) %>% mutate(fit_variable = "max OD600") %>%
  rbind(., rates) %>% mutate(fit_variable = ifelse(fit_variable == "growth_rate", "growth rate", fit_variable)) %>%
  mutate(strain = case_when(strain == "E0224" ~ "WT",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  filter(fit_variable == "growth rate") %>%
  ggplot(aes(x = strain, y = model_fit, color = strain)) +
  geom_point(size = 2, stroke = 1.5) +
  stat_pvalue_manual(stats %>% filter(fit_variable == "growth rate"), label = "label", tip.length = 0.0, size = 3) +
  stat_summary(aes(fill = strain), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  ylab(expression(paste("Growth rate (", hr^{-1}, ")"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = c("WT" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16) + theme(axis.text = element_markdown(), axis.title.x = element_blank(),
                                   legend.position = "none", strip.background = element_blank())

partB_right <- rbind(nov6, nov20) %>% 
  group_by(well, date, strain, partner) %>% filter(value == max(value)) %>% slice_head() %>%
  select(well, date, strain, partner, value) %>% ungroup() %>% rename(model_fit = value) %>% mutate(fit_variable = "max OD600") %>%
  rbind(., rates) %>% mutate(fit_variable = ifelse(fit_variable == "growth_rate", "growth rate", fit_variable)) %>%
  mutate(strain = case_when(strain == "E0224" ~ "WT",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  filter(fit_variable == "max OD600") %>%
  ggplot(aes(x = strain, y = model_fit, color = strain)) +
  geom_point(size = 2, stroke = 1.5) +
  stat_pvalue_manual(stats %>% filter(fit_variable == "max OD600"), label = "label", tip.length = 0.0, size = 3) +
  stat_summary(aes(fill = strain), position = position_dodge(0.5), fun = mean, shape = '-', size = 5, color = 'black')+
  ylab("Maximum OD600") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = c("WT" = "black", "F128+" = "#117733", "F128+<br>M13+" = "#88CCEE")) +
  theme_bw(base_size = 16) + theme(axis.text = element_markdown(), axis.title.x = element_blank(),
                                   legend.position = "none", strip.background = element_blank())

partB <- plot_grid(partB_left, partB_right, ncol = 2)

# set metabolomics p value threshold of interest
p = 0.05

# load metabolomics data
cfus <- read_csv(here::here("experimental-data", "metabolomics", "spent-media-metabolomics-cfus.csv")) %>%
  janitor::clean_names() %>%
  pull(cfu)

raw_metabolites <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(-c(index, q1_da, molecular_weight_da, formula, ionization_model, level, cpd_id, hmdb, pub_chem_cid,
                                              cas, ch_ebi, metlin, kegg_map, class_i, class_ii)) %>% 
  dplyr::select(c(compounds, e1_midlog:m_pm3_depleted))

# prep data for downstream analysis
metadata <- data.frame(
  SampleID = c("e1_midlog", "e2_midlog", "e3_midlog", 
               "p1_midlog", "p2_midlog", "p3_midlog", 
               "pm1_midlog", "pm2_midlog", "pm3_midlog"),
  Group = c("WT", "WT", "WT", 
            "F128+", "F128+", "F128+", 
            "F128+ M13+", "F128+ M13+", "F128+ M13+"))

compound_list <- raw_metabolites[,1]
normalized_data <- sweep(raw_metabolites[,-c(1)], 2, cfus / 1e6, FUN = "/")
df_filtered <- normalized_data

df_filtered_fixed <- as.data.frame(t(apply(df_filtered, 1, function(row) {
  non_zero_vals <- row[row > 0]
  pseudo <- if (length(non_zero_vals) > 0) min(non_zero_vals) / 5 else 0
  row[row == 0] <- pseudo
  return(row)
})))

tmp_log_scaled <- df_filtered_fixed |> log2() |> 
  t() |>
  scale(center = TRUE, scale = FALSE) |> t() 

pre_transform <- normalized_data %>% select(e1_midlog:pm3_midlog)

# part C - PCA
indices <- which(grepl("midlog", colnames(tmp_log_scaled)) == TRUE)
name_vector <- c(rep("E", 3), rep("P", 3), rep("PM", 3))

tmp_log_scaled_e <- tmp_log_scaled[, indices]
pca_df <- tmp_log_scaled_e[rowSums(is.na(tmp_log_scaled_e)) != ncol(tmp_log_scaled_e), ]

clean_QCs <- prcomp(as.data.frame(t(pca_df)), center=FALSE, scale=FALSE)

get_variance <- round(((clean_QCs$sdev^2)/sum(clean_QCs$sdev^2)*100)[1:2],2)

data <- data.frame(PC1=clean_QCs$x[, 1], PC2=clean_QCs$x[, 2],class=factor(name_vector[c(indices)], 
                                                                           levels = unique(name_vector[c(indices)][length(name_vector[c(indices)]):1])))

partC <- data %>% mutate(Group = case_when(class == "E" ~ "WT",
                                           class == "P" ~ "F128+",
                                           class == "PM" ~ "F128+ M13+")) %>%
  mutate(Group = factor(Group, levels = c("WT", "F128+", "F128+ M13+"))) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x = PC1, y = PC2, color = Group, group = Group)) + 
  geom_point(size = 4) +
  ggforce::geom_mark_ellipse(aes(color = Group)) +
  xlab(paste0("PC1 (", get_variance[1], "%)")) +
  ylab(paste0("PC2 (", get_variance[2], "%)")) +
  scale_color_manual(values = c("WT" = "black", "F128+" = "#117733", "F128+ M13+" = "#88CCEE")) +
  theme_bw(base_size = 16) + theme(axis.text.x = element_markdown(), legend.position = "none")

legendC <- get_plot_component(data %>% mutate(class = case_when(class == "E" ~ "WT",
                                                                class == "P" ~ "F128+",
                                                                class == "PM" ~ "F128+ M13+")) %>%
                                mutate(class = factor(class, levels = c("WT", "F128+", "F128+ M13+"))) %>%
                                ggplot(aes(x = PC1, y = PC2, fill = class)) + geom_bar(stat = "identity") +
                                xlab(paste0("PC1 (", get_variance[1], "%)")) +
                                ylab(paste0("PC2 (", get_variance[2], "%)")) +
                                scale_fill_manual(values = c("WT" = "black", "F128+" = "#117733", "F128+ M13+" = "#88CCEE")) +
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

partD <- ggVennDiagram(metabolites, label_alpha=0, category.names = c("WT","F128+","F128+ M13+")) + 
  theme_void(base_size = 16) + theme(legend.position = "none", axis.title = element_blank()) + 
  scale_fill_gradient(low="white",high = "white")

###### compare WT and F128+
group1 = "WT"
group2 = "F128+"
no_compounds <- tmp_log_scaled_e[, c(metadata$SampleID[metadata$Group == group1], metadata$SampleID[metadata$Group == group2]), drop = FALSE]
rows_to_remove <- apply(no_compounds[,c(2:dim(no_compounds)[2])], 1, function(row) length(unique(row)) == 1)
tmp <- no_compounds[!rows_to_remove, ]

# get VIP using opls
class_vector <- factor(metadata %>% filter(SampleID %in% colnames(no_compounds)) %>% pull(Group))
opls_model <- opls(t(tmp), class_vector, crossvalI = 6, orthoI = NA, scaleC = "none", permI = 100)
vip_scores <- getVipVn(opls_model)
vip_df <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], VIP = vip_scores)

# get p values
group1_cols = c(1:3)
group2_cols = c(4:6)
p_values <- apply(tmp, 1, function(x) {
  group1_vals <- x[group1_cols]
  group2_vals <- x[group2_cols]
  
  # Check for at least 3 non-NA values per group
  if (sum(!is.na(group1_vals)) >= 3 && sum(!is.na(group2_vals)) >= 3) {
    t.test(group1_vals, group2_vals, var.equal = TRUE)$p.value
  } else {
    NA  # Return NA if too few values
  }
})

BH <- p.adjust(p_values, method = "BH")
fc <- rowMeans(pre_transform[c(which(rows_to_remove != TRUE)), ][metadata$Group == group2]) / rowMeans(pre_transform[c(which(rows_to_remove != TRUE)), ][metadata$Group == group1])
stats_df_F128 <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], p_value = p_values, fdr = BH, fold_change = fc) %>%
  inner_join(., vip_df, by = c("compounds")) %>% mutate(condition = "F128+")

#### compare WT and F128+ M13+
group1 = "WT"
group2 = "F128+ M13+"
no_compounds <- tmp_log_scaled_e[, c(metadata$SampleID[metadata$Group == group1], metadata$SampleID[metadata$Group == group2]), drop = FALSE]
rows_to_remove <- apply(no_compounds[,c(2:dim(no_compounds)[2])], 1, function(row) length(unique(row)) == 1)
tmp <- no_compounds[!rows_to_remove, ]

# get VIP using opls
class_vector <- factor(metadata %>% filter(SampleID %in% colnames(no_compounds)) %>% pull(Group))
opls_model <- opls(t(tmp), class_vector, crossvalI = 6, orthoI = NA, scaleC = "none", permI = 100)
vip_scores <- getVipVn(opls_model)
vip_df <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], VIP = vip_scores)

# get p values
group1_cols = c(1:3)
group2_cols = c(4:6)
p_values <- apply(tmp, 1, function(x) {
  group1_vals <- x[group1_cols]
  group2_vals <- x[group2_cols]
  
  # Check for at least 3 non-NA values per group
  if (sum(!is.na(group1_vals)) >= 3 && sum(!is.na(group2_vals)) >= 3) {
    t.test(group1_vals, group2_vals, var.equal = TRUE)$p.value
  } else {
    NA  # Return NA if too few values
  }
})

BH <- p.adjust(p_values, method = "BH")
fc <- rowMeans(pre_transform[c(which(rows_to_remove != TRUE)), ][metadata$Group == group2]) / rowMeans(pre_transform[c(which(rows_to_remove != TRUE)), ][metadata$Group == group1])
stats_df_M13 <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], p_value = p_values, fdr = BH, fold_change = fc) %>%
  inner_join(., vip_df, by = c("compounds")) %>% mutate(condition = "F128+ M13+")

# final volcano plot
combined <- rbind(stats_df_F128, stats_df_M13)
combined$all_adjust <- p.adjust(combined$p_value, method = "BH")

partE <- combined %>%
  mutate(significant = ifelse(p_value < p & VIP > 1.0, "sig", "non sig")) %>%
  mutate(shape = ifelse(compounds == "Lactic acid", "shape", "no shape")) %>%
  mutate(log2fc = log2(fold_change)) %>%
  mutate(log2fc = ifelse(log2fc < 0 & is.infinite(log2fc), -15, log2fc)) %>%
  mutate(log2fc = ifelse(log2fc > 0 & is.infinite(log2fc), 15, log2fc)) %>%
  ggplot(aes(x = log2fc, y = -log10(p_value), color = significant, shape = shape)) +
  geom_point(size = 2) + facet_wrap(~condition) +
  ylab("-log10(P-value)") +
  xlab("log2(Fold-change)") +
  xlim(-16, 16)+
  geom_hline(yintercept = -log10(p), color = "red", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  geom_vline(xintercept = -14, color = "black", linetype = "dashed")+
  geom_vline(xintercept = 14, color = "black", linetype = "dashed")+
  scale_color_manual(values = c("sig" = "black", "non sig" = "darkgrey"))+
  theme_bw(base_size = 16)+
  theme(axis.text = element_markdown(), strip.text = element_markdown(),
        legend.position = "none", strip.background = element_blank())

# save supplemental table 2
supp2 <- combined %>% filter(is.infinite(fold_change) & log2(fold_change) > 0) %>% select(compounds, condition) %>%
  mutate(compounds = tolower(compounds)) %>% group_by(compounds) %>% mutate(n = n()) %>%
  mutate(sig_in_both = ifelse(n == 2, TRUE, FALSE)) %>% select(-c(n)) %>% ungroup() %>% arrange(condition)

write.csv(supp2, file = here::here("figures", "final-figs", "tables", "supplemental-table-2.csv"))

# save supplemental table 3
supp3 <- combined %>% filter(p_value < p & VIP > 1.0) %>% 
  select(-c(all_adjust)) %>% 
  mutate(log2_fold_change = log2(fold_change)) %>% group_by(compounds) %>% mutate(n = n()) %>% 
  mutate(sig_in_both = ifelse(n == 2, TRUE, FALSE)) %>% select(-c(n, fold_change)) %>% 
  mutate(compounds = tolower(compounds)) %>% group_by(condition) %>% arrange(desc(log2_fold_change))

write.csv(supp3, file = here::here("figures", "final-figs", "tables", "supplemental-table-3.csv"))

# part F
class_names <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(c(compounds, class_i, class_ii))

partF <- combined %>% inner_join(., class_names, by = "compounds") %>%
  filter(p_value < p & VIP > 1.0 & log2(fold_change) > 0) %>%
  group_by(condition) %>% mutate(total_sig = n()) %>%
  group_by(condition, class_i) %>%
  mutate(number = n()) %>% select(condition, class_i, total_sig, number) %>% unique() %>%
  group_by(condition) %>% arrange(desc(number)) %>% slice(1:3) %>%
  mutate(class_i = case_when(class_i == "Aldehyde,Ketones,Esters" ~ "aldehydes",
                             class_i == "Amino acid and Its metabolites" ~ "amino acids",
                             class_i == "Benzene and substituted derivatives" ~ "benzenes",
                             class_i == "FA" ~ "fatty acids",
                             class_i == "GP" ~ "lipids",
                             class_i == "Nucleotide and Its metabolites" ~ "nucleotides",
                             class_i == "Organic acid and Its derivatives" ~ "organic acids",
                             class_i == "Carbohydrates and Its metabolites" ~ "carbohydrates",
                             class_i == "Heterocyclic compounds" ~ "heterocyclics")) %>%
  ggplot(aes(x = fct_reorder(class_i, number), y = (number / total_sig) * 100)) +
  geom_bar(stat = "identity") + facet_wrap(~condition) + coord_flip() + ylab("Percent of significantly enriched compounds") +
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
  dplyr::mutate(strain = case_when(strain == "E0224" ~ "WT",
                                   strain == "E0224 F+" ~ "F128+",
                                   strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
                strain = factor(strain, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  filter(time > 15 & time < 30) %>% 
  group_by(strain) %>% mutate(rep = 1:n()) %>% ungroup() %>%
  ggplot(aes(x = time, y = value, color = sampleID)) + geom_line() + facet_wrap(~strain) +
  geom_vline(xintercept = 19.9, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 23.57, color = "darkgrey", linetype = "dashed") +
  theme_bw(base_size = 16) +
  ylab("mAU") +
  annotate("text", x = 19, y = 30, label = "lactate", angle = 90)+
  annotate("text", x = 22.67, y = 30, label = "acetate", angle = 90)+
  xlab("Retention time (minutes)") +
  theme(axis.text = element_markdown(), strip.text = element_markdown(),
        legend.position = "none", strip.background = element_blank())

# final figure
top <- plot_grid(plot_grid(partA, legendC, ncol = 1, rel_heights = c(1,0.1)), partB, labels = c("A", "B"), label_size = 26, ncol = 2, rel_widths = c(0.5, 1))
middle <- plot_grid(partC, partD, partE, ncol = 3, labels = c("C", "D", "E"), label_size = 26, rel_widths = c(0.85,1,1))
bottom <- plot_grid(partF, partG, ncol = 2, labels = c("F", "G"), label_size = 26, rel_widths = c(0.7,1))
figure2 <- plot_grid(top, middle, bottom, ncol = 1, rel_heights = c(0.9,0.9,1))

png(here::here("figures", "final-figs", "imgs", "figure-2.png"), res = 300, width = 4250, height = 3000)
figure2
dev.off()




