# ATB
# figure 4
# partner metabolomics

# library 
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

# set metabolomics p value threshold of interest and # of permutations for opls-da
p = 0.05
permutations = 1000

# load metabolomics data for analysis
cfus <- read_csv(here::here("experimental-data", "metabolomics", "spent-media-metabolomics-cfus.csv")) %>%
  janitor::clean_names() %>% 
  pull(cfu) # load cfus at the time of filter sterilzation

raw_metabolites <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(-c(index, q1_da, molecular_weight_da, formula, ionization_model, level, cpd_id, hmdb, pub_chem_cid,
                                              cas, ch_ebi, metlin, kegg_map, class_i, class_ii)) %>% 
  dplyr::select(c(compounds, e1_midlog:m_pm3_depleted)) 

compound_list <- raw_metabolites[,1] # compound list for ease of use downstream
normalized_data <- sweep(raw_metabolites[,-c(1)], 2, cfus / 1e6, FUN = "/") # biomass normalization
df_filtered <- normalized_data

df_filtered_fixed <- as.data.frame(t(apply(df_filtered, 1, function(row) {
  non_zero_vals <- row[row > 0]
  pseudo <- if (length(non_zero_vals) > 0) min(non_zero_vals) / 5 else 0
  row[row == 0] <- pseudo
  return(row)
}))) # missing values imputation with a limit of detection thresholds

tmp_log_scaled <- df_filtered_fixed |> log2() |> 
  t() |>
  scale(center = TRUE, scale = TRUE) |> t() 

# set up metadata and data for analysis
metadata <- data.frame(
  SampleID = c("e1_midlog", "e2_midlog", "e3_midlog", 
               "p1_midlog", "p2_midlog", "p3_midlog", 
               "pm1_midlog", "pm2_midlog", "pm3_midlog",
               "s_e1_depleted", "s_e2_depleted", "s_e3_depleted", 
               "s_p1_depleted", "s_p2_depleted", "s_p3_depleted", 
               "s_pm1_depleted", "s_pm2_depleted", "s_pm3_depleted",
               "m_pm1_depleted", "m_pm2_depleted", "m_pm3_depleted"),
  Group = c("WT", "WT", "WT", 
            "F128+", "F128+", "F128+", 
            "F128+ M13+", "F128+ M13+", "F128+ M13+",
            "WT", "WT", "WT", 
            "F128+", "F128+", "F128+", 
            "F128+ M13+", "F128+ M13+", "F128+ M13+",
            "F128+ M13+", "F128+ M13+", "F128+ M13+"),
  Partner = c(rep(NA, 9), rep("S", 9), rep("M", 3)))

###### compare WT
strain = "WT"
partner = "S"
no_compounds <- tmp_log_scaled[, c(metadata$SampleID[metadata$Group == strain]), drop = FALSE]
rows_to_remove <- apply(no_compounds[,c(2:dim(no_compounds)[2])], 1, function(row) length(unique(row)) == 1)
tmp <- no_compounds[!rows_to_remove, ]

# get VIP using opls
class_vector <- factor(c(rep(strain, 3), rep(partner, 3)), levels = c(strain, partner))
opls_model <- opls(t(tmp), class_vector, crossvalI = 6, orthoI = NA, predI = 1, scaleC = "none", permI = permutations)
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
partner_cols <- metadata %>% filter(Partner == partner & Group == strain) %>% pull(SampleID)
baseline_cols <- metadata %>% filter(is.na(Partner) & Group == strain) %>% pull(SampleID)
fc <- rowMeans(df_filtered[c(which(rows_to_remove != TRUE)), ][partner_cols]) / rowMeans(df_filtered[c(which(rows_to_remove != TRUE)), ][baseline_cols])
stats_df_WT <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], p_value = p_values, fdr = BH, fold_change = fc) %>%
  inner_join(., vip_df, by = c("compounds")) %>% mutate(condition = strain, partner = partner)


###### compare F128+
strain = "F128+"
partner = "S"
no_compounds <- tmp_log_scaled[, c(metadata$SampleID[metadata$Group == strain]), drop = FALSE]
rows_to_remove <- apply(no_compounds[,c(2:dim(no_compounds)[2])], 1, function(row) length(unique(row)) == 1)
tmp <- no_compounds[!rows_to_remove, ]

# get VIP using opls
class_vector <- factor(c(rep(strain, 3), rep(partner, 3)), levels = c(strain, partner))
opls_model <- opls(t(tmp), class_vector, crossvalI = 6, orthoI = NA, predI = 1, scaleC = "none", permI = permutations)
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
partner_cols <- metadata %>% filter(Partner == partner & Group == strain) %>% pull(SampleID)
baseline_cols <- metadata %>% filter(is.na(Partner) & Group == strain) %>% pull(SampleID)
fc <- rowMeans(df_filtered[c(which(rows_to_remove != TRUE)), ][partner_cols]) / rowMeans(df_filtered[c(which(rows_to_remove != TRUE)), ][baseline_cols])
stats_df_F128 <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], p_value = p_values, fdr = BH, fold_change = fc) %>%
  inner_join(., vip_df, by = c("compounds")) %>% mutate(condition = strain, partner = partner)


###### compare F128+ M13+
strain = "F128+ M13+"
partner = "S"
s_phage <- metadata %>% filter(Partner == partner & Group == strain) %>% pull(SampleID)
base <- metadata %>% filter(is.na(Partner) & Group == strain) %>% pull(SampleID)
no_compounds <- cbind(tmp_log_scaled[,base], tmp_log_scaled[,s_phage])
rows_to_remove <- apply(no_compounds[,c(2:dim(no_compounds)[2])], 1, function(row) length(unique(row)) == 1)
tmp <- no_compounds[!rows_to_remove, ]

# get VIP using opls
class_vector <- factor(c(rep(strain, 3), rep(partner, 3)), levels = c(strain, partner))
opls_model <- opls(t(tmp), class_vector, crossvalI = 6, orthoI = NA, predI = 1, scaleC = "none", permI = permutations)
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
partner_cols <- metadata %>% filter(Partner == partner & Group == strain) %>% pull(SampleID)
baseline_cols <- metadata %>% filter(is.na(Partner) & Group == strain) %>% pull(SampleID)
fc <- rowMeans(df_filtered[c(which(rows_to_remove != TRUE)), ][partner_cols]) / rowMeans(df_filtered[c(which(rows_to_remove != TRUE)), ][baseline_cols])
stats_df_phage_s <- data.frame(compounds = compound_list[which(rows_to_remove != TRUE), ], p_value = p_values, fdr = BH, fold_change = fc) %>%
  inner_join(., vip_df, by = c("compounds")) %>% mutate(condition = strain, partner = partner)


###### compare F128+ M13+
strain = "F128+ M13+"
partner = "M"
s_phage <- metadata %>% filter(Partner == partner & Group == strain) %>% pull(SampleID)
base <- metadata %>% filter(is.na(Partner) & Group == strain) %>% pull(SampleID)
no_compounds <- cbind(tmp_log_scaled[,base], tmp_log_scaled[,s_phage])
rows_to_remove <- apply(no_compounds[,c(2:dim(no_compounds)[2])], 1, function(row) length(unique(row)) == 1)
tmp <- no_compounds[!rows_to_remove, ]

# get VIP using opls
class_vector <- factor(c(rep(strain, 3), rep(partner, 3)), levels = c(strain, partner))
opls_model <- opls(t(tmp), class_vector, crossvalI = 6, predI = 1, orthoI = NA, scaleC = "none", permI = permutations)
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
partner_cols <- metadata %>% filter(Partner == partner & Group == strain) %>% pull(SampleID)
baseline_cols <- metadata %>% filter(is.na(Partner) & Group == strain) %>% pull(SampleID)
fc <- rowMeans(df_filtered[,][partner_cols]) / rowMeans(df_filtered[, ][baseline_cols])
stats_df_phage_m <- data.frame(compounds = compound_list, p_value = p_values, fdr = BH, fold_change = fc) %>%
  inner_join(., vip_df, by = c("compounds")) %>% mutate(condition = strain, partner = partner)

combined <- rbind(stats_df_WT, stats_df_F128, stats_df_phage_s, stats_df_phage_m)
combined$all_adjust <- p.adjust(combined$p_value, method = "BH")

# part A
partA <- combined %>%
  mutate(significant = ifelse(p_value < p & VIP > 1, "sig", "non sig")) %>%
  mutate(shape = ifelse(compounds == "Lactic acid", "shape", "no shape")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("WT", "F128+", "F128+ M13+"))) %>%
  mutate(partner = ifelse(partner == "S", "*S. enterica*", "*M. extorquens*"),
         partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  mutate(log2fc = log2(fold_change)) %>%
  mutate(log2fc = ifelse(log2fc < 0 & is.infinite(log2fc), -13, log2fc)) %>%
  mutate(log2fc = ifelse(log2fc > 0 & is.infinite(log2fc), 13, log2fc)) %>%
  ggplot(aes(x = log2fc, y = -log10(p_value), color = significant, shape = shape)) +
  geom_point(size = 2) + facet_grid(condition~partner) +
  ylab("-log10(P-value)") +
  xlab("log2(Fold-change)") +
  xlim(-14, 14)+
  geom_hline(yintercept = -log10(p), color = "red", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  geom_vline(xintercept = -12, color = "black", linetype = "dashed")+
  geom_vline(xintercept = 12, color = "black", linetype = "dashed")+
  scale_color_manual(values = c("sig" = "black", "non sig" = "darkgrey"))+
  theme_bw(base_size = 16)+
  theme(axis.text = element_markdown(), strip.text = element_markdown(),
        legend.position = "none", strip.background = element_blank())

# supp table 5
supp5 <- combined %>% filter(p_value < p & VIP > 1 & log2(fold_change) < 0) %>% mutate(log2_fold_change = log2(fold_change)) %>%
  mutate(compounds = tolower(compounds)) %>% select(compounds, fdr, p_value, log2_fold_change, condition, partner, VIP) %>%
  mutate(partner = ifelse(partner == "S", "S. enterica", "M. extorquens")) %>% group_by(condition, partner) %>% arrange(log2_fold_change)

write.csv(supp5, file = here::here("figures", "final-figs", "tables", "supplemental-table-5.csv"))

# part B
class_names <- read_excel(here::here("experimental-data", "metabolomics", "ALL_sample_data.xlsx")) %>%
  janitor::clean_names() %>% dplyr::select(c(compounds, class_i, class_ii))

partB <- combined %>% inner_join(., class_names, by = "compounds") %>%
  filter(p_value < p & VIP > 1 & log2(fold_change) < 0) %>%
  group_by(condition, partner) %>% mutate(total_sig = n()) %>%
  group_by(condition, class_i, partner) %>%
  mutate(number = n()) %>% select(condition, class_i, total_sig, number, partner) %>% unique() %>%
  group_by(condition, partner) %>% arrange(desc(number)) %>% slice(1:3) %>%
  mutate(class_i = case_when(class_i == "Aldehyde,Ketones,Esters" ~ "aldehydes",
                             class_i == "Amino acid and Its metabolites" ~ "amino acids",
                             class_i == "Benzene and substituted derivatives" ~ "benzenes",
                             class_i == "FA" ~ "fatty acids",
                             class_i == "GP" ~ "lipids",
                             class_i == "Nucleotide and Its metabolites" ~ "nucleotides",
                             class_i == "Organic acid and Its derivatives" ~ "organic acids",
                             class_i == "Carbohydrates and Its metabolites" ~ "carbohydrates",
                             class_i == "Heterocyclic compounds" ~ "heterocyclics")) %>%
  mutate(partner = ifelse(partner == "S", "*S. enterica*", "*M. extorquens*"),
         partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  dplyr::mutate(condition = factor(condition, levels = c("WT", "F128+", "F128+ M13+"))) %>%
  ggplot(aes(x = fct_reorder(class_i, number), y = (number / total_sig) * 100)) +
  geom_bar(stat = "identity") + facet_grid(condition~partner, scales = "free") + coord_flip() + ylab("Percent of significantly depleted compounds") +
  theme_bw(base_size = 16) + theme(strip.background = element_blank(), axis.text.x = element_markdown(), axis.title.y = element_blank(),
                                   strip.text = element_markdown())

# part C
lactic <- as.data.frame(tmp_log_scaled) %>%
  cbind(., compound_list) %>%
  filter(compounds == "Lactic acid")

full_data <- lactic %>% pivot_longer(cols = c(e1_midlog:m_pm3_depleted)) %>% select(value, name) %>%
  mutate(partner = case_when(name %in% c("e1_midlog", "e2_midlog", "e3_midlog",
                                         "p1_midlog", "p2_midlog", "p3_midlog",
                                         "pm1_midlog", "pm2_midlog", "pm3_midlog") ~ NA,
                             name %in% c("m_pm1_depleted", "m_pm2_depleted", "m_pm3_depleted") ~ "*M. extorquens*",
                             TRUE ~ "*S. enterica*")) %>%
  mutate(name = case_when(name %in% c("e1_midlog", "e2_midlog", "e3_midlog") ~ "WT",
                          name %in% c("p1_midlog", "p2_midlog", "p3_midlog") ~ "F128+",
                          name %in% c("pm1_midlog", "pm2_midlog", "pm3_midlog") ~ "F128+ M13+",
                          name %in% c("s_e1_depleted", "s_e2_depleted", "s_e3_depleted") ~ "WT",
                          name %in% c("s_p1_depleted", "s_p2_depleted", "s_p3_depleted") ~ "F128+",
                          name %in% c("s_pm1_depleted", "s_pm2_depleted", "s_pm3_depleted") ~ "F128+ M13+",
                          name %in% c("m_pm1_depleted", "m_pm2_depleted", "m_pm3_depleted") ~ "F128+ M13+")) %>%
  mutate(condition = ifelse(partner %in% c("*S. enterica*", "*M. extorquens*"), "depleted", "midlog"))

stats <- full_data %>%
  rbind(., data.frame(value = full_data %>% filter(is.na(partner) & name == "F128+ M13+") %>% pull(value),
                      name = rep("F128+ M13+", 3), partner = rep("*M. extorquens*", 3), condition = rep("midlog", 3))) %>%
  mutate(partner = ifelse(is.na(partner), "*S. enterica*", partner)) %>%
  mutate(name = ifelse(name  == "F128+ M13+", "F128+<br>M13+", name),
         name = factor(name, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  mutate(condition = factor(condition, levels = c("midlog", "depleted"))) %>%
  group_by(name, partner) %>%
  t_test(value ~ condition, var.equal = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "name") %>%
  mutate(label = paste0(scientific(p.adj,digits = 2), " (", p.adj.signif, ")")) %>%
  mutate(y.position = y.position - 0.1)

partC <- full_data %>%
  rbind(., data.frame(value = full_data %>% filter(is.na(partner) & name == "F128+ M13+") %>% pull(value),
                      name = rep("F128+ M13+", 3), partner = rep("*M. extorquens*", 3), condition = rep("midlog", 3))) %>%
  mutate(partner = ifelse(is.na(partner), "*S. enterica*", partner)) %>%
  mutate(name = ifelse(name  == "F128+ M13+", "F128+<br>M13+", name),
         name = factor(name, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
  mutate(condition = factor(condition, levels = c("midlog", "depleted"))) %>%
  ggplot(aes(x = name, y = value, color = condition)) +
  geom_point(size = 2, stroke = 1.5, position = position_dodge(0.9)) +
  ylab("Normalized lactate intensity") +
  facet_wrap(~partner) +
  stat_summary(aes(fill = condition), position = position_dodge(0.9), fun = mean, shape = '-', size = 5, color = 'black')+
  stat_pvalue_manual(stats, label = "label", tip.length = 0.0, size = 3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown(), text = element_text(color = "black"))

legendC <- get_plot_component(full_data %>%
                                rbind(., data.frame(value = full_data %>% filter(is.na(partner) & name == "F128+ M13+") %>% pull(value),
                                                    name = rep("F128+ M13+", 3), partner = rep("*M. extorquens*", 3), condition = rep("midlog", 3))) %>%
                                mutate(partner = ifelse(is.na(partner), "*S. enterica*", partner)) %>%
                                mutate(name = ifelse(name  == "F128+ M13+", "F128+<br>M13+", name),
                                       name = factor(name, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
                                mutate(partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*"))) %>%
                                mutate(condition = factor(condition, levels = c("midlog", "depleted"))) %>%
                                ggplot(aes(x = name, y = value, fill = condition)) +
                                geom_bar(stat = "identity") +
                                ylab("Normalized lactate intensity") +
                                facet_wrap(~partner) +
                                ylim(-1.75, 2.5)+
                                theme_bw(base_size = 16) +labs(fill = "") +
                                theme(axis.text = element_markdown(),
                                      legend.position = "bottom", strip.background = element_blank(), axis.title.x = element_blank(),
                                      strip.text = element_markdown(), text = element_text(color = "black")),
                              'guide-box-bottom', return_all = TRUE)

full_partC <- plot_grid(partC, legendC, rel_heights = c(1, 0.1), ncol = 1)

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
  mutate(strain = case_when(strain == "E0224" ~ "WT",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 + M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(supplement = ifelse(is.na(supplement), "std", "lact+"),
         supplement = factor(supplement, levels = c("std", "lact+"))) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>% mutate(percent_partner = percent_M + percent_S) %>%
  group_by(partner) %>%
  t_test(percent_partner ~ supplement) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "supplement") %>%
  mutate(label = paste0(scientific(p.adj,digits = 2), " (", p.adj.signif, ")")) %>%
  mutate(y.position = (y.position * 100)) 

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
  mutate(strain = case_when(strain == "E0224" ~ "WT",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("WT", "F128+", "F128+<br>M13+"))) %>%
  mutate(partner = case_when(partner == "S0240" ~ "*S. enterica*",
                             partner == "M0104" ~ "*M. extorquens*",
                             partner == "S0240 + M0104" ~ "*S. enterica*<br>+<br>*M. extorquens*"),
         partner = factor(partner, levels = c("*S. enterica*", "*M. extorquens*", "*S. enterica*<br>+<br>*M. extorquens*"))) %>%
  mutate(supplement = ifelse(is.na(supplement), "std", "lact+"),
         supplement = factor(supplement, levels = c("std", "lact+"))) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>% mutate(percent_partner = percent_M + percent_S) %>%
  ggplot(aes(x = supplement, y = percent_partner * 100)) + 
  stat_summary(fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5) +
  facet_wrap(~partner, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 16) +
  ylab("Percent of co-culture\n composed of partner species") +
  stat_pvalue_manual(stats, label = "label", tip.length = 0.0, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text = element_markdown(),
        legend.position = "none", strip.background = element_blank(), axis.title.x = element_blank(),
        strip.text = element_markdown())

# final figure
top <- plot_grid(partA, partB, ncol = 2, rel_widths = c(1,1), label_size = 26, labels = c("A", "B"))
bottom <- plot_grid(full_partC, partD, ncol = 2, rel_widths = c(0.75, 1), label_size = 26, labels = c("C", "D"))
figure4 <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1,1))

png(here::here("figures", "final-figs", "imgs", "figure-4.png"), res = 300, width = 3700, height = 2700)
figure4
dev.off()


