#clean plasmid data

library("tidyverse")
library("readxl")
library("patchwork")
library("ggtext")
library("ggpubr")
library("cowplot")
library("broom")
library("BBmisc")

source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

OD_path <- generate_paths("od", date, "csv")
JFP_path <- generate_paths("jfp", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
RFP_path <- generate_paths("rfp", date, "csv")
JFP_path <- generate_paths("rfp", date, "csv")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, OD_path), "OD") %>% filter(!is.na(cycle))
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, CFP_path), "CFP") %>% filter(!is.na(cycle))
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, YFP_path), "YFP") %>% filter(!is.na(cycle))
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, RFP_path), "RFP") %>% filter(!is.na(cycle))
JFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, JFP_path), "JFP") %>% filter(!is.na(cycle))
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "plasmid", date, "plate_layout.csv")) %>%
  janitor::clean_names()

if (date == "7July2024"){
  all_data <- OD %>%
    inner_join(., CFP %>% dplyr::select(cycle, well, CFP), by = c("well", "cycle")) %>%
    inner_join(., RFP %>% dplyr::select(cycle, well, RFP), by = c("well", "cycle")) %>%
    inner_join(., YFP %>% dplyr::select(cycle, well, YFP), by = c("well", "cycle")) %>%
    inner_join(., JFP %>% dplyr::select(cycle, well, JFP), by = c("well", "cycle")) %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(strain != "none") %>%
    mutate(partner = ifelse(partner == "none", NA, partner)) %>% filter(well != "C2")
} else if (date == "31July2024") {
  all_data <- OD %>%
    inner_join(., CFP %>% dplyr::select(cycle, well, CFP), by = c("well", "cycle")) %>%
    inner_join(., RFP %>% dplyr::select(cycle, well, RFP), by = c("well", "cycle")) %>%
    inner_join(., YFP %>% dplyr::select(cycle, well, YFP), by = c("well", "cycle")) %>%
    inner_join(., JFP %>% dplyr::select(cycle, well, JFP), by = c("well", "cycle")) %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(strain != "none") %>%
    mutate(partner = ifelse(partner == "none", NA, partner)) %>% filter(!strain %in% c("E0224 F+ M13+", "E0224 O+ helper"))
} else {
  all_data <- OD %>%
    inner_join(., CFP %>% dplyr::select(cycle, well, CFP), by = c("well", "cycle")) %>%
    inner_join(., RFP %>% dplyr::select(cycle, well, RFP), by = c("well", "cycle")) %>%
    inner_join(., YFP %>% dplyr::select(cycle, well, YFP), by = c("well", "cycle")) %>%
    inner_join(., JFP %>% dplyr::select(cycle, well, JFP), by = c("well", "cycle")) %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(strain != "none") %>%
    mutate(partner = ifelse(partner == "none", NA, partner)) %>%
    mutate(well = ifelse(strain == "S0240", gsub("H", "I", well), well))
}


# CFP YFP corrections
strains <- all_data %>% dplyr::select(strain) %>% filter(grepl("E0224", strain) == TRUE) %>% pull() %>% unique()
subsetted <- vector(mode = "list", length = length(strains))

for (i in strains){
  index <- which(strains == i)
  temp <- all_data %>%
    filter(strain %in% c(i, "S0240")) %>%
    mutate(phage = "none", 
           interaction = case_when(strain == i & is.na(partner) ~ "Emono",
                                   strain == "S0240" ~ "Smono",
                                   TRUE ~ strain)) 
  subsetted[[index]] <- temp
}

all_plate_data_corrected <- data.frame()
for (i in 1:length(subsetted)){
  YFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Emono", YFP, CFP, subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% first())
  CFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Smono", CFP, YFP, subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% first())
  all_tecan_adjusted_FP <- adjust_FP_values(subsetted[[i]], YFP_bleed_value = YFP_bleed, CFP_bleed_value = CFP_bleed)
  all_tecan_adjusted_OD <- all_tecan_adjusted_FP %>%
    mutate(E_corrected_OD = adjusted_CFP * subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% last()) %>%
    mutate(S_corrected_OD = adjusted_YFP * subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% last()) %>%
    inner_join(., subsetted[[i]] %>% dplyr::select(cycle, well, OD, hour, RFP, CFP, YFP, strain, partner), by = c("well", "cycle"))
  all_plate_data_corrected <- rbind(all_plate_data_corrected, all_tecan_adjusted_OD)
}

# final datasets
monoculture_dataset <- all_plate_data_corrected %>% 
  filter(is.na(partner) & strain %in% strains) %>% 
  group_by(well) %>%
  summarize(model_fit = fit_baranyi(hour, log(E_corrected_OD)),
            fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>%
  filter(fit_variable != "y0") %>%
  inner_join(., plate_layout, by = "well") %>%
  ungroup() %>% dplyr::select(-c(partner)) %>% 
  pivot_wider(names_from = fit_variable, values_from = model_fit) %>% rename(mono_growth_rate = growth_rate, mono_lag = lag, mono_ymax = ymax) %>% 
  dplyr::select(-well) %>%
  group_by(strain) %>% mutate(rep = 1:n()) %>% ungroup()

if (has_interactions == TRUE){
  coculture_dataset <- all_plate_data_corrected %>% 
    filter(!is.na(partner) & strain %in% strains) %>% 
    pivot_longer(cols = E_corrected_OD:S_corrected_OD) %>%
    group_by(well, name) %>%
    summarize(model_fit = fit_baranyi(hour, log(value)),
              fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>%
    filter(fit_variable != "y0") %>%
    ungroup() %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(fit_variable == "growth_rate") %>%
    pivot_wider(names_from = name, values_from = model_fit) %>%
    mutate(growth_rate_difference = E_corrected_OD - S_corrected_OD) %>% dplyr::select(-c(fit_variable, partner,E_corrected_OD, S_corrected_OD)) %>%
    inner_join(., all_plate_data_corrected, by = c("well", "strain")) %>% group_by(well) %>% filter(OD == max(OD)) %>% ungroup() %>%
    mutate(percent_E = E_corrected_OD / (E_corrected_OD + S_corrected_OD),
           s_to_e = S_corrected_OD / E_corrected_OD,
           total_coculture_density = OD) %>% dplyr::select(well, strain, interaction, growth_rate_difference, percent_E, s_to_e, total_coculture_density) %>%
    dplyr::select(-well) %>%
    group_by(strain, interaction) %>% mutate(rep = 1:n()) %>% ungroup()
} else{
  coculture_dataset <- all_plate_data_corrected %>% 
    filter(!is.na(partner) & strain %in% strains) %>% 
    pivot_longer(cols = E_corrected_OD:S_corrected_OD) %>%
    group_by(well, name) %>%
    summarize(model_fit = fit_baranyi(hour, log(value)),
              fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>%
    filter(fit_variable != "y0") %>%
    ungroup() %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(fit_variable == "growth_rate") %>%
    pivot_wider(names_from = name, values_from = model_fit) %>%
    mutate(growth_rate_difference = E_corrected_OD - S_corrected_OD) %>% dplyr::select(-c(fit_variable, partner,E_corrected_OD, S_corrected_OD)) %>%
    inner_join(., all_plate_data_corrected, by = c("well", "strain")) %>% group_by(well) %>% filter(OD == max(OD)) %>% ungroup() %>%
    mutate(percent_E = E_corrected_OD / (E_corrected_OD + S_corrected_OD),
           s_to_e = S_corrected_OD / E_corrected_OD,
           total_coculture_density = OD) %>% dplyr::select(well, strain, growth_rate_difference, percent_E, s_to_e, total_coculture_density) %>%
    dplyr::select(-well) %>%
    group_by(strain) %>% mutate(rep = 1:n()) %>% ungroup()
}


aggregated_dataset <- monoculture_dataset %>% inner_join(., coculture_dataset, by = c("strain", "rep")) %>%
  mutate(f_presence = ifelse(grepl("F+", strain) == TRUE, TRUE, FALSE),
         p3_presence = ifelse(grepl("p3", strain) == TRUE, TRUE, FALSE),
         delta_M13_presence = ifelse(grepl("helper", strain) == TRUE, TRUE, FALSE),
         m13_presence = ifelse(grepl("M13+", strain) == TRUE, TRUE, FALSE),
         p4_presence = ifelse(grepl("p4", strain) == TRUE, TRUE, FALSE),
         o_presence = ifelse(grepl("pOX38", strain) == TRUE, TRUE, FALSE)) %>%
  mutate(date = date) %>%
  filter(!if_any(everything(), is.na))

all_plate_data_corrected <- all_plate_data_corrected %>% mutate(date = date)


# upload cfus

if (date %in% c("29July2024", "31July2024", "14August2024")){
  cfus <- generate_paths("cfu", date, "xlsx")
  CFU <- load_pfu_data(here::here("experimental-data", "tecan-data","plasmid", date, cfus)) %>%
    filter(!is.na(cfu))
  full_data <- CFU %>% inner_join(., plate_layout, by = "well")
  percent_retained <- full_data %>% dplyr::select(well, cfu, plate, strain, partner) %>% filter(plate != "S") %>% 
    pivot_wider(values_from = "cfu", names_from = "plate") %>% mutate(percent_retained = abx / E) %>% 
    filter(!is.nan(percent_retained)) %>% mutate(date = date)
}

