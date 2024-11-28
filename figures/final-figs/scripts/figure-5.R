# ATB
# figure 5
# reproducibility

# load packages
library("tidyverse")
library("readxl")
library("patchwork")
library("ggtext")
library("ggpubr")
library("cowplot")
library("broom")
library("BBmisc")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

# part A - O and F fac/coop/comp
# competition - all comparisons significant except pOX M13 and uninf, F128 and F128 M13, and pOX and M13F128
# facilitation - all comparisons significant except pOX to uninf and pOX M13 to F128
# mutualism - all comparisons significant except pOX to uninf
dates <- c("29July2024", "31July2024", "14August2024")
has_interactions <- FALSE

aggregated <- data.frame()
time_series <- data.frame()

for (date in dates){
  source(here::here("data-generation", "experimental", "clean-plasmid-tecan-data.R"))
  temp <- aggregated_dataset
  temp1 <- all_plate_data_corrected
  aggregated <- rbind(aggregated, temp)
  time_series <- rbind(time_series, temp1)
}

aggregated_mutualism <- aggregated %>% filter(strain %in% c("E0224", "E0224 F+", "E0224 F+ M13+", "E0224 O+", "E0224 O+ M13+"))

date <- "21March2024"
source(here::here("data-generation", "experimental", "clean-facilitation-tecan-data.R"))
march21 <- all_data_corrected %>% mutate(date = date) %>%
  filter(partner == "S0240") %>% select(well, S_corrected_OD, E_corrected_OD, OD, hour, strain, partner, date)

date <- "27May2024"

OD_path <- generate_paths("od", date, "csv")
JFP_path <- generate_paths("jfp", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
RFP_path <- generate_paths("rfp", date, "csv")
cfus <- generate_paths("cfu", date, "xlsx")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, OD_path), "OD") %>% filter(!is.na(cycle))
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, CFP_path), "CFP") %>% filter(!is.na(cycle))
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, YFP_path), "YFP") %>% filter(!is.na(cycle))
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, RFP_path), "RFP") %>% filter(!is.na(cycle))
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "plasmid", date, "plate_layout.csv")) %>%
  janitor::clean_names()
CFU <- load_pfu_data(here::here("experimental-data", "tecan-data","plasmid", date, cfus)) %>%
  filter(!is.na(cfu))

all_data_plate1 <- OD %>%
  inner_join(., CFP %>% select(cycle, well, CFP), by = c("well", "cycle")) %>%
  inner_join(., RFP %>% select(cycle, well, RFP), by = c("well", "cycle")) %>%
  inner_join(., YFP %>% select(cycle, well, YFP), by = c("well", "cycle")) %>%
  inner_join(., plate_layout %>% filter(plate_number == 1), by = "well") %>%
  filter(strain != "none") %>%
  filter(!well %in% c("B9", "C5", "C8", "G5")) %>%
  mutate(partner = ifelse(partner == "none", NA, partner))

all_plates <- all_data_plate1 %>% select(cycle, well, OD, hour, CFP, RFP, YFP, strain, partner, interaction, plate_number) %>%
  rename(nutrients = interaction)

uninfected <- all_plates %>%
  filter(strain %in% c("E0224", "S0240")) %>%
  mutate(phage = "none", 
         interaction = case_when(strain == "E0224" & is.na(partner) ~ "Emono",
                                 strain == "S0240" ~ "Smono",
                                 TRUE ~ strain)) 

fplasmid <- all_plates %>%
  filter(strain %in% c("E0224 F+", "S0240")) %>%
  mutate(phage = "none", 
         interaction = case_when(strain == "E0224 F+" & is.na(partner) ~ "Emono",
                                 strain == "S0240" ~ "Smono",
                                 TRUE ~ strain)) 

infected <- all_plates %>%
  filter(strain %in% c("E0224 F+ M13+", "S0240")) %>%
  mutate(phage = "none", 
         interaction = case_when(strain == "E0224 F+ M13+" & is.na(partner) ~ "Emono",
                                 strain == "S0240" ~ "Smono",
                                 TRUE ~ strain)) 


subsetted <- list(uninfected = uninfected, fplasmid = fplasmid, infected = infected)
all_plate_data_corrected <- data.frame()
for (i in 1:length(subsetted)){
  YFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Emono", YFP, CFP, subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% first())
  CFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Smono", CFP, YFP, subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% first())
  all_tecan_adjusted_FP <- adjust_FP_values(subsetted[[i]], YFP_bleed_value = YFP_bleed, CFP_bleed_value = CFP_bleed)
  all_tecan_adjusted_OD <- all_tecan_adjusted_FP %>%
    mutate(E_corrected_OD = adjusted_CFP * subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% last()) %>%
    mutate(S_corrected_OD = adjusted_YFP * subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% last()) %>%
    inner_join(., subsetted[[i]] %>% select(cycle, well, OD, hour, RFP, CFP, YFP, strain, partner, nutrients), by = c("well", "cycle"))
  all_plate_data_corrected <- rbind(all_plate_data_corrected, all_tecan_adjusted_OD)
}

may27 <- all_plate_data_corrected %>% mutate(date = date) %>%
  filter(partner == "S0240") %>% select(well, S_corrected_OD, E_corrected_OD, OD, hour, strain, partner, date)

facilitation_only <- rbind(may27, march21) %>% mutate(percent_S = S_corrected_OD / (E_corrected_OD + S_corrected_OD)) %>%
  group_by(well) %>% filter(OD == max(OD)) %>% slice_head(n = 1) %>% ungroup() %>% dplyr::select(strain, partner, date, percent_S) %>% mutate(interaction = "facilitation")

dates <- c("16September2024", "18September2024", "21September2024")
all_dates <- data.frame()
for (date in dates){
  OD_path <- generate_paths("od", date, "csv")
  CFP_path <- generate_paths("cfp", date, "csv")
  YFP_path <- generate_paths("yfp", date, "csv")
  
  OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, OD_path), "OD") %>% filter(!is.na(cycle))
  CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, CFP_path), "CFP") %>% filter(!is.na(cycle))
  YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, YFP_path), "YFP") %>% filter(!is.na(cycle))
  plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "plasmid", date, "plate_layout.csv")) %>%
    janitor::clean_names()
  
  if (date == "16September2024"){
    all_data <- OD %>%
    inner_join(., CFP %>% dplyr::select(cycle, well, CFP), by = c("well", "cycle")) %>%
    inner_join(., YFP %>% dplyr::select(cycle, well, YFP), by = c("well", "cycle")) %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(strain != "none") %>%
    mutate(partner = ifelse(partner == "none", NA, partner)) %>% filter(!strain %in% c("E0224", "E0224 F+ M13+"))
  } else if (date == "18September2024"){
    all_data <- OD %>%
      inner_join(., CFP %>% dplyr::select(cycle, well, CFP), by = c("well", "cycle")) %>%
      inner_join(., YFP %>% dplyr::select(cycle, well, YFP), by = c("well", "cycle")) %>%
      inner_join(., plate_layout, by = "well") %>%
      filter(strain != "none") %>%
      mutate(partner = ifelse(partner == "none", NA, partner)) %>% filter(interaction %in% c("competition", "facilitation") | is.na(interaction)) %>%
      filter(!(interaction == "facilitation" & strain == "E0224 F+")) %>%
      filter(!(is.na(interaction) & strain == "E0224 F+ M13+"))
  } else {
    all_data <- OD %>%
      inner_join(., CFP %>% dplyr::select(cycle, well, CFP), by = c("well", "cycle")) %>%
      inner_join(., YFP %>% dplyr::select(cycle, well, YFP), by = c("well", "cycle")) %>%
      inner_join(., plate_layout, by = "well") %>%
      filter(strain != "none") %>%
      mutate(partner = ifelse(partner == "none", NA, partner)) 
  }
  all_dates <- rbind(all_dates, all_data %>% mutate(date = date))
}

generate_letter_number_combinations <- function() {
  letters <- LETTERS
  numbers <- 1:9
  combinations <- as.vector(outer(letters, numbers, paste0))
  return(combinations)
}
new_wells <- all_dates %>% select(well, strain, partner, interaction, date) %>% 
  unique() %>% mutate(new_well = generate_letter_number_combinations()[1:n()])
all_dates <- all_dates %>% inner_join(., new_wells, by = c("well", "strain", "partner", "interaction", "date"))
all_dates <- all_dates %>% dplyr::select(-c(well)) %>% rename(well = new_well)

strains <- all_dates %>% dplyr::select(strain) %>% filter(grepl("E0224", strain) == TRUE) %>% pull() %>% unique()
subsetted <- vector(mode = "list", length = length(strains))
  
for (i in strains){
    index <- which(strains == i)
    temp <- all_dates %>%
      filter(strain %in% c(i, "S0240")) %>%
      mutate(phage = "none", 
             interaction = case_when(strain == i & is.na(partner) ~ "Emono",
                                     strain == "S0240" ~ "Smono",
                                     TRUE ~ strain)) 
    subsetted[[index]] <- temp
  }
names(subsetted) <- strains
all_plate_data_corrected <- data.frame()
for (i in 1:length(subsetted)){
  if (names(subsetted[i]) %in% c("E0224 O+", "E0224 O+ M13+")){
      YFP_bleed <- mean(c(0.02266644, 0.03120145, 0.02163692))
      CFP_bleed <- mean(c(0.01738057, 0.01738057, 0.01738057))
      all_tecan_adjusted_FP <- adjust_FP_values(subsetted[[i]], YFP_bleed_value = YFP_bleed, CFP_bleed_value = CFP_bleed)
      all_tecan_adjusted_OD <- all_tecan_adjusted_FP %>%
        mutate(E_corrected_OD = adjusted_CFP * subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% last()) %>%
        mutate(S_corrected_OD = adjusted_YFP * subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% last()) %>%
        inner_join(., subsetted[[i]] %>% dplyr::select(cycle, well, OD, hour, CFP, YFP, strain, partner), by = c("well", "cycle"))
      all_plate_data_corrected <- rbind(all_plate_data_corrected, all_tecan_adjusted_OD)
    } else{
      YFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Emono", YFP, CFP, subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% first())
      CFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Smono", CFP, YFP, subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% first())
      all_tecan_adjusted_FP <- adjust_FP_values(subsetted[[i]], YFP_bleed_value = YFP_bleed, CFP_bleed_value = CFP_bleed)
      all_tecan_adjusted_OD <- all_tecan_adjusted_FP %>%
        mutate(E_corrected_OD = adjusted_CFP * subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% last()) %>%
        mutate(S_corrected_OD = adjusted_YFP * subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% last()) %>%
        inner_join(., subsetted[[i]] %>% dplyr::select(cycle, well, OD, hour, CFP, YFP, strain, partner), by = c("well", "cycle"))
      all_plate_data_corrected <- rbind(all_plate_data_corrected, all_tecan_adjusted_OD)
    }
}

all_interactions <- all_plate_data_corrected %>% group_by(well) %>% filter(OD == max(OD)) %>% ungroup() %>% mutate(percent_S = S_corrected_OD / (S_corrected_OD + E_corrected_OD)) %>% 
  dplyr::select(well, strain, partner, percent_S) %>% inner_join(., all_dates %>% dplyr::select(well, interaction, partner, strain, date) %>% unique, by = c("well", "partner", "strain")) %>%
  dplyr::select(-c(well)) %>%
  rbind(., aggregated_mutualism %>% dplyr::select(strain, percent_E, date) %>% mutate(partner = "S0240") %>% mutate(percent_S = 1 - percent_E) %>% 
          mutate(interaction = "mutualism") %>% dplyr::select(-c(percent_E))) %>%
  rbind(., facilitation_only)

partA <- all_interactions %>%
  filter(!is.na(interaction)) %>%
  mutate(strain = case_when(strain %in% c("E0224 F+") ~ "F128+",
                            strain %in% c("E0224 F+ M13+") ~ "F128+<br>M13+",
                            strain %in% c("E0224") ~ "uninf",
                            strain %in% c("E0224 O+") ~ "pOX38+",
                            strain %in% c("E0224 O+ M13+") ~ "pOX38+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+", "pOX38+", "pOX38+<br>M13+"))) %>%
  ggplot(aes(x = strain, y = percent_S, color = strain)) +
  facet_wrap(~interaction, scales = "free") +
  stat_summary(aes(fill = strain), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5)+
  ylab("percent *S. enterica*") +
  theme_bw(base_size = 16)+
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "pOX38+" = "#A95AA1", "F128+<br>M13+" = "#0F2080", "pOX38+<br>M13+" = "#0F2080")) +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown())

# part B - MG1655
# only F128 M13 to uninf is significant
date <- "21May2024"
may21 <- load_pfu_data(here::here("experimental-data", "non-tecan-plates","facilitation", date, "cfu_21May2024.xlsx")) %>%
  filter(!is.na(cfu)) %>% inner_join(., read_csv(here::here("experimental-data", "non-tecan-plates", "facilitation", date, "plate_layout.csv")) %>%
                                       janitor::clean_names(), by = "well") %>% filter(partner == "S0240") %>% select(well, cfu, plate, strain, partner) %>% mutate(date = date)

date <- "20May2024"
may20 <- load_pfu_data(here::here("experimental-data", "non-tecan-plates","facilitation", date, "cfu_20May2024.xlsx")) %>%
  filter(!is.na(cfu)) %>% inner_join(., read_csv(here::here("experimental-data", "non-tecan-plates", "facilitation", date, "plate_layout.csv")) %>%
                                       janitor::clean_names(), by = "well") %>% filter(partner == "S0240") %>% select(cfu, well, plate, strain, partner) %>% mutate(date = date)

partB <- rbind(may21, may20) %>% 
  mutate(strain = ifelse(strain == "MG1655 WT", "MG1655", strain)) %>%
  mutate(strain = ifelse(strain == "MG1655 WT F+", "MG1655 F+", strain)) %>%
  mutate(strain = ifelse(strain == "MG1655 WT F+ M13+", "MG1655 F+ M13+", strain)) %>%
  filter(grepl("old", strain) == FALSE) %>%
  pivot_wider(names_from = "plate", values_from = "cfu") %>% mutate(percent_S = S / (E + S)) %>% 
  mutate(strain = case_when(strain == "MG1655" ~ "uninf",
                            strain == "MG1655 F+" ~ "F128+",
                            strain == "MG1655 F+ M13+" ~ "F128+<br>M13+"),
         strain = factor(strain, levels = c("uninf", "F128+", "F128+<br>M13+"))) %>%
  ggplot(aes(x = strain, y = percent_S, color = strain)) + 
  stat_summary(aes(fill = strain), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5)+
  ylab("percent *S. enterica*") +
  theme_bw(base_size = 16)+
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "F128+<br>M13+" = "#0F2080")) +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
                                   legend.position = "none", strip.background = element_blank(),
                                   axis.title.y = element_markdown())

# part C - biomass-normalized lactate promega
# only pOX M13 to uninf significant for metB
# none significant for MG1655
promega <- read_csv(here::here("experimental-data", "hplc", "4October2024", "lactate_kit.csv"))
for_model <- promega %>% filter(type == "standard") %>% dplyr::select(condition, LU_avg)
for_model$condition <- as.numeric(for_model$condition)
for_model <- for_model %>% filter(condition < 5)

tested_conditions <- promega %>% filter(type != "standard") %>% group_by(condition, bio_rep, midlog_biomass) %>%
  summarize(LU_avg = mean(LU_avg)) %>% ungroup()

linear_model <- lm(condition ~ LU_avg, data = for_model)

tested_conditions$prediction <- predict(linear_model, newdata = data.frame(tested_conditions %>% dplyr::select(LU_avg)))

unused <- tested_conditions %>% mutate(conc_per_cell = prediction / midlog_biomass) %>%
  mutate(conc_per_cell = ifelse(conc_per_cell < 0, 0, conc_per_cell)) %>%
  mutate(background = case_when(grepl("MG1655", condition) ~ "MG1655",
                                grepl("E0224", condition) ~ "*∆metB*"),
         parasite = case_when(condition %in% c("MG1655", "E0224") ~ "uninf",
                              condition %in% c("MG1655 F+", "E0224 F+") ~ "F128+",
                              condition %in% c("MG1655 F+ M13+", "E0224 F+ M13+") ~ "F128+<br>M13+",
                              condition %in% c("E0224 O+") ~ "pOX38+",
                              condition %in% c("E0224 O+ M13+") ~ "pOX38+<br>M13+"),
         parasite = factor(parasite, levels = c("uninf", "F128+", "F128+<br>M13+", "pOX38+", "pOX38+<br>M13+"))) %>%
  ggplot(aes(x = parasite, y = conc_per_cell, color = parasite)) +
  stat_summary(aes(fill = parasite), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5)+
  theme_bw(base_size = 16)+
  ylab("[lactate] per cell") +
  facet_wrap(~background, scales = "free") +
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "F128+<br>M13+" = "#0F2080", "pOX38+" = "#A95AA1", "pOX38+<br>M13+" = "#0F2080")) +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown(), strip.text = element_markdown())

# part C - plasmid retention rates
# only pOX M13 to F128, pOX to F128M13 and pOXM13 to pOX significant
dates <- c("29July2024", "31July2024", "14August2024")

aggregated <- data.frame()
time_series <- data.frame()
plasmid_loss <- data.frame()

for (date in dates){
  source(here::here("data-generation", "experimental", "clean-plasmid-tecan-data.R"))
  temp <- aggregated_dataset
  temp1 <- all_plate_data_corrected
  temp2 <- percent_retained
  aggregated <- rbind(aggregated, temp)
  time_series <- rbind(time_series, temp1)
  plasmid_loss <- rbind(plasmid_loss, temp2)
}

partC <- plasmid_loss %>% filter(strain %in% c("E0224", "E0224 F+", "E0224 F+ M13+", "E0224 O+", "E0224 O+ M13+")) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 F+" ~ "F128+",
                            strain == "E0224 F+ M13+" ~ "F128+<br>M13+",
                            strain == "E0224 O+ M13+" ~ "pOX38+<br>M13+",
                            strain == "E0224 O+" ~ "pOX38+")) %>%
  filter(partner == "S0240") %>% filter(!well %in% c("E3", "E4")) %>%
  ggplot(aes(x = fct_reorder(strain, percent_retained, .desc = TRUE), y = percent_retained, color = strain)) +
  stat_summary(aes(fill = strain), fun = mean, shape = '-', size = 5, color = 'black')+
  geom_point(size = 2, stroke = 1.5)+
  ylab("percent plasmid +") +
  theme_bw(base_size = 16)+
  scale_color_manual(values = c("uninf" = "#F5793A", "F128+" = "#A95AA1", "F128+<br>M13+" = "#0F2080", "pOX38+" = "#A95AA1", "pOX38+<br>M13+" = "#0F2080")) +
  theme(axis.text = element_markdown(), axis.title.x = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_markdown())

# final figure
figure5 <- plot_grid(partA, plot_grid(partB, partC, ncol = 2, label_size = 26, labels = c("B", "C")), label_size = 26, labels = c("A", ""), ncol = 1)

png(here::here("figures", "final-figs", "imgs", "figure-5.png"), res = 300, width = 4000, height = 2000)
figure5
dev.off()




