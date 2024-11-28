#ATB
#Generate data for faciliation tecans

#load packages
library("tidyverse")
library("readxl")
library("patchwork")
library("ggtext")
library("ggpubr")
library("cowplot")
library("abdiv")

#load function
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

#import and clean in vitro data 
date <- date

OD_path <- generate_paths("od", date, "csv")
JFP_path <- generate_paths("jfp", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
RFP_path <- generate_paths("rfp", date, "csv")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "facilitation", date, OD_path), "OD")
JFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "facilitation", date, JFP_path), "JFP")
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "facilitation", date, CFP_path), "CFP")
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "facilitation", date, YFP_path), "YFP")
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "facilitation", date, "plate_layout.csv")) %>%
  janitor::clean_names()
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "facilitation", date, RFP_path), "RFP")

all_data <- OD %>%
  inner_join(., JFP %>% select(cycle, well, JFP), by = c("well", "cycle")) %>%
  inner_join(., CFP %>% select(cycle, well, CFP), by = c("well", "cycle")) %>%
  inner_join(., RFP %>% select(cycle, well, RFP), by = c("well", "cycle")) %>%
  inner_join(., YFP %>% select(cycle, well, YFP), by = c("well", "cycle")) %>%
  inner_join(., plate_layout, by = "well") %>%
  filter(strain != "none") %>%
  mutate(partner = ifelse(partner == "none", NA, partner)) 

#OD traces
uninfected <- all_data %>%
  filter(strain %in% c("E0224", "S0240")) %>%
  mutate(phage = "none", 
         interaction = case_when(strain == "E0224" & is.na(partner) ~ "Emono",
                                 strain == "S0240" ~ "Smono",
                                 TRUE ~ strain)) 

fplasmid <- all_data %>%
  filter(strain %in% c("E0224 F+", "S0240")) %>%
  mutate(phage = "none", 
         interaction = case_when(strain == "E0224 F+" & is.na(partner) ~ "Emono",
                                 strain == "S0240" ~ "Smono",
                                 TRUE ~ strain)) 

infected <- all_data %>%
  filter(strain %in% c("E0224 F+ M13+", "S0240")) %>%
  mutate(phage = "none", 
         interaction = case_when(strain == "E0224 F+ M13+" & is.na(partner) ~ "Emono",
                                 strain == "S0240" ~ "Smono",
                                 TRUE ~ strain)) 

subsetted <- list(uninfected = uninfected, fplasmid = fplasmid, infected = infected)
all_data_corrected <- data.frame()
for (i in 1:length(subsetted)){
  YFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Emono", YFP, CFP, subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% first())
  CFP_bleed <- subsetted[[i]] %>% adjust_FP_vs_FP("Smono", CFP, YFP, subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% first())
  all_tecan_adjusted_FP <- adjust_FP_values(subsetted[[i]], YFP_bleed_value = YFP_bleed, CFP_bleed_value = CFP_bleed)
  all_tecan_adjusted_OD <- all_tecan_adjusted_FP %>%
    mutate(E_corrected_OD = adjusted_CFP * subsetted[[i]] %>% adjust_OD_vs_FP("Emono", CFP) %>% last()) %>%
    mutate(S_corrected_OD = adjusted_YFP * subsetted[[i]] %>% adjust_OD_vs_FP("Smono", YFP) %>% last()) %>%
    inner_join(., subsetted[[i]] %>% select(cycle, well, OD, hour, JFP, RFP, CFP, YFP, strain, partner), by = c("well", "cycle"))
  all_data_corrected <- rbind(all_data_corrected, all_tecan_adjusted_OD)
}
