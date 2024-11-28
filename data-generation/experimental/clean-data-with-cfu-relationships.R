#ATB
#Generate cfu predictions from cfu fp relationships

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
#cfus <- generate_paths("cfu", date, "xlsx")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "mutualism", date, OD_path), "OD")
JFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "mutualism", date, JFP_path), "JFP")
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "mutualism", date, CFP_path), "CFP")
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "mutualism", date, YFP_path), "YFP")
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "mutualism", date, RFP_path), "RFP")
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "mutualism", date, "plate_layout.csv")) %>%
  janitor::clean_names()
#CFU <- load_pfu_data(here::here("experimental-data", "tecan-data","mutualism", date, cfus)) %>%
#filter(!is.na(cfu))

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

all_data_with_cfu <- all_data_corrected %>%
  rbind(., all_data %>% filter(strain == "M0104") %>% 
          mutate(adjusted_CFP = 0, adjusted_YFP = 0, E_corrected_OD = 0, S_corrected_OD = 0) %>%
          select(-c(temp, seconds))) %>%
  inner_join(., CFU %>% select(well, plate, cfu), by = "well")

CFP_split <-  all_data_with_cfu %>% filter(cycle == max(cycle)) %>% filter(plate == "E") %>% select(-plate) %>% rename(CFU = cfu)
YFP_split <-  all_data_with_cfu %>% filter(cycle == max(cycle)) %>% filter(plate == "S") %>% select(-plate) %>% rename(CFU = cfu)
RFP_split <-  all_data_with_cfu %>% filter(cycle == max(cycle)) %>% filter(plate == "M") %>% select(-plate) %>% rename(CFU = cfu)

model_details_cfp <- CFP_split %>%
  group_by(strain) %>%
  do(mod = lm(CFU~CFP, data = .)) %>% ungroup() %>%
  mutate(tidy_mod = map(mod, broom::tidy), 
         glance_mod = map(mod, broom::glance),
         augment_mod = map(mod, broom::augment),
         rsq_mod = glance_mod %>% map_dbl('r.squared'),
         slope_mod = tidy_mod %>% map_dbl(function(x) x$estimate[2]))

model_details_rfp <- RFP_split %>%
  do(mod = lm(CFU~RFP, data = .)) %>% ungroup() %>%
  mutate(tidy_mod = map(mod, broom::tidy), 
         glance_mod = map(mod, broom::glance),
         augment_mod = map(mod, broom::augment),
         rsq_mod = glance_mod %>% map_dbl('r.squared'),
         slope_mod = tidy_mod %>% map_dbl(function(x) x$estimate[2])) %>%
  mutate(strain = "M0104")

model_details_yfp <- YFP_split %>%
  do(mod = lm(CFU~YFP, data = .)) %>% ungroup() %>%
  mutate(tidy_mod = map(mod, broom::tidy), 
         glance_mod = map(mod, broom::glance),
         augment_mod = map(mod, broom::augment),
         rsq_mod = glance_mod %>% map_dbl('r.squared'),
         slope_mod = tidy_mod %>% map_dbl(function(x) x$estimate[2]))%>%
  mutate(strain = "S0240")

nested_data_E <- all_data_with_cfu %>% filter(plate == "E") %>% group_by(strain, well) %>% nest()

predictions_E <- nested_data_E %>% 
  full_join(., model_details_cfp %>% select(strain, mod), by = c("strain")) %>% 
  mutate(model_pred = map2(mod, data, predict))%>% 
  select(strain, well, data, model_pred) %>% 
  unnest(cols = c(data, model_pred)) %>% 
  rename(prediction = model_pred)

nested_data_S <- all_data_with_cfu %>% filter(plate == "S") %>% group_by(well) %>% nest()

predictions_S <- nested_data_S %>% 
  mutate(model_details_yfp %>% select(mod)) %>%
  mutate(model_pred = map2(mod, data, predict))%>% 
  select(well, data, model_pred) %>% 
  unnest(cols = c(data, model_pred)) %>% 
  rename(prediction = model_pred)

nested_data_M <- all_data_with_cfu %>% filter(plate == "M") %>% group_by(well) %>% nest()

predictions_M <- nested_data_M %>% 
  mutate(model_details_rfp %>% select(mod)) %>%
  mutate(model_pred = map2(mod, data, predict))%>% 
  select(well, data, model_pred) %>% 
  unnest(cols = c(data, model_pred)) %>% 
  rename(prediction = model_pred)