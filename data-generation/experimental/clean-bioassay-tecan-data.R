#ATB
#Generate data for KO assays

#load packages
library("broom")
library("tidyverse")
library("readxl")
library("patchwork")
library("ggtext")
library("ggpubr")
library("cowplot")

#load functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

#load data
OD_path <- generate_paths("od", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
cfus <- generate_paths("cfu", date, "xlsx")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, OD_path), "OD")
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, CFP_path), "CFP")
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, YFP_path), "YFP")
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "auxotroph-bioassays", date, "plate_layout.csv")) %>%
  janitor::clean_names()
CFU <- load_pfu_data(here::here("experimental-data", "tecan-data","auxotroph-bioassays", date, cfus)) %>%
  filter(!is.na(cfu))

all_data <- OD %>%
  inner_join(., plate_layout, by = "well") %>%
  filter(strain != "none") %>%
  inner_join(., CFU %>% select(well, cfu), by = "well")
