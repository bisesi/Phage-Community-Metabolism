## ATB
# supplemental tables 1, 4, 6, 
# supplemental table 1 - comets metabolites per cell
# supplemental table 4 - in vitro metabolites per cell
# supplemental table 6 - unique gene products in plasmids
# supplemental table 7 - o and om growth and yield in mono

# load packages
library("chromatographR")
library("tidyverse")
library("readxl")
library("ggtext")
library("cowplot")
library("ggpubr")
library("rstatix")
library("broom")
library("genbankr")

#load helper functions
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

# supplemental table 1
setwd(here::here("fba-data", "figure-1", "comets"))
e_fluxes <- read_csv("E_fluxes_range.csv") %>% select(-c("...1", x, y)) %>% select(., c(lower_bound, optimized, cycle, contains("EX_"))) %>%
  pivot_longer(cols = -c("lower_bound", "optimized", "cycle"), values_to = "flux") 
e_met <- read_csv("E_metabolites_range.csv") %>% select(-c("...1")) %>% pivot_longer(cols = -c("lower_bound", "optimized", "cycle"), values_to = "concentration")
e <- read_csv("E_biomass_range.csv") %>% select(-c("...1")) %>% rename(E0 = iJO1366)

bounds <- e %>% mutate(optimized = ifelse(optimized == "virus" & lower_bound == 0, "host", optimized)) %>% pull(lower_bound) %>% unique()
plasmid_bound <- bounds[10]
virus_bound <- bounds[18]

yield <- e %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+<br>M13+", "F128+"))) %>%
  group_by(optimized) %>% filter(E0 == max(E0)) %>% slice_head(n = 1) %>% ungroup()

growth_rate <- e %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+<br>M13+", "F128+"))) %>%
  group_by(optimized) %>% summarize(model_fit = fit_baranyi(cycle, log(E0), 5),
                                   fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>% filter(fit_variable == "growth_rate") %>%
  ungroup()

met <- e_met %>% mutate(optimized = ifelse(optimized == "plasmid" & lower_bound == 0, "host", optimized)) %>% 
  filter(lower_bound %in% c(0, plasmid_bound, virus_bound)) %>%
  filter(!(lower_bound == 0 & optimized == "virus with plasmid")) %>%
  filter(name %in% c("ac_e", "glyclt_e", "hxa_e", "for_e", "co2_e", "etoh_e")) %>%
  group_by(optimized, name) %>% filter(concentration == max(concentration)) %>% slice_head(n = 1) %>% ungroup() %>%
  mutate(optimized = ifelse(optimized == "plasmid", "F128+", optimized),
         optimized = ifelse(optimized == "virus with plasmid", "F128+<br>M13+", optimized),
         optimized = ifelse(optimized == "host", "uninf", optimized),
         optimized = factor(optimized, levels = c("uninf", "F128+<br>M13+", "F128+")))

supp1 <- met %>% select(optimized, name, concentration) %>% pivot_wider(names_from = name, values_from = concentration) %>%
      inner_join(., yield %>% select(E0, optimized) %>% rename(biomass = E0), by = "optimized") %>%
  inner_join(., growth_rate %>% select(model_fit, optimized) %>% rename(growth_rate = model_fit), by = "optimized")

write.csv(supp1, file = here::here("figures", "final-figs", "tables", "supplemental-table-1.csv"))

# supplemental table 2
hplc <- read_csv(here::here("experimental-data", "hplc", "24April2024", "aggregated_signals_24April2024.csv")) %>% dplyr::select(-c("...1")) %>%
  pivot_longer(cols = c(-time), names_to = "sampleID") %>%
  dplyr::mutate(class = case_when(grepl("std", sampleID) == TRUE ~ "standard",
                                  grepl("blank", sampleID) == TRUE ~ "blank",
                                  sampleID %in% c("e1", "e2", "e3", "f1", "f2", "f4", "fm1", "fm3", "fm4") ~ "sample",
                                  TRUE ~ sampleID))
cfus <- read_csv(here::here("experimental-data", "hplc", "24April2024", "final-CFUs_spent-media_samples17April2024.csv")) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  dplyr::select(sampleID, CFUnotet) %>% rename(CFU = CFUnotet)

peak_integration <- data.frame()
for (i in unique(hplc$sampleID)){
  dat <- as.matrix(hplc %>% filter(sampleID == i) %>% dplyr::select(time, value) %>% rename(`210` = value))
  rownames(dat) <- dat[,1]
  peaks <- get_peaks(dat, lambdas = c(210))
  output <- as.data.frame(peaks$`1`$`210`) %>% dplyr::select(rt, height, area) %>%
    mutate(sampleID = i)
  peak_integration <- rbind(peak_integration, output)
}

lactate_acetate <- peak_integration %>% filter((rt > 19.5 & rt < 20.5) | (rt > 23 & rt < 24)) %>% 
  filter(!sampleID %in% c("blank11", "standard_media")) %>% mutate(type = ifelse(rt < 21, "lactate", "acetate")) %>%
  filter(!(type == "lactate" & sampleID %in% c("e1", "e2", "e3")))

concentration_model <- lactate_acetate %>% filter(grepl("_std", sampleID) == TRUE) %>%
  separate(sampleID, c("mM", "sampleType"), sep = "_") %>%
  mutate(mM = as.numeric(str_remove(mM, "mM"))) %>% dplyr::select(area, mM, type)

acetate <- lm(mM ~ area, data = concentration_model %>% filter(type == "acetate"))
lactate <- lm(mM ~ area, data = concentration_model %>% filter(type == "lactate"))

acetate_predictions <- lactate_acetate %>% filter(grepl("_std", sampleID) == FALSE) %>%
  filter(type == "acetate")
acetate_predictions$predictedmM <- predict(acetate, newdata = data.frame(area = acetate_predictions$area))

lactate_predictions <- lactate_acetate %>% filter(grepl("_std", sampleID) == FALSE) %>%
  filter(type == "lactate")
lactate_predictions$predictedmM <- predict(lactate, newdata = data.frame(area = lactate_predictions$area))

supp4 <- rbind(lactate_predictions %>% mutate(compound = "lactate"), acetate_predictions %>% mutate(compound = "acetate")) %>%
  inner_join(., cfus %>% mutate(sampleID = tolower(sampleID)), by = "sampleID") %>%
  mutate(per_cell = predictedmM / CFU) %>% dplyr::select(-c(height, type)) %>%
  mutate(condition = case_when(sampleID %in% c("f1", "f2", "f4") ~ "F128+",
                              sampleID %in% c("fm1", "fm4", "fm3") ~ "F128+ M13+",
                              sampleID %in% c("e1", "e2", "e3") ~ "uninf")) %>%
  group_by(condition) %>% mutate(biological_replicate = 1:n()) %>%
  ungroup() %>% dplyr::select(-c(sampleID))

write.csv(supp4, file = here::here("figures", "final-figs", "tables", "supplemental-table-4.csv"))

# supplemental table 6 - unique gene products in plasmids
f128 <- readGenBank(here::here("experimental-data", "sequencing", "f-plasmid-plasmidsaurus", "f128-full-sequence.gb"))
pox38 <- readGenBank(here::here("experimental-data", "sequencing", "f-plasmid-plasmidsaurus", "pox38-full-sequence.gb"))

cleaned_pox38 <- exons(pox38) %>% data.frame() %>% select(c("product", "gene", "function."))
products_o <- cleaned_pox38 %>% pull(product)

cleaned_f128 <- exons(f128) %>% data.frame() %>% select(c("product", "gene", "function."))
products_f <- cleaned_f128 %>% pull(product)

f_not_o <- products_f[which(!products_f %in% products_o)]
o_not_f <- products_o[which(!products_o %in% products_f)][-c(1,2)]

supp6 <- cleaned_f128 %>% filter(product %in% f_not_o) %>% select(product, gene) %>% mutate(genotype = "F128") %>%
  rbind(., cleaned_pox38 %>% filter(product %in% o_not_f) %>% select(product, gene) %>% mutate(genotype = "pOX38")) %>%
  rename(unique_products = product, gene_name = gene)

write.csv(supp6, file = here::here("figures", "final-figs", "tables", "supplemental-table-6.csv"))

# supplemental table 7
source(here::here("functions", "tecan-data-helper-functions.R"))
source(here::here("functions", "baranyi-helper-functions.R"))

date <- "29July2024"

OD_path <- generate_paths("od", date, "csv")
JFP_path <- generate_paths("jfp", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
RFP_path <- generate_paths("rfp", date, "csv")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, OD_path), "OD") %>% filter(!is.na(cycle))
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, CFP_path), "CFP") %>% filter(!is.na(cycle))
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, YFP_path), "YFP") %>% filter(!is.na(cycle))
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, RFP_path), "RFP") %>% filter(!is.na(cycle))
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "plasmid", date, "plate_layout.csv")) %>%
  janitor::clean_names()

july29 <- OD %>%
  inner_join(., CFP %>% select(cycle, well, CFP), by = c("well", "cycle")) %>%
  inner_join(., RFP %>% select(cycle, well, RFP), by = c("well", "cycle")) %>%
  inner_join(., YFP %>% select(cycle, well, YFP), by = c("well", "cycle")) %>%
  inner_join(., plate_layout, by = "well") %>%
  filter(strain != "none") %>%
  mutate(partner = ifelse(partner == "none", NA, partner))

date <- "7July2024"

OD_path <- generate_paths("od", date, "csv")
JFP_path <- generate_paths("jfp", date, "csv")
CFP_path <- generate_paths("cfp", date, "csv")
YFP_path <- generate_paths("yfp", date, "csv")
RFP_path <- generate_paths("rfp", date, "csv")

OD <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, OD_path), "OD") %>% filter(!is.na(cycle))
CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, CFP_path), "CFP") %>% filter(!is.na(cycle))
YFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, YFP_path), "YFP") %>% filter(!is.na(cycle))
RFP <- import_tecan_data(here::here("experimental-data", "tecan-data", "plasmid", date, RFP_path), "RFP") %>% filter(!is.na(cycle))
plate_layout <- read_csv(here::here("experimental-data", "tecan-data", "plasmid", date, "plate_layout.csv")) %>%
  janitor::clean_names()

july7 <- OD %>%
  inner_join(., CFP %>% select(cycle, well, CFP), by = c("well", "cycle")) %>%
  inner_join(., RFP %>% select(cycle, well, RFP), by = c("well", "cycle")) %>%
  inner_join(., YFP %>% select(cycle, well, YFP), by = c("well", "cycle")) %>%
  inner_join(., plate_layout, by = "well") %>%
  filter(strain != "none") %>%
  mutate(partner = ifelse(partner == "none", NA, partner))

maxyield <- rbind(july7 %>% mutate(date = "7july"), july29%>% mutate(date = "29july")) %>% 
  filter(strain %in% c("E0224", "E0224 pOX38 WT", "E0224 O+", "E0224 O+ M13+")) %>% 
  filter(is.na(partner)) %>%
  group_by(date, well) %>% filter(OD == max(OD)) %>% slice_head() %>%
  ungroup() %>% mutate(strain = ifelse(strain == "E0224 pOX38 WT", "E0224 O+", strain)) %>%
  group_by(strain) %>% summarize(mean = mean(OD), sd = sd(OD))

growthrate <- rbind(july7 %>% mutate(date = "7july"), july29 %>% mutate(date = "29july")) %>% 
  filter(strain %in% c("E0224", "E0224 pOX38 WT", "E0224 O+", "E0224 O+ M13+")) %>% 
  filter(is.na(partner)) %>%
  group_by(date, well, strain) %>% summarize(model_fit = fit_baranyi(hour, log(CFP), 5),
                                             fit_variable = c("growth_rate", "lag", "ymax", "y0")) %>%
  ungroup() %>% filter(fit_variable == "growth_rate") %>%
  mutate(strain = ifelse(strain == "E0224 pOX38 WT", "E0224 O+", strain)) %>%
  group_by(strain) %>% summarize(mean = mean(model_fit), sd = sd(model_fit))

supp7 <- rbind(growthrate %>% mutate(type = "growth rate"), maxyield %>% mutate(type = "max OD600")) %>%
  mutate(strain = case_when(strain == "E0224" ~ "uninf",
                            strain == "E0224 O+" ~ "pOX38+",
                            strain == "E0224 O+ M13+" ~ "pOX38+ M13+"))

write.csv(supp7, file = here::here("figures", "final-figs", "tables", "supplemental-table-7.csv"))



