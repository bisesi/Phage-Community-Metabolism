#ATB
#Tecan and plaque assay data analysis
#Helper functions

generate_paths <- function(tecan_file_type, date, ending){
  return(paste0(tecan_file_type, "_", date, ".", ending))
}

import_tecan_data <- function(path, file_type = "OD"){
  data <- read.table(path, sep = ",", as.is = TRUE, row.names = 1)
  data <- t(data)
  data <- data.frame(data)
  names(data)[1:3] = c("cycle","seconds","temp")
  data <- data %>%
    gather(well, OD, -cycle, -seconds, -temp) %>%
    dplyr::mutate(hour = seconds / 60 / 60)
  if (file_type == "CFP"){
    data <- data %>% dplyr::rename(CFP = OD)
  }
  if (file_type == "YFP"){
    data <- data %>% dplyr::rename(YFP = OD)
  }
  if (file_type == "RFP"){
    data <- data %>% dplyr::rename(RFP = OD)
  }
  if (file_type == "JFP"){
    data <- data %>% dplyr::rename(JFP = OD)
  }
  return(data)
}

generate_baranyi_growth_data_from_OD <- function(file){
  output <- file %>%
    group_by(well) %>%
    summarize(model_fit = fit_baranyi(hour, log(OD), 5),
              fit_variable = c("growth_rate", "lag", "ymax", "y0"))
  return(output)
}

generate_baranyi_growth_data_from_CFP <- function(file){
  output <- file %>%
    group_by(well) %>%
    summarize(model_fit = fit_baranyi(hour, log(CFP), 5),
              fit_variable = c("growth_rate", "lag", "ymax", "y0"))
  return(output)
}

adjust_OD_vs_FP <- function(input, itx, fluorescence){
  output <- input %>% 
    dplyr::filter(interaction == itx & phage == "none") %>% 
    group_by(well) %>% 
    dplyr::mutate(adjusted = {{fluorescence}} - min({{fluorescence}})) %>% 
    ungroup() %>% 
    dplyr::select(c(cycle, well, adjusted, OD))
  cycle_max <- input %>% dplyr::select(cycle) %>% max()
  best_fit <- c()
  for (i in seq(5, cycle_max, by = 5)){
    best_fit[i] <- summary(lm(OD ~ adjusted, data = output %>% dplyr::filter(cycle < i)))$r.squared
  }
  return(c(which.max(best_fit), summary(lm(OD ~ adjusted, data = output))$coefficients[2,1]))
}

adjust_FP_vs_FP <- function(input, itx, fluorescence1, fluorescence2, cycle_number){
 output <- input %>% 
    dplyr::filter(interaction == itx & phage == "none") %>% 
    group_by(well) %>% 
    dplyr::mutate(adjusted1 = {{fluorescence1}} - min({{fluorescence1}})) %>%
    ungroup() %>% 
    dplyr::mutate(adjusted2 = {{fluorescence2}} - min({{fluorescence2}})) %>%
    ungroup() %>% 
    dplyr::select(c(cycle, well, adjusted1, adjusted2)) %>%
    dplyr::filter(cycle < cycle_number)
  return(summary(lm(adjusted1 ~ adjusted2, data = output))$coefficients[2,1])
}

adjust_FP_values <- function(input, YFP_bleed_value, CFP_bleed_value){
  wells <- c(input %>% dplyr::filter(interaction != "none") %>% dplyr::select(well) %>% unique())$well
  total_data <- input %>% dplyr::select(c(cycle,well,YFP,CFP)) %>% pivot_wider(names_from = well, values_from = c("CFP","YFP"))
  store <- list()
  for (i in wells){
    j <- which(wells == i)
    new_matrix <- total_data %>% dplyr::select(matches(i))
    dat <- apply(new_matrix, 1, function(x) solve(matrix(c(1,YFP_bleed_value, CFP_bleed_value,1), nrow=2, byrow = TRUE), c(x))) %>%
      as.data.frame() %>% 
      t() 
    colnames(dat) <- c("adjusted_CFP", "adjusted_YFP")
    store[[j]] <- dat
  }
  names(store) <- wells
  output <- map_df(store, ~as.data.frame(.x), .id="well") %>%
    group_by(well) %>%
    dplyr::mutate(cycle = c(1:n())) %>%
    ungroup()
  return(output)
}

load_pfu_data <- function(path){
  sheetnames <- excel_sheets(path)
  datalist <- lapply(sheetnames, read_excel, path = path)
  output <- datalist[[1]]
  return(output)
}

clean_pfu_data <- function(input){
  wide_data <- input %>% 
    dplyr::select(condition, interaction, phage, plate, pfu, well) %>%
    pivot_wider(names_from = plate, values_from = pfu) %>%
    dplyr::rename(pfu_S = S, pfu_E = E)
  starting_phage <- wide_data %>%
    dplyr::filter(interaction == "None" & condition == "Start") %>%
    dplyr::mutate(ratio = pfu_E / (pfu_E + pfu_S)) 
  output <- wide_data %>%
    dplyr::filter(condition != "Start") %>%
    dplyr::mutate(phi_doublings = case_when(phage == "Phi" ~ log(pfu_E / starting_phage %>% filter(phage == "Phi") %>% dplyr::select(pfu_E) %>% pull),
                                     phage == "Phi + P22" ~ log(pfu_E / starting_phage %>% filter(phage == "Phi + P22") %>% dplyr::select(pfu_E) %>% pull),
                                     phage == "P22" ~ 0,
                                     TRUE ~ 0)) %>%    
    dplyr::mutate(p22_doublings = case_when(phage == "P22" ~ log(pfu_S / starting_phage %>% filter(phage == "P22") %>% dplyr::select(pfu_S) %>% pull),
                                     phage == "Phi + P22" ~ log(pfu_S / starting_phage %>% filter(phage == "Phi + P22") %>% dplyr::select(pfu_S) %>% pull),
                                     phage == "Phi" ~ 0,
                                     TRUE ~ 0)) %>%
    dplyr::mutate(generalist_fitness = (((pfu_E) / (pfu_E + pfu_S)) / starting_phage %>% filter(phage == "Phi + P22") %>% dplyr::select(ratio) %>% pull)) %>% 
    dplyr::mutate(change_in_percent_generalist = (((pfu_E) / (pfu_E + pfu_S)) - starting_phage %>% filter(phage == "Phi + P22") %>% dplyr::select(ratio) %>% pull)) %>%
    dplyr::mutate_all(~replace(., is.nan(.), 0)) %>%
    dplyr::mutate_all(~replace(., is.infinite(.), log(1 / starting_phage %>% filter(phage == "Phi + P22") %>% dplyr::select(pfu_E) %>% pull))) %>%
    dplyr::mutate(interaction = case_when(interaction == "Smono" ~ "S Monoculture",
                                   interaction == "Emono" ~ "E Monoculture",
                                   interaction == "E Fac" ~ "E Facilitation",
                                   interaction == "Coop" ~ "Mutualism",
                                   interaction == "Comp" ~ "Competition",
                                   interaction == "Fac" ~ "Facilitation",
                                   interaction == "None" ~ "No cells")) %>%
    pivot_longer(cols = ends_with("doublings"),
                 names_to = "doubling_type",
                 values_to = "doublings") %>%
    dplyr::mutate(phage_type = case_when(doubling_type == "p22_doublings" ~ "Specialist phage",
                                  doubling_type == "phi_doublings" ~ "Generalist phage",
                                  TRUE ~ "NA")) %>%
    dplyr::mutate(phage_interaction = case_when(phage == "P22" ~ "Specialist phage",
                             phage == "Phi" ~ "Generalist phage",
                             phage == "Phi + P22" ~ "Phage Competition"))
  return(output)
}

adjust_curves_by_CFU <- function(data, model_frame, model) {
  
  nested_data <- data %>% group_by(strain, well) %>% nest()
  
  predictions <- data %>% 
    inner_join(., nested_data, by = c("well", "strain")) %>% 
    inner_join(., model_frame %>% select(strain, {{model}}), by = c("strain")) %>% 
    dplyr::mutate(model_pred = map2({{model}}, data, predict))
  
  output <- predictions %>% 
    dplyr::select(strain, well, data, model_pred) %>% 
    unnest(cols = c(data, model_pred)) %>% 
    dplyr::rename(prediction = model_pred)
  
  return(output)
}

import_spent_media <- function(date, replicate, path_extension){
  
  date <- date
  
  CFP_path <- generate_paths("cfp", date, "csv")
  CFP <- import_tecan_data(here::here("experimental-data", "tecan-data", path_extension, date, CFP_path), "CFP")
  plate_layout <- read_csv(here::here("experimental-data", "tecan-data", path_extension, date, "plate_layout.csv")) %>%
    janitor::clean_names()
  
  average_cfp <- CFP %>%
    inner_join(., plate_layout, by = "well") %>%
    filter(media %in% c("E0224", "E0224 F+", "E0224 F+ M13+")) %>%
    filter(cycle == min(cycle)) %>%
    dplyr::select(CFP, media) %>%
    dplyr::mutate(replicate = replicate)
  
  return(average_cfp)
}

change_letter <- function(strings, new_letter) {
  # Use gsub() to replace the first character with the new letter
  new_strings <- gsub("^[A-Za-z]", new_letter, strings)
  
  return(new_strings)
}

