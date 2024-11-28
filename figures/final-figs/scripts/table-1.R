# ATB
# table 1

#load packages
library("tidyverse")
library("readxl")
library("patchwork")
library("ggtext")
library("ggpubr")
library("cowplot")
library("chromatographR")
library("broom")
library("genbankr")

# Table 1 - gene products on F but not on O and on O but not on F
f128 <- readGenBank(here::here("experimental-data", "sequencing", "f-plasmid-plasmidsaurus", "f128-full-sequence.gb"))
pox38 <- readGenBank(here::here("experimental-data", "sequencing", "f-plasmid-plasmidsaurus", "pox38-full-sequence.gb"))

cleaned_pox38 <- exons(pox38) %>% data.frame() %>% select(c("product", "gene", "function."))
products_o <- cleaned_pox38 %>% pull(product)

cleaned_f128 <- exons(f128) %>% data.frame() %>% select(c("product", "gene", "function."))
products_f <- cleaned_f128 %>% pull(product)

f_not_o <- products_f[which(!products_f %in% products_o)]
o_not_f <- products_o[which(!products_o %in% products_f)][-c(1,2)]

table3 <- cleaned_f128 %>% filter(product %in% f_not_o) %>% select(product, gene) %>% mutate(genotype = "F128") %>%
  rbind(., cleaned_pox38 %>% filter(product %in% o_not_f) %>% select(product, gene) %>% mutate(genotype = "pOX38")) %>%
  rename(unique_products = product, gene_name = gene)

write.csv(table3, file = here::here("figures", "final-figs", "tables", "table-3.csv"))





