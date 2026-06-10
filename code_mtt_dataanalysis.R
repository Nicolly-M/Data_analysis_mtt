#### Packages ####
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


##### Data preparation
dados <- read_excel("mtt.xlsx")

# Blank data 
blank_vals <- c(0.043, 0.044, 0.045)
blank_mean <- mean(blank_vals)

###### Long format ######
dados_long <- dados %>%
  pivot_longer(cols = starts_with("Rep"),
               names_to = "Replicate",
               values_to = "Abs")

##### Median and Standard Deviation #####
dados_resumo <- dados_long %>%
  group_by(Compound, Cell_line, Time, Concentration) %>%
  summarise(
    mean_abs = mean(Abs, na.rm = TRUE),
    sd_abs = sd(Abs, na.rm = TRUE),
    .groups = "drop"
  )

###### Median per group ######
controle_resumo <- dados_resumo %>%
  filter(Concentration == "Control") %>%
  select(Compound, Cell_line, Time, mean_abs) %>%
  rename(mean_control = mean_abs)

dados_norm <- dados_resumo %>%
  left_join(controle_resumo,
            by = c("Compound", "Cell_line", "Time"))

###### Normalization #####
dados_norm <- dados_norm %>%
  mutate(
    viab = (mean_abs - blank_mean) / (mean_control - blank_mean) * 100
  )

##### Viability Stardard deviation 
dados_norm <- dados_norm %>%
  mutate(
    sd_viab = sd_abs / (mean_control - blank_mean) * 100
  )

###### T-Test Cell line x Control
ttest_resultados <- dados_long %>%
  group_by(Compound, Cell_line, Time) %>%
  group_modify(~ {
    
    df <- .x
    
    controle_vals <- df %>%
      filter(Concentration == "Control") %>%
      pull(Abs)
    
    df %>%
      filter(Concentration != "Control") %>%
      group_by(Concentration) %>%
      summarise(
        p_value = t.test(Abs, controle_vals)$p.value,
        .groups = "drop"
      )
  }) %>%
  ungroup()

# Tuns control to zero
dados_norm <- dados_norm %>%
  mutate(
    Conc_num = ifelse(Concentration == "Control", 0, as.numeric(Concentration))
  )

###### Preparation to excel data [graph format] #######
sele_compound <- "STE"
sele_cell <- "U-251"


table <- dados_norm %>%
  filter(Compound == sele_compound, Cell_line == sele_cell) %>%
  
  left_join(
    ttest_resultados %>%
      select(Compound, Cell_line, Time, Concentration, p_value),
    by = c("Compound", "Cell_line", "Time", "Concentration")
  ) %>%
  
  select(Concentration, Time, viab, sd_viab, p_value) %>%
  
  pivot_wider(
    names_from = Time,
    values_from = c(viab, p_value),
    names_glue = "{.value}_{Time}h"
  )

write_xlsx(table, "C:/Users/seu_usuario/Desktop/table.xlsx")

