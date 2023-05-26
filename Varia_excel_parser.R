#!/usr/bin/Rscript

rm(list = ls())

library(readxl)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_prompt <- args[1]

spreadsheet <- read_excel(input_prompt, sheet = 2)

start_index_1 <- 3
end_index_1 <- which(spreadsheet$...1 == 'Interpreded Results') %>%
  as.numeric()

start_index_2 <- which(spreadsheet$...1 == "Suggested Composition") %>%
  as.numeric()

table_1 <- spreadsheet[start_index_1:(end_index_1-6),] %>%
  setNames(.[1,]) %>%
  .[-1, ] %>%
  .[,1:6]
  
table_1 <- table_1[!is.na(table_1$Read_count),]

table_2 <- spreadsheet[(start_index_2+1):ncol(spreadsheet),] %>%
  setNames(.[1,]) %>%
  .[-1,] %>%
  .[,1:12]

table_2 <- table_2[!is.na(table_2$Read_count),] %>%
  .[,-2]

table_output <- left_join(table_1, table_2, by = "Sample_ID") %>%
  rename("Best_Blast_Hit" = "Best Blast hit",
         "VarDb_hits" = "# hits in varDB",
         "Domain1" = "Position D1 (91% expected accuracy)",
         "Domain2" = "Position D2 (99% expected accuracy)",
         "Domain3" = "Position D3 (96% expected accuracy)",
         "Domain4" = "Position D4 (95% expected accuracy)",
         "Domain5" = "Position D5 (83% expected accuracy)",
         "Domain6" = "Position D6 (78% expected accuracy)",
         "Domain7" = "Position D7 (77% expected accuracy)",
         "Domain8" = "Position D8 (69% expected accuracy)",
         "Domain9" = "Position D9",
         "Domain10" = "Position D10")

write.csv(table_output, file = "Sequences_table.csv", row.names = FALSE)