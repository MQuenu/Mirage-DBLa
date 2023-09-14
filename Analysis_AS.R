library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(entropy)
library(gridExtra)

setwd("~/DBLa/R_analyses")

rm(list = ls())

## format the input table, only do it once

ASprocessed_df <- read.csv("ASformatted_table.csv") %>%
  select(-X) %>%
  mutate(NTS = substr(Domain1, start = 1, stop = 4),
         Read_count = as.numeric(Read_count),
         infection_group = substr(sample, start = 1, stop = 2),
         host = substr(sample, start = 1, stop = 4)) %>% 
  group_by(sample) %>%
  mutate(total_reads = sum(Read_count)) %>%
  ungroup() %>%
  mutate(normalized_reads = Read_count / total_reads,
         Dominant_tags = ifelse(normalized_reads > 0.1, "dominant", "low_expression")) %>%
  group_by(sample) %>%
  mutate(entropy = entropy(normalized_reads)) %>%
  ungroup() %>%
  group_by(host) %>%
  group_by(Dominant_tags, .add = TRUE) %>%
  group_by(DBL_tag, .add = TRUE) %>%
  mutate(recurrence_step1 = ifelse(n() > 1, "recurrent", "single")) %>%
  mutate(recurrence = ifelse(recurrence_step1 == 'recurrent' & Dominant_tags == 'dominant', 'recurrent', 'single')) %>%
  ungroup() %>%
  select(-recurrence_step1)
  
write.csv(ASprocessed_df, "ASformatted_table.csv")

#### Plots of read counts per sample

## DC01 

DC01_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC01_m1" | sample == "DC01_m3" |  sample == "DC01_m4c" | 
           sample == "DC01_m5" | sample == "DC01_m6")

plot_DC01 = ggplot(DC01_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC01

## DC02

DC02_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC02_m1" | sample == "DC02_m2" |  sample == "DC02_m3" |
           sample == "DC02_m4" | sample == "DC02_m5" | sample == "DC02_m6")

plot_DC02 = ggplot(DC02_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC02

## DC03

DC03_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC03_m2" | sample == "DC03_m4" | sample == "DC03_m5" | 
         sample == "DC03_m6")

plot_DC03 = ggplot(DC03_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC03

## DC04

DC04_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC04_m2" | sample == "DC04_m3" |  sample == "DC04_m4c" |
           sample == "DC04_m5" | sample == "DC04_m6")

plot_DC04 = ggplot(DC04_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC04

## DC05

DC05_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC05_m1" | sample == "DC05_m2" |  sample == "DC05_m3" |
           sample == "DC05_m4" | sample == "DC05_m5" | sample == "DC05_m6" )

plot_DC05 = ggplot(DC05_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC05

## DC06

DC06_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC06_m1" | sample == "DC06_m6")

plot_DC06 = ggplot(DC06_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC06

## DC07

DC07_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC07_m1" | sample == "DC07_m2" |  sample == "DC07_m3" |
         sample == "DC07_m5" | sample == "DC07_m6")

plot_DC07 = ggplot(DC07_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC07

## DC08

DC08_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC08_m2" | sample == "DC08_m3" |  sample == "DC08_m4" |
           sample == "DC08_m5" | sample == "DC08_m6")

plot_DC08 = ggplot(DC08_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC08

## DC09

DC09_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC09_m2" | sample == "DC09_m3" |  sample == "DC09_m4" |
           sample == "DC09_m5 ")

plot_DC09 = ggplot(DC09_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC09

## DC10

DC10_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC10_m2" | sample == "DC10_m3")

plot_DC10 = ggplot(DC10_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC10

## DC11

DC11_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC11_m1" | sample == "DC11_m2" |  sample == "DC11_m4" |
           sample == "DC11_m5" | sample == "DC11_m6")

plot_DC11 = ggplot(DC11_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC11

## DC12

DC12_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC12_m1" | sample == "DC12_m2" |  sample == "DC12_m4" |
           sample == "DC12_m5" | sample == "DC12_m6")

plot_DC12 = ggplot(DC12_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC12

## DC13

DC13_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC13_m1" | sample == "DC13_m2" |  sample == "DC13_m3" |
           sample == "DC13_m4" | sample == "DC13_m5")

plot_DC13 = ggplot(DC13_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC13

## DC14

DC14_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC14_m1" | sample == "DC14_m2" |  sample == "DC14_m3")

plot_DC14 = ggplot(DC14_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC14

### DC15

DC15_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC15_m2")

plot_DC15 = ggplot(DC15_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC15

### DC15

DC16_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC16_m1")

plot_DC16 = ggplot(DC16_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC16


## WS

WS1_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "WS05" | sample == "WS06" |  sample == "WS07" |
           sample == "WS08"| sample == "WS09"| sample == "WS10")

plot_WS1 = ggplot(WS1_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_WS1


WS2_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "WS11" |  sample == "WS12" |
          sample == "WS13"| sample == "WS14"| sample == "WS15")

plot_WS2 = ggplot(WS2_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_WS2

## WA

WA1_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "WA05" |  sample == "WA06" |
           sample == "WA07"| sample == "WA08"| sample == "WA09")

plot_WA1 = ggplot(WA1_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_WA1

WA2_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "WA10" |  sample == "WA11" |
           sample == "WA12"| sample == "WA13"| sample == "WA18")

plot_WA2 = ggplot(WA2_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_WA2

WA3_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "WA36" |  sample == "WA37" |
           sample == "WA38"| sample == "WA40"| sample == "WA46")

plot_WA3 = ggplot(WA3_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_WA3

### Histogram of DBLa expressions

df <- read.csv("ASformatted_table.csv", header = TRUE) 

hist_normalized <- ggplot(df, aes(x = normalized_reads)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = 0.1, linetype="dotted", color = "red", size=1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hist_normalized

### Look for shared DBLa

rm(list = ls())

raw_dataframe <- read.csv("ASformatted_table.csv", header = TRUE)

count_dataframe <- raw_dataframe %>%
  group_by(DBL_tag) %>%
  summarize(count = n())

# look through count_dataframe interactively

## Entropy 

rm(list = ls())

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) 

entropy_dataframe <- aggregate(entropy ~ sample + infection_group, data = dataframe, FUN = mean)

plot_entropy <- ggplot(entropy_dataframe, aes(x = infection_group, y = entropy, fill = infection_group)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_brewer() +
  theme_classic()

plot_entropy

### glm between the three infection groups

null_model <- glm(entropy ~ 1, data = entropy_dataframe)
linear_model <- glm(entropy ~ infection_group, data = entropy_dataframe)

lrtest <- anova(null_model, linear_model)
teststat <- lrtest[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

### Looking for tags that would be expressed twice in the same sample (false clustering)

setwd("~/DBLa/R_analyses/false_clustering")

doubled_dataframe <- read.csv("ASformatted_table.txt",header = TRUE) %>%
  mutate(sample = substr(Sample_ID, start = 1, stop = 7)) %>%
  group_by(sample) %>%
  group_by(DBL_tag, add = TRUE) %>%
  mutate(double_tag = ifelse(n() > 1, "doubled", "single")) %>%
  ungroup() %>%
  filter(Read_count > 1)

#count the number of 'doubled' genes

count_double_dataframe <- doubled_dataframe %>%
  count(double_tag)

### Looking at Domain composition of expressed var genes

### functions

rm(list = ls())

Extract_expression_domain <- function(entry_dataframe){
for (column in colnames(entry_dataframe)) {
  domain <- substr(column, start = 1, stop = 6)
  if (domain == "Domain"){
    summary_df <- entry_dataframe %>%
      group_by(.data[[column]]) %>%
      summarize(total_quantity = sum(Read_count)) 
    colnames(summary_df)[1] <- "Expressed_Domain"
    if (exists("combined_data")) {
            combined_data <- full_join(combined_data, summary_df, by = "Expressed_Domain") %>%
        mutate(Total_read_count = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
        select(Expressed_Domain,Total_read_count)
        }
    else {
      combined_data <- summary_df
    }
    }
}
  combined_data <- combined_data %>%
    mutate(first_character = substr(Expressed_Domain,1,1),
           Domain_subclass = ifelse(first_character == "D" | first_character == "N",
                                 substr(Expressed_Domain, start = 1, stop = 4),
                                 substr(Expressed_Domain, start = 1, stop = 5)),
           Domain_class = ifelse(first_character == "D" | first_character == "N",
                                 substr(Expressed_Domain, start = 1, stop = 3),
                                 substr(Expressed_Domain, start = 1, stop = 4))) %>%
    select(Expressed_Domain, Total_read_count, Domain_class, Domain_subclass) %>%
    filter(complete.cases(.)) %>%
    filter(Domain_class != "NTS")
return(combined_data)
}

plot_domains <- function(entry_domain_dataframe) {
  palette <- brewer.pal(11, "Spectral")
  plot <- ggplot(entry_domain_dataframe, aes(x = ID, y = Total_read_count, fill = Domain_subclass)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = palette) +
    theme_classic()
  return(plot)
}

### DC01

DC01_m1_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "DC01_m1")

DC01_m3_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "DC01_m3")

DC01_m5_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "DC01_m5")

DC01_m6_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "DC01_m6")

DC01_m1_domains <- Extract_expression_domain(DC01_m1_dataframe) %>%
  mutate(ID = "DC01_m1")

DC01_m3_domains <- Extract_expression_domain(DC01_m3_dataframe) %>%
  mutate(ID = "DC01_m3")

DC01_m5_domains <- Extract_expression_domain(DC01_m5_dataframe) %>%
  mutate(ID = "DC01_m5")

DC01_m6_domains <- Extract_expression_domain(DC01_m6_dataframe) %>%
  mutate(ID = "DC01_m6")

plot_m1 <- plot_domains(DC01_m1_domains)

plot_m3 <- plot_domains(DC01_m3_domains)

plot_m5 <- plot_domains(DC01_m5_domains)

plot_m6 <- plot_domains(DC01_m6_domains)

grid.arrange(plot_m1, plot_m3, plot_m5, plot_m6, ncol=4, nrow = 1)
