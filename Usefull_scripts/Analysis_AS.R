library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(entropy)
library(gridExtra)
library(emmeans)
library(lme4)
library(stringr)
library(lmtest)

setwd("~/DBLa/R_analyses")

rm(list = ls())

#### format the input table, only do it once ####

ASprocessed_df <- read.csv("ASformatted_table.csv") %>%
  mutate(NTS = substr(Domain1, start = 1, stop = 4),
         Read_count = as.numeric(Read_count),
         infection_group = substr(sample, start = 1, stop = 2),
         host = substr(sample, start = 1, stop = 4),
         timepoint = ifelse(nchar(sample) == 7, substr(sample, start = 6, stop = 7), NA)) %>% 
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
  select(-recurrence_step1) #%>%
  #select(-X)

### Tables to left-join

parasitaemia_df <- read.csv("parasitaemia.csv") %>%
  rename(sample = Alias_updated,
         parasitaemia = Venous.parasitaemia) %>%
  select(-seq)

Age_df <- read.csv("host_age.csv") %>%
  rename(sample = Alias_updated)

ASprocessed_df <- read.csv("ASformatted_table.csv")

Fws_df <- read.csv("fws_table.csv") %>%
  rename(sample = Alias_updated)

updated_df <- left_join(ASprocessed_df, Fws_df, by = "sample")
updated_df <- left_join(ASprocessed_df, Age_df, by = "sample")

## Remove the contaminants from the list of AS formatted dataframe, also DC06

contaminants_tags <- read.table("list_false_tags.txt", header = FALSE)
list_contaminants <- contaminants_tags$V1

ASdataframe <- read.csv("ASformatted_table.csv") %>%
  filter(!DBL_tag %in% list_contaminants) %>%
  filter(host != 'DC06')

## Add the cUPS classifiaction

cUPS_df <- read.csv("cUPS.csv") %>%
  mutate(
    ups = case_when(
      A > 0.9 ~ "A",
      B > 0.9 ~ "B",
      C > 0.9 ~ "C",
      TRUE ~ "NA"
    )
  ) %>%
  rename(DBL_tag = X)

ASdataframe <- read.csv("ASformatted_table.csv") %>%
  select(-X, -X.1, -X.2) 

ASdataframe <- left_join(ASdataframe, cUPS_df, by = "DBL_tag") 

ASdataframe <- ASdataframe %>% 
  rename(PupsA = A,
         PupsB = B,
         PupsC = C)

ASdataframe <- read.csv("ASformatted_table.csv")%>%
  mutate(ups = ifelse(PupsA > 0.5, yes = "A", no = "BC"))

### Adding the tag_host_combination and linking

df_pacbio <- read.csv("Tag_matches_Pacbio.csv", header = TRUE)
ASdataframe <- read.csv("ASformatted_table.csv") %>%
  mutate(tag_host_combination = paste0(host, '_', DBL_tag))

dataframe <- left_join(ASdataframe, df_pacbio, by = "tag_host_combination") %>%
  select(-X, -Host, -tag) %>%
  filter(host == 'DC01' | host == 'DC02'| host == 'DC04'| host =='DC07'| 
        host == 'DC08' | host == 'WA10' | host == 'WA11' | host == 'WS05' |
        host == 'WS06' | host == 'WS07' | host == 'WS08' | host == 'WS11' |
        host == 'WS12' | host == 'WS13' | host == 'WS15'| host == 'DC03') %>%
  group_by(sample) %>%
  mutate(prop_matches = 1- mean(is.na(pacbio_id))) 


write.csv(dataframe, "pacbio_dataframe.csv")
         
## Careful ! the next line will override data

###write.csv(ASdataframe, "ASformatted_table.csv", col.names = FALSE)

#### Barplots of read counts per sample ####

#### function to plot the graph, takes a dataframe with all timepoint from a host as entry

plot_multiple_timepoints <- function(entry_dataframe) {
  
  entry_dataframe <- entry_dataframe %>%
    arrange(timepoint) %>%
    mutate(DBL_tag = factor(DBL_tag, levels = unique(DBL_tag))) %>%
    arrange(timepoint, desc(normalized_reads))
  
  plot_timepoints = ggplot(entry_dataframe, aes(x = DBL_tag, y = normalized_reads)) +
    geom_col(width=0.7) +
    facet_wrap(~ sample, nrow = 10) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Relative var gene expression")
  
  return(plot_timepoints)
}

plot_multiple_timepoints_group <- function(entry_dataframe) {
  
  entry_dataframe <- entry_dataframe %>%
    arrange(timepoint) %>%
    mutate(DBL_tag = factor(DBL_tag, levels = unique(DBL_tag))) %>%
    arrange(timepoint, desc(normalized_reads))
  
  plot_timepoints = ggplot(entry_dataframe, aes(x = DBL_tag, y = normalized_reads, fill = NTS)) +
    geom_col(width=0.7) +
    facet_wrap(~ sample, nrow = 10) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Relative var gene expression")
  
  return(plot_timepoints)
}

#### All graphs per sample:

## DC01 

DC01_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC01_m1" | sample == "DC01_m3" |  sample == "DC01_m4c" | 
           sample == "DC01_m5" | sample == "DC01_m6") %>%
  arrange(timepoint)

plot_DC01 <- plot_multiple_timepoints_group(DC01_dataframe)

plot_DC01

## DC02

DC02_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC02_m1" | sample == "DC02_m2" |  sample == "DC02_m3" |
           sample == "DC02_m4" | sample == "DC02_m5" | sample == "DC02_m6")

plot_DC02 = plot_multiple_timepoints_group(DC02_dataframe)

plot_DC02

## DC03

DC03_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC03_m2" | sample == "DC03_m4" | sample == "DC03_m5" | 
         sample == "DC03_m6")

plot_DC03 = plot_multiple_timepoints(DC03_dataframe)
  
plot_DC03

## DC04

DC04_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC04_m2" | sample == "DC04_m3" |  sample == "DC04_m4c" |
           sample == "DC04_m5" | sample == "DC04_m6")

plot_DC04 = plot_multiple_timepoints_group(DC04_dataframe)

plot_DC04

## DC01 + DC04, same frame

DC01DC04_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC01_m1" | sample == "DC01_m3" |  sample == "DC01_m4c" | 
           sample == "DC01_m5" | sample == "DC01_m6" |sample == "DC04_m2" | sample == "DC04_m3" |  sample == "DC04_m4c" |
           sample == "DC04_m5" | sample == "DC04_m6")

plot_DC01DC04 <- plot_multiple_timepoints_group(DC01DC04_dataframe)

plot_DC01DC04

## DC05

DC05_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC05_m1" | sample == "DC05_m2" |  sample == "DC05_m3" |
           sample == "DC05_m4" | sample == "DC05_m5" | sample == "DC05_m6" )

plot_DC05 = plot_multiple_timepoints_group(DC05_dataframe)

plot_DC05

## DC06

DC06_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC06_m1" | sample == "DC06_m6")

plot_DC06 = plot_multiple_timepoints(DC06_dataframe)

plot_DC06

## DC07

DC07_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC07_m1" | sample == "DC07_m2" |  sample == "DC07_m3" |
         sample == "DC07_m5" | sample == "DC07_m6")

plot_DC07 = plot_multiple_timepoints(DC07_dataframe)

plot_DC07

## DC08

DC08_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC08_m2" | sample == "DC08_m3" |  sample == "DC08_m4" |
           sample == "DC08_m5" | sample == "DC08_m6")

plot_DC08 = plot_multiple_timepoints(DC08_dataframe)

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

plot_DC11 = plot_multiple_timepoints(DC11_dataframe)

plot_DC11

## DC12

DC12_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC12_m1" | sample == "DC12_m2" |  sample == "DC12_m4" |
           sample == "DC12_m5" | sample == "DC12_m6")

plot_DC12 = plot_multiple_timepoints(DC12_dataframe)

plot_DC12

## DC13

DC13_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC13_m1" | sample == "DC13_m2" |  sample == "DC13_m3" |
           sample == "DC13_m4" | sample == "DC13_m5")

plot_DC13 = plot_multiple_timepoints(DC13_dataframe)

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

#### Piecharts plots ####

AS_dataframe <- read.csv("ASformatted_table.csv") 

my_palette <- c("#FFC20A", "#0C7BDC","#3498DB")

plot_piecharts = ggplot(AS_dataframe, aes(x = "", y = normalized_reads, fill = ups)) +
  geom_col() +
  coord_polar(theta = "y") +
  labs(title = NULL)+
  facet_wrap(~ sample) +
  scale_fill_manual(values = my_palette) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#FFFFFF"),
        plot.background = element_rect(fill = "#FFFFFF"),
        legend.background = element_rect(fill = "#FFFFFF")) 

plot_piecharts


#### Tests frequencies of A vs BC #####

rm(list = ls())

df_ups_pacbio <- read.table("ups_varIDs.txt", header = FALSE)

df_ups_pacbio <- read.csv("ups_varIDs.txt", header = FALSE) %>%
  rename(gene_id = V1, PacBioUTS = V2) %>%
  filter(PacBioUTS != 'other') %>%
  mutate(genome = str_split(gene_id, "_", simplify = TRUE)[, 1]) %>%
  mutate(boolean_group = ifelse(PacBioUTS == 'A', yes = 'A', no = 'BC')) %>%
  group_by(genome) %>%
  count(boolean_group) %>%
  mutate(frequency = n / sum(n)) 

write.csv(df_ups_pacbio,'frequency_upsA_pacbio.csv')

df_DBLa_upsA <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  select(sample, normalized_reads, ups) %>%
  group_by(sample, ups) %>%
  mutate(cumulative_frequency = sum(normalized_reads)) %>%
  distinct(ups, .keep_all = TRUE)

write.csv(df_DBLa_upsA,'frequency_upsA_DBLa.csv')

## Stats

df_glm <- read.csv("frequency_upsBC.csv") %>%
  rename(freqB = frequency) %>%
  mutate(freqA = 1 - freqB)

### WA group

df_glm_WA <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaWA') 
  
full_model <- glm(freqA ~ Group, data = df_glm_WA)
null_model <- glm(freqA ~ 1, data = df_glm_WA)
summary(full_model)
  
lr_test <- anova(null_model,full_model)
print(lr_test)

teststat <- lr_test[2, "Deviance"]
teststat
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)
pvalue

wald_test <- waldtest(full_model)

em <- emmeans(full_model, "Group")
contrast(em, "pairwise", adjust = "Tukey")

## WS group

df_glm_WS<- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaWS')

full_model <- glm(freqA ~ Group, data = df_glm_WS)
null_model <- glm(freqA ~ 1, data = df_glm_WS)

lr_test <- anova(null_model,full_model)
print(lr_test)

teststat <- lr_test[2, "Deviance"]
teststat
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

wald_test <- waldtest(full_model)
print(wald_test)

summary(full_model)

## DC group

df_glm_DC<- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC')

full_model <- glm(freqA ~ Group, data = df_glm_DC)
null_model <- glm(freqA ~ 1, data = df_glm_DC)

lr_test <- anova(null_model,full_model)
print(lr_test)

teststat <- lr_test[2, "Deviance"]
teststat
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue
## Boxplot

boxplot_upsA <- ggplot(df_glm, aes(x = Group, y = freqA)) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()

boxplot_upsA

### Within each DC host

df_glm_DC01 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC01')

m1 <- glm(freqA ~ host, data = df_glm_DC01)
summary(m1)

df_glm_DC02 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC02')

m2 <- glm(freqA ~ host, data = df_glm_DC02)
summary(m2)

df_glm_DC04 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC04')

m4 <- glm(freqA ~ host, data = df_glm_DC04)
summary(m4)

df_glm_DC05 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC05')

m5 <- glm(freqA ~ host, data = df_glm_DC05)
summary(m5)

df_glm_DC07 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC07')

m7 <- glm(freqA ~ host, data = df_glm_DC07)
summary(m7)

df_glm_DC08 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC08')

m8 <- glm(freqA ~ host, data = df_glm_DC08)
summary(m8)

df_glm_DC11 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC11')

m11 <- glm(freqA ~ host, data = df_glm_DC11)
summary(m11)

df_glm_DC12 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC12')

m12 <- glm(freqA ~ host, data = df_glm_DC12)
summary(m12)

df_glm_DC13 <- df_glm %>%
  filter(Group == 'PacBio' | Group == 'dblaDC') %>%
  mutate(host = substr(sample_id, start = 1, stop = 4)) %>%
  mutate(host = ifelse(Group == "dblaDC", yes = host, no = "PacBio")) %>%
  filter(host == 'PacBio' | host == 'DC13')

m13 <- glm(freqA ~ host, data = df_glm_DC13)
summary(m13)


### Histogram of DBLa expressions ####

df <- read.csv("ASformatted_table.csv", header = TRUE) 

hist_normalized <- ggplot(df, aes(x = normalized_reads)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = 0.1, linetype="dotted", color = "red", size=1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hist_normalized

### Count the number of tags without varia prediction
df <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  distinct(DBL_tag, .keep_all = TRUE)

na_count <- sum(is.na(df$Domain1))

### Count the number of recurrent tags
df <- read.csv("ASformatted_table.csv", header = TRUE) #%>%
  #filter(infection_group == "DC") %>%
  distinct(DBL_tag, .keep_all = TRUE)

count_reccurent <- df %>%
  count(recurrence) #%>%
  filter(recurrence == "reccurent")

### Look for shared DBLa

rm(list = ls())

raw_dataframe <- read.csv("ASformatted_table.csv", header = TRUE)

count_dataframe <- raw_dataframe %>%
  group_by(DBL_tag) %>%
  summarize(count = n()) 

count_of_ones <- count_dataframe %>%
  filter(count == 1) %>%
  nrow()

# look through count_dataframe interactively

##### Parasitaemia plots ######

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(infection_group == 'DC') %>%
  mutate(timepoint_continuous = gsub("m", "", timepoint)) %>%
  mutate(timepoint_continuous = as.numeric(timepoint_continuous)) %>%
  mutate(log_parasitaemia = log(parasitaemia)) %>%
  filter(!is.infinite(log_parasitaemia)) %>%
  group_by(timepoint) 

parasitaemia_dataframe_timepoints_data <- aggregate(parasitaemia ~ sample + infection_group + host + timepoint, data = dataframe, FUN = mean) %>%
  group_by(timepoint) %>%
  mutate(median_parasitaemia = median(parasitaemia)) %>%
  ungroup() %>%
  mutate(timepoint_continuous = gsub("m", "", timepoint)) %>%
  mutate(timepoint_continuous = as.numeric(timepoint_continuous)) %>%
  arrange(sample)
    
timepoints_parasitaemia <- ggplot(parasitaemia_dataframe_timepoints_data) +
  geom_point(aes(x = timepoint, y = parasitaemia, group = host, color = host), size = 2) + 
  geom_line(aes(x = timepoint, y = parasitaemia, group = host, color = host), linetype = 2) +
  scale_y_continuous(trans='log10') +
  theme_classic()

timepoints_parasitaemia <- ggplot(parasitaemia_dataframe_timepoints_data) +
  geom_point(aes(x = timepoint, y = parasitaemia, group = host, color = host), alpha = 0.5, size = 0.5) + 
  geom_line(aes(x = timepoint, y = parasitaemia, group = host, color = host), linetype = 2, alpha = 0.5) +
  geom_point(aes(x = timepoint, y = median_parasitaemia), size = 2.5) + 
  geom_path(aes(x = timepoint, y = median_parasitaemia), color = "black", size = 1) +
  scale_y_continuous(trans='log10') +
  ylab("parasitaemia (p/ul)") +
  theme_classic()

timepoints_parasitaemia

### glm

# with timepoints as a factor

null_model <- lmer(log_parasitaemia ~  (1 | host), data = dataframe)
linear_model <- lmer(log_parasitaemia ~ timepoint + (1 | host), data = dataframe)

lrtest <- anova(null_model, linear_model)

print(lrtest)

summary(linear_model)

# with timepoint as a continuous variable

null_model <- lmer(log_parasitaemia ~  (1 | host), data = dataframe)
linear_model <- lmer(log_parasitaemia ~ timepoint_continuous + (1 | host), data = dataframe)

lrtest <- anova(null_model, linear_model)

print(lrtest)

summary(linear_model)

##### comparing with quantity of DBLa expressed #####

rm(list = ls())

df <- read.csv('ASformatted_table.csv')

dataframe_qDBLA <- aggregate(qDBLa ~ sample + entropy + infection_group, data = df, FUN = mean) %>%
  arrange(sample)

cor_qDBLa_entropy <- ggplot(data = dataframe_qDBLA, aes(x = qDBLa, y = entropy, color = infection_group))+
  geom_point() + 
  scale_x_log10() +
  theme_classic()

cor_qDBLa_entropy

null_model <- glm(entropy ~ 1, data = dataframe_qDBLA)
linear_model <- glm(entropy ~ qDBLa, data = dataframe_qDBLA)

lr_test <- anova(linear_model, null_model)

teststat <- lr_test[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

##### Entropy analyses ##### 

rm(list = ls())

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) 

entropy_dataframe_all_data <- aggregate(entropy ~ sample + infection_group, data = dataframe, FUN = mean) %>%
  arrange(sample)

entropy_dataframe_clonality <- aggregate(entropy ~ sample + infection_group + host + clonality, data = dataframe, FUN = mean) %>%
  arrange(sample)

hist_entropy <- ggplot(entropy_dataframe_all_data, aes(x = entropy)) +
  geom_histogram(binwidth = 0.25, color="black", fill="gray") + 
  theme_classic()

hist_entropy

boxplot_entropy <- ggplot(entropy_dataframe_all_data, aes(x = infection_group, y = entropy)) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()

boxplot_entropy

barplot_entropy <- ggplot(entropy_dataframe_all_data, aes(x = sample, y = entropy)) +
  geom_bar(stat = "identity", na.rm = FALSE) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_entropy

boxplot_entropy <- ggplot(entropy_dataframe_clonality, aes(x = clonality, y = entropy)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_brewer() +
  theme_classic()

boxplot_entropy

#### Timepoints data

entropy_dataframe_timepoint_data <- aggregate(entropy ~ sample + infection_group + host + timepoint, data = dataframe, FUN = mean) %>%
  group_by(timepoint) %>%
  mutate(mean_entropy = mean(entropy))%>%
  ungroup %>%
  arrange(sample) 

timepoints_entropy <- ggplot(entropy_dataframe_timepoint_data, aes(x = timepoint, y = entropy, group = host, color = host)) +
  geom_point(size = 2) + 
  geom_line(linetype = 2) +
  scale_color_viridis_b( palette = 'Cividis') +
  theme_classic()

timepoints_entropy <- ggplot(entropy_dataframe_timepoint_data) +
  geom_point(aes(x = timepoint, y = entropy, group = host, color = host), size = 1, alpha = 0.5) + 
  geom_line(aes(x = timepoint, y = entropy, group = host, color = host),linetype = 2, alpha = 0.5) +
  geom_point(aes(x = timepoint, y = mean_entropy), size = 2, color = "black") + 
  geom_line(aes(x = timepoint, y = mean_entropy),linetype = 2, color = "black") +
  theme_classic()

timepoints_entropy

## correlation with parasitaemia

plot_corParasitaemia <- ggplot(data = entropy_dataframe_all_data, aes(x = entropy, y = parasitaemia)) +
  geom_point(size = 2) + 
  scale_y_log10() +
  ylab("parasitaemia (p/ul)") +
  theme_classic()

plot_corParasitaemia

# fitting glm

initial_values <- c(coef1 = 0, coef2 = 0)

glm_model <- glm(formula = entropy ~ parasitaemia, data = entropy_dataframe_all_data)
glm_model_bis <- glm(formula = entropy ~ 1, data = entropy_dataframe_all_data)

lr_test <- anova(glm_model, glm_model_bis)

teststat <- lr_test[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

## Correlation with Fws

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) 

entropy_dataframe_fws <- aggregate(entropy ~ sample + host + Fws, data = dataframe, FUN = mean) %>%
  arrange(sample)

plot_corFws <- ggplot(data = entropy_dataframe_fws, aes(x = Fws, y = entropy)) +
  geom_point(size = 2) + 
  theme_classic()

plot_corFws

## fitting simple linear model

linear_model <- glm(data = entropy_dataframe_fws, 
                    formula = Fws ~ entropy)

summary(linear_model)

### Entropy and antibody levels

entropy_dataframe_antibody <- aggregate(entropy ~ sample + infection_group + host + Antibody_level + timepoint, data = dataframe, FUN = mean) %>%
  arrange(sample) %>%
  mutate(Antibody_level = as.numeric(Antibody_level)) %>%
  filter(timepoint == 'm1')

plot_antibody_ent <- ggplot(data = entropy_dataframe_antibody, aes(x = Antibody_level, y = entropy)) +
  geom_point() +
  theme_classic()

plot_antibody_ent

####### glms and mixed effect models #####

rm(list = ls())

## per infection group

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) 

entropy_dataframe_all_data <- aggregate(entropy ~ sample + infection_group + host, data = dataframe, FUN = mean) %>%
  arrange(sample) 

#simple glm

null_model <- lmer(entropy ~ (1 | host), data = entropy_dataframe_all_data)
full_model <- lmer(entropy ~ infection_group + (1 | host), data = entropy_dataframe_all_data)

lrtest <- anova(null_model, full_model)

print(lrtest)

em <- emmeans(full_model, "infection_group")
contrast(em, "pairwise", adjust = "Tukey")

## per clonality

entropy_dataframe_clonality <- aggregate(entropy ~ sample + host + clonality, data = dataframe, FUN = mean) %>%
  arrange(sample)

null_model <- glm(entropy ~ 1, data = entropy_dataframe_clonality)
full_model <- glm(entropy ~ clonality, data = entropy_dataframe_clonality)

lrtest <- anova(null_model, full_model)

print(lrtest)
teststat <- lrtest[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

## per timepoint

entropy_dataframe_timepoint_data <- aggregate(entropy ~ sample + infection_group + host + timepoint, data = dataframe, FUN = mean) %>%
  group_by(timepoint) %>%
  mutate(mean_entropy = mean(entropy))%>%
  ungroup %>%
  arrange(sample) 

### using mixed models
null_model <- lmer(entropy ~  (1 | host), data = entropy_dataframe_timepoint_data)
linear_model <- lmer(entropy ~ timepoint + (1 | host), data = entropy_dataframe_timepoint_data)

lrtest <- anova(null_model, linear_model)

print(lrtest)

### Looking for tags that would be expressed twice in the same sample (false clustering)

doubled_dataframe <- read.csv("ASformatted_table.csv",header = TRUE) %>%
  mutate(sample = substr(sample, start = 1, stop = 7)) %>%
  group_by(sample) %>%
  group_by(DBL_tag, .add = TRUE) %>%
  mutate(double_tag = ifelse(n() > 1, "doubled", "single")) %>%
  ungroup() %>%
  filter(Read_count > 1)

#count the number of 'doubled' genes

count_double_dataframe <- doubled_dataframe %>%
  count(double_tag)
#### Counting the domains and comparative analysis ####

dataframe <- read.csv("ASformatted_table.csv")

dataframe_NTS <- aggregate(data = dataframe, . ~ NTS, FUN = length)

### Plotting success matches to Pacbio genomes ####

dataframe <- read.csv("pacbio_dataframe.csv") %>%
  distinct(sample, .keep_all = TRUE)

barplot_success <- ggplot(dataframe, aes(x = sample, y = prop_matches)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_success

dataframe_bis <- read.csv("pacbio_dataframe.csv") %>%
  filter(!is.na(pacbio_id)) %>%
  group_by(sample) %>%
  mutate(prop_matches_reads = sum(normalized_reads)) %>%
  distinct(sample, .keep_all = TRUE)

barplot_reads_success <- ggplot(dataframe_bis, aes(x = sample, y = prop_matches_reads)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_reads_success

#### Comparing reads percentage associated with different domains ####

rm(list = ls())

tags_dataframe = read.csv("pacbio_dataframe.csv") %>%
  select('sample', 'normalized_reads', 'DBL_tag', 'ups', 'pacbio_id', 
         'infection_group', 'host') %>%
  filter(host =='DC01' | host == 'DC04'| host == 'DC07'| host == 'DC08'|
         host == 'WS05'| host == 'WS06'| host =='WS07' | host == 'WS08' | 
           host == 'WS11' | host == 'WS13')

domain_dataframe <- read.csv("domains_only.csv") %>%
  select('gene_id', 'Domain1', 'Domain2', 'Domain3', 'Domain4', 'Domain5',
         'Domain6', 'Domain7', 'Domain8', 'Domain9', 'Domain10', 'Domain11',
         'Domain12', 'Domain13', 'Domain14', 'Domain15', 'Domain16', 'Domain17',
         'Domain18') %>%
  mutate(pacbio_id = str_replace(gene_id, "\\.1$", ""))

dataframe <- left_join(tags_dataframe, domain_dataframe, by = "pacbio_id") %>%
  filter(!is.na(gene_id)) %>%
  rowwise %>%
  mutate(DBLb = any(str_detect(c_across(where(is.character)), "DBLb"))) %>%
  mutate(DBLd = any(str_detect(c_across(where(is.character)), "DBLd"))) %>%
  mutate(DBLe = any(str_detect(c_across(where(is.character)), "DBLe"))) %>%
  mutate(DBLz = any(str_detect(c_across(where(is.character)), "DBLz"))) %>%
  mutate(CIDRb = any(str_detect(c_across(where(is.character)), "CIDRb"))) %>%
  mutate(CIDRg = any(str_detect(c_across(where(is.character)), "CIDRg"))) %>%
  mutate(CIDRd = any(str_detect(c_across(where(is.character)), "CIDRd"))) %>%
  ungroup() %>%
  mutate(reads_DBLb = ifelse(DBLb == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLd = ifelse(DBLd == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLe = ifelse(DBLe == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLz = ifelse(DBLz == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_CIDRb = ifelse(CIDRb == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_CIDRg = ifelse(CIDRg == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_CIDRd = ifelse(CIDRd == 'TRUE', normalized_reads, 0)) %>%
  group_by(sample) %>%
  mutate(sum_reads_DBLb = sum(reads_DBLb)) %>%
  mutate(sum_reads_DBLd = sum(reads_DBLd)) %>%
  mutate(sum_reads_DBLe = sum(reads_DBLe)) %>%
  mutate(sum_reads_DBLz = sum(reads_DBLz)) %>%
  mutate(sum_reads_CIDRb = sum(reads_CIDRb)) %>%
  mutate(sum_reads_CIDRg = sum(reads_CIDRg)) %>%
  mutate(sum_reads_CIDRd = sum(reads_CIDRd)) %>%
  ungroup() %>%
  distinct(sample, .keep_all = TRUE)

### Plotting everything

plot_DBLb <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_DBLb)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

plot_DBLd <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_DBLd)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

plot_DBLe <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_DBLe)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

plot_DBLz <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_DBLz)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

plot_CIDRb <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_CIDRb)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

plot_CIDRd <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_CIDRd)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

plot_CIDRg <- ggplot(data = dataframe, aes(x = infection_group, y = sum_reads_CIDRg)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

grid.arrange(plot_DBLb, plot_DBLd, plot_DBLe, plot_DBLz, plot_CIDRb, plot_CIDRd,
             plot_CIDRg, ncol = 2)

### statistical tests

m_DBLb <- lmer(data = dataframe, formula = sum_reads_DBLb ~ infection_group  + (1|host))
summary
  
### Looking at Domain composition of expressed var genes (archives, keeping it here for now) ####

### functions

rm(list = ls())

Extract_expression_domain <- function(entry_dataframe){
for (column in colnames(entry_dataframe)) {
  domain <- substr(column, start = 1, stop = 6)
  if (domain == "Domain"){
    summary_df <- entry_dataframe %>%
      group_by(.data[[column]]) %>%
      summarize(total_quantity = sum(normalized_reads)) 
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


### differences in the domain expression of recurrent vs dominant var genes (archives) 

recurrent_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(recurrence == "recurrent")

recurrent_tags_domains <- Extract_expression_domain(recurrent_dataframe) %>%
  mutate(ID = "recurrent")

plot_recurrent <- plot_domains(recurrent_tags_domains)

plot_recurrent

dominant_single_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(recurrence == "single" & Dominant_tags == "dominant")

dominant_single_domains <- Extract_expression_domain(dominant_single_dataframe) %>%
  mutate(ID = "single_dominant")

plot_dominant_single <- plot_domains(dominant_single_domains)

plot_dominant_single
