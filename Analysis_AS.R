library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(entropy)
library(gridExtra)
library(emmeans)
library(lme4)
library(stringr)
library(lmtest)
library(blandr)
library(ggbeeswarm)

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

## Adding qDBLa for DBla, ATS1 and 2

ASdataframe <- read.csv("ASformatted_table.csv") %>%
  select(-qRT_ATS1, -qRT_ATS2, -qRT_DBLa)

qPCR_dataframe <- read.csv("qRT_pcr_var.csv") %>%
  rename(qRT_ATS1 = varATS1,
         qRT_ATS2 = varATS2,
         qRT_DBLa = varDBLa) %>%
  select(sample, qRT_ATS1, qRT_ATS2, qRT_DBLa)

dataframe_fromatted <- left_join(ASdataframe, qPCR_dataframe, by = 'sample') %>%
  mutate(groupWvsD = ifelse(groupWvsD == "dry", "DRY", "WET")) 


write.csv(dataframe_fromatted, "ASformatted_table.csv")

### Adding the tag_host_combination and linking

# Pacbio

df_pacbio <- read.csv("Tag_matches_Pacbio.csv", header = TRUE) %>%
  mutate(pacbio_id = str_remove_all(pacbio_id, "\\.1:pep"))
ASdataframe <- read.csv("ASformatted_table.csv") %>%
  mutate(tag_host_combination = paste0(host, '_', DBL_tag))
new_IDs <- read.csv("new_var_IDs.csv")

dataframe <- left_join(ASdataframe, df_pacbio, by = "tag_host_combination") %>%
  select(-X, -Host, -tag) %>%
  filter(host == 'DC01' | host == 'DC02'| host == 'DC04'| host =='DC07'| 
        host == 'DC08' | host == 'DC09' |host == 'WA10' | host == 'WA11' | 
        host == 'WS05' | host == 'WS06' | host == 'WS07' | host == 'WS08' | 
        host == 'WS11' | host == 'WS12' | host == 'WS13' | host == 'WS15'| 
        host == 'DC03' | host == 'WS10'| host == 'DC05') %>%
  group_by(sample) %>%
  mutate(prop_matches = 1- mean(is.na(pacbio_id)))%>%
  mutate(original_ID = pacbio_id) 

dataframe_step2 <- left_join(dataframe, new_IDs, by = 'original_ID') %>%
  mutate(tagID = ifelse(is.na(var_ID), yes = DBL_tag, no = var_ID))

#write.csv(dataframe_step2, "pacbio_dataframe.csv")

# SWGA

df_swga <- read.csv("Swga_dataframe.csv")
ASdataframe <- read.csv("ASformatted_table.csv") %>%
  mutate(tag_host_combination = paste0(host, '_', DBL_tag))

dataframe <- left_join(ASdataframe, df_swga, by = "tag_host_combination") %>%
  select(-host.y, -Domain1, -Domain2, -Domain3, -Domain4, -Domain5, -Domain6,
         -Domain7, -Domain8, -Domain9, -Domain10, - cohortID, -E_value, -VarDb_hits) %>%
  rename(host = host.x) %>%
  group_by(sample) %>%
  mutate(prop_matches_swga = 1- mean(is.na(full_id))) 

write.csv(dataframe, "SWGA_dataframe.csv")

### 'Category 3'

df_all_matches <- read.table('blast_matches_alldata.txt', header = FALSE) %>%
  rename(DBL_tag = V1, gene_id = V2) %>%
  select(DBL_tag, gene_id) %>%
  distinct(DBL_tag, .keep_all = TRUE)
ASdataframe <- read.csv("ASformatted_table.csv")

dataframe <- left_join(ASdataframe, df_all_matches, by = "DBL_tag") %>% 
  group_by(sample) %>%
  mutate(prop_matches_cat3 = 1- mean(is.na(gene_id))) 

write.csv(dataframe, "category3_dataframe.csv")

#### Barplots of read counts per sample ####

#### function to plot the graph, takes a dataframe with all timepoint from a host as entry

plot_multiple_timepoints <- function(entry_dataframe) {
  
  entry_dataframe <- entry_dataframe %>%
    arrange(timepoint) %>%
    mutate(DBL_tag = factor(tagID, levels = unique(tagID))) %>%
    arrange(timepoint, desc(normalized_reads))
  
  plot_timepoints = ggplot(entry_dataframe, aes(x = tagID, y = normalized_reads)) +
    geom_col(width=0.7) +
    facet_wrap(~ sample, nrow = 10) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Relative var gene expression")
  
  return(plot_timepoints)
}

plot_multiple_timepoints_noPB <- function(entry_dataframe) {
  
  entry_dataframe <- entry_dataframe #%>%
    #arrange(timepoint) %>%
    #mutate(DBL_tag = factor(DBL_tag, levels = unique(DBL_tag))) %>%
    #arrange(timepoint, desc(normalized_reads))
  
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
  
  plot_timepoints = ggplot(entry_dataframe, aes(x = tagID, y = normalized_reads, fill = NTS)) +
    geom_col(width=0.7) +
    facet_wrap(~ sample, nrow = 10) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Relative var gene expression")
  
  return(plot_timepoints)
}

#### All graphs per sample:

## DC01 

DC01_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  filter(sample == "DC01_m1" | sample == "DC01_m3" |  sample == "DC01_m4c" | 
           sample == "DC01_m5" | sample == "DC01_m6") %>%
  arrange(timepoint)

plot_DC01 <- plot_multiple_timepoints(DC01_dataframe)

plot_DC01

## DC02

DC02_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  filter(sample == "DC02_m1" | sample == "DC02_m2" |  sample == "DC02_m3" |
           sample == "DC02_m4" | sample == "DC02_m5" | sample == "DC02_m6")

plot_DC02 = plot_multiple_timepoints(DC02_dataframe)

plot_DC02

## DC03

DC03_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(sample == "DC03_m2" | sample == "DC03_m4" | sample == "DC03_m5" | 
         sample == "DC03_m6")

plot_DC03 = plot_multiple_timepoints(DC03_dataframe)
  
plot_DC03

## DC04

DC04_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC04_m2" | sample == "DC04_m3" |  sample == "DC04_m4c" | 
           sample == "DC04_m5" |  sample == "DC04_m6")

plot_DC04 = plot_multiple_timepoints(DC04_dataframe)

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
  filter(sample == "DC05_m1" | sample == "DC05_m2" | sample == "DC05_m3" |
           sample == "DC05_m4" | sample == "DC05_m5" | sample == "DC05_m6" )

plot_DC05 = plot_multiple_timepoints_noPB(DC05_dataframe)

plot_DC05

## DC07

DC07_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
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
  filter(sample == "DC09_m2" | sample == "DC09_m3" |  sample == "DC09_m4" |
           sample == "DC09_m5")

plot_DC09 = plot_multiple_timepoints(DC09_dataframe)

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

plot_DC11 = plot_multiple_timepoints_noPB(DC11_dataframe)

plot_DC11

## DC12

DC12_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC12_m1" | sample == "DC12_m2" |  sample == "DC12_m4" |
           sample == "DC12_m5" | sample == "DC12_m6")

plot_DC12 = plot_multiple_timepoints_noPB(DC12_dataframe)

plot_DC12

## DC13

DC13_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC13_m1" | sample == "DC13_m2" |  sample == "DC13_m3" |
           sample == "DC13_m4" | sample == "DC13_m5")

plot_DC13 = plot_multiple_timepoints_noPB(DC13_dataframe)

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


#### Finding recurrent var genes ####

## DC01

DC01_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  filter(sample == "DC01_m1" | sample == "DC01_m3" |  sample == "DC01_m4c" | 
           sample == "DC01_m5" | sample == "DC01_m6") %>%
  arrange(timepoint) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, tagID)

## DC02

DC02_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  filter(sample == "DC02_m1" | sample == "DC02_m2" |  sample == "DC02_m3" |
           sample == "DC02_m4" | sample == "DC02_m5" | sample == "DC02_m6") %>%
  arrange(tagID) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, tagID)

## DC04

DC04_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC04_m2" | sample == "DC04_m3" |  sample == "DC04_m4c" | 
           sample == "DC04_m5" |  sample == "DC04_m6") #%>%
  arrange(tagID) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, tagID)

## DC05

DC05_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC05_m1" | sample == "DC05_m2" | sample == "DC05_m3" |
           sample == "DC05_m4" | sample == "DC05_m5" | sample == "DC05_m6" ) %>%
  arrange(DBL_tag) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, DBL_tag)

##DC07

DC07_dataframe <- read.csv("pacbio_dataframe.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC07_m1" | sample == "DC07_m2" |  sample == "DC07_m3" |
           sample == "DC07_m5" | sample == "DC07_m6") %>%
  arrange(tagID) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, tagID)

## DC11

DC11_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC11_m1" | sample == "DC11_m2" |  sample == "DC11_m4" |
           sample == "DC11_m5" | sample == "DC11_m6")%>%
  arrange(DBL_tag) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, DBL_tag)

### DC12

DC12_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC12_m1" | sample == "DC12_m2" |  sample == "DC12_m4" |
           sample == "DC12_m5" | sample == "DC12_m6") %>%
  arrange(DBL_tag) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, DBL_tag) 

## DC13

DC13_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC13_m1" | sample == "DC13_m2" |  sample == "DC13_m3" |
           sample == "DC13_m4" | sample == "DC13_m5") %>%
  arrange(DBL_tag) %>%
  filter(normalized_reads > 0.05) %>%
  select(sample,normalized_reads, DBL_tag) 
  

#### Piecharts plots ####

## annotated via cUPS

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

## pacbio annotation

pacbio_dataframe <- read.csv("pacbio_dataframe.csv") 
df_ups_pacbio <- read.csv("ups_varIDs.txt", header = FALSE) %>%
  rename(pacbio_id = V1)

p_dataframe <- left_join(x = pacbio_dataframe, y = df_ups_pacbio, by = 'pacbio_id')
  
my_palette <- c("#FFC20A", "#0C7BDC","#ed4d87", "#a2a2a2")

plot_piecharts = ggplot(p_dataframe, aes(x = "", y = normalized_reads, fill = V2)) +
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


#### Summary statistics DBLa ####

AS_dataframe <- read.csv("ASformatted_table.csv") %>%
  group_by(sample) %>%
  mutate(nbDBLa = n_distinct(DBL_tag)) %>%
  ungroup() %>%
  distinct(sample, .keep_all = TRUE)

median(AS_dataframe$nbDBLa)
min(AS_dataframe$nbDBLa)
max(AS_dataframe$nbDBLa)



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

# comparing WS vs DC vs WA A / BC

dataframe <- read.csv("frequency_upsA_DBLa.csv") %>%
  mutate(infection_group = substr(sample, start = 1, stop = 2)) %>%
  mutate(host = substr(sample, start = 1, stop = 4)) %>%
  filter(ups == 'A')

m_groupA <- lmer(data = dataframe, formula = cumulative_frequency ~ infection_group  + (1|host))
null_groupA <- lmer(data = dataframe, formula = cumulative_frequency ~ (1|host))
anova(m_groupA, null_groupA)
summary(m_groupA)

em <- emmeans(m_groupA, "infection_group")
contrast(em, "pairwise", adjust = "Tukey")


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

##### qRT-DBLa  #####

rm(list = ls())

df <- read.csv('ASformatted_table.csv') %>%
  mutate(qRT_ATS1 = as.numeric(qRT_ATS1),
         qRT_ATS2 = as.numeric(qRT_ATS2),
         qRT_DBLa = as.numeric(qRT_DBLa))  

## DBLa

dataframe_qDBLA <- aggregate(qRT_DBLa ~ sample + groupWvsD + entropy, data = df, FUN = mean) %>%
  arrange(sample)

qDBLa_per_group <- ggplot(data = dataframe_qDBLA, aes(x = groupWvsD, y = qRT_DBLa))+
  geom_violin(width = .6) + 
  geom_beeswarm() +
  theme_classic()

qDBLa_per_group

null_model <- glm(qRT_DBLa ~ 1, data = dataframe_qDBLA)
linear_model <- glm(qRT_DBLa ~ groupWvsD, data = dataframe_qDBLA)

lr_test <- anova(linear_model, null_model)

teststat <- lr_test[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

## ATS1

dataframe_qATS1 <- aggregate(qRT_ATS1 ~ sample + groupWvsD, data = df, FUN = mean) %>%
  arrange(sample)

qATS1_per_group <- ggplot(data = dataframe_qATS1, aes(x = groupWvsD, y = qRT_ATS1))+
  geom_jitter() + 
  geom_boxplot() +
  theme_classic()

qATS1_per_group

null_model <- glm(qRT_ATS1 ~ 1, data = dataframe_qATS1)
linear_model <- glm(qRT_ATS1 ~ groupWvsD, data = dataframe_qATS1)

lr_test <- anova(linear_model, null_model)

teststat <- lr_test[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

## Correlation with entropy

cor_qDBLa_entropy <- ggplot(data = dataframe_qDBLA, aes(x = qRT_DBLa, y = entropy, color = groupWvsD))+
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

## Correlation qDBLa parasitaemia

dataframe_qDBLA <- aggregate(qRT_DBLa ~ groupWvsD + parasitaemia, data = df, FUN = mean) 

cor_qDBLa_par <- ggplot(data = dataframe_qDBLA, aes(x = qRT_DBLa, y = parasitaemia))+
  geom_point() + 
  scale_x_log10() +
  theme_classic()

cor_qDBLa_par

null_model <- glm(entropy ~ 1, data = dataframe_qDBLA)
linear_model <- glm(entropy ~ qDBLa, data = dataframe_qDBLA)

lr_test <- anova(linear_model, null_model)

teststat <- lr_test[2, "Deviance"]
pvalue <- pchisq(teststat, df = 1, lower.tail = FALSE)

pvalue

## qDBLa over time in chronic infections

df_qDBLa_chronic <- aggregate(data = df, qRT_DBLa ~ sample + host + timepoint, FUN = mean) %>%
  mutate(month = as.numeric(substr(timepoint, 2, 2))) %>%
  group_by(timepoint) %>%
  mutate(mean_qRT_DBLa = mean(qRT_DBLa)) %>%
  ungroup() %>%
  arrange(sample)

timepoints_qRT_DBLa <- ggplot(df_qDBLa_chronic) +
  geom_point(aes(x = timepoint, y = qRT_DBLa, group = host, color = host), size = 1, alpha = 0.5) + 
  geom_line(aes(x = timepoint, y = qRT_DBLa, group = host, color = host),linetype = 2, alpha = 0.5) +
  geom_point(aes(x = timepoint, y = mean_qRT_DBLa), size = 2, color = "black") +
  theme_classic()

timepoints_qRT_DBLa 

model_time <- lmer(data = df_qDBLa_chronic, formula = qRT_DBLa ~ month + (1|host))
null_model <- lmer(data = df_qDBLa_chronic, formula = qRT_DBLa ~ (1|host))

lrtest <- anova(model_time, null_model)

print(lrtest)

##### Entropy analyses ##### 

rm(list = ls())

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) 

entropy_dataframe_all_data <- aggregate(entropy ~ sample + groupWvsD, data = dataframe, FUN = mean) %>%
  arrange(sample)

entropy_dataframe_qDBLa <- aggregate(entropy ~ sample + infection_group + qRT_DBLa, data = dataframe, FUN = mean) %>%
  arrange(sample)

## infection group, summary of all data

hist_entropy <- ggplot(entropy_dataframe_all_data, aes(x = entropy)) +
  geom_histogram(binwidth = 0.25, color="black", fill="gray") + 
  theme_classic()

hist_entropy

violinplot_entropy_group <- ggplot(entropy_dataframe_all_data, aes(x = groupWvsD, y = entropy)) +
  geom_violin(width=.6) +
  geom_beeswarm() +
  theme_classic()

violinplot_entropy_group

barplot_entropy <- ggplot(entropy_dataframe_all_data, aes(x = sample, y = entropy)) +
  geom_bar(stat = "identity", na.rm = FALSE) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_entropy



### CLonality

entropy_dataframe_clonality <- aggregate(entropy ~ sample + groupWvsD+ host + clonality, data = dataframe, FUN = mean) %>%
  arrange(sample)


violinplot_entropy_clonality <- ggplot(entropy_dataframe_clonality, aes(x = clonality, y = entropy)) +
  geom_violin(width=.6) +
  geom_beeswarm() +
  scale_fill_brewer() +
  theme_classic()

violinplot_entropy_clonality

corplot_parasitaemia <- ggplot(entropy_dataframe_all_data, aes(x = parasitaemia, y = entropy)) +
  geom_point() +
  scale_x_log10() +
  theme_classic()

corplot_parasitaemia

corplot_qDBLa <- ggplot(entropy_dataframe_qDBLa, aes(x = qRT_DBLa, y = entropy)) +
  geom_point() +
  theme_classic()

corplot_qDBLa

grid.arrange(corplot_parasitaemia, corplot_qDBLa, ncol = 2)

#### Timepoints data

entropy_dataframe_timepoint_data <- aggregate(entropy ~ sample + infection_group + host + timepoint, data = dataframe, FUN = mean) %>%
  group_by(timepoint) %>%
  mutate(mean_entropy = mean(entropy))%>%
  ungroup %>%
  arrange(sample) 


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
summary(glm_model)

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

##### DBLa vs transcriptomics ####

## Bland altman method

df_exon1 <- read.csv("transcripto_entropy_exon1.csv") %>%
  mutate(normalized_expression = ifelse(var_expression > 0.01, var_expression, 0)) %>%
  rename(normalized_exon1 = normalized_expression) %>%
  mutate(gene_id = str_replace(gene_id, '.1-exon1', '')) %>%
  mutate(gene_sample_combination = paste0(sample, '_', gene_id)) %>%
  distinct(gene_sample_combination, .keep_all = TRUE)

df_dbla <- read.csv("pacbio_dataframe.csv") %>%
  select(sample, pacbio_id, normalized_reads, DBL_tag) %>%
  rename(gene_id = pacbio_id, normalized_dbla = normalized_reads) %>%
  mutate(gene_id = str_replace(gene_id, ':pep', '')) %>%
  mutate(tag_host_c = paste0(sample, '_', DBL_tag)) %>%
  mutate(gene_sample_combination = ifelse(is.na(gene_id), tag_host_c, paste0(sample, '_', gene_id))) %>%
  mutate(normalized_dbla = ifelse(row_number() == 40, 0.1622, normalized_dbla)) %>%
  distinct(gene_sample_combination, .keep_all = TRUE)

df_bland_altmann <- full_join(df_exon1, df_dbla, by = 'gene_sample_combination') %>%
  select(gene_sample_combination, sample.x, sample.y, normalized_dbla, normalized_exon1) %>%
  mutate(sample = ifelse(sample.x == 'NA', sample.y, sample.x)) %>%
  select(-sample.x, -sample.y) %>%
  mutate(normalized_dbla = ifelse(is.na(normalized_dbla), 0, normalized_dbla)) %>%
  mutate(normalized_exon1 = ifelse(is.na(normalized_exon1), 0, normalized_exon1)) %>%
  filter(!(normalized_exon1 == 0 & normalized_dbla == 0)) %>%
  filter(sample == 'DC01_m3' | sample == 'DC01_m6' | sample == 'DC07_m1' |
        sample == 'DC07_m2' | sample == 'WA11' | sample == 'WS06')

blandr.statistics (df_bland_altmann$normalized_dbla, df_bland_altmann$normalized_exon1, sig.level=0.95 )

stats_blandAltman <- blandr.statistics (df_bland_altmann$normalized_dbla, df_bland_altmann$normalized_exon1, sig.level=0.95 )

## plot with the built-in function

blandr.draw(df_bland_altmann$normalized_dbla, df_bland_altmann$normalized_exon1)

## Using ggplot 

## Loading data 

rm(list = ls())

DBLa_dataframe <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  select(sample, entropy) %>%
  rename(entropy_dbla = entropy) %>%
  distinct(sample, .keep_all = TRUE)

Transcripto_dataframe_complete <- read.csv("transcripto_entropy_complete.csv", header = TRUE) %>%
  mutate(var_expression_normalized = ifelse(var_expression > 0.01, var_expression, 0)) %>%
  group_by(sample) %>% 
  mutate(entropy_transcripto = entropy(var_expression_normalized)) %>%
  distinct(sample, .keep_all = TRUE)

### complete == mapping on all var gene, not just exon 1 or 2

Transcripto_dataframe_exon1 <- read.csv("transcripto_entropy_exon1.csv", header = TRUE) %>%
  mutate(var_expression_normalized = ifelse(var_expression > 0.01, var_expression, 0)) %>%
  group_by(sample) %>% 
  mutate(entropy_transcripto = entropy(var_expression_normalized)) %>%
  distinct(sample, .keep_all = TRUE)

quality_dataframe <- read.csv("quality_transcripto.csv", header = TRUE)

## 'complete / fullvar' dataframe

dataframe_step1_complete <- left_join(DBLa_dataframe, Transcripto_dataframe_complete, by = "sample") %>%
  filter(!is.na(entropy_transcripto))

dataframe_complete <- left_join(dataframe_step1_complete, quality_dataframe, by = "sample")

scatterplot_complete_entropies <- ggplot(dataframe_complete, aes(x = entropy_dbla, 
                                               y = entropy_transcripto,
                                               color = Detected_gene)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "red")

scatterplot_complete_entropies

lm_entropies <- lm(data = dataframe_complete, formula = entropy_dbla ~ entropy_transcripto)

summary(lm_entropies)

### Exon 1 dataframe

dataframe_step1_exon1 <- left_join(DBLa_dataframe, Transcripto_dataframe_exon1, by = "sample") %>%
  filter(!is.na(entropy_transcripto))

dataframe_exon1 <- left_join(dataframe_step1_exon1, quality_dataframe, by = "sample") %>%
  rename(entropy_exon1 = entropy_transcripto) %>%
  filter(sample == 'DC01_m3' | sample == 'DC01_m6' | sample == 'DC02_m1' |
           sample == 'DC04_m5' | sample == 'DC03_m4' | sample == 'DC07_m1' |
           sample == 'WA11' | sample == 'WS06')

scatterplot_exon1_entropies <- ggplot(dataframe_exon1, aes(x = entropy_dbla, 
                                                                 y = entropy_exon1,
                                                                 color = Detected_gene)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "red") + 
  ggtitle("correlation entropy exon 1 / DBLa")

scatterplot_exon1_entropies

lm_entropies <- lm(data = dataframe_exon1, formula = entropy_dbla ~ entropy_exon1)

summary(lm_entropies)

### correlate entropies between exon 1 / fullvar

dataframe_fullvar_simplified <- dataframe_complete %>%
  rename(entropy_fullvar = entropy_transcripto) %>%
  select(sample, entropy_fullvar) 

dataframe_exon1_simplifies <- dataframe_exon1 %>%
  rename(entropy_exon1 = entropy_exon1) %>%
  select(sample, entropy_exon1, Detected_gene) 

dataframe_correlation <- left_join(dataframe_fullvar_simplified, dataframe_exon1_simplifies,
                                   by = 'sample')

plot_two_entropies <- ggplot(data = dataframe_correlation,
                             aes(x = entropy_fullvar, y = entropy_exon1, color = Detected_gene)) +
  geom_point() + 
  ggtitle("Correlation entropies full var gene / exon 1") +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "red")

plot_two_entropies

grid.arrange(plot_two_entropies, scatterplot_exon1_entropies, ncol = 2)

lm_entropies <- lm(data = dataframe_correlation, formula = entropy_fullvar ~ entropy_exon1)

summary(lm_entropies)

##### DBLa vs transcripto, samples barplots

### DC02_m2 (outlier)

rm(list = ls())

plot_dbla_vs_exon1 <- function(sample_id) {
  
  sample_id_exon1 <- str_c(sample_id, "_exon1")
  sample_id_dbla <-  str_c(sample_id, "_dbla")
  
  df_exon1 <- read.csv("transcripto_entropy_exon1.csv") %>%
    mutate(normalized_expression = ifelse(var_expression > 0.01, var_expression, 0)) %>%
    filter(sample == sample_id) %>%
    mutate(gene_id = str_replace(gene_id, '.1-exon1', '')) %>%
    mutate(sample = str_replace(sample, sample_id, sample_id_exon1))
  
  df_dbla <- read.csv("pacbio_dataframe.csv") %>%
    select(sample, pacbio_id, normalized_reads) %>%
    filter(sample == sample_id) %>%
    rename(gene_id = pacbio_id, normalized_expression = normalized_reads) %>%
    mutate(gene_id = str_replace(gene_id, ':pep', '')) %>%
    mutate(sample = str_replace(sample, sample_id, sample_id_dbla))
  
  df_bind <- bind_rows(df_dbla, df_exon1)
  
  plot_timepoints = ggplot(df_bind, aes(x = gene_id, y = normalized_expression)) +
    geom_col(width=0.7) +
    theme_bw() +
    facet_wrap(~ sample, nrow = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Relative var gene expression")
  
  plot_timepoints
  
}

plot_dbla_vs_exon1('DC01_m3')
plot_dbla_vs_exon1('DC01_m6')
plot_dbla_vs_exon1('DC03_m4')
plot_dbla_vs_exon1('DC07_m1')
plot_dbla_vs_exon1('DC07_m2')
plot_dbla_vs_exon1('WA11')
plot_dbla_vs_exon1('WS06')

####### glms and mixed effect models of samples dbla entropies#####

rm(list = ls())

## per infection group

dataframe <- read.csv("ASformatted_table.csv", header = TRUE) 

entropy_dataframe_all_data<- aggregate(entropy ~ sample + groupWvsD + host, data = dataframe, FUN = mean) %>%
  arrange(sample)

#simple glm

null_model <- lmer(entropy ~ (1 | host), data = entropy_dataframe_all_data)
full_model <- lmer(entropy ~ groupWvsD + (1 | host), data = entropy_dataframe_all_data)
lrtest <- anova(null_model, full_model)
print(lrtest)
summary(full_model)

## looking for difference between WS and WA

dataframe_WAvsWS <- read.csv("ASformatted_table.csv", header = TRUE) %>%
  filter(infection_group == "WS" | infection_group == "WA")

entropy_dataframe_WAvsWS <- aggregate(entropy ~ sample + infection_group + host, data = dataframe_WAvsWS, FUN = mean) %>%
  arrange(sample) 

null_model <- glm(entropy ~ 1, data = entropy_dataframe_WAvsWS)
full_model <- glm(entropy ~ infection_group, data = entropy_dataframe_WAvsWS)
lrtest <- anova(null_model, full_model)
print(lrtest)
summary(full_model)

entropy_dataframe_WAvsWS <- aggregate(entropy ~ sample + infection_group + host, data = dataframe_WAvsWS, FUN = mean) %>%
  filter(infection_group == "WS") %>%
  count()

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

entropy_dataframe_clonalitycount <- entropy_dataframe_clonality %>%
  filter(clonality == 'polyclonal') %>%
  count()
  

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
### Plotting success matches to Pacbio genomes / SWGA ####

## Pacbio

dataframe <- read.csv("pacbio_dataframe.csv") %>%
  group_by(sample) %>%
  mutate(prop_matches = 1- mean(is.na(pacbio_id))) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(host == 'DC05')

barplot_success <- ggplot(dataframe, aes(x = sample, y = prop_matches)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.8,           
             linetype = "dashed",        
             color = "red",              
             size = 0.5)   

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

# SWGA

dataframe <- read.csv("SWGA_dataframe.csv") %>%
  distinct(sample, .keep_all = TRUE)

barplot_success <- ggplot(dataframe, aes(x = sample, y = prop_matches_swga)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_success

dataframe_bis <- read.csv("SWGA_dataframe.csv") %>%
  distinct(cluster_id, .keep_all = TRUE) %>%
  filter(!is.na(full_id)) %>%
  group_by(sample) %>%
  mutate(prop_matches_reads = sum(normalized_reads)) %>%
  distinct(sample, .keep_all = TRUE)

barplot_reads_success <- ggplot(dataframe_bis, aes(x = sample, y = prop_matches_reads)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_reads_success

## 'Category_3'

dataframe <- read.csv("category3_dataframe.csv") %>%
  distinct(sample, .keep_all = TRUE)

barplot_success <- ggplot(dataframe, aes(x = sample, y = prop_matches_cat3)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

barplot_success

dataframe_bis <- read.csv("category3_dataframe.csv") %>%
  filter(!is.na(gene_id)) %>%
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
         'infection_group', 'host',) %>%
  mutate(groupWvsD = ifelse(infection_group == 'DC', yes = 'DRY', no = 'WET')) %>%
  filter(sample =='DC01_m1' | sample =='DC01_m4c' | sample =='DC01_m5' | 
         sample =='DC01_m6' | sample =='DC02_m1' | sample =='DC02_m2' |
         sample =='DC02_m3' | sample == 'DC02_m4' | sample == 'DC02_m5' |
         sample =='DC02_m6' | sample == 'DC04_m2' | sample == 'DC04_m4c' |
         sample =='DC04_m5' | sample == 'DC04_m6' | sample == 'DC07_m1' |
         sample =='DC07_m3' | sample == 'DC07_m5' | sample == 'DC07_m6' |
         sample =='DC08_m4' | sample =='DC08_m5' | sample =='DC08_m2' |
         sample =='WS05' | sample =='WS06' | sample =='WS07' | 
         sample =='WS08' | sample =='WS11' | sample =='WS13' |
         sample == 'WS10' | sample == 'WA11' | sample == 'DC05_m1' | 
         sample == 'DC05_m2' | sample == 'DC05_m3' | sample == 'DC05_m4' |
         sample == 'DC05_m5' | sample == 'DC05_m6') %>%
  mutate(pacbio_id = str_replace(pacbio_id, ".1:pep", ""))

domain_dataframe <- read.csv("domains_only.csv") %>%
  select('gene_id', 'Domain1', 'Domain2', 'Domain3', 'Domain4', 'Domain5',
         'Domain6', 'Domain7', 'Domain8', 'Domain9', 'Domain10', 'Domain11',
         'Domain12', 'Domain13', 'Domain14', 'Domain15', 'Domain16', 'Domain17',
         'Domain18') %>%
  mutate(pacbio_id = str_replace(gene_id, "\\.1$", ""))
  
ups_dataframe <- read.csv("ups_varIDs.txt", header = FALSE) %>%
  rename(pacbio_id = V1,
         ups_group = V2)

CIDRa_dataframe <- read.csv("CIDRa_classification.csv") %>%
  rename(pacbio_id = gene_id)

dataframe_stepA <- left_join(tags_dataframe, domain_dataframe, by = "pacbio_id") %>%
  filter(!is.na(pacbio_id)) %>%
  rowwise %>%
  mutate(DBLb = any(str_detect(c_across(where(is.character)), "DBLb"))) %>%
  mutate(DBLd = any(str_detect(c_across(where(is.character)), "DBLd"))) %>%
  mutate(DBLg = any(str_detect(c_across(where(is.character)), "DBLg"))) %>%
  mutate(DBLe = any(str_detect(c_across(where(is.character)), "DBLe"))) %>%
  mutate(DBLz = any(str_detect(c_across(where(is.character)), "DBLz"))) %>%
  mutate(CIDRb = any(str_detect(c_across(where(is.character)), "CIDRb"))) %>%
  mutate(CIDRg = any(str_detect(c_across(where(is.character)), "CIDRg"))) %>%
  mutate(CIDRd = any(str_detect(c_across(where(is.character)), "CIDRd"))) %>%
  ungroup() %>%
  mutate(reads_DBLb = ifelse(DBLb == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLd = ifelse(DBLd == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLg = ifelse(DBLg == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLe = ifelse(DBLe == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_DBLz = ifelse(DBLz == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_CIDRb = ifelse(CIDRb == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_CIDRg = ifelse(CIDRg == 'TRUE', normalized_reads, 0)) %>%
  mutate(reads_CIDRd = ifelse(CIDRd == 'TRUE', normalized_reads, 0)) %>%
  group_by(sample) %>%
  mutate(sum_reads_DBLb = sum(reads_DBLb)) %>%
  mutate(sum_reads_DBLd = sum(reads_DBLd)) %>%
  mutate(sum_reads_DBLg = sum(reads_DBLg)) %>%
  mutate(sum_reads_DBLe = sum(reads_DBLe)) %>%
  mutate(sum_reads_DBLz = sum(reads_DBLz)) %>%
  mutate(sum_reads_CIDRb = sum(reads_CIDRb)) %>%
  mutate(sum_reads_CIDRg = sum(reads_CIDRg)) %>%
  mutate(sum_reads_CIDRd = sum(reads_CIDRd)) %>%
  ungroup() 

dataframe_stepB <- left_join(dataframe_stepA, ups_dataframe, by = "pacbio_id") %>%
  mutate(reads_upsA = ifelse(ups_group == 'A', normalized_reads, 0),
         reads_upsB = ifelse(ups_group == 'B', normalized_reads, 0),
         reads_upsC = ifelse(ups_group == 'C', normalized_reads, 0)) %>%
  group_by(sample) %>%
  mutate(sum_reads_upsA = sum(reads_upsA, na.rm = TRUE),
         sum_reads_upsB = sum(reads_upsB, na.rm = TRUE),
         sum_reads_upsC = sum(reads_upsC, na.rm = TRUE)) %>%
  ungroup()

dataframe <- left_join(dataframe_stepB, CIDRa_dataframe, by = "pacbio_id") %>%
  mutate(reads_CIDRa1 = ifelse(CIDR == 'CIDRa1', normalized_reads, 0)) %>%
  group_by(sample) %>%
  mutate(sum_reads_CIDRa1 = sum(reads_CIDRa1, na.rm = TRUE)) %>%
  ungroup #%>%
  distinct(sample, .keep_all = TRUE)

### Plotting everything

plot_upsA <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_upsA)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'upstream A') +
  theme_bw()

plot_upsB <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_upsB)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'upstream B') +
  theme_bw()

plot_upsC <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_upsC)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'upstream C') +
  theme_bw()

plot_DBLb <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_DBLb)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'DBLb') +
  theme_bw()

plot_DBLg <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_DBLg)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'DBLg') +
  theme_bw()

plot_DBLd <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_DBLd)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'DBLd') +
  theme_bw()

plot_DBLe <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_DBLe)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'DBLe') +
  theme_bw()

plot_DBLz <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_DBLz)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'DBLz') +
  theme_bw()

plot_CIDRb <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_CIDRb)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'CIDRb') +
  theme_bw()

plot_CIDRd <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_CIDRd)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'CIDRd') +
  theme_bw()

plot_CIDRg <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_CIDRg)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  labs(x = 'Infection group', y = 'CIDRg') +
  theme_bw()

plot_CIDRa <- ggplot(data = dataframe, aes(x = groupWvsD, y = sum_reads_CIDRa1)) +
  geom_violin() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) +
  ylim(0,1) +
  theme_bw()

grid.arrange(plot_upsA, plot_upsB, plot_upsC, plot_DBLb, plot_DBLg, plot_DBLd, plot_DBLe, plot_DBLz, plot_CIDRb, plot_CIDRd,
             plot_CIDRg, ncol = 2)

### statistical tests

m_upsA <- lmer(data = dataframe, formula = sum_reads_upsA ~ groupWvsD  + (1|host))
null_upsA <- lmer(data = dataframe, formula = sum_reads_upsA ~ (1|host))
anova(m_upsA, null_upsA)
summary(m_upsA)

m_upsB <- lmer(data = dataframe, formula = sum_reads_upsB ~ groupWvsD  + (1|host))
null_upsB <- lmer(data = dataframe, formula = sum_reads_upsB ~ (1|host))
anova(m_upsB, null_upsB)
summary(m_upsB)

m_upsC <- lmer(data = dataframe, formula = sum_reads_upsC ~ groupWvsD  + (1|host))
null_upsC <- lmer(data = dataframe, formula = sum_reads_upsC ~ (1|host))
anova(m_upsC, null_upsC)
summary(m_upsC)

m_DBLb <- lmer(data = dataframe, formula = sum_reads_DBLb ~ groupWvsD  + (1|host))
null_DBLb <- lmer(data = dataframe, formula = sum_reads_DBLb ~ (1|host))
anova(m_DBLb, null_DBLb)
summary(m_DBLb)

m_DBLg <- lmer(data = dataframe, formula = sum_reads_DBLg ~ groupWvsD  + (1|host))
null_DBLg <- lmer(data = dataframe, formula = sum_reads_DBLg ~ (1|host))
anova(m_DBLg, null_DBLg)
summary(m_DBLg)

m_DBLd <- lmer(data = dataframe, formula = sum_reads_DBLd ~ groupWvsD  + (1|host))
null_DBLd <- lmer(data = dataframe, formula = sum_reads_DBLd ~ (1|host))
anova(m_DBLd, null_DBLd)
summary(m_DBLd)

m_DBLe <- lmer(data = dataframe, formula = sum_reads_DBLe ~ groupWvsD  + (1|host))
null_DBLe <- lmer(data = dataframe, formula = sum_reads_DBLe ~ (1|host))
anova(m_DBLe, null_DBLe)
summary(m_DBLe)

m_DBLz <- lmer(data = dataframe, formula = sum_reads_DBLz ~ groupWvsD  + (1|host))
null_DBLz <- lmer(data = dataframe, formula = sum_reads_DBLz ~ (1|host))
anova(m_DBLz, null_DBLz)
summary(m_DBLz)

m_CIDRb <- lmer(data = dataframe, formula = sum_reads_CIDRb ~ groupWvsD  + (1|host))
null_CIDRb <- lmer(data = dataframe, formula = sum_reads_CIDRb ~ (1|host))
anova(m_CIDRb, null_CIDRb)
summary(m_CIDRb)

m_CIDRd <- lmer(data = dataframe, formula = sum_reads_CIDRd ~ groupWvsD  + (1|host))
null_CIDRd <- lmer(data = dataframe, formula = sum_reads_CIDRd ~ (1|host))
anova(m_CIDRd, null_CIDRd)
summary(m_CIDRd)

m_CIDRg <- lmer(data = dataframe, formula = sum_reads_CIDRg ~ groupWvsD  + (1|host))
null_CIDRg <- lmer(data = dataframe, formula = sum_reads_CIDRg ~ (1|host))
anova(m_CIDRg, null_CIDRg)
summary(m_CIDRg)

m_CIDRa <- lmer(data = dataframe, formula = sum_reads_CIDRa1 ~ groupWvsD  + (1|host))
null_CIDRa <- lmer(data = dataframe, formula = sum_reads_CIDRa1 ~ (1|host))
anova(m_CIDRa, null_CIDRa)

#### Length of var genes in different groups ####

rm(list = ls())

tags_dataframe = read.csv("pacbio_dataframe.csv") %>%
  select('sample', 'normalized_reads', 'DBL_tag', 'ups', 'pacbio_id', 
         'infection_group', 'host', 'Dominant_tags', 'var_ID') %>%
  mutate(groupWvsD = ifelse(infection_group == 'DC', yes = 'DRY', no = 'WET')) %>%
  filter(sample =='DC01_m1' | sample =='DC01_m4c' | sample =='DC01_m5' | 
           sample =='DC01_m6' | sample =='DC02_m1' | sample =='DC02_m2' |
           sample =='DC02_m3' | sample == 'DC02_m4' | sample == 'DC02_m5' |
           sample =='DC02_m6' | sample == 'DC04_m2' | sample == 'DC04_m4c' |
           sample =='DC04_m5' | sample == 'DC04_m6' | sample == 'DC07_m1' |
           sample =='DC07_m3' | sample == 'DC07_m5' | sample == 'DC07_m6' |
           sample =='DC08_m4' | sample =='DC08_m5' | sample =='DC08_m2' |
           sample =='WS05' | sample =='WS06' | sample =='WS07' | 
           sample =='WS08' | sample =='WS11' | sample =='WS13' |
           sample == 'WS10' | sample == 'WA11') %>%
  mutate(pacbio_id = str_replace(pacbio_id, ".1:pep", ""))

domain_dataframe <- read.csv("domains_only.csv") %>%
  select('gene_id', 'Domain1', 'Domain2', 'Domain3', 'Domain4', 'Domain5',
         'Domain6', 'Domain7', 'Domain8', 'Domain9', 'Domain10', 'Domain11',
         'Domain12', 'Domain13', 'Domain14', 'Domain15', 'Domain16', 'Domain17',
         'Domain18') %>%
  mutate(pacbio_id = str_replace(gene_id, "\\.1$", ""))

dataframe <- left_join(tags_dataframe, domain_dataframe, by = 'pacbio_id') %>%
  rowwise() %>%
  mutate(length_domain = sum(!is.na(c_across(Domain1: Domain18)) & c_across(Domain1:Domain18) != "")) %>%
  ungroup() %>%
  filter(normalized_reads > 0.05) %>%
  filter(length_domain != 0)


violin_varL <- ggplot(data = dataframe, aes(x = groupWvsD, y = length_domain)) +
  geom_violin() + 
  geom_beeswarm() +
  labs(x = 'Infection group', y = 'var gene length') +
  theme_bw()

violin_varL

model_varL <- lmer(data = dataframe, length_domain ~ groupWvsD + (1|host))
null_varL <- lmer(data = dataframe, length_domain ~ (1|sample))
anova(model_varL, null_varL)
summary(model_varL)

export_df_recurrent <- dataframe %>%
  select(var_ID, pacbio_id, Domain1, Domain2, Domain3, Domain4, Domain5, Domain6, Domain7,
         Domain8, Domain9, Domain10, Domain11) %>%
  filter(var_ID == "DC01m1c_var14" | var_ID == "DC01m1c_var30" | var_ID == "DC01m1c_var50" |
         var_ID == "DC01m1c_var54" | var_ID == "DC02m1c_var22" | var_ID == "DC02m1c_var32" |
         var_ID == "DC02m1c_var5" | var_ID == "DC04m2c_var11" | var_ID == "DC04m2c_var37")

#write.csv(export_df_recurrent, "table_domains_recurrents.csv")

#### Internal vs external var genes ####

dataframe <- read.csv("localisation_var.csv") %>%
  filter(!is.na(localisation)) %>%
  group_by(sample) %>%
  mutate(count_external = sum(localisation == 'external'),
         total_count = n()) %>%
  distinct(sample, .keep_all = TRUE) %>%
  mutate(prop_external = count_external / total_count) %>%
  filter(!sample %in% c('DC01m5', 'WS10', 'WS13'))

violin_localisation <- ggplot(data = dataframe, aes(x = Group, y = prop_external)) +
  geom_violin() + 
  geom_beeswarm() +
  labs(x = 'Infection group', y = 'proportion of var genes expressed from external cluster') +
  theme_bw() 

violin_localisation
  
model_loca <- lmer(data = dataframe, prop_external ~ Group + (1|host))
null_loca <- lmer(data = dataframe, prop_external ~ (1|host))
anova(model_loca, null_loca)
summary(model_loca)

grid.arrange(violin_varL, violin_localisation, ncol = 2)
