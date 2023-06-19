library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(entropy)
library(gridExtra)

setwd("~/DBLa/R_analyses")

rm(list = ls())

#### Correlation between the number of detected dbla and total read count / parasitaemia

raw_df <- read.csv("ASformatted_table.txt") 

write.csv(raw_df, "ASformatted_table.txt")

parasitaemia_df <- read.csv("parasitaemia.csv", header= TRUE) %>%
  mutate(Venous_parasitaemia = as.numeric(Venous_parasitaemia), Microscopy = as.numeric(Microscopy))

nbdbl_df <- aggregate(DBL_tag ~ sample, data = raw_df, FUN = length)

sample_df <- aggregate(Read_count ~ sample, data = raw_df, FUN = sum) 

group_df <- aggregate(raw_df$group, by = list(raw_df$sample, raw_df$group), FUN = identity) %>%
  select(-x) %>%
  rename(sample = Group.1, group = Group.2)

merged_df <- merge(nbdbl_df, sample_df) %>%
  merge(group_df) %>%
  left_join(parasitaemia_df, by = "sample")

plot1 = ggplot(merged_df, aes(x = DBL_tag, y = read_count, color = group)) +
  geom_point() +
  xlab(label = "Number of DBLas") +
  ylab(label = "Total read count") +
  theme_classic()

plot1

plot2 = ggplot(merged_df, aes(x = Venous_parasitaemia, y = read_count, color = group)) +
  geom_point() +
  xlab(label = "Venous parasitaemia") +
  ylab(label = "Total read count") +
  #scale_x_log10() +
  theme_classic()

plot2

#### Histogram of read_counts per DBLa tag

rm(list = ls())

DBL_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  group_by(DBL_tag) %>%
  summarise(total_reads = sum(Read_count)) %>%
  arrange(desc(total_reads)) 

barplot1 <- ggplot(data = DBL_dataframe, aes(x = reorder(DBL_tag,total_reads), y = total_reads)) + 
  geom_col() + 
  scale_y_log10() +
  scale_x_discrete(labels = NULL) +
  theme_classic()

barplot1

#### Plots of read counts per sample

rm(list = ls())

## DC01 and DC04

DC01_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "DC01_m1" | sample == "DC01_m3" |  sample == "DC01_m5" |
           sample == "DC01_m6")

plot_DC01 = ggplot(DC01_dataframe, aes(x = DBL_tag, y = Read_count, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC01

DC04_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC04_m2" | sample == "DC04_m3" |  sample == "DC04_m4" |
           sample == "DC04_m5" | sample == "DC04_m6")

plot_DC04 = ggplot(DC04_dataframe, aes(x = DBL_tag, y = Read_count, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC04

DC05_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC05_m1" | sample == "DC05_m3" |  sample == "DC05_m5" |
           sample == "DC05_m6")

plot_DC05 = ggplot(DC05_dataframe, aes(x = DBL_tag, y = Read_count, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC05

DC08_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  mutate(NTS = substr(Domain1, 1, 4)) %>%
  filter(sample == "DC08_m2" | sample == "DC08_m3" |  sample == "DC08_m4" |
           sample == "DC08_m5" | sample == "DC08_m6")

plot_DC08 = ggplot(DC08_dataframe, aes(x = DBL_tag, y = Read_count, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC08

DC07_dataframe <- read.csv("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "DC07_m1" | sample == "DC07_m2" |  sample == "DC07_m3" |
           sample == "DC07_m5" | sample == "DC07_m6")

plot_DC07 = ggplot(DC07_dataframe, aes(x = DBL_tag, y = Read_count, fill = NTS)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_DC07

## WS

WS1_dataframe <- read.table("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "WS05" | sample == "WS06" |  sample == "WS07" |
           sample == "WS08"| sample == "WS09"| sample == "WS10")

plot_WS1 = ggplot(WS1_dataframe, aes(x = DBL_tag, y = read_count)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 6) +
  theme_minimal()

plot_WS1

WS2_dataframe <- read.table("ASformatted_table.txt", header = TRUE) %>%
  filter(sample == "WS10" | sample == "WS11" |  sample == "WS12" |
           sample == "WS13"| sample == "WS14"| sample == "WS15")

plot_WS2 = ggplot(WS2_dataframe, aes(x = DBL_tag, y = read_count)) +
  geom_col(width=0.7) +
  facet_wrap(~ sample, nrow = 6) +
  theme_minimal()


### comparison between batches

rm(list = ls())

setwd("~/Varia/DBLa/Results_processed/Batch_comparisons")

batch_DC01m1_dataframe <- read.table("formatted_table.txt", header = TRUE)

plot_batch_DC01m1 <- ggplot(data = batch_DC01m1_dataframe, aes(x = DBL_tag, y = read_count)) +
  geom_col(width=0.7) +
  facet_wrap(~ batch, nrow = 6) +
  theme_minimal()

plot_batch_DC01m1

### Look for shared DBLa

rm(list = ls())

raw_dataframe <- read.csv("ASformatted_table.txt", header = TRUE)

count_dataframe <- raw_dataframe %>%
  group_by(DBL_tag) %>%
  summarize(count = n())

# look through count_dataframe interactively

## Entropy by group

rm(list = ls())

raw_dataframe <- read.table("ASformatted_table.txt", header = TRUE) %>% 
  group_by(sample) %>%
  mutate(entropy = entropy(read_count))

entropy_dataframe <- aggregate(entropy ~ sample + group, data = raw_dataframe, FUN = mean)

plot_entropy <- ggplot(entropy_dataframe, aes(x = group, y = entropy, fill = group)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_brewer() +
  theme_classic()

plot_entropy

null_model <- glm(entropy ~ 1, data = entropy_dataframe)
linear_model <- glm(entropy ~ group, data = entropy_dataframe)

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

### Recurring DBLa genes analysis

rm(list = ls())

### Create a dataframe with a column for 'recurring' genes

chronic_dataframe <- read.csv("ASformatted_table.txt",header = TRUE) %>%
  mutate(infection_group = substr(sample, start = 1, stop = 2)) %>%
  mutate(host = substr(sample, start = 1, stop = 4)) %>% 
  subset(infection_group == "DC") %>%
  group_by(host) %>%
  group_by(DBL_tag, add = TRUE) %>%
  mutate(recurrence = ifelse(n() > 1, "recurrent", "single")) %>%
  mutate(expression = ifelse(Read_count < 9, "background", "expressed")) %>%
  ungroup() %>%
  filter(Read_count > 1)

### Proportions of recurring and single genes over the total number of genes

counts_recurrent_expressed <- chronic_dataframe %>%
  filter(expression == "background") #%>%
  filter(recurrence == "recurrent")
