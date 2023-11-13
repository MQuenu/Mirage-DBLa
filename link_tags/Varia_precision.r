library(dplyr)
library(ggplot2)
library(stringr)

setwd("~/DBLa/link_tag_genome/R_analyses")

rm(list = ls())

#### data prep and getting the table, only do it once ####

link_dataframe <- read.csv("tags_varIDs_for_R.csv", header = FALSE) %>%
  select('V1', 'V2') %>%
  rename('DBL_tag' = 'V1',
         'gene_id' = 'V2') %>%
  mutate(gene_id = sub(":pep", "", gene_id))

tags_dataframe <- read.csv("ASformatted_table.csv") %>%
  select('DBL_tag', 'Domain1', 'Domain2', 'Domain3', 'Domain4', 'Domain5',
         'Domain6', 'Domain7', 'Domain8', 'Domain9', 'Domain10') %>%
  rename('Domain_tag1' = 'Domain1',
         'Domain_tag2' = 'Domain2',
         'Domain_tag3' = 'Domain3',
         'Domain_tag4' = 'Domain4',
         'Domain_tag5' = 'Domain5',
         'Domain_tag6' = 'Domain6',
         'Domain_tag7' = 'Domain7',
         'Domain_tag8' = 'Domain8',
         'Domain_tag9' = 'Domain9',
         'Domain_tag10' = 'Domain10')

domain_dataframe <- read.csv("domains_only.csv") %>%
  select('gene_id', 'Domain1', 'Domain2', 'Domain3', 'Domain4', 'Domain5',
         'Domain6', 'Domain7', 'Domain8', 'Domain9', 'Domain10') %>%
  rename('Domain_genome1' = 'Domain1',
         'Domain_genome2' = 'Domain2',
         'Domain_genome3' = 'Domain3',
         'Domain_genome4' = 'Domain4',
         'Domain_genome5' = 'Domain5',
         'Domain_genome6' = 'Domain6',
         'Domain_genome7' = 'Domain7',
         'Domain_genome8' = 'Domain8',
         'Domain_genome9' = 'Domain9',
         'Domain_genome10' = 'Domain10')

comparison_dataframe <- link_dataframe %>%
  left_join(tags_dataframe, by = "DBL_tag") %>%
  left_join(domain_dataframe, by = "gene_id") %>%
  distinct(DBL_tag, .keep_all = TRUE)

write.csv(comparison_dataframe, file = "domain_comparison.csv")

#### Domain analysis ####

## format the dataframe for comparisons

simplify_domain <- function(domain){
  first_character = substr(domain, start = 1, stop = 1)
  domain_simplified = case_when(
    first_character == 'A' ~ substr(domain, start = 1, stop = 3),
    first_character == 'N' ~ substr(domain, start = 1, stop = 3),
    first_character == 'D' ~ substr(domain, start = 1, stop = 4),
    first_character == 'C' ~ substr(domain, start = 1, stop = 5)
  )
  return(domain_simplified)
}

comparison_dataframe <- read.csv("domain_comparison.csv") %>%
  filter(!is.na(Domain_tag1) & !is.na(Domain_tag2)) %>%
  mutate(Domain_tag1 = simplify_domain(Domain_tag1),
         Domain_tag2 = simplify_domain(Domain_tag2),
         Domain_tag3 = simplify_domain(Domain_tag3),
         Domain_tag4 = simplify_domain(Domain_tag4),
         Domain_tag5 = simplify_domain(Domain_tag5),
         Domain_tag6 = simplify_domain(Domain_tag6),
         Domain_tag7 = simplify_domain(Domain_tag7),
         Domain_tag8 = simplify_domain(Domain_tag8),
         Domain_tag9 = simplify_domain(Domain_tag9),
         Domain_tag10 = simplify_domain(Domain_tag10)) %>%
  mutate(Domain_match1 = ifelse(Domain_tag1 == Domain_genome1, 1, 0),
         Domain_match2 = ifelse(Domain_tag2 == Domain_genome2, 1, 0),
         Domain_match3 = ifelse(Domain_tag3 == Domain_genome3, 1, 0),
         Domain_match4 = ifelse(Domain_tag4 == Domain_genome4, 1, 0),
         Domain_match5 = ifelse(Domain_tag5 == Domain_genome5, 1, 0),
         Domain_match6 = ifelse(Domain_tag6 == Domain_genome6, 1, 0),
         Domain_match7 = ifelse(Domain_tag7 == Domain_genome7, 1, 0),
         Domain_match8 = ifelse(Domain_tag8 == Domain_genome8, 1, 0),
         Domain_match9 = ifelse(Domain_tag9 == Domain_genome9, 1, 0),
         Domain_match10 = ifelse(Domain_tag10 == Domain_genome10, 1, 0))

## Proportions of correct assignment for each domain

### Get a function to extract the proportion of matches 

Extract_proportion_matches <- function(dataframe, domain) {
  if (!domain %in% colnames(dataframe)) {
    stop(paste("Column", domain, "not found in the dataframe"))
  }
  column_domain <- dataframe[[domain]]
  proportion_of_match <- dataframe %>%
    summarise(Proportion = mean(column_domain == 1, na.rm = TRUE))
  result <- data.frame(Domain = domain, Proportion = as.numeric(proportion_of_match))
  return(result)
}

### Loop the function over all domains 

domains <- c("Domain_match1", "Domain_match2", "Domain_match3", "Domain_match4",
             "Domain_match5", "Domain_match6", "Domain_match7", "Domain_match8",
             "Domain_match9", "Domain_match10")

proportions_df <- data.frame(Domain = character(0), Proportion = numeric(0))

for (domain in domains) {
  proportion <- Extract_proportion_matches(comparison_dataframe, domain)
  proportions_df <- rbind(proportions_df, proportion)
}

proportions_df$Domain <- factor(proportions_df$Domain, levels = domains)

### Use gg

plot_domain_accuracy <- ggplot(proportions_df, aes(y = Proportion, x = Domain, group = Domain)) +
  geom_point() + 
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

plot_domain_accuracy
