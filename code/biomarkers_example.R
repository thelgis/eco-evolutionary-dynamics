library(tidyverse)
library(here)

data <- read.csv(here("data/GSE53987_combined.csv")) %>%  as_tibble() %>% 
  dplyr::rename(sample=1)

### data 

data %>% 
  summarise(
    unique_sample = n_distinct(sample),
    unique_Patient = n_distinct(Patient),
    unique_disp = n_distinct(Source.name),
    unique_tissue=n_distinct(Tissue),
    unique_phenotype=n_distinct(Disease.state)
  )

data %>% 
  ggplot() +
  geom_histogram(aes(x=Age)) +
  geom_density(aes(x=Age))

head(data)

library("GGally")
data %>%
  select(Disease.state, Age, Gender, Race, Pmi) %>% 
  GGally::ggpairs(aes(color = Disease.state)) 