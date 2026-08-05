library(tidyverse)

population_data_2 <- read_csv(file = "population-by-state.csv")

summarize(.data = population_data_2, mean_population = mean(Pop))

population_data_2 %>%
  summarize(mean_population = mean(Pop))

population_data_2 %>%
  filter(rank <= 5) %>%
  summarize(mean_population = mean(Pop))

# Calculate the mean population of the five largest states
population_data_2 %>%
  filter(rank <= 5) %>%
  summarize(mean_population = mean(Pop))