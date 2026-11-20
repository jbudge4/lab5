# Part 1. Configure your project and examine the data

## ----part1-packages
library(tidyverse) # Includes dplyr, ggplot2, and forcats
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(emmeans)
library(broom.mixed)
library(gt)
library(AER)

## ----part1-data
# Load NMES1988 explicitly from AER before creating the analysis dataset.
data("NMES1988", package = "AER")
data <- NMES1988

# Rebuild the factors to clear any contrasts stored in the original dataset.
data <- data |>
  select(visits, chronic, health, gender, school, insurance) |>
  mutate(
    health = fct_relevel(factor(health), "average"),
    gender = fct_relevel(factor(gender), "male"),
    insurance = fct_relevel(factor(insurance), "no")
  )

## ----part1-describe
# Summarize physician office visits.
visits_summary <- data |>
  summarise(
    n      = n(),
    mean   = mean(visits),
    var    = var(visits),
    p_zero = mean(visits == 0),
    max    = max(visits)
  )

## ----part1-table
visits_table <- visits_summary |>
  gt() |>
  tab_header(title = "Descriptive statistics for physician office visits") |>
  cols_label(
    n = "N",
    mean = "Mean",
    var = "Variance",
    p_zero = "Proportion with zero visits",
    max = "Maximum"
  ) |>
  fmt_integer(columns = c(n, max)) |>
  fmt_number(columns = c(mean, var), decimals = 2) |>
  fmt_percent(columns = p_zero, decimals = 2)

## ----part1-plot
descriptive_plot <- ggplot(data, aes(x = visits)) +
  geom_bar() +
  labs(x = "Physician office visits", y = "Number of participants")

## ----part1-console
# Run lab5.R in RStudio to inspect all six variables.
print(head(data))
glimpse(data)
print(summary(data))
print(visits_summary)

if (interactive()) {
  print(descriptive_plot)
}
