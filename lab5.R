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

# Part 2. Specify competing count models

## ----part2-gaussian
m_gauss <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = data,
  family = gaussian(link = "identity")
)

# Use the demonstration's seed for this and later DHARMa checks.
dharma_seed <- 9041
res_gauss <- simulateResiduals(m_gauss, n = 500, plot = FALSE, seed = dharma_seed)

## ----part2-poisson
# M0: Intercept only.
m0_null <- glmmTMB(
  visits ~ 1,
  data = data,
  family = poisson(link = "log")
)

# M1: Number of chronic conditions.
m1_chronic <- glmmTMB(
  visits ~ chronic,
  data = data,
  family = poisson(link = "log")
)

# M2: All five predictors.
m2_full <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = data,
  family = poisson(link = "log")
)

## ----part2-anova
tests <- anova(m0_null, m1_chronic, m2_full)

comparison_table <- tests |>
  as.data.frame() |>
  mutate(Model = c("M0: Null", "M1: Chronic", "M2: Full")) |>
  select(Model, AIC, BIC, logLik, Chisq, `Chi Df`, `Pr(>Chisq)`) |>
  rename(
    Nested.df = `Chi Df`,
    Nested.Chisq.stat = Chisq,
    Nested.p.value = `Pr(>Chisq)`
  )

## ----part2-table
poisson_table <- comparison_table |>
  gt() |>
  tab_header(title = "Likelihood-ratio tests for the nested Poisson models") |>
  cols_label(
    logLik = "Log-likelihood",
    Nested.df = "df",
    Nested.Chisq.stat = "χ²",
    Nested.p.value = "p-value"
  ) |>
  fmt_number(
    columns = c(AIC, BIC, logLik, Nested.Chisq.stat),
    decimals = 2
  ) |>
  fmt(
    columns = Nested.p.value,
    fns = function(x) format.pval(x, digits = 3, eps = .001)
  ) |>
  sub_missing(
    columns = c(Nested.df, Nested.Chisq.stat, Nested.p.value),
    missing_text = "—"
  ) |>
  tab_source_note("Each test compares the model in that row with the model above it.")

## ----part2-console
print(tests)

if (interactive()) {
  plot(res_gauss)
}

# Part 3. Evaluate the Poisson model

## ----part3-dharma
# Use the same seed and number of simulations as for the Gaussian model.
res_pois <- simulateResiduals(m2_full, n = 500, plot = FALSE, seed = dharma_seed)

## ----part3-tests
disp_pois <- testDispersion(res_pois, plot = FALSE)
zi_pois <- testZeroInflation(res_pois, plot = FALSE)

## ----part3-console
print(disp_pois)
print(zi_pois)

if (interactive()) {
  plot(res_pois)
}

# Part 4. Fit and evaluate a negative binomial model

## ----part4-nb
m_nb <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = data,
  family = nbinom2(link = "log")
)

## ----part4-anova
lrt_nb <- anova(m2_full, m_nb)

## ----part4-dharma
# Set the shared seed before calling DHARMa so clean R and Quarto runs agree.
set.seed(dharma_seed)
res_nb <- simulateResiduals(m_nb, n = 500, plot = FALSE, seed = NULL)

## ----part4-tests
disp_nb <- testDispersion(res_nb, plot = FALSE)
zi_nb <- testZeroInflation(res_nb, plot = FALSE)

## ----part4-predictions
# Predict mean visits at every observed number of chronic conditions.
pred_nb <- predict_response(
  m_nb,
  terms = "chronic [all]",
  condition = c(school = mean(data$school))
)

pred_nb_plot <- plot(pred_nb) +
  labs(
    title = "Predicted physician office visits by chronic conditions",
    x = "Number of chronic conditions",
    y = "Predicted mean visits"
  )

## ----part4-console
print(summary(m_nb))
print(lrt_nb)
print(disp_nb)
print(zi_nb)
print(pred_nb, n = Inf)

if (interactive()) {
  plot(res_nb)
  print(pred_nb_plot)
}

# Part 5. Decide whether the zeros need their own model

## ----part5-zero-models
# Zero-inflated NB with the same predictors in the count and zero parts.
m_zinb_full <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  ziformula = ~ chronic + health + gender + school + insurance,
  data = data,
  family = nbinom2(link = "log")
)

# Hurdle NB uses the zero part for all zeros and a truncated model for positive counts.
m_hurd <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  ziformula = ~ chronic + health + gender + school + insurance,
  data = data,
  family = truncated_nbinom2(link = "log")
)

## ----part5-compare
lrt_zinb <- anova(m_nb, m_zinb_full)

aic_zero <- AIC(m_nb, m_zinb_full, m_hurd)
bic_zero <- BIC(m_nb, m_zinb_full, m_hurd)

zero_comparison <- tibble(
  Model = c("NB2", "ZINB", "Hurdle NB"),
  AIC = aic_zero$AIC,
  BIC = bic_zero$BIC
)

zero_comparison_table <- zero_comparison |>
  gt() |>
  tab_header(title = "Fit comparison for the negative binomial models") |>
  fmt_number(columns = c(AIC, BIC), decimals = 2)

## ----part5-dharma
# Re-check the hurdle model, which has the lowest AIC and BIC.
set.seed(dharma_seed)
res_hurd <- simulateResiduals(m_hurd, n = 500, plot = FALSE, seed = NULL)
disp_hurd <- testDispersion(res_hurd, plot = FALSE)
zi_hurd <- testZeroInflation(res_hurd, plot = FALSE)

## ----part5-console
print(lrt_zinb)
print(zero_comparison)
print(disp_hurd)
print(zi_hurd)

if (interactive()) {
  plot(res_hurd)
}
