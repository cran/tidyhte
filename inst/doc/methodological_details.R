## ----setup, echo=FALSE--------------------------------------------------------
suppressMessages({
library(palmerpenguins)
library(dplyr)
library(tidyhte)
library(magrittr)
library(SuperLearner)
library(ggplot2)

theme_set(theme_minimal())
})

set.seed(100)
n <- nrow(penguins)

## ----dgp----------------------------------------------------------------------
penguins <- within(penguins, {
  id <- 1:n
  propensity_score <- 0.5
  treatment <- rbinom(n, 1, propensity_score)
  tau <- 0.2 * (species == "Gentoo") - 0.2 * (species == "Adelie") + rnorm(n, sd = 0.05)
  food_consumed_g <- rnorm(n, 500, 5) * (1 + tau * treatment)
  body_mass_g <- body_mass_g * (1 + tau * treatment)
})

## ----cfg----------------------------------------------------------------------
cfg <- basic_config() %>%
    add_known_propensity_score("propensity_score") %>%
    add_outcome_model("SL.glmnet", alpha = c(0.0, 1.0)) %>%
    add_moderator("Stratified", species, island, sex, year) %>%
    add_moderator("KernelSmooth", bill_length_mm, bill_depth_mm, flipper_length_mm) %>%
    add_vimp(sample_splitting = FALSE)

## ----attach-------------------------------------------------------------------
penguins %<>% attach_config(cfg)

## ----split--------------------------------------------------------------------
penguins %<>% make_splits(id, species, sex, flipper_length_mm, .num_splits = 3)

## ----sl, eval = FALSE---------------------------------------------------------
#  learners <- create.Learner(
#      "SL.glmnet",
#      tune = list(
#          alpha = c(0.05, 0.15, 0.2, 0.25, 0.5, 0.75)
#      ),
#      detailed_names = TRUE,
#      name_prefix = paste0("SLglmnet")
#  )
#  
#  CV.SuperLearner(label, covariates, SL.library = learners$names)

## ----eval = FALSE-------------------------------------------------------------
#  add_outcome_model(cfg, "SL.glmnet", alpha = c(0.0, 0.25, 0.5, 0.75, 1.0))

## ----plugins, results='hide', message=FALSE, warning=FALSE--------------------
penguins %<>% produce_plugin_estimates(
  # outcome
  food_consumed_g,
  # treatment
  treatment,
  # covariates
  species, island, sex, year, bill_length_mm, bill_depth_mm, flipper_length_mm
)

## ----psi, results='hide', message=FALSE, warning=FALSE------------------------
penguins %<>% construct_pseudo_outcomes(food_consumed_g, treatment)

## ----qoi, results='hide', message=FALSE, warning=FALSE------------------------
penguins %>%
  estimate_QoI(species, island, sex, year, bill_length_mm, bill_depth_mm, flipper_length_mm) ->
  results

## ----show_results-------------------------------------------------------------
results

## ----discrete_mcate, fig.height = 4, fig.width = 8----------------------------
filter(results, estimand == "MCATE", is.na(value)) %>%
ggplot(aes(level, estimate)) +
geom_point() +
geom_linerange(aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error)) +
geom_hline(yintercept = 0, linetype = "dashed") +
coord_flip() +
facet_wrap(~term, scales = "free_y")

## ----cts_mcate, fig.height = 6, fig.width = 8---------------------------------
filter(results, estimand == "MCATE", is.na(level)) %>%
ggplot(aes(value, estimate)) +
geom_line() +
geom_ribbon(
  aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
  alpha = 0.5
) +
geom_hline(yintercept = 0, linetype = "dashed") +
scale_x_continuous("Covariate value") +
scale_y_continuous("CATE") +
coord_flip() +
facet_wrap(~term, scales = "free_y")

## ----risk, fig.height = 2, fig.width = 8--------------------------------------
filter(results, estimand == "SL risk") %>%
ggplot(aes(reorder(term, estimate), estimate)) +
geom_point() +
geom_linerange(
  aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
  alpha = 0.5
) +
geom_hline(yintercept = 0, linetype = "dashed") +
scale_x_discrete("") +
scale_y_continuous("Risk") +
facet_wrap(~level) +
coord_flip()

## ----coef, fig.height = 2, fig.width = 8--------------------------------------
filter(results, estimand == "SL coefficient") %>%
ggplot(aes(reorder(term, estimate), estimate)) +
geom_point() +
geom_linerange(
  aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
  alpha = 0.5
) +
geom_hline(yintercept = 0, linetype = "dashed") +
scale_x_discrete("") +
scale_y_continuous("Coefficient") +
facet_wrap(~level) +
coord_flip()

## ----vimp, fig.height = 4, fig.width = 8--------------------------------------
filter(results, estimand == "VIMP") %>%
ggplot(aes(reorder(term, estimate), estimate)) +
geom_point() +
geom_linerange(
  aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
  alpha = 0.5
) +
geom_hline(yintercept = 0, linetype = "dashed") +
scale_x_discrete("") +
scale_y_continuous("Reduction in R²") +
coord_flip()

## ----repeat, message=FALSE,warning=FALSE--------------------------------------
penguins %<>% produce_plugin_estimates(
  # outcome
  body_mass_g,
  # treatment
  treatment,
  # covariates
  species, island, sex, year, bill_length_mm, bill_depth_mm, flipper_length_mm
) %>%
construct_pseudo_outcomes(body_mass_g, treatment) %>%
estimate_QoI(
  species, island, sex, year, bill_length_mm, bill_depth_mm, flipper_length_mm
) -> results_mass

## ----combine------------------------------------------------------------------
results_all <- bind_rows(
  results %>% mutate(outcome = "food_consumed_g"),
  results_mass %>% mutate(outcome = "body_mass_g")
)

