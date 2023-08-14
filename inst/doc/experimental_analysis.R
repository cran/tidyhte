## ----setup, print=FALSE, message=FALSE----------------------------------------
library(tidyhte)
library(ggplot2)
library(dplyr)

## ----sim_data-----------------------------------------------------------------
set.seed(100)
n <- 500
data <- tibble(
    uid = 1:n
) %>%
    mutate(
        a = rbinom(n, 1, 0.3),
        ps = rep(0.3, n),
        x1 = rnorm(n),
        x2 = factor(sample(1:4, n, prob = c(1 / 100, 39 / 100, 1 / 5, 2 / 5), replace = TRUE)),
        x3 = factor(sample(1:3, n, prob = c(1 / 5, 1 / 5, 3 / 5), replace = TRUE)),
        x4 = (x1 + rnorm(n)) / 2,
        x5 = rnorm(n),
        y = (
            a + x1 - a * (x1 - mean(x1)) + (4 * rbinom(n, 1, 0.5) - 1) * a * (x2 == 2) +
            a * (x2 == 3) + 0.5 * a * (x2 == 4) +
            0.25 * rnorm(n)
        ),
        w = 0.1 + rexp(n, 1 / 0.9)
    )

## ----recipe-------------------------------------------------------------------
basic_config() %>%
    add_known_propensity_score("ps") %>%
    add_outcome_model("SL.glm.interaction") %>%
    add_outcome_model("SL.glmnet", alpha = c(0, 1)) %>%
    add_outcome_model("SL.glmnet.interaction", alpha = c(0, 1)) %>%
    add_outcome_diagnostic("RROC") %>%
    add_effect_model("SL.glm.interaction") %>%
    add_effect_model("SL.glmnet", alpha = c(0, 1)) %>%
    add_effect_model("SL.glmnet.interaction", alpha = c(0, 1)) %>%
    add_effect_diagnostic("RROC") %>%
    add_moderator("Stratified", x2, x3) %>%
    add_moderator("KernelSmooth", x1, x4, x5) %>%
    add_vimp(sample_splitting = FALSE) ->
    hte_cfg

## ----estimate, message=FALSE--------------------------------------------------
data %>%
    attach_config(hte_cfg) %>%
    make_splits(uid, .num_splits = 3) %>%
    produce_plugin_estimates(
        y,
        a,
        x1, x2, x3, x4, x5,
    ) %>%
    construct_pseudo_outcomes(y, a) -> prepped_data

prepped_data %>%
    estimate_QoI(x1, x2, x3, x4, x5) -> results

## ----show_qoi, message=FALSE--------------------------------------------------
results

## ----ates---------------------------------------------------------------------
filter(results, grepl("SATE|PATE", estimand))

## ----sl_coef------------------------------------------------------------------
filter(results, grepl("SL coefficient", estimand)) %>%
mutate(level = factor(level, levels = c("Control Response", "Treatment Response"))) %>%
ggplot(aes(
            x = reorder(term, estimate),
            y = estimate,
            ymin = estimate - 1.96 * std_error,
            ymax = estimate + 1.96 * std_error
    )) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    geom_pointrange() +
    expand_limits(y = 0) +
    scale_x_discrete("Model name") +
    scale_y_continuous("Coefficient in SuperLearner Ensemble") +
    facet_wrap(~level) +
    coord_flip() +
    ggtitle("SuperLearner Ensemble") +
    theme_minimal()

## ----sl_risk------------------------------------------------------------------
filter(results, grepl("SL risk", estimand)) %>%
mutate(
    level = factor(level, levels = c("Control Response", "Treatment Response", "Effect Surface"))
) %>%
ggplot() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    geom_pointrange(
        aes(
            x = reorder(term, -estimate),
            y = estimate,
            ymin = estimate - 1.96 * std_error,
            ymax = estimate + 1.96 * std_error)
        ) +
    expand_limits(y = 0) +
    scale_x_discrete("Model name") +
    scale_y_continuous("CV Risk in SuperLearner Ensemble") +
    facet_wrap(~level, scales = "free_x") +
    coord_flip() +
    ggtitle("Submodel Risk Estimates") +
    theme_minimal()

## ----rroc---------------------------------------------------------------------
filter(results, grepl("RROC", estimand)) %>%
mutate(
    level = factor(level, levels = c("Control Response", "Treatment Response", "Effect Surface"))
) %>%
ggplot() +
    geom_line(
        aes(
            x = value,
            y = estimate
        )
    ) +
    geom_point(
        aes(x = value, y = estimate),
        data = filter(results, grepl("RROC", estimand)) %>% group_by(level) %>% slice_head(n = 1)
    ) +
    expand_limits(y = 0) +
    scale_x_continuous("Over-estimation") +
    scale_y_continuous("Under-estimation") +
    facet_wrap(~level, scales = "free_x") +
    coord_flip() +
    ggtitle("Regression ROC Curves") +
    theme_minimal()

## ----vimp---------------------------------------------------------------------
ggplot(filter(results, estimand == "VIMP")) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    geom_pointrange(
        aes(
            x = term,
            y = estimate,
            ymin = estimate - 1.96 * std_error,
            ymax = estimate + 1.96 * std_error
        )
    ) +
    expand_limits(y = 0) +
    scale_x_discrete("Covariate") +
    scale_y_continuous("Reduction in R² from full model") +
    coord_flip() +
    ggtitle("Covariate Importance") +
    theme_minimal()

## ----cts_mcate_plot, message=FALSE--------------------------------------------
for (cov in c("x1", "x4", "x5")) {
    ggplot(filter(results, estimand == "MCATE", term == cov)) +
        geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
        geom_ribbon(
            aes(
                x = value,
                ymin = estimate - 1.96 * std_error,
                ymax = estimate + 1.96 * std_error
            ),
            alpha = 0.75
        ) +
        geom_line(
            aes(x = value, y = estimate)
        ) +
        expand_limits(y = 0) +
        scale_x_continuous("Covariate level") +
        scale_y_continuous("CATE") +
        ggtitle(paste("Marginal effects across", cov)) +
        theme_minimal() -> gp
    print(gp)
}

## ----discrete_mcate_plot------------------------------------------------------
for (cov in c("x2", "x3")) {
    ggplot(filter(results, estimand == "MCATE", term == cov)) +
        geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
        geom_pointrange(
            aes(
                x = level,
                y = estimate,
                ymin = estimate - 1.96 * std_error,
                ymax = estimate + 1.96 * std_error
            )
        ) +
        expand_limits(y = 0) +
        scale_x_discrete("Covariate level") +
        scale_y_continuous("CATE") +
        ggtitle(paste("Marginal effects across", cov)) +
        theme_minimal() -> gp
    print(gp)
}

## ----session_info-------------------------------------------------------------
print(sessionInfo())

