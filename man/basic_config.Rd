% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recipe-api.R
\name{basic_config}
\alias{basic_config}
\title{Create a basic config for HTE estimation}
\usage{
basic_config()
}
\value{
\code{HTE_cfg} object
}
\description{
This provides a basic recipe for HTE estimation that can
be extended by providing additional information about models
to be estimated and what quantities of interest should be
returned based on those models. This basic model includes
only linear models for nuisance function estimation, and
basic diagnostics.
}
\details{
Additional models, diagnostics and quantities of interest should
be added using their respective helper functions provided as part
of the Recipe API.

To see an example analysis, read \code{vignette("experimental_analysis")} in the context
of an experiment, \code{vignette("experimental_analysis")} for an observational study, or
\code{vignette("methodological_details")} for a deeper dive under the hood.
}
\examples{
library("dplyr")
basic_config() \%>\%
   add_known_propensity_score("ps") \%>\%
   add_outcome_model("SL.glm.interaction") \%>\%
   add_outcome_model("SL.glmnet", alpha = c(0.05, 0.15, 0.2, 0.25, 0.5, 0.75)) \%>\%
   add_outcome_model("SL.glmnet.interaction", alpha = c(0.05, 0.15, 0.2, 0.25, 0.5, 0.75)) \%>\%
   add_outcome_diagnostic("RROC") \%>\%
   add_effect_model("SL.glm.interaction") \%>\%
   add_effect_model("SL.glmnet", alpha = c(0.05, 0.15, 0.2, 0.25, 0.5, 0.75)) \%>\%
   add_effect_model("SL.glmnet.interaction", alpha = c(0.05, 0.15, 0.2, 0.25, 0.5, 0.75)) \%>\%
   add_effect_diagnostic("RROC") \%>\%
   add_moderator("Stratified", x2, x3) \%>\%
   add_moderator("KernelSmooth", x1, x4, x5) \%>\%
   add_vimp(sample_splitting = FALSE) -> hte_cfg
}
\seealso{
\code{\link[=add_propensity_score_model]{add_propensity_score_model()}}, \code{\link[=add_known_propensity_score]{add_known_propensity_score()}},
\code{\link[=add_propensity_diagnostic]{add_propensity_diagnostic()}}, \code{\link[=add_outcome_model]{add_outcome_model()}}, \code{\link[=add_outcome_diagnostic]{add_outcome_diagnostic()}},
\code{\link[=add_effect_model]{add_effect_model()}}, \code{\link[=add_effect_diagnostic]{add_effect_diagnostic()}}, \code{\link[=add_moderator]{add_moderator()}}, \code{\link[=add_vimp]{add_vimp()}}
}
