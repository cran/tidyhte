% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_cfg.R
\name{PCATE_cfg}
\alias{PCATE_cfg}
\title{Configuration of Partial CATEs}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
\code{PCATE_cfg} is a configuration class for estimating marginal
response surfaces based on heterogeneous treatment effect estimates.
"Partial" in this context is used similarly to the use in partial
dependence plots or in partial regression. In essence, a PCATE
attempts to partial out the contribution to the CATE from all other
covariates. Two highly correlated variables may have very different
PCATE surfaces.
}
\examples{
PCATE_cfg$new(
   cfgs = list(x1 = KernelSmooth_cfg$new(neval = 100)),
   model_covariates = c("x1", "x2", "x3"),
   num_mc_samples = list(x1 = 100)
)

## ------------------------------------------------
## Method `PCATE_cfg$new`
## ------------------------------------------------

PCATE_cfg$new(
   cfgs = list(x1 = KernelSmooth_cfg$new(neval = 100)),
   model_covariates = c("x1", "x2", "x3"),
   num_mc_samples = list(x1 = 100)
)

## ------------------------------------------------
## Method `PCATE_cfg$add_moderator`
## ------------------------------------------------

cfg <- PCATE_cfg$new(
   cfgs = list(x1 = KernelSmooth_cfg$new(neval = 100)),
   model_covariates = c("x1", "x2", "x3"),
   num_mc_samples = list(x1 = 100)
)
cfg <- cfg$add_moderator("x2", KernelSmooth_cfg$new(neval = 100))
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{cfgs}}{Named list of covariates names to a \code{Model_cfg} object defining
how to present that covariate's CATE surface.}

\item{\code{model_covariates}}{A character vector of all the covariates
to be included in the second-level effect regression.}

\item{\code{num_mc_samples}}{A named list from covariate name to the number
of Monte Carlo samples to take to calculate the double integral (See Details).}

\item{\code{estimand}}{String indicating the estimand to target.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-PCATE_cfg-new}{\code{PCATE_cfg$new()}}
\item \href{#method-PCATE_cfg-add_moderator}{\code{PCATE_cfg$add_moderator()}}
\item \href{#method-PCATE_cfg-clone}{\code{PCATE_cfg$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PCATE_cfg-new"></a>}}
\if{latex}{\out{\hypertarget{method-PCATE_cfg-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PCATE_cfg} object with specified model name and hyperparameters.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PCATE_cfg$new(model_covariates, cfgs, num_mc_samples = 100)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model_covariates}}{A character vector of all the covariates to be
included in the second-level effect regression.}

\item{\code{cfgs}}{Named list from moderator name to a \code{Model_cfg} object defining how to
present that covariate's CATE surface.}

\item{\code{num_mc_samples}}{A named list from covariate name to the number of Monte Carlo
samples to take to calculate the double integral (See Details). If all covariates
should use the same number of samples, simply pass the (integer) number of samples.}

\item{\code{effect_cfg}}{A \code{Model_cfg} object indicating how to fit the second level effect
regression (joint across all selected covariates).}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PCATE_cfg} object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{PCATE_cfg$new(
   cfgs = list(x1 = KernelSmooth_cfg$new(neval = 100)),
   model_covariates = c("x1", "x2", "x3"),
   num_mc_samples = list(x1 = 100)
)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PCATE_cfg-add_moderator"></a>}}
\if{latex}{\out{\hypertarget{method-PCATE_cfg-add_moderator}{}}}
\subsection{Method \code{add_moderator()}}{
Add a moderator to the \code{PCATE_cfg} object. This entails adding it to the joint
model of effects and defines a configuration for displaying the effect surface
for that moderator.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PCATE_cfg$add_moderator(var_name, cfg)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{var_name}}{The name of the moderator to add (and the name of the column in
the dataset).}

\item{\code{cfg}}{A \code{Model_cfg} defining how to display the selected moderator's effect
surface.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An updated \code{PCATE_cfg} object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{cfg <- PCATE_cfg$new(
   cfgs = list(x1 = KernelSmooth_cfg$new(neval = 100)),
   model_covariates = c("x1", "x2", "x3"),
   num_mc_samples = list(x1 = 100)
)
cfg <- cfg$add_moderator("x2", KernelSmooth_cfg$new(neval = 100))
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PCATE_cfg-clone"></a>}}
\if{latex}{\out{\hypertarget{method-PCATE_cfg-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PCATE_cfg$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
