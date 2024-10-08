% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_sup.R
\name{plot.probsup}
\alias{plot.probsup}
\title{Plots for the \code{probsup} object}
\usage{
\method{plot}{probsup}(x, ..., category = "perfo", level = "across")
}
\arguments{
\item{x}{An object of class \code{probsup}.}

\item{...}{currently not used}

\item{category}{A string indicating which plot to build. See options in the Details section.}

\item{level}{A string indicating the information level to be used for building
the plots. Options are \code{"across"} for focusing on the probabilities across environments,
or \code{"within"} to focus on the within-environment effects. Defaults to \code{"across"}.}
}
\description{
Build plots using the outputs stored in the \code{probsup} object.
}
\details{
The available options are:
\itemize{
\item \code{hpd} : a caterpillar plot representing the marginal genotypic value of
each genotype, and their respective highest posterior density interval (95\% represented by the
thick line, and 97.5\% represented by the thin line). Available only if \code{level = "across"}.
\item \code{perfo} : if \code{level = "across"}, a lollipop plot illustrating the probabilities of superior performance.
If \code{level = "within"}, a heatmap with the probabilities of superior performance within
environments. If a model with \code{reg} and/or \code{year} is fitted, multiple plots are produced.
\item \code{stabi}: a lollipop plot with the probabilities of superior stability.
If a model with \code{reg} and/or \code{year} is fitted, multiple plots are produced.
Available only if \code{level = "across"}. Unavailable if an entry-mean model was used in \code{bayes_met}.
\item \code{pair_perfo} : if \code{level = "across"}, a heatmap representing the pairwise probability of superior
performance (the probability of genotypes at the \emph{x}-axis being superior.
to those on the \emph{y}-axis). If \code{level = "within"}, a list of heatmaps representing the pairwise probability of superior
performance within environments.  If a model with \code{reg} and/or \code{year} is fitted, multiple plots (and multiple lists) are produced.
Should this option is set, it is mandatory to store the outputs in an object.
(e.g., \code{pl <- plot(obj, category = "pair_perfo", level = "within")}) so they can be visualized one at a time.
The option \code{level = "within"} is unavailable if an entry-mean model was used in \code{bayes_met}.
\item \code{pair_stabi}: a heatmap with the pairwise probabilities of superior stability
(the probability of genotypes at the \emph{x}-axis being more stable than those on the \emph{y}-axis).
If a model with \code{reg} and/or \code{year} is fitted, multiple plots are produced. Available only if \code{level = "across"}.
Unavailable if an entry-mean model was used in \code{bayes_met}.
\item \code{joint}: a lollipop plot with the joint probabilities of superior performance and stability. Unavailable if an entry-mean model was used in \code{bayes_met}.
}
}
\examples{
\donttest{
mod = bayes_met(data = maize,
                gen = "Hybrid",
                loc = "Location",
                repl = c("Rep","Block"),
                trait = "GY",
                reg = "Region",
                year = NULL,
                res.het = TRUE,
                iter = 2000, cores = 2, chain = 4)

outs = extr_outs(model = mod,
                 probs = c(0.05, 0.95),
                 verbose = TRUE)

results = prob_sup(extr = outs,
                   int = .2,
                   increase = TRUE,
                   save.df = FALSE,
                   verbose = FALSE)

plot(results, category = "hpd")
plot(results, category = "perfo", level = "across")
plot(results, category = "perfo", level = "within")
plot(results, category = "stabi")
plot(results, category = "pair_perfo", level = "across")
plwithin = plot(results, category = "pair_perfo", level = "within")
plot(results, category = "pair_stabi")
plot(results, category = "joint")
}

}
\seealso{
\link{prob_sup}
}
