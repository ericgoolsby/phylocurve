\name{compare.models}
\alias{compare.models}
\title{
Compare evolutionary hypotheses
}
\description{
This function compares two evolutionary hypotheses via phylogenetic simulation under the null model (\code{model1}) and alternative model (\code{model2}) (Boettiger et al. 2012; Goolsby 2016). Typically the alternative model is nested within the null model, but is not a requirement of this testing procedure. By default, data are randomly simulated under model1 and model2, and model1 and model2 are fit to generate a distribution of likelihood ratios for data simulated under the null hypothesis and a distribution of likelihood ratios from data simulated under the alternative hypothesis. Comparison of the observed likelihood ratio to the null distribution provides an estimate of the P-value, whereas comparison of the alternative distribution to the critical value of the test statistic provides an estimate of statistical power.
}
\usage{
compare.models(model1, model2, nsim = 1000,
plot = TRUE, estimate_power = TRUE, parallel = TRUE,
conf.int = TRUE)
}
\arguments{
  \item{model1}{
The null model (an object of class \code{rate.mult}).
}
  \item{model2}{
The alternative model (an object of class \code{rate.mult}).
}
  \item{nsim}{
The number of iterations for phylogenetic simulation.
}
  \item{plot}{
Whether or not to plot the null and alternative distributions of likelihood ratios.
}
  \item{estimate_power}{
Whether to simulate the alternative distribution (default=TRUE).
}
  \item{parallel}{
Whether to use parallel processing to speed up computations (default=TRUE).
}
  \item{conf.int}{
Whether to estimate confidence inervals for tree transformation parameters in the alternative model (e.g., OU, EB, lambda, kappa, delta). If TRUE and the alternative evolutionary model is not 'BM', estimate_power is automatically set to TRUE.
}
}
\value{
An object of class \code{compare.models}.
}
\references{
Boettiger C., Coop G., Ralph P. 2012. Is your phylogeny informative? Measuring the power of comparative methods. Evolution. 66:2240-2251.

Golsby E.W. 2015. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link[phylocurve]{evo.model}}
}
\examples{
rand.data <- sim.traits()
X <- rowMeans(rand.data$trait_data)
null.model <- evo.model(tree = rand.data$tree,
  Y = rand.data$trait_data,method = "Pairwise ML")
alt.model <- evo.model(tree = rand.data$tree,
  Y = rand.data$trait_data,fixed.effects = X,method = "Pairwise ML")

compare.models(model1 = null.model,model2 = alt.model,
  nsim = 100,parallel = FALSE)
}
