\name{phylocurve-package}
\alias{phylocurve-package}
\docType{package}
\title{
Phylogenetic Comparative Methods for Function-Valued and Other High-Dimensional Traits
}
\description{
Tools for studying the evolution of high-dimensional traits (in particular, function-valued traits), including ancestral state reconstruction, estimating phylogenetic signal, and assessing correlated trait evolution.
}
\details{
\tabular{ll}{
Package: \tab phylocurve\cr
Type: \tab Package\cr
Version: \tab 2.0.4\cr
Date: \tab 2016-05-12\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Eric W. Goolsby <eric.goolsby.evolution@gmail.com>
}
\references{
Goolsby, E.W. 2015. "Phylogenetic comparative methods for evaluating the evolutionary history of function-valued traits." Systematic Biology. In press.

Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional multivariate data. Evolution. 68:2675-2688.

Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Systematic Biology. 63:685-697.

Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
}
\keyword{ phylogenetic comparative method }
\keyword{ ancestral state reconstruction }
\keyword{ function-valued traits }
\keyword{ infnite-dimensional traits }
\keyword{ evolution of tolerance }
\keyword{ reaction norms }
\keyword{ phenotypic plasticity }
\keyword{ dose-response curves }
\keyword{ phylogenetic signal }
\keyword{ PGLS }
\examples{
require(phytools)

# simulate evolution of a function-valued trait (glm with logit link)
sim_data <- sim.curves()

# perform ancestral curve reconstruction
anc_recon <- phylocurve(y~x,tree = sim_data$tree,data = sim_data$data)

# get tip coefficients and aligned function data
tip_coefficients <- get.tip.coefficients(formula = y~x,tree = sim_data$tree,data = sim_data$data)
data <- get.aligned.function.data(tip_coefficients)

# estimate evolutionary rates
evo.model.fitted <- evo.model(sim_data$tree,data)
}
