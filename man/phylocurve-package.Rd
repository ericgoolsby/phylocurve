\name{phylocurve-package}
\alias{phylocurve-package}
\docType{package}
\title{
phylocurve: Phylogenetic comparative methods for function-valued traits
}
\description{
phylocurve offers a variety of phylogenetic comparative tools for studying the evolution of function-valued traits, including ancestral state reconstruction, estimating phylogenetic signal, and assessing correlated function-valued trait evolution.
}
\details{
\tabular{ll}{
Package: \tab phylocurve\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-26\cr
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
sim_data <- simcurves()

# ancestral curve reconstruction
anc_recon <- phylocurve(formula = y~x,tree = sim_data$tree,data = sim_data$data)

# estimate phylogenetic signal
phylocurve.signal(tip_coefficients = anc_recon$tip_coefficients,tree = sim_data$tree)

# assess correlated trait evolution with a univariate trait
X <- fastBM(sim_data$tree)
phylocurve.pgls(tip_coefficients = anc_recon$tip_coefficients,univariate_trait = X,tree = sim_data$tree)
}
