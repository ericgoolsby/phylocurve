\name{phylocurve.signal}
\alias{phylocurve.signal}
\title{
Phylogenetic signal of function-valued traits
}
\description{
This function wraps physignal from the geomorph package (Adams 2014) to calculate a multivariate version of Blomberg's K (Kmult).
}
\usage{
phylocurve.signal(tip_coefficients, tree, ymin = 0.01, ymax = 0.99, ylength = 30, iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tip_coefficients}{
A matrix of regression coefficients can be provided instead of raw data. Row names should be species names, the first column should be the glm logit intercept parameter, and the second column should be the glm logit slope paramter.
}
  \item{tree}{
An phylogenetic tree of class "phylo"
}
  \item{ymin}{
Because 0 and 1 are undefined for inverse logit functions, the minimum and maximum values are defaulted to .01 and .99. Can be adjusted as needed.
}
  \item{ymax}{
Because 0 and 1 are undefined for inverse logit functions, the minimum and maximum values are defaulted to .01 and .99. Can be adjusted as needed.
}
  \item{ylength}{
How many landmarks (points on the curve) to evaluate.
}
  \item{iter}{
Number of iterations for significance testing
}
}
\details{
A multivariate generalization of Blomberg's K. See \code{\link[geomorph]{physignal}} for additional details.
}
\value{
\item{phy.signal }{The estimate of phylogenetic signal}
\item{pvalue }{The significance level of the observed signal}
}
\references{
Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Systematic Biology. 63:685-697.

Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution, 57:717-745.

Goolsby, E.W. 2015. "Phylogenetic comparative methods for evaluating the evolutionary history of function-valued traits." Systematic Biology. In press.
}\examples{
require(phytools)

# simulate evolution of a function-valued trait (glm with logit link)
sim_data <- simcurves()

# ancestral curve reconstruction
anc_recon <- phylocurve(formula = y~x,tree = sim_data$tree,data = sim_data$data)

# estimate phylogenetic signal
phylocurve.signal(tip_coefficients = anc_recon$tip_coefficients,tree = sim_data$tree)
}