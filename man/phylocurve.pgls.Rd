\name{phylocurve.pgls}
\alias{phylocurve.pgls}
\title{
Phylogenetic ANOVA/regression for function-valued data
}
\description{
This function wraps procD.pgls from the geomorph package (Adams 2014) to assess the significance of correlated evolution between a function-valued and a univariate trait.
}
\usage{
phylocurve.pgls(tip_coefficients, univariate_trait, tree, ymin = 0.01, ymax = 0.99, ylength = 30, iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tip_coefficients}{
A matrix of regression coefficients can be provided instead of raw data. Row names should be species names, the first column should be the glm logit intercept parameter, and the second column should be the glm logit slope paramter.
}
  \item{univariate_trait}{
A univariate trait (with species names) to assess for correlated evolution with a function-valued trait.
}
  \item{tree}{
A phylogenetic tree of class "phylo"
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
See \code{\link[geomorph]{procD.pgls}} for additional details.
}
\value{
Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS, F ratio, Prand, and Rsquare.
}
\references{
Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional multivariate data. Evolution. 68:2675-2688.

Goolsby, E.W. 2015. "Phylogenetic comparative methods for evaluating the evolutionary history of function-valued traits." Systematic Biology. In press.
}
\examples{
# simulate evolution of a function-valued trait (glm with logit link)
sim_data <- simcurves()

# ancestral curve reconstruction
anc_recon <- phylocurve(formula = y~x,tree = sim_data$tree,data = sim_data$data)

# assess correlated trait evolution with a univariate trait
X <- fastBM(sim_data$tree)
multivar_pgls(tip_coefficients = anc_recon$tip_coefficients,univariate_trait = X,tree = sim_data$tree)
}
