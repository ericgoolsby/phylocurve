\name{phylocurve}
\alias{phylocurve}
\title{
Ancestral curve reconstruction for logistic regression (glm with logit link)
}
\description{
This function performs ancestral reconstruction of function-valued traits (assuming Brownian motion).
}
\usage{
phylocurve(formula, tree, data, ymin = 0.01, ymax = 0.99, ylength = 30,
tip.coefficients, species.identifier = "species", verbose = FALSE)
}
\arguments{
  \item{formula}{
Formula for function-valued trait (currently only supports models of the form Y~X)
}
  \item{tree}{
An phylogenetic tree of class "phylo"
}
  \item{data}{
A data frame with data for tip curve estimation, where each row contains a single data point. A column named "species" has the species names corresponding to each data point, a predictor (X) variable and the response (Y) variable, which must be scaled between 0 and 1.
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
  \item{tip.coefficients}{
A matrix of regression coefficients can be provided instead of raw data. Row names should be species names, the first column should be the glm logit intercept parameter, and the second column should be the glm logit slope paramter.
}
  \item{species.identifier}{
Default is "species". Can be changed if the column in data has a different species identifier name.
}
  \item{verbose}{
Whether to print progress during ancestral curve reconstruction.
}
}
\details{
This function uses a PGLS-based method described in Goolsby (2015) to perform ancestral curve reconstruction. This function uses a fast tree transversal method via the phylolm package (Ho and Ane, 2014)
}
\value{
\item{node_coefficients }{Estimated regression coefficients of internal nodes}
\item{fitted_x}{Reconstructed x-values for each internal node. Curve coordinates for internal node i are (fitted_x[,i],y_vals)}
\item{lower_CI_x}{Lower 95 percent confidence interval for fitted x-values for each internal node}
\item{upper_CI_x}{Upper 95 percent confidence interval for fitted x-values for each internal node}
\item{y_vals}{Vector of y-values used as inputs for inverse function evaluation. Curve coordinates for internal node i are (fitted_x[,i],y_vals)}
\item{tip_coefficients}{Estimated regression coefficients of tip curves.}
}
\references{
Goolsby, E.W. 2015. "Phylogenetic comparative methods for evaluating the evolutionary history of function-valued traits." Systematic Biology. In press.

Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
}
\author{
Eric W. Goolsby
}
\examples{
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
