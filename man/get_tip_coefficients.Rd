\name{get_tip_coefficients}
\alias{get_tip_coefficients}
\title{
Estimate regression coefficients for tip species curves
}
\description{
Estiamtes regression coefficients for tip species curves from raw data. May be useful if wanting to perform methdos without first performing ancestral curve reconstruction via the phylocurve() function.
}
\usage{
get_tip_coefficients(formula, tree, data, ymin = 0.01, ymax = 0.99, ylength = 30, species_identifier = "species", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula for function-valued trait (currently only supports models of the form Y~X)
}
  \item{tree}{
A phylogenetic tree of class "phylo"
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
  \item{species_identifier}{
Default is "species". Can be changed if the column in data has a different species identifier name.
}
  \item{verbose}{
ether to print progress during tip curve coefficient estimation.
}
}
\value{
Estimated regression coefficients of tip curves.
}
\examples{
# simulate evolution of a function-valued trait (glm with logit link)
sim_data <- simcurves()

tip_coefficients <- get_tip_coefficients(formula = y~x,tree = sim_data$tree,data = sim_data$data)

# ancestral curve reconstruction
anc_recon <- phylocurve(tree = sim_data$tree,tip_coefficients = tip_coefficients)

# estimate phylogenetic signal
multivar_phylosig(tip_coefficients = tip_coefficients,tree = sim_data$tree)

# assess correlated trait evolution with a univariate trait
X <- fastBM(sim_data$tree)
multivar_pgls(tip_coefficients = tip_coefficients,univariate_trait = X,tree = sim_data$tree)
}
