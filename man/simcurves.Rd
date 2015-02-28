\name{simcurves}
\alias{simcurves}
\title{
Simulate function-valued curve evolution
}
\description{
This function simulates evolution of a function-valued trait (glm with logit link) under Brownian motion.
}
\usage{
simcurves(nspecies = 30, x_length = 20, nreps = 20, startree = FALSE, lambda = 1, seed)
}
\arguments{
  \item{nspecies}{
Number of species
}
  \item{x_length}{
Number of points on the curve to simulate
}
  \item{startree}{
Simulate using a star tree
}
  \item{lambda}{
Pagel's lambda tree transformation parameter (lambda = 0 for star tree; lambda = 1 for original tree)
}
  \item{seed}{
Optional seed to set for reproducibility
}
}
\value{
\item{data }{Data frame with simulated curves formatted appropriately for phylocurve()}
\item{tree }{Phylogenetic tree of class "phylo"}
\item{true_coefs }{True tip and ancestral coefficients for glm logit link}
}
\references{
Pagel, M. 1999. "Inferring the historical patterns of biological evolution". Nature 401:877-884.
}
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
