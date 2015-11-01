\name{sim.curves}
\alias{sim.curves}
\title{
Simulate function-valued curve evolution
}
\description{
This function simulates evolution of a function-valued trait (glm with logit link) under Brownian motion.
}
\usage{
sim.curves(nspecies = 30, x_length = 20, startree = FALSE, lambda = 1, seed)
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
# simulate evolution of a function-valued trait (glm with logit link)
sim_data <- sim.curves()

# perform ancestral curve reconstruction
anc_recon <- phylocurve(y~x,tree = sim_data$tree,data = sim_data$data)

# get tip coefficients and aligned function data
tip.coefficients <- get.tip.coefficients(formula = y~x,tree = sim_data$tree,data = sim_data$data)
data <- get.aligned.function_data(tip.coefficients)

# estimate evolutionary rates
rate.mult.fitted <- rate.mult(sim_data$tree,data)
}
