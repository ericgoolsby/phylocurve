\name{sim.mult}
\alias{sim.mult}
\title{
Simulate multivariate trait evolution
}
\description{
Simulates multivariate trait evolution
}
\usage{
sim.mult(nspecies, R, error, nreps = 1, nmissing = 0,
model, parameters, anc, tree, seed, nsims = 1)
}
\arguments{
  \item{nspecies}{
Number of species to simulate
}
  \item{R}{
Covariance matrix for traits to simulate
}
  \item{error}{
Vector of variances to simulate within-species variation for each trait
}
  \item{nreps}{
Number of replicates to simulate for each trait per species
}
  \item{nmissing}{
Number of data points to randomly remove
}
  \item{model}{
Model of evolution (if not "BM", can be "OUrandomRoot", "OUfixedRoot", "EB", "lambda", "kappa", or "delta")
}
  \item{parameters}{
List of evolutionary model parameters. For either OU model, use alpha. For EB, use rate. For lambda, kappa, or delta, use lambda, kappa, or delta.
}
  \item{anc}{
Vector of ancestral values for each trait (default is 0).
}
  \item{tree}{
Can provide phylogenetic tree for simulations here; otherwise it is randomly generated based on nspecies
}
  \item{seed}{
Whether or not to set a seed for reproducible results
}
  \item{nsims}{
Number of simulations (default is 1)
}
}
\value{
\item{Y_means }{Means of simulated data}
\item{Y_raw }{Raw simulated data (only different from means if nreps>1)}
\item{Y_array }{Means of simulated data in array format}
\item{Y_original }{Original simulated species means (only different if nerps>1)}
\item{tree }{An object of class phylo}
\item{tree_sim }{The tree on which data was simulated (only different if model!="BM")}
}
\author{
Eric W. Goolsby
}
\examples{
# Simulate random data for a multivariate trait with 50 species and 10 trait dimensions
s <- sim.mult(nspecies = 50,R = diag(10))

# Estimate evolutionary rates
rate.mult.fitted <- rate.mult(tree = s$tree,Y = s$Y_means)

# Estimate multivariate phylogenetic signal
K.mult(rate.mult.fitted = rate.mult.fitted)
}
