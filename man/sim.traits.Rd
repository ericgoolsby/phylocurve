\name{sim.traits}
\alias{sim.traits}
\title{
Simulate multivariate trait evolution
}
\description{
This function wraps the \code{\link[geiger]{sim.char}} in geiger function to simulate evolutionary models under a given evolutionary rate matrix, with the option of simulating multiple within-species observations and missing data.
}
\usage{
sim.traits(ntaxa = 15, ntraits = 4, nreps = 1,
nmissing = 0, tree, v, anc, intraspecific,
model = "BM", parameters, nsim = 1, return.type = "matrix")
}
\arguments{
  \item{ntaxa}{
Number of species to randomly generate on a pure-birth phylogeny (used if \code{tree} is not provided).
}
  \item{ntraits}{
Number of traits to randomly simulate according to a random evolutionary rate matrix (used if \code{v} is not provided).
}
  \item{nreps}{
Number of within-species observations to simulate (default is 1).
}
  \item{nmissing}{
Number of observations to randomly remove from the simulated data (default is 0).
}
  \item{tree}{
Optional: an object of class \code{phylo} (if none provided, a random pure-birth phylogeny of size ntraxa is generated).
}
  \item{v}{
The evolutionary covariance matrix. If none provided, a random covariance matrix of size ntraits-by-ntraits is generated.
}
  \item{anc}{
An optional vector of ancestral values.
}
  \item{intraspecific}{
An optional vector of within-species variance.
}
  \item{model}{
The evolutionary model to simulate trait evolution (BM, OUrandomRoot, OUfixedRoot, EB, lambda, kappa, or delta).
}
  \item{parameters}{
The evolutionary model parameter in a named list. Example: \code{list(lambda = 0.5)}
}
  \item{nsim}{
Number of datasets to simulate.
}
  \item{return.type}{
"data.frame" or "matrix"
}
}
\value{
\item{trait_data }{Simulated trait data}
\item{tree }{The phylogenetic tree}
\item{sim_tree }{The phylogenetic tree with specified branch length transformations}
}
\references{
Harmon Luke J, Jason T Weir, Chad D Brock, Richard E Glor, and Wendell
  Challenger. 2008. GEIGER: investigating evolutionary radiations. Bioinformatics
  24:129-131.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link[geiger]{sim.char}}
}
