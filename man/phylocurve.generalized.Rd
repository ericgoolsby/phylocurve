\name{phylocurve.generalized}
\alias{phylocurve.generalized}
\title{
Align curves of any shape
}
\description{
Aligns curves by x,y coordinates so that species curves can be analyzed as a single high-dimensional trait.
}
\usage{
phylocurve.generalized(tree, X, Y)
}
\arguments{
  \item{tree}{
An object of class phylo
}
  \item{X}{
A numeric vector of length M, used for each species
}
  \item{Y}{
An N x M marix where each row corresponds to a species and each column corresponds to index M_i of X
}
}
\value{
\item{aligned_data }{x and y coordinates arranged in a single row vector for each species}
\item{aligned_coordinates }{x and y coordinates for each species (one data point per row)}
\item{aligned_X}{Aligned X coordinates in N x M matrix}
\item{aligned_Y}{Aligned Y coordinates in N x M matrix}
\item{nr}{Number of aligned landmarks}
\item{anc_X}{Vector of root X values (assuming Brownian motion)}
\item{anc_Y}{Vector of root Y values (assuming Brownian motion)}
\item{tree}{Phylogenetic tree supplied into the function}
}
\references{
Goolsby, E.W. 2015. "Phylogenetic comparative methods for evaluating the evolutionary history of function-valued traits." Systematic Biology. In press.
}
\author{
Eric W. Goolsby
}
\examples{
# NOT RUN (takes several minutes)
# x_length <- 50
# nspecies <- 25
#sim_data <- sim.curves(nspecies = nspecies,x_length = x_length)
#
# Perform ancestral curve reconstruction without any assumptions about function shape
# P <- phylocurve.generalized(tree = sim_data$tree,X = sim_data$data$x[1:x_length],
# Y = t(matrix(sim_data$data$y,nrow = x_length,
# dimnames = list(1:x_length,sim_data$tree$tip.label))))
#
# Data might need trimming for values where there is no variation
# plot(P$anc_X,P$anc_Y)
#
# P.trimmed <- phylocurve.trim(phylocurve.generalized = P,min_X = 4,max_X = 8)
}
