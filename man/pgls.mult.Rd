\name{pgls.mult}
\alias{pgls.mult}
\title{
Multivariate phylogenetic generalized least squares
}
\description{
Performs phylogenetic generalized least squares using the method described in Adams (2014) using a covariance-based approach (Goolsby 2015)
}
\usage{
pgls.mult(rate.mult.fitted, X)
}
\arguments{
  \item{rate.mult.fitted}{
A fitted rate.mult object
}
  \item{X}{
A matrix of predictor variables X for PGLS. Note, an intercept is estimated automatically (do not add a column of ones to X). Rows should be named according to species names.
}
}
\value{
An object of class "anova"
}
\references{
Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional multivariate data. Evolution. 68:2675-2688.

Golsby E.W. 2015. Covariance-based maximum likelihood estimation is equivalent to and more flexible than distance-based phylogenetic comparative methods for high-dimensional multivariate traits. In review.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link[geomorph]{procD.pgls}}
}
\examples{
# Simulate random data for a multivariate trait with 50 species and 10 trait dimensions
s <- sim.mult(nspecies = 50,R = diag(10))

# Simulate a univariate trait
require(phytools)
X <- as.matrix(fastBM(s$tree))

# Estimate evolutionary rates
rate.mult.fitted <- rate.mult(tree = s$tree,Y = s$Y_means)

# Perform multivariate phylogenetic regression
pgls.mult(rate.mult.fitted = rate.mult.fitted,X = X)
}
