\name{K.mult}
\alias{K.mult}
\title{
Multivariate phylogenetic signal
}
\description{
Estimates multivariate phylogenetic signal (Adams 2014) using a covariance-based impelmentation (Goolsby 2015).
}
\usage{
K.mult(rate.mult.fitted, iter = 1000)
}
\arguments{
  \item{rate.mult.fitted}{
A fitted rate.mult object
}
  \item{iter}{
Number of phylogenetic permutations for significance testing
}
}
\value{
\item{phy.signal }{Multivariate phylogenetic signal (Kmult)}
\item{pvalue }{Significance of Kmult}
}

\references{
Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Systematic Biology. 63:685-697.

Golsby E.W. 2015. Covariance-based maximum likelihood estimation is equivalent to and more flexible than distance-based phylogenetic comparative methods for high-dimensional multivariate traits. In review.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link[geomorph]{physignal}}
}
\examples{
# Simulate random data for a multivariate trait with 50 species and 10 trait dimensions
s <- sim.mult(nspecies = 50,R = diag(10))

# Estimate evolutionary rates
rate.mult.fitted <- rate.mult(tree = s$tree,Y = s$Y_means)

# Estimate multivariate phylogenetic signal
K.mult(rate.mult.fitted = rate.mult.fitted)
}