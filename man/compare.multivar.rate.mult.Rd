\name{compare.multivar.rate.mult}
\alias{compare.multivar.rate.mult}
\title{
Compare multiple multivariate evolutionary rates
}
\description{
Compares hypotheses about evolutionary rates for multivariate traits combining methods from Adams (2013) and Adams (2014), using a covariance-based approach as described in Goolsby (2015).
}
\usage{
compare.multivar.rate.mult(null_model, alt_model_list)
}
\arguments{
  \item{null_model}{
A fitted rate.mult object in which all data is fit under the null hypothesis (e.g., all traits fit under the same evolutionary rate).
}
  \item{alt_model_list}{
A fitted rate.mult object in which all data is fit under the alternative hypothesis (e.g., traits are fit under different evolutionary rates).
}
}
\value{
\item{null.logL }{Log-likelihood of the null model}
\item{null.pars }{Number of parameters used to fit the null model}
\item{alt.logL }{Log-likelihood of the alternative model}
\item{alt.pars }{Number of parameters used to fit the alternative model}
\item{df }{Degrees of freedom for likelihood ratio test}
\item{chi_sq }{Chi-square value for likelihood ratio test}
\item{p }{p-value for likelihood ratio test}
}
\references{
Adams D.C. 2013. Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. Systematic Biology. 62:181-192.

Adams D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for shape and other high-dimensional phenotypic data. Systematic Biology. 63:166-177.

Golsby E.W. 2015. Covariance-based maximum likelihood estimation is equivalent to and more flexible than distance-based phylogenetic comparative methods for high-dimensional multivariate traits. In review.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link{rate.mult}}, \code{\link{compare.rate.mult}}, \code{\link[geomorph]{compare.evol.rates}}
}
\examples{
# Simulate random data for 2 multivariate traits
# (5 dimensions each) with different evolutionary rates
s <- sim.mult(nspecies = 50,R = diag(c(rep(1,5),c(rep(2,5)))))

# Estimate the null model (a single rate parameter)
null <- rate.mult(tree = s$tree,Y = s$Y_means)

# Estimate components of the alternative hypothesis and store in a list
alt1 <- rate.mult(s$tree,Y = s$Y_means[,c(1,2:6)])
alt2 <- rate.mult(s$tree,Y = s$Y_means[,c(1,7:11)])
alt_list <- list(alt1,alt2)

# Compare the null to the alternative to hypothesis
comparison <- compare.multivar.rate.mult(null_model = null,alt_model_list = alt_list)
}
