\name{compare.rate.mult}
\alias{compare.rate.mult}
\title{
Compare evolutionary rates
}
\description{
This function tests for differences in evolutionary rates (and optionally, evolutionary models) using the methods described in Adams (2014) and Goolsby (2015) in a covariance-based approach.
}
\usage{
compare.rate.mult(rate.mult.fitted, groups, fit_individual = FALSE)
}
\arguments{
  \item{rate.mult.fitted}{
A fitted rate.mult object in which all data is fit under the null hypothesis (e.g., all trait dimensions fit under the same evolutionary rate and/or evolutionary model).
}
  \item{groups}{
A named factor by species specifying which group to assign individual species. For example, given a tree with 16 species, the command setNames(object = factor(c(rep("a",8),rep("b",8)),levels=c("a","b")),nm = tree$tip.label) assigns 8 species to group "a" and 8 species to group "b".
}
  \item{fit_individual}{
Whether or not to fit individual evolutionary models to each trait dimension (if TRUE), or if a single evolutionary model is fit to all trait dimensions (if FALSE, the default). Ignored if model="BM".
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
\item{null_model_list }{List of fitted rate.mult objects for the null model corresponding to each group}
\item{alt_model_list }{List of fitted rate.mult objects for the alternative model corresponding to each group}
}
\references{
Adams D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for shape and other high-dimensional phenotypic data. Systematic Biology. 63:166-177.

Golsby E.W. 2015. Covariance-based maximum likelihood estimation is equivalent to and more flexible than distance-based phylogenetic comparative methods for high-dimensional multivariate traits. In review.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link{rate.mult}}, \code{\link{compare.multivar.rate.mult}}, \code{\link[geomorph]{compare.evol.rates}}
}
\examples{
# Simulate random data for a multivariate trait with 50 species and 10 trait dimensions
s <- sim.mult(nspecies = 50,R = diag(10))

# Estimate the null model (a single rate parameter for the entire tree)
null <- rate.mult(tree = s$tree,Y = s$Y_means)

# Set 2 groups of species to test for different evolutionary rates
groups <- setNames(factor(c(rep("a",25),rep("b",25)),levels=c("a","b")),s$tree$tip.label)

# Test the hypothesis of different evolutionary rates for each group
compare.rate.mult(rate.mult.fitted = null,groups = groups)
}