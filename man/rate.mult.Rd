\name{rate.mult}
\alias{rate.mult}
\title{
Estimate multivariate evolutionary rates and models
}
\description{
This function estimates multivariate evolutionary rates (Adams 2014) and alternative evolutionary models using a covariance-based approach (Goolsby 2015). The function must be run prior to running other multivariate functions, including K.mult, pgls.mult, compare.rate.mult, and compare.multivar.rate.mult.
}
\usage{
rate.mult(tree = tree, Y = Y, type = c("mult", "diag", "all"),
method = c("REML", "ML"), error = c("none", "estimate", "supply"),
error_n = 20, error_supply, model = "BM", fixed_sigma2, fixed_model_pars)
}
\arguments{
  \item{tree}{
An object of class phylo
}
  \item{Y}{
A data frame in which the first column is labeled "species", and is filled with species names. The remaining columns contain values for the multivariate trait Y. Multiple observations per species (using multiple rows) and missing data are allowed.
}
  \item{type}{
The form of evolutionary rate estimation to use. If type="mult", a single evolutionary rate is estimated for all trait dimensions. If type="diag", a separate evolutionary rate is estimated for all trait dimensions. NOTE: the option type="all" is also available, in which an entire covariance matrix R is estimated for traits (not compatbile with other functions!), in which case, only BM evolution with no missing data and no within-species error is allowed.
}
  \item{method}{
Whether to use restricted maximum likelihood (REML, the default) or maximum likelihood.
}
  \item{error}{
Whether to incorporate within-species measurement error. Can be set to "none" (the default) if only a single observation is available per species. If raw data with multiple observations per species is supplied, select "estimate" to estimate pooled within-species variance for species with less than error_n observations (for species with >error_n observations, variance is used). Alternatively, a matrix of the square of standard errors (variance divided by the number of within-species observations) for each species and trait dimension can be supplied if error is set to "supply".
}
  \item{error_n}{
Only relevant if error="estimate". For trait dimensions with less than error_n within-species observations, pooled within-species variance is used. If the number of within-species observations for a trait dimension is greater than or equal to error_n, variance is estimated directly.
}
  \item{error_supply}{
Only relevent if error="supply". A matrix of species and trait dimension squared standard errors (or variance divided by the number of within-species observations). Row names should be species names.
}
  \item{model}{
Evolutionary model specification -- default is "BM". Other options include "OUrandomRoot", "OUfixedRoot", "EB", "lambda", "kappa", and "delta". More than one model may be specified (e.g., c("lambda","EB")), in which case the branch length transformations are applied in order of the specified models. NOTE: if incorporating alternative models, DO NOTE include "BM".
}
  \item{fixed_sigma2}{
If the evolutionary rate (single value if type="mult") or rates (if type="diag") is desired to be fixed, entire it here.
}
  \item{fixed_model_pars}{
If the evolutionary model parameter or parameters are desired to be fixed, entire them here in a vector named according to the appropraite model (should be in order of the models specified by the model argument).
}
}
\value{
\item{sigma2 }{Estimated evolutionary rate}
\item{pars }{Estimated evolutionary model parameters (NULL If model="BM")}
\item{logL }{Log-likelihood of the estimated parameters}
\item{method }{REML or ML}
\item{model }{Evoluionary model}
\item{rate.mult.args }{List of arguments to be called with do.call(rate.mult,rate.mult$rate.mult.args)}
}
\references{
Adams D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for shape and other high-dimensional phenotypic data. Systematic Biology. 63:166-177.

Golsby E.W. 2015. Covariance-based maximum likelihood estimation is equivalent to and more flexible than distance-based phylogenetic comparative methods for high-dimensional multivariate traits. In review.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link{compare.multivar.rate.mult}}, \code{\link{compare.rate.mult}}, \code{\link[geomorph]{compare.evol.rates}}
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