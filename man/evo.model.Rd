\name{evo.model}
\alias{evo.model}
\title{
Fit an evolutionary model
}
\description{
This function fits the evolutionary model as specified (e.g., trait groups, BM/EB/OU model, fixed effects, etc.).
}
\usage{
evo.model(tree, Y, fixed.effects = NA, species.groups, trait.groups,
model = "BM", diag.phylocov = FALSE, method = "Pairwise REML",
force.zero.phylocov = character(), species.id = "species",
max.combn = 10000, painted.edges, ret.level = 2, plot.LL.surface = TRUE,
par.init.iters = 50, fixed.par = numeric(), multirate = FALSE, subset = TRUE)
}
\arguments{
  \item{tree}{
An object of class 'phylo'
}
  \item{Y}{
An nspecies-by-ntraits matrix with row names corresponding to species names. Alternatively, a data frame with a column for species names (with name species.id, which by default is 'species').
}
  \item{fixed.effects}{
Optional. An matrix of fixed effects (predictor variables for phylogenetic regression), of the same format as \code{Y}. 
}
  \item{species.groups}{
Optional. A named factor assigning species to regimes. Names correspond to species names.
}
  \item{trait.groups}{
Optional. A named factor assigning species to trait groups,for each a generalized evolutionary rate is fit (see Adams and Collyer 2015, Evolution). NOTE: if the argument is called, \code{multirate} is automatically set to \code{TRUE}.
}
  \item{model}{
Default = 'BM': The model for branch length transformations (default is 'BM'). Other options include 'OU', 'EB', 'lambda', 'kappa', and 'delta'.
}
  \item{diag.phylocov}{
Default = FALSE. Whether to assume trait independence (a diagonal evolutionary rate matrix). Not recommended unless specifically testing for trait independence or if trait independence is deliberately assumed a priori.
}
  \item{method}{
Default = 'Pairwise REML'. The method for calculating log-likelihood. Other options inlcude 'Pairwise ML', 'Full ML', and 'Full REML'. In general, it is recommended to always use 'Pairwise REML', but when comparing fixed effects, use an 'ML' method should be used. NOTE: 'Full ML' and 'Full REML' are numerically unstable as the number of traits approaches the number of species, and is undefined if ntraits >= nspecies.
}
  \item{force.zero.phylocov}{
Optional. The names of traits to be assigned to group A, in which group A is assumed to be uncorrelated with the remaining traits in group B. This procedure is used to test for evolutionary covariance between two groups of multivariate traits, similar to phylogeneic partial least squares (Adams and Felice 2014).
}
  \item{species.id}{
Default = 'species'. The name of the column in 'Y' corresponding to species names. Only applicable if 'Y' is a \code{data.frame}.
}
  \item{max.combn}{
Default = 10000. The maximum number of pairwise combinations allowed for pairwise log-likelihood calculations. If ncol(combn(ncol(Y))) > max.combn, a Monte Carlo approach with max.combn random samples is used to approximate the composite log-likelihood.
}
  \item{painted.edges}{
Optional. An nedge-by-nspecies.groups matrix corresponding to the proportion of each tree edge to assign to an evolutioanry regime. This is fit automatically if not supplied, or can be fit using the \link{paint.edges} function. Edges MUST be in postorder.
}
  \item{ret.level}{
Default = 2. The level of information to return. ret.level = 1 only returns the log-likelihood; ret.level = 2 returns most information but not the evolutionary rate matrix. ret.level = 3 returns most information including the evolutionary matrix.
}
  \item{plot.LL.surface}{
Default = TRUE. Whether to plot the log-likelihood surface when estimating tree transformation parameters (i.e. OU, EB, lambda, delta, or kappa).
}
  \item{par.init.iters}{
Default = 50. The number of values to evaluate along the feasible range of tree transformations parameters (not applicable if model='BM').
}
  \item{fixed.par}{
Optional. If a tree transformation is used (i.e. OU, EB, lambda, delta, or kappa), this argument can be used to fix the tree transformation parameter to a fixed value.
}
  \item{multirate}{
Default = FALSE. Whether or not to use a single value for evolutionary rates instead of an unconstrained evolutionary rate matrix. For compatibility with \code{\link[geomorph]{compare.evol.rates}} and \code{\link[geomorph]{compare.multi.evol.rates}}.
}
  \item{subset}{
Default = TRUE. Whether or not traits are subsets of a single landmark configuration (implemented for compatability with the \code{\link[geomorph]{compare.multi.evol.rates}} function.
}
}
\details{
Coming soon.
}
\value{
An object of class \code{evo.model}.
}
\references{
Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.

Adams D.C. 2014a. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Syst. Biol. 63:685-697.

Adams D.C. 2014b. A method for assessing phylogenetic least squares models for shape and other high-dimensional multivariate data. Evolution 68:2675-2688.

Adams D.C. 2014c. Quantifying and comparing phylogenetic evolutionary rates for shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.

Adams D.C., Felice R. 2014. Assessing phylogenetic morphological integration and trait covariation in morphometric data using evolutionary covariance matrices. PLOS ONE 9(4):e94335.

Adams D.C., Collyer M. 2015. Permutation tests for phylogenetic comparative analyses of high-dimensional shape data: What you shuffle matters. Evolution 69:823-829.
}
\author{
Eric W. Goolsby
}
\seealso{
\code{\link[phylocurve]{compare.models}}
\code{\link[phylocurve]{K.mult}}
\code{\link[geomorph]{compare.evol.rates}}
\code{\link[geomorph]{compare.multi.evol.rates}}
\code{\link[geomorph]{physignal}}
\code{\link[geomorph]{procD.pgls}}
\code{\link[geomorph]{phylo.pls}}
}
\examples{
rand.data <- sim.traits()
X <- rowMeans(rand.data$trait_data)
null.model <- evo.model(tree = rand.data$tree,
  Y = rand.data$trait_data,method = "Pairwise ML")
alt.model <- evo.model(tree = rand.data$tree,
  Y = rand.data$trait_data,fixed.effects = X,method = "Pairwise ML")

compare.models(model1 = null.model,model2 = alt.model,
  nsim = 100,parallel = FALSE)

K.mult(model = alt.model,nsim = 100)
}
