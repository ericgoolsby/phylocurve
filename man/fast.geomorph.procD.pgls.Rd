\name{fast.geomorph.procD.pgls}
\alias{fast.geomorph.procD.pgls}
\title{
Fast covariance-based implementations of distance-based methods
}
\description{
The functions \code{fast.geomorph.compare.evol.rates}, \code{fast.geomorph.compare.multi.evol.rates}, \code{fast.geomorph.phylo.pls}, \code{fast.geomorph.procD.pgls}
, and \code{fast.geomorph.physignal} are covariance-baesd implementations of the geomorph functions \link[geomorph]{compare.evol.rates}, \link[geomorph]{compare.multi.evol.rates}, \link[geomorph]{phylo.pls}, \link[geomorph]{procD.pgls}, and \link[geomorph]{physignal} using a fast linear-time algorithm. Code is directly modified from the original geomorph code for direct comparison between distance-based and covariance-based methods.
}
\usage{
fast.geomorph.procD.pgls(f1, phy, iter = 999,
int.first = FALSE, verbose = FALSE)
}
\arguments{
  \item{f1}{
From geomorph: A formula for the linear model (e.g., y~x1+x2)
}
  \item{phy}{
From geomorph: A phylogenetic tree of class phylo
}
  \item{iter}{
From geomorph: Number of iterations for significance testing
}
  \item{int.first}{
From geomorph: A logical value to indicate if interactions of first main effects should precede subsequent main effects
}
  \item{verbose}{
From geomorph: A logical value specifying whether additional output should be displayed
}
}
\details{
See \link[geomorph]{procD.pgls}
}
\value{
See \link[geomorph]{procD.pgls}
}
\references{
Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.

Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional multivariate data. Evolution. 68:2675-2688.
}
\author{
Eric W. Goolsby
}
\seealso{
\link[geomorph]{procD.pgls}
}
\examples{
### NOTE: this example is identical 
### to the example code for the
### analogous geomorph function
### for direct comparisons with
### 'fast.geomorph' phylocurve functions
require(geomorph)

### Example of D-PGLS for high-dimensional data
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
procD.pgls(Y.gpa$coords ~ Y.gpa$Csize,plethspecies$phy,iter=999)
fast.geomorph.procD.pgls(Y.gpa$coords ~ Y.gpa$Csize,plethspecies$phy)
}