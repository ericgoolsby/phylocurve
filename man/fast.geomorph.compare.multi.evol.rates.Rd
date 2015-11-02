\name{fast.geomorph.compare.multi.evol.rates}
\alias{fast.geomorph.compare.multi.evol.rates}
\title{
Fast covariance-based implementations of distance-based methods
}
\description{
The functions \code{fast.geomorph.compare.evol.rates}, \code{fast.geomorph.compare.multi.evol.rates}, \code{fast.geomorph.phylo.pls}, \code{fast.geomorph.procD.pgls}
, and \code{fast.geomorph.physignal} are covariance-baesd implementations of the geomorph functions \link[geomorph]{compare.evol.rates}, \link[geomorph]{compare.multi.evol.rates}, \link[geomorph]{phylo.pls}, \link[geomorph]{procD.pgls}, and \link[geomorph]{physignal} using a fast linear-time algorithm. Code is directly modified from the original geomorph code for direct comparison between distance-based and covariance-based methods.
}
\usage{
fast.geomorph.compare.multi.evol.rates(A, gp,
phy, Subset = TRUE, method = "ML", ShowPlot = TRUE, iter = 999)
}
\arguments{
  \item{A}{
From geomorph: A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
}
  \item{gp}{
From geomorph: A factor array designating group membership
}
  \item{phy}{
From geomorph: A phylogenetic tree of class phylo
}
  \item{Subset}{
From geomorph: A logical value indicating whether or not the traits are subsets from a single landmark configuration (default is TRUE)
}
  \item{method}{
Maximum likelihood "ML" or restricted maximum likelihood "REML"
}
  \item{ShowPlot}{
From geomorph: A logical value indicating whether or not the plot should be returned
}
  \item{iter}{
From geomorph: Number of iterations for significance testing
}
}
\details{
See \link[geomorph]{compare.multi.evol.rates}
}
\value{
See \link[geomorph]{compare.multi.evol.rates}
}
\references{
Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.

Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.

Denton, J.S.S., and D.C. Adams. 2015. A new phylogenetic test for comparing multiple high-dimensional evolutionary rates suggests interplay of evolutionary rates and modularity in lanternfishes (Myctophiformes; Myctophidae). Evolution. 69: doi:10.1111/evo.12743
}
\author{
Eric W. Goolsby
}
\seealso{
\link[geomorph]{compare.multi.evol.rates}
}
\examples{
### NOTE: this example is identical 
### to the example code for the
### analogous geomorph function
### for direct comparisons with
### 'fast.geomorph' phylocurve functions
require(geomorph)
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
land.gp<-c("A","A","A","A","A","B","B","B","B","B","B")  #mandible and cranium subsets

compare.multi.evol.rates(Y.gpa$coords,land.gp,plethspecies$phy,iter=99)
fast.geomorph.compare.multi.evol.rates(Y.gpa$coords,land.gp,plethspecies$phy)
}
