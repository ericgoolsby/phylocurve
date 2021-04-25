\name{fast.geomorph.physignal}
\alias{fast.geomorph.physignal}
\title{
DEPRECATED: Fast covariance-based implementations of distance-based methods
}
\description{
All fast.geomorph functions are deprecated and will be removed from phylocurve soon. If you are trying to perform distanced-based phylogenetic comparative methods on high-dimensional data, please use and cite the geomorph package.

The functions \code{fast.geomorph.compare.evol.rates}, \code{fast.geomorph.compare.multi.evol.rates}, \code{fast.geomorph.phylo.integration}, \code{fast.geomorph.procD.pgls}
, and \code{fast.geomorph.physignal} are covariance-baesd implementations of the geomorph functions \link[geomorph]{compare.evol.rates}, \link[geomorph]{compare.multi.evol.rates}, \link[geomorph]{phylo.integration}, \link[geomorph]{procD.pgls}, and \link[geomorph]{physignal} using a fast linear-time algorithm. Code is directly modified from the original geomorph code for direct comparison between distance-based and covariance-based methods.
}
\usage{
fast.geomorph.physignal(phy, A, iter = 1000,
ShowPlot = TRUE, method = c("Kmult", "SSC"))
}
\arguments{
  \item{phy}{
From geomorph: A phylogenetic tree of class phylo
}
  \item{A}{
From geomorph: A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
}
  \item{iter}{
From geomorph: Number of iterations for significance testing
}
  \item{ShowPlot}{
From geomorph: A logical value indicating whether or not the plot should be returned
}
  \item{method}{
From geomorph: Method for estimating phylogenetic signal (Kmult or SSC)
}
}
\details{
See \link[geomorph]{physignal}
}
\value{
See \link[geomorph]{physignal}
}
\references{
Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.

Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution, 57:717-745.

Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.

Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Systematic Biology. 63:685-697.
}
\author{
Eric W. Goolsby
}
\seealso{
\link[geomorph]{physignal}
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

#Test for phylogenetic signal in shape
physignal(phy=plethspecies$phy,Y.gpa$coords,iter=99)
fast.geomorph.physignal(plethspecies$phy,Y.gpa$coords,method="Kmult")

#Test for phylogenetic signal in size
physignal(phy=plethspecies$phy,Y.gpa$Csize,iter=99)
fast.geomorph.physignal(plethspecies$phy,Y.gpa$Csize,method="Kmult")
}
