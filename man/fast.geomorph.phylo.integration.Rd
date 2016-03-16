\name{fast.geomorph.phylo.integration}
\alias{fast.geomorph.phylo.integration}
\title{
Fast covariance-based implementations of distance-based methods
}
\description{
The functions \code{fast.geomorph.compare.evol.rates}, \code{fast.geomorph.compare.multi.evol.rates}, \code{fast.geomorph.phylo.integration}, \code{fast.geomorph.procD.pgls}
, and \code{fast.geomorph.physignal} are covariance-baesd implementations of the geomorph functions \link[geomorph]{compare.evol.rates}, \link[geomorph]{compare.multi.evol.rates}, \link[geomorph]{phylo.integration}, \link[geomorph]{procD.pgls}, and \link[geomorph]{physignal} using a fast linear-time algorithm. Code is directly modified from the original geomorph code for direct comparison between distance-based and covariance-based methods.
}
\usage{
fast.geomorph.phylo.integration(A1, A2, phy, iter = 999,
label = NULL, verbose = FALSE, ShowPlot = TRUE)
}
\arguments{
  \item{A1}{
From geomorph: A 2D array (n x [p1 x k1]) or 3D array (p1 x k1 x n) containing landmark coordinates for the first block
}
  \item{A2}{
From geomorph: A 2D array (n x [p2 x k2]) or 3D array (p2 x k2 x n) containing landmark coordinates for the second block
}
  \item{phy}{
From geomorph: A phylogenetic tree of class phylo 
}
  \item{iter}{
From geomorph: Number of iterations for significance testing
}
  \item{label}{
From geomorph: An optional vector indicating labels for each specimen that are to be displayed
}
  \item{verbose}{
From geomorph: A logical value indicating whether the output is basic or verbose
}
  \item{ShowPlot}{
From geomorph: A logical value indicating whether or not the plot should be returned
}
}
\details{
See \code{phylo.integration}
}
\value{
See \code{phylo.integration}
}
\references{
Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.

Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological integration and trait covariation in morphometric data using evolutionary covariance matrices. PLOS ONE. 9(4):e94335.
}
\author{
Eric W. Goolsby
}
\seealso{
\link[geomorph]{phylo.integration}
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

# NOTE: phylo.integration is phylo.integration in newer versions of geomorph
# phylo.integration(Y.gpa$coords[1:5,,],Y.gpa$coords[6:11,,],plethspecies$phy,iter=99)
# phylo.integration(Y.gpa$coords[1:5,,],Y.gpa$coords[6:11,,],plethspecies$phy,iter=99)
fast.geomorph.phylo.integration(Y.gpa$coords[1:5,,],Y.gpa$coords[6:11,,],plethspecies$phy)
}
