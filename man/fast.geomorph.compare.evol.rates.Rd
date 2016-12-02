\name{fast.geomorph.compare.evol.rates}
\alias{fast.geomorph.compare.evol.rates}
\title{
Fast covariance-based implementations of distance-based methods
}
\description{
The functions \code{fast.geomorph.compare.evol.rates}, \code{fast.geomorph.compare.multi.evol.rates}, \code{fast.geomorph.phylo.pls}, \code{fast.geomorph.procD.pgls}
, and \code{fast.geomorph.physignal} are covariance-baesd implementations of the geomorph functions \link[geomorph]{compare.evol.rates}, \link[geomorph]{compare.multi.evol.rates}, \link[geomorph]{phylo.pls}, \link[geomorph]{procD.pgls}, and \link[geomorph]{physignal} using a fast linear-time algorithm. Code is directly modified from the original geomorph code for direct comparison between distance-based and covariance-based methods.
}
\usage{
fast.geomorph.compare.evol.rates(phy, A, gp,
method = "ML", ShowPlot = TRUE, iter = 1000,
censored = FALSE, force.diag = FALSE)
}
\arguments{
  \item{phy}{
An object of class \code{phylo}.
}
  \item{A}{
From geomorph: A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
}
  \item{gp}{
From geomorph: A factor array designating group membership
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
  \item{censored}{
Whether or not to use the censored approach for estimating regime-specific evolutionary rates (O'Meara et al. 2006). The censored approach is more efficient and has identical statistical power to the default method based on phylogenetic transformation of residuals.
}
  \item{force.diag}{
Whether or not to assume trait independence (a diagonal evolutionary rate matrix), as was originally described in Adams (2014c) and implemented in geomorph up to version 2.1.5. Should normally be set to FALSE (the default) unless trait independence is deliberately assumed.
}
}
\details{
See \link[geomorph]{compare.evol.rates}
}
\value{
See \link[geomorph]{compare.evol.rates}
}
\references{
Goolsby E.W. 2016. Likelihood-Based Parameter Estimation for High-Dimensional Phylogenetic Comparative Models: Overcoming the Limitations of 'Distance-Based' Methods. In review.

Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.
}
\author{
Eric W. Goolsby
}
\seealso{
\link[geomorph]{compare.evol.rates}
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

 gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
 names(gp.end)<-plethspecies$phy$tip

#Calculate rates of shape
compare.evol.rates(phy=plethspecies$phy,Y.gpa$coords,gp=gp.end,iter=49)
fast.geomorph.compare.evol.rates(plethspecies$phy,Y.gpa$coords,gp=gp.end)

#Calculate rates of size
compare.evol.rates(phy=plethspecies$phy,Y.gpa$Csize,gp=gp.end,iter=49)
fast.geomorph.compare.evol.rates(plethspecies$phy,Y.gpa$Csize,gp=gp.end,iter=49)
}
