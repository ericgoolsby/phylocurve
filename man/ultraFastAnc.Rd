\name{ultraFastAnc}
\alias{ultraFastAnc}
\title{
Ultra-fast maximum likelihood ancestral state reconstruction
}
\description{
This function performs ancestral state reconstruction using a fast algorithm based on phylogenetically independent contrasts (Felsenstein 1985). Contrary to \link[phytools]{fastAnc}, which reroots the tree at each internal node and recalculates the root state for each tree, this algorithm only requires a single tree traversal (Goolsby, In Review). This function is several hundred to several thousand times faster than the rerooting method used in \link[phytools]{fastAnc}, which is the next-fastest maximum likelihood ancestral state reconstruction method currently described (Revell 2012). As in the \link[phytools]{fastAnc} function, the variance and 95\% intervals of estimates are optionally returned.

NOTE: trees with polytomies are supported but (currently) slow down the algorithm considerably.
}
\usage{
ultraFastAnc(phy, x, vars = FALSE, CI = FALSE)
}
\arguments{
  \item{phy}{
A tree of class 'phylo'
}
  \item{x}{
A named vector of trait values
}
  \item{vars}{
Whether to return the variances of the (restricted) maximum likelihood estimates
}
  \item{CI}{
Whether to return 95\% confidence intervals of the (restricted) maximum likelihood estimates
}
}
\value{
A named vector of maximum likelihood ancestral states (with names corresponding to node number). If vars or CI is set to TRUE, a list is returned with these values included.
}
\references{
Goolsby E.W. (In Review). "Ultra-fast ancestral state reconstruction of continuous characters: a rerooting-free maximum likelihood approach."

Felsenstein, J. (1985) Phylogenies and the comparative method. American Naturalist, 125, 1-15.

Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.
}
\author{
Eric W. Goolsby
}
\seealso{
\link[phytools]{fastAnc}, \link[ape]{ace}, \link[ape]{pic}
}
\examples{
require(ape)
tree <- rtree(1e4) # random tree with 10,000 taxa
x <- setNames(rnorm(1e4),tree$tip.label) # random trait data
recon <- ultraFastAnc(phy=tree,x=x,CI=TRUE)
}
