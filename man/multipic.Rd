\name{multipic}
\alias{multipic}
\title{
Efficient PIC algorithm for multiple traits
}
\description{
Following the function prep_multipic, this function efficiently performs the PIC algorithm for multiple traits simultaneously, allowing for rapid repeated calls using different datasets (e.g., for bootstrap simulations). NOTE: this function is intended for use via the do.call function (see example).
}
\usage{
multipic(ntip, nnode, edge1, edge2, edge_len, phe, contr,
var_contr, scaled, pic_len, pic_recon)
}
\arguments{
  \item{ntip}{
This quantity is computed by the function prep_multipic
}
  \item{nnode}{
This quantity is computed by the function prep_multipic
}
  \item{edge1}{
This quantity is computed by the function prep_multipic
}
  \item{edge2}{
This quantity is computed by the function prep_multipic
}
  \item{edge_len}{
This quantity is computed by the function prep_multipic
}
  \item{phe}{
This quantity is computed by the function prep_multipic
}
  \item{contr}{
This quantity is computed by the function prep_multipic
}
  \item{var_contr}{
This quantity is computed by the function prep_multipic
}
  \item{scaled}{
This quantity is computed by the function prep_multipic
}
  \item{pic_len}{
This quantity is computed by the function prep_multipic
}
  \item{pic_recon}{
This quantity is computed by the function prep_multipic
}
}
\value{
\item{contrasts }{A matrix of PICs, one column for each trait (rows correspond to node numbers on a postorder tree)}
\item{sum_invV }{Sum of the inverse of the phylogenetic covariance matrix: sum(solve(vcv(tree)))}
\item{log_detV }{The natural log of the determinant of the phylogenetic covariance matrix}
\item{root }{Vector of estimated root values for each trait}
\item{pic_len }{Vector of PIC-rescaled branch lengths}
\item{pic_recon }{Matrix of PIC-reconstructed ancestral values. NOTE: this is not equivalent to the ML ancestral estimates except for the root value!}
}
\references{
Felsenstein, J. (1985) Phylogenies and the comparative method. American Naturalist, 125, 1-15.

Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289-290.

Paradis, E. (2012) Analysis of Phylogenetics and Evolution with R (Second Edition). New York: Springer.
}
\seealso{
\link[ape]{pic}
}
\examples{
require(phylocurve)
require(ape)

nspecies <- 1000
ntraits <- 500
tree <- rtree(n = nspecies) # Simulate a random 1000-species phylogeny
Y <- matrix(rnorm(nspecies*ntraits),ncol=ntraits) # Generate random data
rownames(Y) <- tree$tip.label

# Call prep_multipic
prep.Y <- prep_multipic(x = Y,phy = tree)

# Calculate phylogenetically independent contrasts using pic (ape package)
Y.pics.ape <- apply(X = Y,MARGIN = 2,FUN = pic,phy = tree)

# Calculate PICs using multipic
# returns a list with contrasts, sum_invV, log_detV, root
# 'contrast' is a matrix of PICs
# 'sum_invV' is equal to sum(solve(vcv(tree)))
# 'log_detV' is equal to log(det(vcv(tree)))
# 'root' is the maximum likelihood phenotypic value at the root
Y.pics.phylocurve <- do.call(multipic,prep.Y)

# Verify that results are identical
range(Y.pics.ape - Y.pics.phylocurve$contrasts)

# Generate 50 random datasets (NOT RUN)
#niter <- 50
#randomY <- vector(mode = "list",length = niter)
#for(i in 1:niter)
#{
#  randomY[[i]] <- matrix(rnorm(nspecies*ntraits),ncol=ntraits)
#  rownames(randomY[[i]]) <- tree$tip.label
#}


###########################################################################
# Compare time to calculate PICs on 50 random datasets using pic vs multipic
#
##### pic function (NOT RUN)
#system.time(for(i in 1:niter) apply(X = randomY[[i]],MARGIN = 2,FUN = pic,phy = tree))
##### user  system elapsed 
##### 18.35    0.23   18.61 

##### multipic function (NOT RUN)
#system.time(for(i in 1:niter)
#{
#  prep.Y$phe[1:nspecies,] <- randomY[[i]] # update prep.Y$phe with new data
#  do.call(multipic,prep.Y)
#})
##### user  system elapsed 
##### 1.38    0.14    1.52 
#
###########################################################################
}