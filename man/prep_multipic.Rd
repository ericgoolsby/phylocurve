\name{prep_multipic}
\alias{prep_multipic}
\title{
Prep multipic
}
\description{
Prepares data and phylogeny for fast repeated phylogenetically independent contrast calculations in the function multipic.
}
\usage{
prep_multipic(x, phy, scaled = TRUE, var.contrasts = FALSE,
rescaled.tree = FALSE, pic_recon = FALSE)
}
\arguments{
  \item{x}{
A matrix with rownames corresponding to species in the phylogeny (one observation per species and no missing data), and columns corresponding to individual traits.
}
  \item{phy}{
An object of class "phylo"
}
  \item{scaled}{
logical, indicates whether the contrasts should be scaled with their expected variances (default to TRUE).
}
  \item{var.contrasts}{
logical, indicates whether the expected variances of the contrasts should be returned (default to FALSE).
}
  \item{rescaled.tree}{
logical, if TRUE the rescaled tree is returned together with the main results.
}
  \item{pic_recon}{
logical, whether to return PIC ancestral reconstruction values (NOTE: this is not equivalent to the ML ancestral estimates except for the root value!)
}
}
\value{
A list consisting of quantities used by the function multipic
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