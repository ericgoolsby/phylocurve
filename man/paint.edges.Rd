\name{paint.edges}
\alias{paint.edges}
\title{
Paint tree edges based on species regimes
}
\description{
This function estimates the proportion of edges to be assigned to regimes based on ancestral regime reconstruction (Yang et al. 1995; Pupko et al. 2000; Paradis et al. 2004; Revell 2012).
}
\usage{
paint.edges(tree, species.groups, average.nodes = TRUE, root.edge = TRUE)
}
\arguments{
  \item{tree}{
An object of class \code{phylo}.
}
  \item{species.groups}{
A named factor assigning species to regimes.
}
  \item{average.nodes}{
Whether or not to average reconstructed probabilities of regime states from parent and daughter nodes for tree edges.
}
  \item{root.edge}{
Whether or not to return a root edge.
}
}
\value{
A matrix of dimension nedge-by-nregimes. NOTE: edges are in postorder.
}
\references{
Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Pupko, T., Pe'er, I, Shamir, R., and Graur, D. (2000) A fast algorithm for joint reconstruction of ancestral amino acid sequences. Molecular Biology and Evolution, 17, 890–896.

Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.

Yang, Z., Kumar, S., Nei, M. (1995) A new method of inference of ancestral nucleotide and amino acid sequences. Genetics, 141, 1641-1650.
}
\author{
Eric W. Goolsby
}
\examples{
require(ape)
tree <- rtree(n=10)
groups <- setNames(factor(c(rep("groupA",5),rep("groupB",5))),tree$tip.label)
painted.edges <- paint.edges(tree = tree,species.groups = groups)
}
