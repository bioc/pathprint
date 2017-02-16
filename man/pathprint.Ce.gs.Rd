\name{pathprint.Ce.gs}
\alias{pathprint.Ce.gs}
\docType{data}
\title{
Pathprint genesets - C. elegans
}
\description{
Pathways and genesets used by pathprint for \emph{C. elegans} arrays, referenced by Entrez Gene ID
}
\usage{pathprint.Ce.gs}
\details{
Gene sets were inferred by homology from the human genesets, \code{\link{pathprint.Hs.gs}},
 using the HomoloGene database, \url{www.ncbi.nlm.nih.gov/homologene}
}
\source{
O. Hofmann
}
data("pathprint.Ce.gs")

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(result)$fingerprint
GEO.metadata.matrix = colData(result)
\references{
Sayers et al. Database resources of the National Center for Biotechnology Information.
 Nucleic Acids Research (2011) vol. 39 (Database issue) pp. D38-51
\cr Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\examples{
pathprint.Ce.gs[grep("ZN175_7728", names(pathprint.Ce.gs))]
}
\keyword{datasets}
