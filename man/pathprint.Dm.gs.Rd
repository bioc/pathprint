\name{pathprint.Dm.gs}
\alias{pathprint.Dm.gs}
\docType{data}
\title{
Pathprint genesets - D. melanogaster
}
\description{
Pathways and genesets used by pathprint for \emph{D. melanogaster} arrays, referenced by Entrez Gene ID
}
\usage{pathprint.Dm.gs}
\details{
Gene sets were inferred by homology from the human genesets, \code{\link{pathprint.Hs.gs}},
 using the HomoloGene database, \url{www.ncbi.nlm.nih.gov/homologene}
}
\source{
O. Hofmann
}

\references{
Sayers et al. Database resources of the National Center for Biotechnology Information.
 Nucleic Acids Research (2011) vol. 39 (Database issue) pp. D38-51
\cr Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\examples{
data("pathprint.Dm.gs")
pathprint.Dm.gs[grep("ZN175_7728", names(pathprint.Dm.gs))]
}
\keyword{datasets}

