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
\cr Altschuler, G. M., O. Hofmann, I. Kalatskaya, R. Payne, S. J. Ho Sui, U. Saxena, A. V. Krivtsov, S. A. Armstrong, T. Cai, L. Stein and W. A. Hide (2013). "Pathprinting: An integrative approach to understand the functional basis of disease." Genome Med 5(7): 68.
}
\examples{
data("pathprint.Dm.gs")
pathprint.Dm.gs[grep("ZN175_7728", names(pathprint.Dm.gs))]
}
\keyword{datasets}

