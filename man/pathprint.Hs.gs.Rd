\name{pathprint.Hs.gs}
\alias{pathprint.Hs.gs}
\docType{data}
\title{
Pathprint genesets - H. sapiens
}
\description{
Pathways and genesets used by pathprint for \emph{H.sapiens} arrays, referenced by Entrez Gene ID
}
\usage{pathprint.Hs.gs}

\source{
O. Hofmann
}

\references{
Altschuler, G. M., O. Hofmann, I. Kalatskaya, R. Payne, S. J. Ho Sui, U. Saxena, A. V. Krivtsov, S. A. Armstrong, T. Cai, L. Stein and W. A. Hide (2013). "Pathprinting: An integrative approach to understand the functional basis of disease." Genome Med 5(7): 68.
}
\examples{
data("pathprint.Hs.gs")
pathprint.Hs.gs[grep("ZN175_7728", names(pathprint.Hs.gs))]
}
\keyword{datasets}

