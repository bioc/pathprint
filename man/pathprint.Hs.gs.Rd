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
Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\examples{
pathprint.Hs.gs[grep("ZN175_7728", names(pathprint.Hs.gs))]
}
\keyword{datasets}

