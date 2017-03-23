\name{platform.thresholds}
\alias{platform.thresholds}
\docType{data}
\title{ Pathway fingerprint threshold values
}
\description{
Ternary threshold values for conversion of continuous geneset enrichment scores to
 discrete Pathway Fingerprint scores - high (1), mid (0), low (-1) for each geneset and platform
 covered by the Pathway Fingerprint.
}
\usage{platform.thresholds}

\references{
Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\examples{
require(pathprintGEOData)

data("platform.thresholds")

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(result)$fingerprint
GEO.metadata.matrix = colData(result)
head(platform.thresholds[[1]])
}

\keyword{datasets}
