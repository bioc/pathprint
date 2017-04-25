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
Altschuler, G. M., O. Hofmann, I. Kalatskaya, R. Payne, S. J. Ho Sui, U. Saxena, A. V. Krivtsov, S. A. Armstrong, T. Cai, L. Stein and W. A. Hide (2013). "Pathprinting: An integrative approach to understand the functional basis of disease." Genome Med 5(7): 68.
}
\examples{
require(pathprintGEOData)
library(SummarizedExperiment)

# load  the data
data(compressed_result)
data("platform.thresholds")

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(result)$fingerprint
GEO.metadata.matrix = colData(result)
head(platform.thresholds[[1]])
}

\keyword{datasets}
