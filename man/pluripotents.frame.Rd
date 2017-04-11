\name{pluripotents.frame}
\alias{pluripotents.frame}
\docType{data}
\title{
Manually curated list of pluripotent arrays
}
\description{
A manually compiled list of pluripotent arrays (induced pluripotent cells and embryonic stem cells)
 together with their GEO IDs and descriptions
}
\usage{pluripotents.frame}
\format{
  A data frame with 278 observations on the following 5 variables.
  \describe{
    \item{\code{GSM}}{GEO sample ID}
    \item{\code{GSE}}{GEO series ID}
    \item{\code{GPL}}{GEO platform ID}
    \item{\code{source}}{GEO description - Source}
    \item{\code{Characteristics}}{GEO description - Characteristic}
  }
}
\source{
\url{http://www.ncbi.nlm.nih.gov/geo/}
}
\references{
Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\seealso{
\code{\link{consensusDistance}}, \code{\link{consensusFingerprint}}
}
\examples{
require(pathprintGEOData)
library(SummarizedExperiment)

# load  the data
data(compressed_result)

ds = c("chipframe", "genesets","pathprint.Hs.gs",
  "platform.thresholds","pluripotents.frame")
data(list = ds)

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(result)$fingerprint
GEO.metadata.matrix = colData(result)

head(pluripotents.frame)

# Use pathway fingerprints to search for
# additional pluripotent arrays across GEO
# create consensus pluripotent fingerprint
pluripotent.consensus<-consensusFingerprint(
	GEO.fingerprint.matrix[,pluripotents.frame$GSM], threshold=0.9)

# calculate distance from the pluripotent consensus
geo.pluripotentDistance<-consensusDistance(
	pluripotent.consensus, GEO.fingerprint.matrix)

# plot histograms
par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
geo.pluripotentDistance.hist<-hist(geo.pluripotentDistance[,"distance"],
	nclass = 50, xlim = c(0,1), main = "Distance from pluripotent consensus")
par(mar = c(7, 4, 4, 2))
hist(geo.pluripotentDistance[pluripotents.frame$GSM, "distance"],
	breaks = geo.pluripotentDistance.hist$breaks, xlim = c(0,1), 
	main = "", xlab = "above: all GEO, below: pluripotent samples")

# annotate top 100 matches not in original seed with metadata
geo.pluripotentDistance.noSeed<-geo.pluripotentDistance[
	!(rownames(geo.pluripotentDistance) \%in\% pluripotents.frame$GSM),
	]

top.noSeed.meta<-GEO.metadata.matrix[
	match(head(rownames(geo.pluripotentDistance.noSeed), 1000),
	GEO.metadata.matrix$GSM),
	]
# head(top.noSeed.meta[,c("GSM", "GPL", "Source")],10)

}
\keyword{datasets}
