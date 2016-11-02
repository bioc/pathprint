\name{GEO.fingerprint.matrix}
\alias{GEO.fingerprint.matrix}
\docType{data}
\title{ Matrix of GEO fingerprints
}
\description{
A matrix containing Pathway Fingerprint vectors have been pre-calculated for ~160,000 publicly
 available arrays from the GEO corpus, spanning 6 species and 31 platforms 
}
\usage{GEO.fingerprint.matrix}

\source{
Primary data was retrieved from \url{http://www.ncbi.nlm.nih.gov/geo/}
}
\references{
Barrett et al. NCBI GEO: mining tens of millions of expression profiles--database and tools update.
 \emph{Nucleic acids research} (2007) vol. 35 (Database issue) pp. D760-5
\cr Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\examples{
# create consensus fingerprint for pluripotent samples
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
head(top.noSeed.meta[,c("GSM", "GPL", "Source")],10)
}
\keyword{datasets}
