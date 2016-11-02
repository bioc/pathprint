\name{pathprint}
\alias{Introduction to pathprint}
\alias{pathprint}
\docType{package}
\title{Pathway fingerprinting for analysis of gene expression arrays}
\description{
Algorithms to convert a gene expression array provided as an expression table to a 'pathway fingerprint'.
 The pathway fingerprint provides an unbiased, consistent annotation of expression data as a molecular phenotype,
 represented by activation status in 633 pathways. This is a vector of
 discrete ternary scores to represent high (1),
 low(-1) or insignificant (0) expression in a suite of pathways.
 Systematic definition of these functional relationships provides a tool for searching
 a pathway activation map of gene expression spanning species and technologies.
}
\details{
\tabular{ll}{
Package: \tab pathprint\cr
Type: \tab Package\cr
Version: \tab 1.2.3\cr
Date: \tab 2012-06-14\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Gabriel Altschuler, Uma Saxena.
\cr Maintainer: Gabriel Altschuler \email{galtschu@hsph.harvard.edu}
}
\references{
Pathway activation mapping for large scale comparison of gene expression, 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\keyword{ package }
\seealso{
\code{\link{exprs2fingerprint}}, \code{\link{consensusFingerprint}}, \code{\link{consensusDistance}}
}
\examples{

# Use fingerprints to analyze the ALL dataset
require(ALL)
data(ALL)
annotation(ALL)

# The chip used was the Affymetrix Human Genome U95 Version 2 Array
# The correspending GEO ID is GPL8300

# Extract portion of the expression matrix
ALL.exprs<-exprs(ALL)
ALL.exprs.sub<-ALL.exprs[,1:5]

# Process fingerprints
ALL.fingerprint<-exprs2fingerprint(exprs = ALL.exprs.sub,
platform = "GPL8300",
species = "human",
progressBar = TRUE
)

head(ALL.fingerprint)


####
# Construct consensus fingerprint based on pluripotent records
# Use this consensus to find similar arrays

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
