\name{thresholdFingerprint}
\alias{thresholdFingerprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Apply threshold values to produce a Pathway Fingerprint
}
\description{
Function to produce ternary threshold values, Pathway Fingerprint scores, from continuous
 geneset enrichment values. Returns ternary scores for each pathway, high (1), mid (0), low (-1)
}
\usage{
thresholdFingerprint(SCE, platform)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SCE}{
Pathway enrichment matrix from \code{\link{single.chip.enrichment}}
}
  \item{platform}{
GEO platform ID for array used
}
}
\details{
The thresholds have been pre-calculated and optimized against a panel of tissue samples (see ref).
}
\value{
Matrix containing ternary scores for each sample in the SCE input matrix.
 Rownames are genesets and colnames are the columns of the SCE matrix.
}
\references{
Altschuler, G. M., O. Hofmann, I. Kalatskaya, R. Payne, S. J. Ho Sui, U. Saxena, A. V. Krivtsov, S. A. Armstrong, T. Cai, L. Stein and W. A. Hide (2013). "Pathprinting: An integrative approach to understand the functional basis of disease." Genome Med 5(7): 68.
}
\author{
Gabriel Altshuler
}

\seealso{
\code{\link{exprs2fingerprint}}, \code{\link{platform.thresholds}}
}
\examples{
require(pathprintGEOData)
library(SummarizedExperiment)

# load  the data
data(compressed_result)

# Comparing workflows

# 1. Pathway Fingerprint scores from exprs2fingerprint

# Use ALL dataset as an example

require(ALL)
data(ALL)
annotation(ALL)

ds = c("chipframe", "genesets","pathprint.Hs.gs",
  "platform.thresholds","pluripotents.frame")
data(list = ds)

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(result)$fingerprint
GEO.metadata.matrix = colData(result)

# The chip used was the Affymetrix Human Genome U95 Version 2 Array
# The correspending GEO ID is GPL8300

# Analyze patients with ALL1/AF4 and BCR/ABL translocations
ALL.eset <- ALL[,1:5]
ALL.exprs<-exprs(ALL.eset)
# Process fingerprints
ALL.fingerprint<-exprs2fingerprint(exprs = ALL.exprs,
	platform = "GPL8300",
	species = "human",
	progressBar = TRUE
	)
									   
# 2. Thresholded pathway enrichment values

# Annotate
ALL.exprs.entrez <- customCDFAnn(ALL.exprs, chipframe$GPL8300$ann)

# Pathway enrichment
ALL.enrichment <- single.chip.enrichment(exprs = ALL.exprs.entrez,
	geneset = pathprint.Hs.gs,
	transformation = "squared.rank",
	statistic = "mean",
	normalizedScore = FALSE,
	progressBar = TRUE
	)

# Threshold
ALL.enrichment.threshold <- thresholdFingerprint(
	ALL.enrichment, "GPL8300")

# Compare 1. and 2.
all.equal(ALL.enrichment.threshold, ALL.fingerprint)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
