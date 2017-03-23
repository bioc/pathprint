\name{single.chip.enrichment}
\alias{single.chip.enrichment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate enrichment of a list of genesets in an array
}
\description{
Function to assess enrichment of gene sets in an array or matrix of arrays using various summary statistics
}
\usage{
single.chip.enrichment(exprs,
	geneset,
	transformation = "rank",
	statistic = "mean",
	normalizedScore = FALSE,
	progressBar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{
An expression matrix, rownames correspond to gene ids used in the list of genesets
}
  \item{geneset}{
list of pathways or genesets over which to assess statistic
}
  \item{transformation}{
Initial transformation applied to each column of exprs, can be one of "rank", "squared.rank" or "log.rank" 
}
  \item{statistic}{
Summary statistic to be applied, either "mean" or "median"
}
  \item{normalizedScore}{
Logical. If statistic = "mean" and normalizedScore = TRUE, option to calculate a parametric
significance score based on the expected distribution of scores.
 Other summary statistics currently not supported
}
  \item{progressBar}{
Logical. Shows progress of script, good to check running okay, set to FALSE for possible faster running}
}
\details{
This is the worker function for exprs2fingerprint, in conjuction with an exprs based on
 Entrez Gene IDs and the standard pathprint genesets e.g. \code{\link{pathprint.Hs.gs}}.
  The (un-normalized) results are passed onto thresholdFingerprint to produce the Pathway Fingerprint scores
}
\value{
Matrix containing pathway enrichment scores for each sample in the exprs input matrix.
 Rownames are genesets and colnames are the columns of the exprs matrix.
}
\references{
Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\author{
Gabriel Altschuler
}

\seealso{
\code{\link{exprs2fingerprint}}
}
\examples{
require(pathprintGEOData)

# Compare continuous pathway enrichment values to Pathway Fingerprint scores

# Use ALL dataset as an example

require(ALL)
data(ALL)
annotation(ALL)

#ds = c("chipframe", "genesets","pathprint.Hs.gs",
#  "platform.thresholds","pathprintSummarizedGEOData", 
#  "pluripotents.frame")
#data(list = ds)

## extract GEO.fingerprint.matrix and GEO.metadata.matrix
#GEO.fingerprint.matrix = assays(result)$fingerprint
#GEO.metadata.matrix = colData(result)

ds = c("chipframe","pathprint.Hs.gs","genesets","platform.thresholds")
data(list = ds)


# The chip used was the Affymetrix Human Genome U95 Version 2 Array
# The correspending GEO ID is GPL8300

# Analyze patients with ALL1/AF4 and BCR/ABL translocations
ALL.eset <- ALL[, ALL$mol.biol \%in\% c("BCR/ABL", "ALL1/AF4")]
ALL.exprs<-exprs(ALL.eset)

patient.type<-as.character(ALL$mol.biol[
	ALL$mol.biol \%in\% c("BCR/ABL", "ALL1/AF4")])

# Process fingerprints
ALL.fingerprint<-exprs2fingerprint(exprs = ALL.exprs,
	platform = "GPL8300",
	species = "human",
	progressBar = TRUE
	)

color.map <- function(mol.biol) {
	if (mol.biol=="ALL1/AF4") "#00FF00" else "#FF00FF"
	}
patientcolors <- sapply(ALL$mol.biol[
	ALL$mol.biol \%in\% c("BCR/ABL", "ALL1/AF4")],
	function(x){
	if (x == "ALL1/AF4") "#00FF00" else "#FF00FF"
	})


# define list of differentially activated pathways between the two groups
signif.pathways<-diffPathways(ALL.fingerprint,
	fac = patient.type,
	threshold = 0.6)

# draw heatmap
heatmap(ALL.fingerprint[signif.pathways,],
	ColSideColors = patientcolors,
	col = c("blue", "white", "red"),
	scale = "none", mar = c(10,20),
	cexRow = 0.75)
title(sub = "Pathways differentially activated in patients
 with ALL1/AF4 (green) and BCR/ABL(purple) translocations",
      cex.sub = 0.75)
      

######
# Compare to continous values
ALL.exprs.entrez <- customCDFAnn(ALL.exprs, chipframe$GPL8300$ann)
ALL.enrichment <- single.chip.enrichment(exprs = ALL.exprs.entrez,
	geneset = pathprint.Hs.gs,
	transformation = "squared.rank",
	statistic = "mean",
	normalizedScore = FALSE,
	progressBar = TRUE
	)

heatmap(ALL.enrichment[signif.pathways,],
	ColSideColors = patientcolors,
	col = colorRampPalette(c("blue", "white", "red"))(100),
	scale = "row", mar = c(10,20),
	cexRow = 0.75)
title(sub = "Continuous pathway enrichment scores for patients
 with ALL1/AF4 (green) and BCR/ABL(purple) translocations",
      cex.sub = 0.75)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
