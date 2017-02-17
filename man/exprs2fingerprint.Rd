\name{exprs2fingerprint}
\alias{exprs2fingerprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a pathway fingerprint from a gene expression table.}
\description{
The function converts the gene expression values to a ternary matrix of pathway expression values,
 (-1,0,1) corresponding to (low, background, high).
 This is based on applying a pre-calculated threshold to pathway enrichment scores.
}
\usage{
exprs2fingerprint(exprs, platform, species, progressBar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{
matrix containing a probe expression table, can be one or more columns
}
  \item{platform}{
microarray platform GEO ID
}
  \item{species}{
character string to define the species of the experiment, see details.
}
  \item{progressBar}{
logical. If TRUE, a progress bar is displayed while the script is running
}
}
\details{
exprs should be a matrix or dataframe of the expression values, with rownames containing probe names and colnames the experiment IDs.
 Platforms should be of the type listed in GEO (e.g. "GPL570"). Species can be full latin names
\cr "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster",
 "Caenorhabditis elegans".
\cr or corresponding common-use names
\cr "human",  "mouse", "rat", "zebrafish", "drosophila", "C.elegans".
\cr The array is first annotated with Entrez Gene IDs using annotations contained in \code{\link{chipframe}}.
 Pathway expression scores are calculated by the mean-squared rank of the gene expression and
 normalized against the appropriate distribution for the given platform in the GEO corpus.
 There is a progressBar to track the script, can be set to FALSE for (possibly) marginally faster running
}
\value{
Returns a dataframe containing the pathway fingerprint for each of column in the expression table.
 Rownames correspond to pathways and colnames to the experiment IDs.
}
\references{
Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\author{
Gabriel Altschuler
}
\seealso{
\code{\link{consensusFingerprint}}, \code{\link{single.chip.enrichment}}, \code{\link{customCDFAnn}},
 \code{\link{thresholdFingerprint}}
}
\examples{

# Use ALL dataset as an example
require(ALL)
data(ALL)
annotation(ALL)

ds = c("chipframe", "genesets","pathprint.Hs.gs",
  "platform.thresholds","pathprintSummarizedGEOData", 
  "pluripotents.frame")
data(list = ds)

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(result)$fingerprint
GEO.metadata.matrix = colData(result)

# The chip used was the Affymetrix Human Genome U95 Version 2 Array
# The correspending GEO ID is GPL8300

# Extract portion of the expression matrix and process fingerprints
ALL.First5.fingerprint<-exprs2fingerprint(exprs = ALL[1:5],
    platform = "GPL8300",
    species = "human",
    progressBar = TRUE
)
head(ALL.First5.fingerprint)								   

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
