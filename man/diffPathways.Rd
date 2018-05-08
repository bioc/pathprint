\name{diffPathways}
\alias{diffPathways}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect differentially activated pathways between fingerprints
}
\description{
A function to return pathways consistently differentially expressed between two 
groups of pathway fingerprints
}
\usage{
diffPathways(fingerprints, fac, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{fingerprints}{matrix of fingerprints, the number of columns should 
    correspond to the length of fac
}
    \item{fac}{
vector of characters or factors, in an order corresponding to the order of 
columns in the fingerprint matrix. Contains two levels, denoting the groups 
to be compared.
}
    \item{threshold}{
numeric, between 0 and 2 - the threshold at which to assign an average 
difference in pathway usage. 
}
}
\details{
The vector of factors must contain only two levels (or two unique values for a 
character vector).
}
\value{
Returns a list of the rownames (i.e. pathways for the pathway fingerprint)
corresponding to the rows for which the difference
in the means between the two groups is greater than the threshold value.
\cr For a ternary fingerprint (-1,0,1), setting the threshold between 0.5 and 1
ensures that rownames are selected that differ across the majority of the 
arrays in the two groups. 
with values closer to 1 representing higher stringency.
This can break down and allow false positives in the case where one group 
contains a significant
but minority number of +1 and the other -1s.
}
\references{
Altschuler, G. M., O. Hofmann, I. Kalatskaya, R. Payne, S. J. Ho Sui, U. Saxena,
A. V. Krivtsov, S. A. Armstrong, T. Cai, L. Stein and W. A. Hide (2013). 
"Pathprinting: An integrative approach to understand the functional basis of 
disease." Genome Med 5(7): 68.
}
\author{
Gabriel Altschuler
}
\seealso{
\code{\link{exprs2fingerprint}}, \code{\link{consensusDistance}}, 
\code{\link{consensusFingerprint}}
}
\examples{
require(pathprintGEOData)

# Use ALL dataset as an example

require(ALL)
data(ALL)
annotation(ALL)
library(SummarizedExperiment)

# load  the data
data(SummarizedExperimentGEO)

ds = c("chipframe", "genesets","pathprint.Hs.gs",
    "platform.thresholds","pluripotents.frame")
data(list = ds)

# extract part of the GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(geo_sum_data[,1:200000])$fingerprint
GEO.metadata.matrix = colData(geo_sum_data[,1:200000])

# free up space by removing the geo_sum_data object
remove(geo_sum_data)

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


# define differentially activated pathways between the two groups
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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
