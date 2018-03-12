\name{consensusDistance}
\alias{consensusDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate a distribution of distances from a consensus fingerprint
}
\description{
Calculates the distance from a consensus for a series of pathway fingerprints,
accounting only for significantly high or low (-1 or 1) pathways in the 
consensus
}
\usage{
consensusDistance(consensus, fingerprintframe)
}
\arguments{
    \item{consensus}{
consensus fingerprint
}
    \item{fingerprintframe}{
dataframe of sample fingerprints from which the distance will be calculated
}
}
\details{
The consensus fingerprint can be calculated using 
\code{\link{consensusFingerprint}} or alternatively can be a single 
fingerprint vector
}
\value{
A dataframe with rows corresponding to each sample contained in the 
fingerprintframe with the following columns
\item{distance}{Manhattan distance of sample from the consensus fingerprint, 
scaled by the maximum possible distance}
\item{pvalue}{p-value representing the probabilty that the samples are not 
phenotypically matched.
N.B. this is only valid when the fingerprint frame represents a sufficiently 
broad coverage of phenotypes, e.g. the GEO corpus.
This p-value is based on an assumption that the distances are normally 
distributed}
}
\references{
Altschuler, G. M., O. Hofmann, I. Kalatskaya, R. Payne, S. J. Ho Sui, U. Saxena,
A. V. Krivtsov, S. A. Armstrong, T. Cai, L. Stein and W. A. Hide (2013). 
"Pathprinting: An integrative approach to understand the functional basis 
of disease." Genome Med 5(7): 68.
}
\author{Gabriel Altschuler}
\seealso{
\code{\link{consensusFingerprint}}
}
\examples{
require(pathprintGEOData)
library(SummarizedExperiment)

# load  the data
data(SummarizedExperimentGEO)

ds = c("chipframe", "genesets",
    "pathprint.Hs.gs","platform.thresholds", "pluripotents.frame")
data(list = ds)

# extract GEO.fingerprint.matrix and GEO.metadata.matrix
GEO.fingerprint.matrix = assays(geo_sum_data)$fingerprint
GEO.metadata.matrix = colData(geo_sum_data)

# search for pluripotent arrays
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
    main = "", xlab = "above: all GEO, below: curated pluripotent samples")
}
