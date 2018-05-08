\name{consensusFingerprint}
\alias{consensusFingerprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a consensus fingerprint
}
\description{
Produces a pathway fingerprint that represents the consensus of a series of 
pathway fingerprints, according to a user-defined threshold
}
\usage{
consensusFingerprint(fingerprintframe, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{fingerprintframe}{
        matrix of fingerprints from which the consensus will be calculated
    }
    \item{threshold}{threshold value (between 0 and 1)
    }
}
\details{
For each pathway the mean fingerprint score, m, is calculated, and the consensus
defined as
\cr +1 if m > threshold
\cr -1 if m < threshold
\cr 0 otherwise
}
\value{
Vector of consensus pathway fingerprint scores with names corresponding to 
pathways
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
\code{\link{consensusDistance}}
}
\examples{
require(pathprintGEOData)
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

# search for pluripotent arrays
# load fingerprint matrix and pluripotent reference

# create consensus fingerprint
pluripotent.consensus<-consensusFingerprint(
    GEO.fingerprint.matrix[,pluripotents.frame$GSM], threshold=0.9)

# calculate distance from the pluripotent consensus
geo.pluripotentDistance<-consensusDistance(pluripotent.consensus,
    GEO.fingerprint.matrix)

# plot histograms
par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
geo.pluripotentDistance.hist<-hist(geo.pluripotentDistance[,"distance"],
    nclass = 50, xlim = c(0,1), main = "Distance from pluripotent consensus")
par(mar = c(7, 4, 4, 2))
hist(geo.pluripotentDistance[pluripotents.frame$GSM, "distance"],
    breaks = geo.pluripotentDistance.hist$breaks, xlim = c(0,1), 
    main = "", xlab = "above: all GEO, below: curated pluripotent samples")

# annotate top 100 matches not in original seed with metadata
geo.pluripotentDistance.noSeed<-geo.pluripotentDistance[
    !(rownames(geo.pluripotentDistance) \%in\% pluripotents.frame$GSM),]

top.noSeed.meta<-GEO.metadata.matrix[
    match(head(rownames(geo.pluripotentDistance.noSeed), 1000),
    GEO.metadata.matrix$GSM),]
    
# head(top.noSeed.meta[,c("GSM", "GPL", "Source")],10)
}

