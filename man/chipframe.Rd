\name{chipframe}
\alias{chipframe}
\docType{data}
\title{
Probe to Entrez Gene ID mappings
}
\description{
Probe to Entrez Gene ID mappings for platforms covered by the pathway 
fingerprint
}
\usage{chipframe}
\format{
    List with an entry for each GEO ID covered by pathprint (35 elements), 
    each containing a list
    \describe{
        \item{\code{ann}}{dataframe containing array probe 
        to Entrez Gene ID mappings}
    \item{\code{title}}{character, array name}
    }
}
\details{
The platform IDs correspond to GEO references
\url{http://www.ncbi.nlm.nih.gov/geo/}
}
\source{
Annotations obtained from the AILUN server \url{http://ailun.stanford.edu/}
}
\references{
Chen R., Li L., and Butte AJ (2007) AILUN: Reannotating Gene
Expression Data Automatically, \emph{Nature Methods}, 4(11), 879
}
\seealso{
\code{\link{customCDFAnn}}
}

\examples{
data(chipframe)
names(chipframe)
chipframe$GPL570$title
head(chipframe$GPL570$ann)
}
\keyword{datasets}
