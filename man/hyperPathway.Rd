\name{hyperPathway}
\alias{hyperPathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	 Produces a list of pathways enrichments calculated using the hypergeometric distribution

	}
\description{
A function that returns a list enrichment statistics for a gene list in a set of pathways or modules
}
\usage{
hyperPathway(genelist, geneset, Nchip)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{
%%     ~~Describe \code{genelist} here~~
}
  \item{geneset}{
list of pathways or genesets over which to assess statistic
}
  \item{Nchip}{
The background number of genes on which to base the hypergeometric distribution
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Pathway activation mapping for large scale comparison of gene expression 
Altschuler G, Hofmann O, Payne R, Kalatskaya I, Stein L, Cai T, Hide W \emph{in preparation}
}
\author{
Gabriel Altschuler
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
