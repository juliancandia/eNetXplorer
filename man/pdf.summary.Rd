\name{pdf.summary}
\alias{pdf.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates pdf report with summary of main results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pdf.summary(x, path = getwd(), filename = "eNetXplorer.summary.pdf", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{eNetXplorer object}
  \item{path}{Directory path for output pdf file. Default is current working directory}
  \item{filename}{Name for output pdf file.}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
%% ~put references to the literature/web site here ~
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{eNetXplorer}
}
\examples{
\dontrun{
data(QuickStartExample)
result = eNetXplorer(x=QuickStartExample$predictor,y=QuickStartExample$response,family="gaussian",n_run=25,n_perm_null=20,seed=111)
suppressWarnings(pdf.summary(result))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
