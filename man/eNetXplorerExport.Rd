\name{eNetXplorerExport}
\alias{eNetXplorerExport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates plain text tables from eNetXplorer object
}
\description{

}
\usage{
eNetXplorerExport(x, path=getwd(), delimiter=c("tab","csv"), input.data=T, output.data.summary=T, output.data=T)
}
\arguments{
  \item{x}{\code{eNetXplorer} object from function eNetXplorer.}
  \item{path}{Destination file path. Default is the current working directory.}
  \item{delimiter}{Delimiter for the generated files, either tab- or comma-separated values. Default is tab.}
  \item{input.data}{Logical variable indicating whether to generate files containing input data (i.e. data fed into the models and model arguments).} 
  \item{output.data.summary}{Logical variable indicating whether to generate a file with summary results from the models.} 
  \item{output.data}{Logical variable indicating whether to generate files with detailed results from the models.}
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

\seealso{\code{eNetXplorer}}

\examples{
\dontrun{
data(QuickStartExample)
result = eNetXplorer(x=QuickStartExample$predictor,y=QuickStartExample$response,family="gaussian",n_run=25,n_perm_null=20,seed=111)
eNetXplorerExport(result)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
