\name{export}
\alias{export}
\title{
generates plain text files from eNetXplorer object
}
\description{
This function enables the extraction of three different levels of data (input, summary, and detailed output results) 
from an \code{eNetXplorer} object. Plain text data files are produced with tab- or comma-separated-value formats.
}
\usage{
export(x, path=getwd(), delim=c("tab","csv"), input.data=T, 
summary.data=T, output.data=T)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{path}{Destination file path. Default is the current working directory.}
  \item{delim}{Delimiter for the generated files, either tab- or comma-separated values. Default is tab.}
  \item{input.data}{Logical variable indicating whether to generate files containing input data (i.e. data fed into the models and model arguments). Default is TRUE.}   \item{summary.data}{Logical variable indicating whether to generate a file with summary results from the models. Default is TRUE.} 
  \item{output.data}{Logical variable indicating whether to generate files with detailed results from the models. Default is TRUE.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{eNetXplorer}}}
\examples{
\dontrun{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor, y=QuickStartEx$response, 
family="gaussian", n_run=20, n_perm_null=10, seed=111)
export(fit)
}
}

