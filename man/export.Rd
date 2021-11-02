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
export(x, dest_dir=getwd(), dest_dir_create=TRUE, delim=c("tsv","csv"), 
input.data=TRUE, summary.data=TRUE, output.data=TRUE)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{dest_dir}{Destination directory. Default is the working directory.}
  \item{dest_dir_create}{Creates destination directory if it does not exist already. Default is \code{TRUE}.}
  \item{delim}{Delimiter for the generated files, either tab-separated (\code{"tsv"}) or comma-separated (\code{"csv"}) values. Default is \code{"tsv"}.}
  \item{input.data}{Logical variable indicating whether to generate files containing input data (i.e. data fed into the models and model arguments). Default is \code{TRUE}.}   
  \item{summary.data}{Logical variable indicating whether to generate a file with summary results from the models. Default is \code{TRUE}.} 
  \item{output.data}{Logical variable indicating whether to generate files with detailed results from the models. Default is \code{TRUE}.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{eNetXplorer}}}
\examples{
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response, 
family="gaussian",n_run=20,n_perm_null=10,seed=111)
export(x=fit,dest_dir=tempdir())
}
}

