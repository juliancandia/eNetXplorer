\name{summaryPDF}
\alias{summaryPDF}
\title{
generates PDF report with summary of main results
}
\description{
This function generates a PDF report that contains a plot of model performance across the \code{alpha} range, 
followed by plots showing detailed results for each value of \code{alpha}.
}
\usage{
summaryPDF(x, dest_dir=getwd(), dest_dir_create=TRUE, dest_file="eNetSummary.pdf")
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{dest_dir}{Destination directory. Default is the working directory.}
\item{dest_dir_create}{Creates destination directory if it does not exist already. Default is \code{TRUE}.}
  \item{dest_file}{Name for output PDF file.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response, 
family="gaussian",n_run=20,n_perm_null=10,seed=111)
summaryPDF(x=fit,dest_dir=tempdir())
}
}
