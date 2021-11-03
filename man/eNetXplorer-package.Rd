\name{eNetXplorer-package}
\alias{eNetXplorer-package}
\docType{package}
\title{explores elastic net families for generalized linear models}
\description{
Provides a quantitative toolkit to explore elastic net families and to uncover correlates contributing to prediction under a cross-validation framework. Fits linear, binomial (logistic), multinomial, and Cox regression models.
}
\details{
\tabular{ll}{
Package: \tab eNetXplorer\cr
Type: \tab Package\cr
Version: \tab 1.1.3\cr
Date: \tab 2021-11-03\cr
License: \tab GPL-3\cr
}
This package provides a full pipeline of analysis: 
\code{\link{eNetXplorer}} takes in \code{x}, \code{y} data to generate a family of elastic net models over a range of \code{alpha} values; 
\code{\link{summary}} generates a summary of results in tabular format; \code{\link{plot}} provides a variety of plots to visualize results; \code{\link{summaryPDF}} generates a 
report in PDF format; \code{\link{export}} creates plain text output files for downstream
processing; and \code{\link{mergeObj}} merges \code{eNetXplorer} objects with different \code{alpha} values.
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\references{Candia J and Tsang JS.
  \emph{eNetXplorer: an R package for the quantitative exploration of elastic net families 
  for generalized linear models}, BMC Bioinformatics (2019) 20:189.
}
\keyword{ package }
\seealso{
\code{\link{eNetXplorer}}, \code{\link{summary}}, \code{\link{plot}}, \code{\link{summaryPDF}}, \code{\link{export}}, \code{\link{mergeObj}}.
}
\examples{
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",n_run=20,n_perm_null=10,seed=111)
summary(fit)
plot(x=fit,plot.type="measuredVsOOB",alpha.index=4)
summaryPDF(x=fit,dest_dir=tempdir())
export(x=fit,dest_dir=tempdir())
}
}
