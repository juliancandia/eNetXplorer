\name{eNetXplorer-internal}
\title{Internal eNetXplorer functions}
\alias{eNetXplorer.default}
\alias{eNetXplorerBinomial}
\alias{eNetXplorerGaussian}
\alias{eNetXplorerMultinomial}
\alias{plot.feature.caterpillar.multinom.eNetXplorer}
\alias{plot.feature.heatmap.multinom.eNetXplorer}
\alias{plot.measured.vs.oob.categ.eNetXplorer}
\alias{plot.measured.vs.oob.numer.eNetXplorer}
\alias{print.eNetXplorer}
\alias{print.summary.eNetXplorer}
\alias{summary.eNetXplorer}
\alias{plot.contingency.eNetXplorer}
\alias{plot.feature.caterpillar.eNetXplorer}
\alias{plot.feature.heatmap.eNetXplorer}
\alias{plot.lambda.vs.QF.eNetXplorer}
\alias{plot.measured.vs.oob.eNetXplorer}
\alias{plot.summary.eNetXplorer}
\description{Internal eNetXplorer functions}
\usage{
fix.lam(lam)
lambda.interp(lambda,s)
nonzeroCoef(beta, bystep=FALSE)
glmnet_softmax(x)
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\details{These are not intended for use by users. \code{lambda.interp}
  does linear interpolation of the lambdas to obtain a prediction at a
  new point s. \code{glmnet_softmax} does the classification for multinomial models.
  \code{nonzeroCoef} determines in an efficient manner which variables
  are nonzero in each fit. \code{jerr} prints out error messages from
  the fortran. \code{plotCoef} is called by the \code{plot} method for
  \code{glmnet} objects. \code{check_dots} is used in \code{coef} and
  \code{predict} with argument \code{exact=TRUE}, to make sure user
  supplies original data used to fit the \code{"glmnet"} object.
}
\keyword{internal}
