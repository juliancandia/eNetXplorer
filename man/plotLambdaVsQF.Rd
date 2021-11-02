\name{plotLambdaVsQF}
\alias{plotLambdaVsQF}
\title{
generates plot of quality function across \code{lambda}
}
\description{
Given \code{alpha}, this function generates a plot of the quality (objective) function across \code{lambda}, which is useful to examine how was the "best \code{lambda}" value selected. 
}
\usage{
plotLambdaVsQF(x, alpha.index=NULL, xlab="lambda", 
ylab="QF (response vs out-of-bag predicted)", cex.lab=0.95, main=NULL, 
col.main="black", cex.main=0.95, log="x", type="b", \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{alpha.index}{Integer indices to select \code{alpha} values. Default is \code{1:length(alpha)}}
  \item{xlab}{Custom x-axis label.}
  \item{ylab}{Custom y-axis label.}
  \item{cex.lab}{Axis label size.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{log}{Log scale axis.}
  \item{type}{Plot type.}
  \item{\dots}{Additional plotting parameters.}
}
\details{By definition, the "best \code{lambda}" value for a given \code{alpha} is the one that maximizes the quality function (QF) over the range of \code{lambda} values considered. Therefore, QF vs \code{lambda} distributions with sharp, narrow, well-defined peaks provide more confidence in the selection of the optimal \code{lambda} value than those with less-defined peaks. Sometimes, and particularly for the ridge (\code{alpha}=0) solutions, QF is observed to increase or decrease monotonically with \code{lambda} over its entire range, causing a boundary \code{lambda} value to be selected;  
we conservatively recommend to disregard \code{alpha}-models generated under such circumstances. If interested in investigating further, we suggest to re-run those \code{alpha}-models by extending the default range of \code{lambda} values (via the argument \code{nlambda.ext}) or its density (via the argument \code{nlambda}). On occasion, the range of \code{lambda} values is effectively limited by convergence issues of the underlying \code{glmnet} model; in such scenario, we recommend to increment the maximum allowed number of iterations (via the argument \code{mxit}, which is passed on to \code{glmnet.control}) or to limit the complexity of the model (e.g. by filtering and reducing the number of features fed into \code{eNetXplorer}).
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\dontshow{
set.seed(123)
fit = eNetXplorer(x=matrix(rnorm(75),ncol=3),y=rnorm(25),family="gaussian",
n_run=5, n_perm_null=4,alpha=c(0.5,1))
plotLambdaVsQF(x=fit,alpha.index=1)
}
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",n_run=20,n_perm_null=10,seed=111)
plot(x=fit,plot.type="lambdaVsQF")
plotLambdaVsQF(x=fit,alpha.index=c(1,3),main="custom title",col.main="red")
}
}
