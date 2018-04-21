\name{plotLambdaVsQF}
\alias{plotLambdaVsQF}
\title{
generates plot of quality function across \code{lambda}
}
\description{
Given \code{alpha}, this function generates a plot of the quality function across \code{lambda}, which is useful to examine the selection of the best penalty parameter.
}
\usage{
plotLambdaVsQF(x, alpha.index=NULL, xlab="lambda", ylab="QF (response vs out-of-bag predicted)", cex.lab=0.95, main=NULL, col.main="black", cex.main=0.95, log="x", type="b", \dots)
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
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\dontrun{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",n_run=20,n_perm_null=10,seed=111)
plot(x=fit,plot.type="lambdaVsQF")
plotLambdaVsQF(x=fit,alpha.index=c(1,3),main="custom title",col.main="red")
}
}
