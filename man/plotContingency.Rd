\name{plotContingency}
\alias{plotContingency}
\title{generates plot of response vs out-of-bag predictions across classes
}
\description{
For categorical models, this function generates a graphical representation of the 
true vs predicted contingency matrix across classes for a given \code{alpha}.
}
\usage{
plotContingency(x, alpha.index=NULL, xlab="class (true)", ylab="class (predicted)", 
cex.lab=0.95, main=NULL, col.main="black", cex.main=0.85, cex.axis=1, 
symbol.size.inches=0.5, bg.color="steelblue2", fg.color=NULL, margin=0.2, 
frequency.label=TRUE, frequency.label.cex=1, frequency.label.offset=0, \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{alpha.index}{Integer indices to select alpha values. Default is \code{1:length(alpha)}}
  \item{xlab}{Custom x-axis label.}
  \item{ylab}{Custom y-axis label.}
  \item{cex.lab}{Axis label size.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{cex.axis}{Axis size.}
  \item{symbol.size.inches}{Symbol size.}
  \item{bg.color}{Symbol color.}
  \item{fg.color}{Color of symbol background.}
  \item{margin}{Margin size  to accomodate symbols.}
  \item{frequency.label}{Logical to display class frequency labels. Default is \code{TRUE}.}
  \item{frequency.label.cex}{Size of class frequency labels.}
  \item{frequency.label.offset}{Offset of class frequency labels.}
  \item{\dots}{Additional plotting parameters.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\dontshow{
set.seed(123)
fit = eNetXplorer(x=matrix(rnorm(120),ncol=4),y=c(rep(0,15),rep(1,15)),
family="binomial",n_run=5, n_perm_null=4,alpha=c(0.5,1))
plot(x=fit,plot.type="contingency")
plotContingency(x=fit,alpha.index=1)
}
\donttest{
data(QuickStartEx)
binarized=rep("low",length(QuickStartEx$response))
binarized[QuickStartEx$response>median(QuickStartEx$response)]="high"
fit = eNetXplorer(x=QuickStartEx$predictor,y=binarized,family="binomial",n_run=20, 
n_perm_null=10,seed=111)
plot(x=fit,plot.type="contingency")
plotContingency(x=fit,alpha.index=6)
}
}
