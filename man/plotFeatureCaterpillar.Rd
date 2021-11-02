\name{plotFeatureCaterpillar}
\alias{plotFeatureCaterpillar}
\title{generates caterpillar plot of feature statistics}
\description{
Given \code{alpha}, this function generates a caterpillar plot of feature statistics compared to permutation null models, which includes statistical significance annotations for individual features. 
By default, features are selected (and ordered top-down) by statistical significance; options are provided to customize feature selection and display.}
\usage{
plotFeatureCaterpillar(x, alpha.index=NULL, stat=c("freq","coef"), 
feature.all=FALSE, feature.pval.thres=NULL, feature.set=NULL, feature.top.n=25, 
signif.code=TRUE, xlab=NULL, ylab=NULL, main=NULL, col.main="black", 
cex.main=0.85, line=1.5, subtitle=NULL, col.subtitle="darkgray", 
line.subtitle=0.5, cex.subtitle=0.55, cexRow=NULL, cex.lab=0.95, legend=TRUE, \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{alpha.index}{Integer indices to select \code{alpha} values. Default is \code{1:length(alpha)}}
  \item{stat}{
  Feature statistic: \code{"freq"} for nonzero frequency, \code{"coef"} for mean nonzero coefficient.}
  \item{feature.all}{(Feature selection argument 1) Logical to show all features. Default is \code{FALSE}.}
  \item{feature.pval.thres}{(Feature selection argument 2) P-value threshold to select features. Default is \code{NULL}.}
  \item{feature.set}{(Feature selection argument 3) Character vector of feature names to display. Default is \code{NULL}.}
  \item{feature.top.n}{(Feature selection argument 4) Number of top features (ordered by p-value) to display. Default is 25.}
  \item{signif.code}{Logical to display significance annotations. Default is \code{TRUE}.}
  \item{xlab}{Label for x axis.}
  \item{ylab}{Label for y axis.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{line}{Title location.}
  \item{subtitle}{Custom subtitle.}
  \item{col.subtitle}{Subtitle color.}
  \item{line.subtitle}{Subtitle location.}
  \item{cex.subtitle}{Subtitle size.}
  \item{cexRow}{Size of row labels.}
  \item{cex.lab}{Axis label size.}
  \item{legend}{Logical to display legend. Default is \code{TRUE}.}
  \item{\dots}{Additional plotting parameters.}
}
\details{Feature selection criteria are hierarchical based on arguments 1 through 4 (see argument description above), with argument 1 at the top of the hierarchy. E.g. if \code{feature.all} is explicitly set to \code{TRUE}, it will take precedence over any other feature selection argument. By default, the top 25 features are displayed, ordered top-down by significance based on the given value of \code{alpha}.}
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
plotFeatureCaterpillar(x=fit,alpha.index=2,stat="coef")
}
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor, y=QuickStartEx$response,
family="gaussian", n_run=20, n_perm_null=10, seed=111)
plot(x=fit,plot.type="featureCaterpillar",stat="coef")
plotFeatureCaterpillar(x=fit,alpha.index=3,stat="coef",main="custom title")
}
}
