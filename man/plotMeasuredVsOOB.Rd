\name{plotMeasuredVsOOB}
\alias{plotMeasuredVsOOB}
\title{generates plot of response vs out-of-bag predictions across instances}
\description{
Given \code{alpha}, this function generates plots of response vs out-of-bag predictions across instances, which can be used to assess individual instances, examine outliers, etc. 
For linear regression models, it generates a response vs out-of-bag prediction scatterplot; 
it also displays the best linear fit and its 95\% confidence level region. 
For categorical models, it generates a boxplot across classes showing the frequency of out-of-bag correct predictions.
}
\usage{
plotMeasuredVsOOB(x, alpha.index=NULL, xlab=NULL, ylab=NULL, 
cex.lab=0.95, main=NULL, col.main="black", cex.main=0.85, instance.label=TRUE,
instance.label.cex=NULL, instance.label.offset=NULL, 
instance.label.added.margin=NULL, col=NULL, box.wex=NULL, box.range=NULL,
box.col=NULL, transparency=NULL, jitter=NULL, cex.pt=NULL, class.color=NULL, \dots)
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
  \item{instance.label}{Logical to display instance labels. Default is \code{TRUE}.}
  \item{instance.label.cex}{Size of instance labels.}
  \item{instance.label.offset}{Offset of instance labels.}
  \item{instance.label.added.margin}{(linear regression only) Margin size  to accomodate instance label display.}
  \item{col}{(linear regression only) Symbol color.}
  \item{box.wex}{(categorical models only) Boxplot \code{boxwex} parameter. Default is 0.5.}
  \item{box.range}{(categorical models only) Boxplot \code{range} parameter. Default is 0.}
  \item{box.col}{(categorical models only) Boxplot \code{col} parameter. Default is \code{white}.}
  \item{transparency}{(categorical models only) Symbol transparency. Default is 70.}
  \item{jitter}{(categorical models only) Symbol jitter. Default is 0.25.}
  \item{cex.pt}{(categorical models only) Symbol size. Default is 1.7}
  \item{class.color}{(categorical models only) Vector of class colors. Default is the default palette.}
  \item{\dots}{Additional plotting parameters.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\donttest{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",n_run=20,n_perm_null=10,seed=111)
plot(x=fit,plot.type="measuredVsOOB")
plotMeasuredVsOOB(x=fit,alpha.index=2)
}
\donttest{
data(QuickStartEx)
binarized=rep("low",length(QuickStartEx$response))
binarized[QuickStartEx$response>median(QuickStartEx$response)]="high"
fit = eNetXplorer(x=QuickStartEx$predictor,y=binarized,family="binomial",n_run=20, 
n_perm_null=10,seed=111)
plot(x=fit,plot.type="measuredVsOOB")
plotMeasuredVsOOB(x=fit,alpha.index=2)
}
}
