\name{plotKaplanMeier}
\alias{plotKaplanMeier}
\title{generates Kaplan-Meier plot of survival probability as a function of time
}
\description{
For Cox regression models, this function generates a Kaplan-Meier plot of survival probability as a function of time for a given \code{alpha}. The default behavior is to partition the cohort in two groups by the predicted risk median, but custom partitions in two or more groups (specified by a vector of predicted risk percentiles) are also possible. In the former case, provided that the \code{eNetXplorer} object was generated with the \code{logrank=TRUE} argument, the corresponding cross-validated log-rank test p-value is displayed in the default title.
}
\usage{
plotKaplanMeier(x, alpha.index=NULL, xlab="Time", ylab="Probability of Survival", 
cex.lab=1, main=NULL, col.main="black", cex.main=0.95, conf.int=TRUE,
breaks_ptiles=NULL, risk.col=NULL, legend=TRUE, legend.cex=0.75, \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object (must be \code{family="cox"}).}
  \item{alpha.index}{Integer indices to select alpha values. Default is \code{1:length(alpha)}}
  \item{xlab}{Custom x-axis label.}
  \item{ylab}{Custom y-axis label.}
  \item{cex.lab}{Axis label size.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{conf.int}{Logical to display 95\% confidence intervals. Default is \code{TRUE}.}
  \item{breaks_ptiles}{Vector of percentiles (in 0-1 range) to partition the cohort based on predicted risk. 
  	Default is 0.5.}
  \item{risk.col}{Vector of colors to display the predicted risk-based subcohorts.}
  \item{legend}{Logical to display legend. Default is \code{TRUE}.}
  \item{legend.cex}{Legend size.}
  \item{\dots}{Additional plotting parameters.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{
\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\donttest{
data(breastCancerSurv)
fit = eNetXplorer(x=breastCancerSurv$predictor, y=breastCancerSurv$response, family="cox", 
n_run=25, n_perm_null=15, seed=111, logrank=TRUE)
plot(x=fit, plot.type="KaplanMeier")
plotKaplanMeier(x=fit, alpha.index=6, conf.int=FALSE, breaks_ptiles=c(0.333,0.667))
}
}
