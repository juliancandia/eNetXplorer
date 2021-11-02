\name{plotSurvROC}
\alias{plotSurvROC}
\title{generates time-dependent ROC plots from Cox predicted risks
}
\description{
For Cox regression models, this function generates time-dependent ROC plot(s) (true positive rate vs false positive rate) for a given \code{alpha} at the timepoint(s) provided based on median predicted risk.
Provided that the \code{eNetXplorer} object was generated with \code{survAUC=T}, the cross-validated
median AUC and 95\% CI are shown in the default title. 
For more details, see Heagerty et al and package \code{survivalROC}.}
\usage{
plotSurvROC(x, alpha.index=NULL, survAUC_time,
xlab="False positive rate (1 - Specificity)",
ylab="True positive rate (Sensitivity)", cex.lab=1, main=NULL, col.main="black", 
cex.main=0.95, status0="censored", status1="events", \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object (must be \code{family="cox"}).}
  \item{alpha.index}{Integer indices to select alpha values. Default is \code{1:length(alpha)}}
  \item{survAUC_time}{Timepoint(s) of interest. Must be in the same time units as the survival time
provided to build the \code{eNetXplorer} object.}
  \item{xlab}{Custom x-axis label.}
  \item{ylab}{Custom y-axis label.}
  \item{cex.lab}{Axis label size.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{status0}{Title label for censoring ("status"=0).}
  \item{status1}{Title label for events ("status"=1).}
  \item{\dots}{Additional parameters.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\references{Blanche P, Dartigues J-F and Jacqmin-Gadda H.
    \emph{Estimating and comparing time-dependent areas under receiver operating
    characteristic curves for censored event times with competing risks},
    Statistics in Medicine (2013) 32:5381-5397.
}

\seealso{
\code{\link{eNetXplorer}}, \code{\link{plot}}
}
\examples{
\donttest{
data(breastCancerSurv)
fit = eNetXplorer(x=breastCancerSurv$predictor, y=breastCancerSurv$response, family="cox", 
n_run=25, n_perm_null=15, seed=111, survAUC=TRUE, survAUC_time=c(1,5)*365)
plot(x=fit, plot.type="survROC", survAUC_time=c(1,5)*365, status0="censored", status1="deaths")
}
}
