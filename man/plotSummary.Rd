\name{plotSummary}
\alias{plotSummary}
\title{
generates summary plots of model performance across alpha
}
\description{
This function generates summary plots to display the performance of all models in the elastic net family. Two measures are used: 
1) mean quality function of response vs out-of-bag predictions, and 2) model vs null p-values. Taken together, these 
plots enable visual assessments of the relative performance among different member models in the 
elastic net family, as well as in relation to permutation null models. 
}
\usage{
plotSummary(x, show.pval.ref = T, main = NULL, col.main = "black", 
cex.main = 0.95, line = 1, \dots)
}
\arguments{
  \item{x}{\code{eNetXplorer} object.}
  \item{show.pval.ref}{Logical to display reference lines of significance (if within range of model vs null p-values). Default is TRUE.}
  \item{main}{Custom title.}
  \item{col.main}{Title color.}
  \item{cex.main}{Title size.}
  \item{line}{Title position.}
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
plot(x=fit, plot.type="summary")
plotSummary(x=fit,show.pval.ref=F,main="custom title",col.main="blue")
}
}
