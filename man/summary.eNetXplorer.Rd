\name{summary}
\alias{summary}
\alias{summary.eNetXplorer}
\title{generates list of model statistics}
\description{
This function generates a standard list of model statistics. For each \code{alpha}, it contains the 
best value of \code{lambda} (obtained by maximizing a quality function over out-of-bag instances), 
the corresponding maximum value of the quality function, and the model significance (p-value 
based on comparison to permutation null models).}
\usage{
\method{summary}{eNetXplorer}(object, \dots)
}
\arguments{
  \item{object}{\code{eNetXplorer} object.}
  \item{\dots}{Additional parameters.}
}
\value{
  \item{alpha}{Vector of \code{alpha} values.}
  \item{best_lambda}{Best \code{lambda} obtained by maximization of the quality function.}
  \item{model_QF_est}{Maximum of the quality function.}
  \item{QF_model_vs_null_pval}{P-value from model vs null comparison to assess statistical significance.}
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{eNetXplorer}}}
\examples{
\dontrun{
data(QuickStartEx)
fit = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",n_run=20,n_perm_null=10,seed=111)
summary(fit)
}
}