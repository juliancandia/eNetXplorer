\name{mergeObj}
\alias{mergeObj}
\title{
merges eNetXplorer objects with different alphas
}
\description{
Upon sequential or parallel execution of two or more \code{eNetXplorer} runs with different values of the mixing parameter \code{alpha}, and assuming the objects from those runs have been saved, this function creates a new \code{eNetXplorer} object that merges the \code{alpha} values. It currently supports linear (gaussian), logistic (binomial), and Cox regression models.
}
\usage{
mergeObj(source_obj, source_dir=getwd(), dest_obj="eNet_merged.Robj",
dest_dir=NULL)
}
\arguments{
  \item{source_obj}{Vector with the names of two or more \code{eNetXplorer} objects.}
  \item{source_dir}{Source directory. Default is the working directory.}
  \item{dest_obj}{Name of the merged \code{eNetXplorer} object.}
  \item{dest_dir}{Destination directory. If not specified, it will use \code{source_dir} as default.}
}
\value{
An object with S3 class \code{"eNetXplorer"}.
}
\author{Julian Candia and John S. Tsang\cr 
Maintainer: Julian Candia \email{julian.candia@nih.gov}}
\seealso{\code{\link{eNetXplorer}}}
\examples{
\donttest{
# we first generate two objects over different alpha values, then merge them.
# we generate summary PDFs to compare the results before and after merging.
data(QuickStartEx)
working_dir = tempdir()

fit1 = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",alpha=seq(0,1,by=0.2),save_obj=TRUE,dest_dir=working_dir,
dest_obj="eNet1.Robj",n_run=20,n_perm_null=10,seed=111)
summaryPDF(fit1, dest_file="eNet1.pdf",dest_dir=working_dir)

fit2 = eNetXplorer(x=QuickStartEx$predictor,y=QuickStartEx$response,
family="gaussian",alpha=seq(0.1,0.9,by=0.2),save_obj=TRUE,dest_dir=working_dir,
dest_obj="eNet2.Robj",n_run=20,n_perm_null=10,seed=111)
summaryPDF(fit2, dest_file="eNet2.pdf",dest_dir=working_dir)

eNet_merged=mergeObj(source_obj=c("eNet1.Robj","eNet2.Robj"),source_dir=working_dir)
summaryPDF(eNet_merged,dest_file="eNet_merged.pdf",dest_dir=working_dir)
}
}

