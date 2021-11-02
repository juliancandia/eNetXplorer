\name{Leuk_miR_filt}
\alias{Leuk_miR_filt}
\docType{data}
\title{microRNA expression of acute leukemia phenotypes (filtered dataset)}
\description{
Data of human microRNA (miR) expression of 370 miRs from 80 acute myeloid (AML) and acute lymphoblastic (ALL) leukemia cell lines and 60 primary (patient) samples. Acute lymphoblastic leukemia samples are further classified by B-cell (B-ALL) and T-cell (T-ALL) subphenotypes.
}

\format{
The filtered dataset \code{Leuk_miR_filt} consists of 
a numerical matrix of predictors (with samples as rows and miRs as columns) and two categorical response vectors provided for binomial (AML, ALL) and multinomial (AML, B-ALL, T-ALL) classification.
}

\details{
The filtered dataset \code{Leuk_miR_filt} is a subset of the full dataset \code{Leuk_miR_full}, which includes only miRs with median expression >1.2 across all samples. Only leukemia samples (cell lines and primary) were kept. Data available at GEO under Accession Number GSE51908.}

\references{
Tan YS et al.
\emph{Regulation of RAB5C is important for the growth inhibitory effects of MiR-509 in human precursor-B acute lymphoblastic leukemia}, PLoS One (2014) 9:e111777.

Candia J et al.
  \emph{Uncovering low-dimensional, miR-based signatures of acute myeloid and lymphoblastic leukemias with a machine-learning-driven network approach},
  Converg Sci Phys Oncol (2015) 1:025002.\cr
}

\examples{
\donttest{
# Multinomial models to predict acute myeloid (AML), acute B-cell lymphoblastic 
# (B-ALL) and acute T-cell lymphoblastic (T-ALL) leukemias:
data(Leukemia_miR)
fit = eNetXplorer(x=Leuk_miR_filt$predictor, y=Leuk_miR_filt$response_multinomial,
family="multinomial", n_run=25, n_perm_null=15, seed=111)
}
}

\keyword{datasets}

