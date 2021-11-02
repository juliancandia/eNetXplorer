\name{Leuk_miR_full}
\alias{Leuk_miR_full}
\docType{data}
\title{microRNA expression of acute leukemia phenotypes (full dataset)}
\description{
Data of human microRNA (miR) expression of 847 miRs from 80 acute myeloid (AML) and acute lymphoblastic (ALL) leukemia cell lines, 60 primary (patient) samples, and 50 normal control samples sorted by cell type (CD34+ HSPC, Granulocytes, Monocytes, T-cells and B-cells). Acute lymphoblastic leukemia samples are further classified by B-cell (B-ALL) and T-cell (T-ALL) subphenotypes.
}

\format{
The full dataset \code{Leuk_miR_full} consists of a numerical matrix of expression (with samples as rows and miRs as columns) and two data frames with sample and miR metadata.
}

\details{
Data available at GEO under Accession Number GSE51908.}

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
predictor = Leuk_miR_full$expression_matrix
rownames(predictor) = Leuk_miR_full$sample_metadata$sample
colnames(predictor) = Leuk_miR_full$miRNA_short
response = Leuk_miR_full$sample_metadata$sample_class
fit = eNetXplorer(x=predictor, y=response,
family="multinomial", n_run=25, n_perm_null=15, seed=111)
}
}

\keyword{datasets}

