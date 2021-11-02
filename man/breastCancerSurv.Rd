\name{breastCancerSurv}
\alias{breastCancerSurv}
\docType{data}
\title{gene signature for breast cancer survival}
\description{
Gene signature proposed by Desmedt et al associated with breast cancer survival. Microarray data from van't Veer et al and van de Vijver et al. Dataset adapted from Schroeder et al.
}
\usage{data(breastCancerSurv)}

\format{
A numerical matrix of predictors is provided with subjects as rows and genes as columns. The response is a two-column matrix with survival time (in days) and status (0=censored, 1=dead).}

\details{
Desmedt et al proposed a gene signature associated with breast cancer clinical outcome that captures different biological processes: AURKA (proliferation), PLAU (tumor invasion/metastasis), STAT1 (immune response), VEGFA (angiogenesis), CASP3 (apoptosis), ESR1 (ER signaling) and ERBB2 (HER2 signaling). Microarray data was obtained from van't Veer et al and van de Vijver et al. The dataset was adapted from Schroeder et al. Missing predictor data was imputed using the \code{missForest} package. Subjects with missing survival data were removed.}

\references{
Desmedt C et al.
\emph{Biological Processes Associated with Breast Cancer Clinical Outcome Depend on the Molecular Subtypes}, Clinical Cancer Research (2008) 14(16):5158-5165.

Schroeder MS et al.
\emph{survcomp: an R/Bioconductor package for performance assessment and comparison of survival models}, Bioinformatics (2011) 27(22):3206-3208.

van de Vijver MJ et al.
\emph{A Gene Expression Signature as a Predictor of Survival in Breast Cancer}, New England Journal of Medicine (2002) 347(25):1999-2009.

van't Veer LJ et al.
\emph{Gene expression profiling predicts clinical outcome of breast cancer}, Nature (2002) 415:530-536.
}

\keyword{datasets}

