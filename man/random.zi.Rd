% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{random.zi}
\alias{random.zi}
\title{Title   random.zi: Simulation data. Genes with random pattern.}
\usage{
random.zi(
  lambda = 0.5,
  spots,
  g,
  type = "ZINB",
  g.p,
  g.size = 0.25,
  g.mu,
  outlier = F
)
}
\arguments{
\item{lambda}{The spots detection efficiency.It should be between 0 and 1.
The default value is 0.5.}

\item{spots}{The number of spots.}

\item{g}{The number of genes.}

\item{type}{A character.'ZINB' or 'ZIP. The default value is ZINB.}

\item{g.p}{A number. The zero proportion of the zero generation process.}

\item{g.size}{The size para in the NB distribution.}

\item{g.mu}{the lambda para in the poisson distribution or the mu para in the NB distribution}

\item{outlier}{A numeric value should be between 0 and 1.}
}
\value{
A data.frame of gene expression counts and coordinates.
spots: the sample size
g : the number of genes
type: corresponding to the 'family' parameter of the simu_zi function
g.p: corresponding to the 'zi.p' parameter of the simu_zi function
}
\description{
Title   random.zi: Simulation data. Genes with random pattern.
}
