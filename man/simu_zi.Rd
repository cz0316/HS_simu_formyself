% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_fucntions.R
\name{simu_zi}
\alias{simu_zi}
\title{Title:simu_zi
The basic function of simulation generation functions.}
\usage{
simu_zi(family, subject.n, zi.p = 0.5, mu = 0.5, size = 0.25)
}
\arguments{
\item{family}{'ZINB' or 'ZIP',zero-inflated negative binomial or zero-inflated Poisson distribution.}

\item{subject.n}{the number of spots,because each gene follow a zero-infalted P/NB distribution}

\item{zi.p}{The zero proportion of the zero generation process.
if  zi.p=0, it generates data which follows Poisson or NB distribution.}

\item{mu}{the lambda para in the poisson distribution or the mu para in the NB distribution}

\item{size}{the size para in the NB distribution}
}
\value{
A numeric vector.
}
\description{
Title:simu_zi
The basic function of simulation generation functions.
}
