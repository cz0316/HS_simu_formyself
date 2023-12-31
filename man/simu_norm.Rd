% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0907_new_basic.R
\name{simu_norm}
\alias{simu_norm}
\title{Title simu_norm}
\usage{
simu_norm(subject.n, zi.p = 0.5, mu = 0.5, sd = 0.25)
}
\arguments{
\item{subject.n}{sample size}

\item{zi.p}{The zero proportion of the zero generation process.
if  zi.p=0, it generates data from normal distribution.}

\item{mu}{the mean of the normal distribution.}

\item{sd}{the sd of the naromal distribution.}
}
\value{
simulation data.
}
\description{
A intermediate step in the new.norm.
}
