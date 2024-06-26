% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Update_alpha_SLCMRT}
\alias{Update_alpha_SLCMRT}
\title{Update Alpha}
\usage{
Update_alpha_SLCMRT(
  N,
  J,
  nClass,
  Y,
  CLASS,
  PY_a,
  pis,
  Atable,
  alpha,
  d0,
  logt,
  Eta,
  Tau,
  H,
  sigma2_t
)
}
\arguments{
\item{N}{Number of respondents.}

\item{J}{Number of items.}

\item{nClass}{Number of classes.}

\item{Y}{Matrix of categorical responses.}

\item{CLASS}{Vector of class labels.}

\item{PY_a}{Cube of conditional probabilities.}

\item{pis}{Vector of latent class proportion parameters.}

\item{Atable}{Matrix representing the attribute profile table.}

\item{alpha}{Matrix of attribute profiles by respondents.}

\item{d0}{Dirichlet distribution parameter.}

\item{logt}{Matrix of log-transformed response times.}

\item{Eta}{Vector of item response time parameters.}

\item{Tau}{Vector of latent speed coefficients.}

\item{H}{Matrix of item-attribute response time coefficients.}

\item{sigma2_t}{Vector of item response time variances.}
}
\value{
A list containing updated values.
}
\description{
Update alpha values using MCMC sampling.

This function updates alpha values based on provided parameters using MCMC sampling.
}
