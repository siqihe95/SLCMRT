% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{permute_attr_order}
\alias{permute_attr_order}
\title{Permute Attribute Order}
\usage{
permute_attr_order(nClass, K, M, perm)
}
\arguments{
\item{nClass}{Number of classes.}

\item{K}{Number of attributes.}

\item{M}{Levels of attribute.}

\item{perm}{Matrix representing the attribute profile table, permat = permuteGeneral(c(0:(K-1))).}
}
\value{
A matrix representing the permuted indices of the attribute profile table.
}
\description{
Permute the order of attributes in the attribute profile table.

This function permutes the order of attributes in the attribute profile table (\code{perm}) based on the provided parameters including
the number of classes (\code{nClass}), the number of attributes (\code{K}), the number of levels of attribute (\code{M}).
}
