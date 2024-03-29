% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{simulate_compressedMatrixGP_TP}
\alias{simulate_compressedMatrixGP_TP}
\title{simulate_compressedMatrixGP_TP}
\usage{
simulate_compressedMatrixGP_TP(N, GPobj, tP)
}
\arguments{
\item{N}{Number of observations to draw.}

\item{GPobj}{The H matrix, assembled and factorized for use.}

\item{tP}{A scale factor (or scale vector??) **}
}
\value{
Given N and GPobj, generate N random N(0, 1)s, scale using tP, then
   use HODLR code to do a quick symmetric factor product. Returns
   the product, which are the draws from the process using compression.
}
\description{
Simulate a draw from a compressed Gaussian Process with **something scaled**.
}
