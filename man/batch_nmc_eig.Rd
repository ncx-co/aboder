% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eig.R
\name{batch_nmc_eig}
\alias{batch_nmc_eig}
\title{Batch Nested Monte Carlo Expected Information Gain}
\usage{
batch_nmc_eig(
  designs,
  model,
  N = 1000,
  M = 1000,
  cores = 1,
  update_model = TRUE,
  ...
)
}
\arguments{
\item{designs}{a list of design dataframes}

\item{model}{a fitted brms object}

\item{N}{int, how many outer loop Monte Carlo samples}

\item{M}{int, how many inner loop Monte Carlo samples}

\item{cores}{int, how many cores and chains to run MCMC with}

\item{update_model}{bool, update model with new samples}

\item{...}{optional parameters passed to \code{update_model}}
}
\value{
expected information gain (float)
}
\description{
Calculate EIG for a batch of designs.
}
