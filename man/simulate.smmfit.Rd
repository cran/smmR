% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smmfit.R
\name{simulate.smmfit}
\alias{simulate.smmfit}
\title{Simulates semi-Markov chains}
\usage{
\method{simulate}{smmfit}(object, nsim = 1, seed = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{smmfit} (inheriting from the S3 classes
\code{smm}, \link{smmnonparametric} or \link{smmparametric}).}

\item{nsim}{An integer or vector of integers (for multiple sequences)
specifying the length of the sequence(s).}

\item{seed}{\code{seed} for the random number generator.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
A list of vectors representing the sequences.
}
\description{
Simulates sequences from a fitted semi-Markov model.
}
\details{
If \code{nsim} is a single integer then a chain of that length is
produced. If \code{nsim} is a vector of integers, then \code{length(nsim)}
sequences are generated with respective lengths.
}
\references{
V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov
Models Toward Applications - Their Use in Reliability and DNA Analysis.
New York: Lecture Notes in Statistics, vol. 191, Springer.
}
\seealso{
\link{smmnonparametric}, \link{smmparametric}, \link{fitsmm}
}
