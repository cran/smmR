% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgenerics.R
\name{mttf}
\alias{mttf}
\title{Mean Time To Failure (MTTF) Function}
\usage{
mttf(x, upstates = x$states, level = 0.95, klim = 10000)
}
\arguments{
\item{x}{An object of S3 class \code{smmfit} or \code{smm}.}

\item{upstates}{Vector giving the subset of operational states \eqn{U}.}

\item{level}{Confidence level of the asymptotic confidence interval. Helpful
for an object \code{x} of class \code{smmfit}.}

\item{klim}{Optional. The time horizon used to approximate the series in the
computation of the mean sojourn times vector \eqn{m} (cf.
\link{meanSojournTimes} function) for the asymptotic variance.}
}
\value{
A matrix with \eqn{\textrm{card}(U) = s_{1}} rows, and with columns
giving values of the mean time to failure for each state \eqn{i \in U},
variances, lower and upper asymptotic confidence limits (if \code{x} is an
object of class \code{smmfit}).
}
\description{
Consider a system \eqn{S_{ystem}} starting to work at time
\eqn{k = 0}. The mean time to failure (MTTF) is defined as the mean
lifetime.
}
\details{
Consider a system (or a component) \eqn{S_{ystem}} whose possible
states during its evolution in time are \eqn{E = \{1,\dots,s\}}.
Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of
the system (the up states) and by \eqn{D = \{s_1 + 1,\dots,s\}} the
subset of failure states (the down states), with \eqn{0 < s_1 < s}
(obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset},
\eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states
of \eqn{U} as different operating modes or performance levels of the
system, whereas the states of \eqn{D} can be seen as failures of the
systems with different modes.

We are interested in investigating the mean time to failure of a
discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
that the evolution in time of the system is governed by an E-state space
semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at
instant \eqn{0} and the state of the system is given at each instant
\eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain
\eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode
\eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain
\eqn{j \in D}, means that the system is not operational at time \eqn{k}
due to the mode of failure \eqn{j} or that the system is under the
repairing mode \eqn{j}.

Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called
the lifetime of the system, i.e.,

\deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}

The mean time to failure (MTTF) is defined as the mean lifetime, i.e., the
expectation of the hitting time to down set \eqn{D},

\deqn{MTTF = E[T_{D}]}
}
\references{
V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov
Models Toward Applications - Their Use in Reliability and DNA Analysis.
New York: Lecture Notes in Statistics, vol. 191, Springer.

I. Votsi & A. Brouste (2019) Confidence interval for the mean time to
failure in semi-Markov models: an application to wind energy production,
Journal of Applied Statistics, 46:10, 1756-1773
}
