% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgenerics.R
\name{availability}
\alias{availability}
\title{Availability Function}
\usage{
availability(x, k, upstates = x$states, level = 0.95, klim = 10000)
}
\arguments{
\item{x}{An object of S3 class \code{smmfit} or \code{smm}.}

\item{k}{A positive integer giving the time at which the availability
should be computed.}

\item{upstates}{Vector giving the subset of operational states \eqn{U}.}

\item{level}{Confidence level of the asymptotic confidence interval. Helpful
for an object \code{x} of class \code{smmfit}.}

\item{klim}{Optional. The time horizon used to approximate the series in the
computation of the mean sojourn times vector \eqn{m} (cf.
\link{meanSojournTimes} function) for the asymptotic variance.}
}
\value{
A matrix with \eqn{k + 1} rows, and with columns giving values of
the availability, variances, lower and upper asymptotic confidence limits
(if \code{x} is an object of class \code{smmfit}).
}
\description{
The pointwise (or instantaneous) availability of a system
\eqn{S_{ystem}} at time \eqn{k \in N} is the probability that the system
is operational at time \eqn{k} (independently of the fact that the system
has failed or not in \eqn{[0, k)}).
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

We are interested in investigating the availability of a discrete-time
semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose that the
evolution in time of the system is governed by an E-state space
semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given
at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}},
for a certain \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in
operating mode \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a
certain \eqn{j \in D}, means that the system is not operational at time
\eqn{k} due to the mode of failure \eqn{j} or that the system is under the
repairing mode \eqn{j}.

The pointwise (or instantaneous) availability of a system \eqn{S_{ystem}}
at time \eqn{k \in N} is the probability that the system is operational
at time \eqn{k} (independently of the fact that the system has failed or
not in \eqn{[0, k)}).

Thus, the pointwise availability of a semi-Markov system at time
\eqn{k \in N} is

\deqn{A(k) = P(Z_k \in U) = \sum_{i \in E} \alpha_i A_i(k),}

where we have denoted by \eqn{A_i(k)} the conditional availability of the
system at time \eqn{k \in N}, given that it starts in state \eqn{i \in E},

\deqn{A_i(k) = P(Z_k \in U | Z_0 = i).}
}
\references{
V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov
Models Toward Applications - Their Use in Reliability and DNA Analysis.
New York: Lecture Notes in Statistics, vol. 191, Springer.
}
