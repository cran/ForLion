% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nu_log_self.R
\name{nu_log_self}
\alias{nu_log_self}
\title{function to calculate w = nu(eta) given eta for log link}
\usage{
nu_log_self(x)
}
\arguments{
\item{x}{Numeric vector of eta, eta = X*beta.}
}
\value{
A numeric vector representing the diagonal elements of the W matrix (nu(eta)).
}
\description{
function to calculate w = nu(eta) given eta for log link
}
\examples{
eta = c(1,2,3,4)
nu_log_self(eta)

}
