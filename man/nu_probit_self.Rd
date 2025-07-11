% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nu_probit_self.R
\name{nu_probit_self}
\alias{nu_probit_self}
\title{function to calculate w = nu(eta) given eta for probit link}
\usage{
nu_probit_self(x)
}
\arguments{
\item{x}{vector of eta, eta=X*beta}
}
\value{
diagonal element of W matrix which is nu(eta)
}
\description{
function to calculate w = nu(eta) given eta for probit link
}
\examples{
eta = c(1,2,3,4)
nu_probit_self(eta)


}
