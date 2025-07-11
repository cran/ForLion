% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nu1_log_self.R
\name{nu1_log_self}
\alias{nu1_log_self}
\title{function to calculate first derivative of nu function given eta for log link}
\usage{
nu1_log_self(x)
}
\arguments{
\item{x}{vector of eta, eta=X*beta}
}
\value{
the first derivative of nu function given eta for log link
}
\description{
function to calculate first derivative of nu function given eta for log link
}
\examples{
eta = c(1,2,3,4)
nu1_log_self(eta)

}
