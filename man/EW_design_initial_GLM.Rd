% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EW_design_initial_GLM.R
\name{EW_design_initial_GLM}
\alias{EW_design_initial_GLM}
\title{function to generate a initial EW Design for generalized linear models}
\usage{
EW_design_initial_GLM(
  k.continuous,
  factor.level,
  Integral_based,
  b_matrix,
  joint_Func_b,
  Lowerbounds,
  Upperbounds,
  xlist_fix = NULL,
  lvec,
  uvec,
  h.func,
  link = "continuation",
  delta0 = 1e-06,
  epsilon = 1e-12,
  maxit = 1000
)
}
\arguments{
\item{k.continuous}{number of continuous variables}

\item{factor.level}{list of distinct factor levels, “(min, max)” for continuous factors that always come first, finite sets for discrete factors.}

\item{Integral_based}{TRUE or FALSE, whether or not integral-based EW D-optimality is used, FALSE indicates sample-based EW D-optimality is used.}

\item{b_matrix}{matrix of bootstrapped or simulated parameter values.}

\item{joint_Func_b}{prior distribution function of model parameters}

\item{Lowerbounds}{vector of lower ends of ranges of prior distribution for model parameters.}

\item{Upperbounds}{vector of upper ends of ranges of prior distribution for model parameters.}

\item{xlist_fix}{list of discrete factor experimental settings under consideration, default NULL indicating a list of all possible discrete factor experimental settings will be used.}

\item{lvec}{lower limit of continuous variables}

\item{uvec}{upper limit of continuous variables}

\item{h.func}{function, is used to transfer the design point to model matrix (e.g. add interaction term, add intercept)}

\item{link}{link function, default "continuation", other options "baseline", "adjacent" and "cumulative"}

\item{delta0}{tuning parameter, the distance threshold, || x_i(0) - x_j(0) || >= delta0}

\item{epsilon}{determining f.det > 0 numerically, f.det <= epsilon will be considered as f.det <= 0}

\item{maxit}{maximum number of iterations}
}
\value{
X        matrix of initial design point

p0       initial random approximate allocation

f.det    the determinant of the expected Fisher information matrix for the initial design
}
\description{
function to generate a initial EW Design for generalized linear models
}
