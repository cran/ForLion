#' function for calculating X=h(x) and w=nu(beta^T h(x)) given a design point x = (x1,...,xd)^T
#'
#' @param x x=(x1,...,xd) -- design point/experimental setting
#' @param b b=(b1,...,bp) -- assumed parameter values
#' @param link link = "logit"  -- link function, default: "logit", other links: "probit", "cloglog", "loglog", "cauchit", "log"
#' @param h.func function h(x)=(h1(x),...,hp(x)), default (1,x1,...,xd)
#'
#' @return X=h(x)=(h1(x),...,hp(x)) -- a row for design matrix
#' @return w -- nu(b^t h(x))
#' @return link -- link function applied
#' @export
#'
#' @examples
#' # y -> h(y)=(y1,y2,y3,y4,y5,y4*y5,1) in hfunc
#' hfunc.temp = function(y) {c(y,y[4]*y[5],1);};
#' link.temp="logit"
#' x.temp = c(25,1,1,1,1)
#' b.temp = c(-7.533386, 1.746778, -0.1937022, -0.09704664, 0.1077859, 0.2729715, 0.4293171)
#' Xw_maineffects_self(x.temp, b.temp, link=link.temp, h.func=hfunc.temp)
#'


Xw_maineffects_self <- function(x, b, link="logit", h.func=NULL) {
  if(is.null(h.func)) h.func = function(y) {c(1,y);}; # default: main-effects
  xrow = h.func(x);          # X=(h1(x),...,hp(x))
  eta = sum(b*xrow);
  w = NULL ;
  if(link=="probit") w = nu_probit_self(eta);
  if(link=="cloglog") w = nu_loglog_self(eta);
  if(link=="loglog") w = nu_loglog_self(eta);
  if(link=="cauchit") w = nu_cauchit_self(eta);
  if(link=="log") w = nu_log_self(eta);
  if(is.null(w)) { link="logit"; w=nu_logit_self(eta);};
  list(X=xrow, w=w, link=link);
}
