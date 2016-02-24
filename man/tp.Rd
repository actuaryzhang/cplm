\name{tp}
\alias{tp}
\alias{bsp}
\title{Generate basis functions for penalized spline smoothing.}
\usage{
tp(x, degree=1, k=15, by=NULL, allPen=FALSE, varying=NULL, diag=FALSE,
    knots=quantile(x, probs = (1:(k - degree))/(k - degree  + 1)), 
    centerscale=NULL, scaledknots=FALSE)
    
bsp(x, k=15, spline.degree=3, diff.ord=2, knots, by,
    allPen=FALSE, varying, diag=FALSE)
    }
\note{These are from the \code{amer} package that has retired from CRAN.}
\description{\code{tp} generates a truncated power basis and \code{bsp} generates a reparameterized b-spline basis for penalized spline smoothing.}

\details{

\code{tp} generates truncated power bases which have \code{degree} unpenalized basis functions, namely \eqn{x^1,\dots, x^{degree}} and \eqn{k-}\code{degree} penalized basis functions that contain the positive part \eqn{(x-\kappa_j)^{degree}} for knots \eqn{\kappa_j, j=1,dots,k-}\code{degree}.
This function can be used as a reference when implementing other \code{basisGenerators} that can be used for additive models through \code{cpglmm}. 

\code{bsp} generate a b-spline basis with equidistant knots in mixed model reparameterization.}

\value{list with entries:
\code{"X"}: For \code{tp}, it is an \code{n x degree} design matrix for unpenalized part (without intercept) (or a list of those for every level of by if allPen=F); and for \code{bsp}, it is an \code{n x (diff.ord - 1)} design matrix for unpenalized part (without intercept).

\code{"Z"}: For \code{tp}, it is an \code{n x (k-degree)} design matrix for penalized part (or a list of those for every level of by if allPen=F); and for \code{bsp}, it is an \code{n x (k - diff.ord+1)} design matrix for penalized part.}


\arguments{\item{x}{covariate for the smooth function}
\item{degree}{integer: degree of truncated polynomials (0: piecewise constant, 1: piecewise linear etc..)}
\item{k}{integer: dimensionality of the basis (i.e.: number of knots + degree)}
\item{by}{factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.}
\item{allPen}{boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X}
\item{varying}{numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying}
\item{diag}{logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}?}
\item{knots}{vector of knot locations (optional). Defaults to quantile-based knots at the \eqn{i/(k+1-}degree)-quantiles 
		  for \eqn{i=1,\dots,k-}degree.}
\item{centerscale}{numeric(2): center&scale x by these values if not NULL}
\item{scaledknots}{boolean:	are knot locations given for the rescaled x-values?}
\item{spline.degree}{integer: degree of B-splines (defaults to cubic)}
\item{diff.ord}{integer: order of the difference penalty on the un-reparamerized spline coefficients. Defaults to 2, that is, penalized deviations from linearity.}
}

\author{
Fabian Scheipl \email{fabian.scheipl@googlemail.com}
}
