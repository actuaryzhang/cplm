\name{getF}
\alias{getF}
\alias{plotF}
\title{Get and plot the smoothing function values}
\usage{
getF(object, which, n=100, newdata, interval=c("NONE", "MCMC",
    "RW"), addConst=TRUE, varying=1, level=0.9, sims=1000)
    
plotF(object, which, n=100, interval="RW", addConst=TRUE,
    trans=I, level=0.9, sims=1000, auto.layout=TRUE, rug=TRUE,
    legendPos="topright", ...)
    }
\description{Get and plot the estimated smoothing function values}
\note{These are from the \code{amer} package that has retired from CRAN. The formula used for the pointwise bias-adjusted CIs is taken from Ruppert and Wand's  'Semiparametric Regression' (2003), p. 140. 
These leave out the uncertainty associated with the variance component estimates. }
\value{a list with one \code{data.frame} for each function, giving \code{newdata} or the values of the generated grid plus the fitted values (and confidence/HPD intervals).}


\arguments{\item{object}{a fitted \code{cpglmm} object.}
\item{which}{(optional) an integer vector or a character vector of names giving the smooths for which fitted values are desired. Defaults to all.}
\item{n}{if no \code{newdata} is given, fitted values for a regular grid with n values in the range of the respective covariates are returned}
\item{newdata}{An optional data frame in which to look for variables with which to predict}
\item{interval}{what mehod should be used to compute pointwise confidence/HPD intervals: RW= bias-adjusted empirical bayes}
\item{addConst}{boolean should the global intercept and intercepts for the levels of the by-variable be included in the fitted values (and their CIs) can also be a vector of the same length as \code{which}}
\item{varying}{value of the\code{varying}-covariate (see \code{\link{tp}}) to be used if no newdata is supplied. 
Defaults to 1.}
\item{level}{level for the confidence/HPD intervals}
\item{sims}{how many iterates should  be generated for the MCMC-based HPD-intervals}
\item{trans}{a function that should be applied to the fitted values and ci's before plotting (e.g. the inverse link function to get plots on the scale of the reponse)}
\item{auto.layout}{automagically set plot layout via \code{par()$mfrow}}
\item{rug}{add \code{\link{rug}}-plots of the observed covariate locations}
\item{legendPos}{a (vector of) keyword(s) where to put labels of by-variables (see \code{\link[graphics]{legend}}). "none" if you don't want a legend.}
\item{...}{arguments passed on to the low-level plot functions (\code{plot}, \code{matlines}), \code{legend}, and \code{title}}
}
\author{
Fabian Scheipl \email{fabian.scheipl@googlemail.com}
}
\seealso{See the vignette for examples}
