\name{cpglmm}
\alias{cpglmm}
\title{
Compound Poisson Generalized Linear Mixed Models
}
\description{
Laplace approximation and adaptive Gauss-Hermite quadrature methods for compound Poisson  mixed and additive models. 
}
\usage{
cpglmm(formula, link = "log", data, weights, offset, subset, 
    na.action, inits = NULL,  contrasts = NULL, 
    control = list(), basisGenerators = c("tp", "bsp", "sp2d"),
    optimizer = "nlminb", doFit = TRUE, nAGQ = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a two-sided linear formula object describing the  model structure, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. The vertical bar character "|" separates an expression for a model matrix and a grouping factor. The right side can also include basis generators. See \code{lme4} and \code{basisGenerators} below.}

  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link="log"}.
}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
}
 
  \item{subset, weights, na.action, offset, contrasts}{further model specification arguments as in \code{\link{cpglm}}; see there for details.}
  
\item{inits}{a named list with three components 'beta', 'phi', 'p', 'Sigma' that supply the initial values used in the optimization. If not supplied, the function will generate initial values automatically, which are based on a GLM with the supplied model structure.
}
  \item{control}{
a list of parameters for controlling the fitting process. See \code{\link{cpglm}}. The parameter \code{PQL.init} is not used.
}
  \item{basisGenerators}{
a character vector of names of functions that generate spline bases. This is used when smoothing effects are to be included in the model. See \code{\link{tp}} for details. 
}
  \item{optimizer}{
a character string that determines which optimization routine is to be used. Possible choices are \code{"nlminb"} (the default, see \code{\link[stats]{nlminb}}), \code{"bobyqa"} (\code{\link[minqa]{bobyqa}}) and \code{"L-BFGS-B"} (\code{\link[stats]{optim}}).   
}
  \item{doFit}{
if \code{FALSE}, the constructed \code{"cpglmm"} object is returned before the model is fitted. 
}
  \item{nAGQ}{
  a positive integer - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. This defaults to 1, corresponding to the Laplacian approximation. Values greater than 1 produce greater accuracy in the evaluation of the log-likelihood at the expense of speed.
  }
}

\details{
Estimation of compound Poisson mixed models in  existing software has been limited to the Penalized Quasi-Likelihood [PQL] approach (e.g., see \code{\link[MASS]{glmmPQL}}). While straightforward and fast, this method is not equipped to estimate the unknown variance function, i.e., the index parameter. In contrast, the function \code{cpglmm} implements true likelihood-based inferential procedures, i.e., the Laplace approximation and the Adaptive Gauss-Hermite Quadrature (for single grouping factor), so that all parameters in the model can be estimated using maximum likelihood estimation. 
 
This implementation is based on the older \code{lme4} package (the \code{0.9*} version), with changes made on updating of the mean, the variance function and the marginal loglikelihood. For the Laplace method, the contribution of the dispersion parameter to the approximated loglikelihood is explicitly accounted for, which should be more accurate and more consistent with the quadrature estimate. Indeed, both the dispersion parameter and  the index parameter are included as a part of the optimization process. In computing the marginal loglikelihood, the density of the compound Poisson distribution is approximated using numerical methods provided in the \code{tweedie} package. For details of the Laplace approximation and the Gauss-Hermite quadrature method for generalized linear mixed models, see the documentation associated with \code{lme4}. 

In addition, similar to the package \code{amer} (already retired from CRAN), we provide convenient interfaces for fitting additive models using penalized splines.  See the 'example' section for one such application.  



}
\value{
  \code{cpglmm} returns an object of class \code{cpglmm}. See \code{\link{cpglmm-class}} for details of the return values as well as various method available for this class.
}


\references{
\cite{ Zhang Y (2013). Likelihood-based and Bayesian Methods for Tweedie Compound Poisson Linear Mixed Models, \emph{Statistics and Computing}, 23, 743-757. 
} 
\url{https://github.com/actuaryzhang/cplm/files/144051/TweediePaper.pdf}

\cite{Bates D, Maechler M, Bolker B and Walker S (2015). \code{lme4}: Linear mixed-effects models using Eigen and S4..
}
}

\author{
Yanwei (Wayne)) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see \code{\link{cpglm}} for a general introduction to the compound Poisson distribution, \code{lme4} for syntax and usage of mixed-effect models and \code{\link{cpglmm-class}} for detailed explanation of the return value.
}
\examples{
\dontrun{
# use Stock and Spacing as main effects and Plant as random effect
(f1 <- cpglmm(RLD ~ Stock + Spacing +  (1|Plant), data = FineRoot))
            
coef(f1); fixef(f1); ranef(f1)  #coefficients
VarCorr(f1)  #variance components

# add another random effect
(f2 <- update(f1, . ~ . + (1|Zone)))
# test the additional random effect
anova(f1,f2)

# try a different optimizer 
(f3 <- cpglmm(RLD ~  Stock + Spacing +  (1|Plant), 
            data = FineRoot, optimizer = "bobyqa", 
            control = list(trace = 2)))

# adaptive G-H quadrature  
(f4 <- cpglmm(RLD ~  Stock + Spacing +  (1|Plant), 
            data = FineRoot, nAGQ = 3))

# a model with smoothing effects
(f5 <- cpglmm(increLoss ~ tp(lag, k = 4) + (1|year) , 
            data = ClaimTriangle))
}
}

\keyword{models}
