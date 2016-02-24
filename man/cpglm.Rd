\name{cpglm}
\alias{cpglm}

\title{Compound Poisson Generalized Linear Models
}
\description{This function fits compound Poisson generalized linear models. 
}
\usage{
cpglm(formula, link = "log", data, weights, offset, 
          subset, na.action = NULL, contrasts = NULL, 
          control = list(), chunksize = 0, 
          optimizer = "nlminb", ...)                      
}

\arguments{
  \item{formula}{an object of class \code{formula}. See also in \code{\link[stats]{glm}}.
}
  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link = "log"}.
}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
}
  \item{weights}{an optional vector of weights. Should be either \code{NULL} or a numeric vector. When it is numeric, it must be positive. Zero weights are not allowed in \code{cpglm}. 
}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset.  Another possible value is \code{NULL}, no action. Value \code{na.exclude} can be useful.
}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be either \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. 
}
  \item{contrasts}{an optional list. See \code{contrasts.arg}.
}
  \item{control}{a list of parameters for controling the fitting process. See 'Details' below. 
} 
  \item{chunksize}{an integer that indicates the size of chunks for processing the data frame as used in \code{\link[biglm]{bigglm}}. The value of this argument also determines how the model is estimated. When it is \code{0} (the default), regular Fisher's scoring algorithms are used, which may run into memory issues when handling large data sets. In contrast,  a value  greater than \code{0} indicates that the \code{bigglm} is employed to fit the GLMs. The function \code{bigglm} relies on the bounded memory regression technique, and thus is well suited to large data GLMs. 
} 
  \item{optimizer}{a character string that determines which optimization routine is to be used in estimating the index and the dispersion parameters. Possible choices are \code{"nlminb"} (the default, see \code{\link[stats]{nlminb}}), \code{"bobyqa"} (\code{\link[minqa]{bobyqa}}) and \code{"L-BFGS-B"} (\code{\link[stats]{optim}}).   
}
  \item{\dots}{ additional arguments to be passed to \code{bigglm}. Not used when \code{chunksize = 0}. The \code{maxit} argument defaults to \code{50} in \code{cpglm} if not specified. 
}

}

\details{
 
This function implements the profile likelihood approach in Tweedie compound Poisson generalized linear models. First, the index and the dispersion parameters are estimated  by maximizing (numerically) the profile likelihood (profile out the mean parameters as they are determined for a given value of the index parameter). Then the mean parameters are estimated using a GLM with the above-estimated index parameter. To compute the profile likelihood, one must resort to numerical methods provided in the  \code{tweedie} package for approximating the density of the compound Poisson distribution. Indeed, the function  \code{\link[tweedie]{tweedie.profile}} in that package makes available the profile likelihood approach. The \code{cpglm} function  differs from \code{\link[tweedie]{tweedie.profile}} in two aspects. First, the user does not need to specify the grid of possible values the index parameter can take. Rather, the optimization of the profile likelihood is automated.   Second, big data sets can be handled where the \code{bigglm} function from the \code{biglm} package is used in fitting GLMs. The \code{bigglm} is invoked when the argument \code{chunksize} is greater than 0. It is also to be noted that only MLE estimate for the dispersion parameter is included here, while \code{\link[tweedie]{tweedie.profile}} provides several other possibilities.

The package used to implement a second approach using the Monte Carlo EM algorithm, but it is now removed because it does not offer obvious advantages over the profile likelihood approach for this model.


The \code{control} argument is a list that can supply various controlling elements used in the optimization process, and it has the following components:
\describe{
\item{\code{bound.p}}{a vector of lower and upper bounds for the index parameter \eqn{p} used in the optimization. The default is \code{c(1.01, 1.99)}. }
\item{\code{trace}}{if greater than 0, tracing information on the progress of the fitting is produced. For \code{optimizer = "nlminb"} or \code{optimizer = "L-BFGS-B"}, this is the same as the \code{trace} control parameter, and for \code{optimizer = "bobyqa"}, this is the same as the \code{iprint} control parameter. See the corresponding documentation for details. 
}
\item{\code{max.iter}}{maximum number of iterations allowed in the optimization. The default is \code{300}. }
\item{\code{max.fun}}{maximum number of function evaluations allowed in the optimizer. The default is \code{2000}.}
}

}

\value{
  \code{cpglm} returns an object of class \code{"cpglm"}. See \code{\link{cpglm-class}} for details of the return values as well as various methods available for this class. 
}

\references{
\cite{ Dunn, P.K. and Smyth, G.K. (2005). Series evaluation of Tweedie exponential dispersion models densities. \emph{Statistics and Computing}, 15, 267-280.}
}

\author{
Yanwei (Wayne) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{cpglm-class}}, \code{\link[stats]{glm}}, \code{\link[statmod]{tweedie}}, and \code{\link[tweedie]{tweedie.profile}} for related information.
}


\examples{

fit1 <- cpglm(RLD ~ factor(Zone) * factor(Stock),
  data = FineRoot)
     
# residual and qq plot
parold <- par(mfrow = c(2, 2), mar = c(5, 5, 2, 1))
# 1. regular plot
r1 <- resid(fit1) / sqrt(fit1$phi)
plot(r1 ~ fitted(fit1), cex = 0.5)
qqnorm(r1, cex = 0.5)
# 2. quantile residual plot to avoid overlapping
u <- tweedie::ptweedie(fit1$y, fit1$p, fitted(fit1), fit1$phi)
u[fit1$y == 0] <- runif(sum(fit1$y == 0), 0, u[fit1$y == 0])
r2 <- qnorm(u)
plot(r2 ~ fitted(fit1), cex = 0.5)
qqnorm(r2, cex = 0.5)
par(parold)

# use bigglm 
fit2 <- cpglm(RLD ~ factor(Zone), 
  data = FineRoot, chunksize = 250)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}

