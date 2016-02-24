\name{zcpglm}
\alias{zcpglm}
\title{Zero-inflated Compound Poisson Generalized Linear Models
}
\description{This function fits zero-inflated compound Poisson generalized linear models. 
}
\usage{
zcpglm(formula, link = "log", data, weights, offset, 
      subset, na.action = NULL, contrasts = NULL, 
      control = list(), optimizer = "nlminb")                     
}

\arguments{
  \item{formula}{an object of class \code{formula}. See details below. 
}
  \item{link}{a specification for the compound Poisson model link function. See \code{\link{cpglm}}. For the zero-inflation model, the logit link is used. 
}
  \item{data, subset, na.action, contrasts, control}{see \code{\link{cpglm}} for details.
}
  \item{weights, offset}{prior weights and offset. If specified, they will be used in both the zero-inflation and the compound Poisson model. See details below for specifying different offsets for each model. 
}
  \item{optimizer}{
a character string that determines which optimization routine is to be used. Possible choices are \code{"nlminb"} (the default, see \code{\link[stats]{nlminb}}), \code{"bobyqa"} (\code{\link[minqa]{bobyqa}}) and \code{"L-BFGS-B"} (\code{\link[stats]{optim}}).   
}
}

\details{
 
This function implements zero-inflated compound Poisson generalized linear models. This is similar to the zero-inflated Poisson model for count data, with the Poisson distribution replaced by the compound Poisson distribution. Specifically, the observation is allowed to com from a degenerate distribution at zero with a positive probability, in addition to the regular compound Poisson process. This latent zero-inflation part is specified using a logistic regression structure, while the compound Poisson component is modeled the same as in \code{cpglm}. Parameters are estimated by maximizing the marginal likelihood, and the variance-covariance matrix is computed numerically.   

The formula specification is similar to that in \code{zeroinfl} (package \code{pscl}) except that we use \code{||} instead of \code{|} to separate the two parts of the model. For example, the formula \code{y ~ x1 + x2 || 1} indicates that \code{y ~ x1 + x2} is used for the compound Poison part and \code{~ 1} for the zero-inflation part. If the inflation part is omitted, it defaults to a model with an intercept. Offsets can be specified in both components of the model, e.g., \code{y ~ x1 + offset(x2) | z1 + z2 + offset(z3)}. 

More details of this model are available in the package vignettes. 
}

\value{
  \code{zcpglm} returns an object of class \code{"zcpglm"}. See \code{\link{zcpglm-class}} for details of the return values as well as various methods available for this class. 
}


\author{
Yanwei (Wayne) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{zcpglm-class}}, \code{\link{cpglm}} and \code{zeroinfl} for related information.
}

\examples{
da <- subset(AutoClaim, IN_YY == 1) # use data in the Yip and Yau paper
da <- transform(da, CLM_AMT5 = CLM_AMT5/1000,
                 INCOME = INCOME/10000)
(Z0 <- zcpglm(CLM_AMT5 ~ CAR_USE + MARRIED + AREA + MVR_PTS||
             MVR_PTS + INCOME, data = da))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}

