\name{cplm-package}
\alias{cplm-package}
\alias{cplm}
\docType{package}
\title{
  Tweedie compound Poisson linear models 
}
\description{The Tweedie compound Poisson distribution is a mixture of a degenerate distribution at the origin and a continuous distribution on the positive real line. It has been applied in a wide range of fields in which continuous data with exact zeros regularly arise. Nevertheless, statistical inference based on full likelihood and Bayesian methods is not available in most statistical software, largely because the distribution has an intractable density function and numerical methods that allow fast and accurate evaluation of the density did not appear until fairly recently. The \code{cplm} package provides likelihood-based and Bayesian procedures for fitting common Tweedie compound Poisson linear models. In particular, models with hierarchical structures or extra zero inflation can be handled. Further, the package implements the Gini index based on an ordered version of the Lorenz curve as a robust model comparison tool involving zero-inflated and highly skewed distributions.  

The following features of the package may be of special interest to the users:

\enumerate{
\item All methods available in the package enable the index parameter (i.e., the unknown variance function) to be estimated from the data.
\item The compound Poisson generalized linear model handles large data set using the bounded memory regression facility in \code{biglm}.
\item For mixed models, we provide likelihood-based methods using Laplace approximation and adaptive Gauss-Hermit quadrature. 
\item A convenient interface is offered to fit additive models (penalized splines) using the mixed model estimation procedure.
\item Self-tuned Markov chain Monte Carlo procedures are available for both GLM-type and  mixed models.
\item The package also implements a zero-inflated compound Poisson model, in which the observed frequency of zeros can generally be more adequately modeled. 
\item We provide the Gini index based on an ordered Lorenz curve, which is better suited for model comparison involving the compound Poisson distribution.
}
}

\author{
Yanwei (Wayne) Zhang <actuary_zhang@hotmail.com>
}
\references{
\cite{Dunn, P.K. and Smyth, G.K. (2005). Series evaluation of Tweedie exponential dispersion models densities. \emph{Statistics and Computing}, 15, 267-280.}

\cite{Frees, E. W., Meyers, G. and Cummings, D. A. (2011). Summarizing Insurance Scores Using
a Gini Index. \emph{Journal of the American Statistical Association}, 495, 1085 - 1098.
}

\cite{ Zhang, Y (2013). Likelihood-based and Bayesian Methods for Tweedie Compound Poisson Linear Mixed Models, \emph{Statistics and Computing}, 23, 743-757.
}
}


\keyword{ package }


