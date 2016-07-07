\name{bcplm}
\alias{bcplm}


\title{
Bayesian Compound Poisson Linear Models
}
\description{
This function fits Tweedie compound Poisson linear models using Markov Chain Monte Carlo methods.
}
\usage{

bcplm(formula, link = "log", data, inits = NULL,
  weights, offset, subset, na.action, contrasts = NULL, 
  n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter / 2),
  n.thin = max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
  n.sims = 1000, n.report = 2, prior.beta.mean = NULL, 
  prior.beta.var = NULL, bound.phi = 100, bound.p = c(1.01, 1.99), 
  tune.iter = 5000, n.tune = floor(tune.iter/100),  
  basisGenerators = c("tp", "bsp", "sp2d"), doFit = TRUE, ...)  
    
}

\arguments{
 \item{formula}{an object of class \code{formula}. See \code{\link[stats]{glm}} and \code{cpglmm} for details.
}
  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link = "log"}.
}
 \item{inits}{a list of initial values to be used for each chain. It must be of length \code{n.chains}. Each element is a named list with the following components: 'beta' (fixed effects), 'phi' (dispersion), and 'p' (index parameter). If the formula indicates a mixed model, it must also contain two additional members 'u' (random effects) and 'Sigma' (variance components). 'Sigma' must be a list of the same format as the \code{ST} slot in \code{cpglmm}. If not supplied, the function will generate initial values automatically.
}

  \item{data, subset, weights, na.action, offset, contrasts}{further model specification arguments as in \code{\link{cpglm}}; see there for details.}
    
  \item{n.chains}{an integer indicating the number of Markov chains (default: \code{3}).  
}
  \item{n.iter}{
the number of total iterations per chain (including burn in; default: \code{2000})
}
  \item{n.burnin}{
the length of burn in, i.e. number of iterations to discard at the beginning. Default
is \code{n.iter/2}, that is, discarding the first half of the simulations.
}
  \item{n.thin}{
thinning rate. Must be a positive integer. Set \code{n.thin > 1} to save memory and
computation time if \code{n.iter} is large. Default is \code{max(1, floor(n.chains * (n.iter - n.burnin) / 1000))} which will only thin if there are at
least \code{2000} simulations.
}
  \item{n.sims}{
The approximate number of simulations to keep after thinning (all chains combined).
}
  \item{n.report}{if greater than zero, fitting information will be printed out  \code{n.report} times for each chain. 
}
  \item{prior.beta.mean}{a vector of prior means for the fixed effects. Default is a vector of zeros. 
}
  \item{prior.beta.var}{
a vector of prior variances for the fixed effects. Default is a vector of \code{10000}'s.
}
  \item{bound.phi}{a numeric value indicating the upper bound of the uniform prior for the dispersion parameter. The default is \code{100}. The lower bound is set to be \code{0} in the function.
}
  \item{bound.p}{a vector of lower and upper bounds for the index parameter \eqn{p}. The default is \code{c(1.01, 1.99)}.
}

\item{tune.iter}{the number of iterations used for tuning the proposal variances used in the Metropolis-Hastings updates. These iterations will not be included in the final output. Default is \code{5000}. Set it to be zero if the tuning process is not desired.  
}
\item{n.tune}{a positive integer (default: \code{20}). The \code{tune.iter} iterations is divided into \code{n.tune} loops. Proposal variances are updated at the end of each loop if acceptance rates are outside the desired interval.
}

 \item{basisGenerators}{
a character vector of names of functions that generate spline bases. See \code{\link{tp}} for details. 
}

  \item{doFit}{
if \code{FALSE}, the constructed \code{"bcplm_input"} object is returned before the model is fitted. 
}

  \item{\dots}{ not used.
}
}
\details{
This function provides Markov chain Monte Carlo [MCMC] methods for fitting Tweedie compound Poisson linear models within the Bayesian framework. Both  generalized linear models and mixed models can be handled. In computing the posterior distribution, the series evaluation method (see, e.g., \code{\link[tweedie]{dtweedie}}) is employed to evaluate the compound Poisson density.    

In the Bayesian model, prior distributions have to be specified for all parameters in the model. Here, Normal distributions are used for the fixed effects (\eqn{\beta}), a Uniform distribution for the  dispersion parameter (\eqn{\phi}),  a Uniform distribution for the index parameter (\eqn{p}). If a mixed model is specified,  prior distributions must be specified for the  variance component. If there is one random effect in a group, the inverse Gamma (scale = \code{0.001}, shape = \code{0.001}) is specified as the prior. If there is more than one random effects in a group, the inverse Wishart (identity matrix as the scale and the dimension of the covariance matrix as the shape)  is specified as the prior.

Prior means and variances of the fixed effects can be supplied using the argument \code{prior.beta.mean} and \code{prior.beta.var}, respectively. The prior distribution of \eqn{\phi} is uniform on (0, \code{bound.phi}). And the bounds of the Uniform for \eqn{p} can be specified in the argument \code{bound.p}. See details in section 'Arguments'. 

In implementing the MCMC, a Gibbs sampler is constructed in which parameters are updated one at a time given the current values of all the other parameters. Specifically, we use the random-walk Metropolis-Hastings algorithm in updating each parameter except for the variance components, which can be simulated directly due to conjugacy.

Before the MCMC, there is a tuning process where the proposal variances of the (truncated) Normal proposal distributions are updated according to the sample variances computed from the simulations in each tuning loop. The goal is to make the acceptance rate roughly between 40\% and 60\% for univariate M-H updates. The argument \code{tune.iter} determines how many iterations are used for the tuning process, and \code{n.tune} determines how many loops these iterations should be divided into. These iterations will not be used in the final output. 

The simulated values of all model parameters are stored in the  \code{sims.list} slot of the returned \code{bcplm} object. It is a list of \code{n.chains} matrices and each matrix has approximately \code{n.sims} rows. The \code{sims.list} slot is further coerced to be of class \code{"mcmc.list"} so that various methods from the  \code{coda} package can be directly applied  to get Markov chain diagnostics, posterior summary and plots. See \code{coda} for available methods.    
}

\value{
  \code{bcplm} returns an object of class \code{"bcplm"}. See \code{\link{bcplm-class}} for details of the return values as well as various methods available for this class.
}

\references{
\cite{ Zhang, Y (2013). Likelihood-based and Bayesian Methods for Tweedie Compound Poisson Linear Mixed Models, \emph{Statistics and Computing}, 23, 743-757.
}
}

\author{
Yanwei (Wayne) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{bcplm-class}}, \code{\link{cpglm}}, \code{\link{cpglmm}}, \code{\link[coda]{mcmc}}, and \code{\link[statmod]{tweedie}} for related information.
}

\examples{

\dontrun{

# fit the FineRoot data with Bayesian models
# Bayesian cpglm
set.seed(10)
fit1 <- bcplm(RLD ~ factor(Zone) * factor(Stock), 
            data = FineRoot, tune.iter = 2000, 
            n.iter = 6000, n.burnin = 1000, n.thin = 5)

gelman.diag(fit1$sims.list)
# diagnostic plots                             
acfplot(fit1$sims.list, lag.max = 20)
xyplot(fit1$sims.list)                              
densityplot(fit1$sims.list)               
summary(fit1)
plot(fit1)
                          

# now fit the Bayesian model to an insurance loss triangle 
# (see Peters et al. 2009)
fit2 <- bcplm(increLoss ~ factor(year) + factor(lag), 
            data = ClaimTriangle, n.iter = 12000, 
            n.burnin = 2000, n.thin = 10, bound.p = c(1.1, 1.95))
gelman.diag(fit2$sims.list)                   
summary(fit2)

# mixed models 
set.seed(10)
fit3 <- bcplm(RLD ~ Stock * Zone + (1|Plant), 
            data = FineRoot, n.iter = 15000, 
            n.burnin = 5000, n.thin = 10)
gelman.diag(fit3$sims.list)                   
summary(fit3)



}

}

\keyword{models}

