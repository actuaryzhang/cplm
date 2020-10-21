\name{class-methods}
\docType{class}
\alias{cplm-class}
\alias{$,cplm-method}
\alias{[,cplm,character,missing,missing-method}
\alias{[,cplm,numeric,missing,missing-method}
\alias{[[,cplm,character,missing-method}
\alias{[[,cplm,numeric,missing-method}
\alias{model.matrix,cplm-method}
\alias{names,cplm-method}
\alias{show,cplm-method}
\alias{terms,cplm-method}
\alias{formula,cplm-method}
\alias{vcov,cplm-method}

\alias{cpglm-class}
\alias{coef,cpglm-method}
\alias{fitted,cpglm-method}
\alias{residuals,cpglm-method}
\alias{resid,cpglm-method}
\alias{AIC,cpglm,missing-method}
\alias{deviance,cpglm-method}
\alias{summary,cpglm-method}
\alias{predict,cpglm-method}

\alias{cpglmm-class}
\alias{coef,cpglmm-method}
\alias{fixef,cpglmm-method}
\alias{ranef,cpglmm-method}
\alias{logLik,cpglmm-method}
\alias{anova,cpglmm-method}
\alias{fitted,cpglmm-method}
\alias{residuals,cpglmm-method}
\alias{resid,cpglmm-method}
\alias{print,cpglmm-method}
\alias{summary,cpglmm-method}
\alias{show,cpglmm-method}
\alias{VarCorr,cpglmm-method}
\alias{vcov,cpglmm-method}
\alias{predict,cpglmm-method}



\alias{summary.cpglmm-class}

\alias{bcplm-class}
\alias{plot,bcplm,missing-method}
\alias{summary,bcplm-method}
\alias{show,bcplm-method}
\alias{VarCorr,bcplm-method}
\alias{fixef,bcplm-method}

\alias{NullNum-class}
\alias{NullList-class}
\alias{NullFunc-class}
\alias{ListFrame-class}

\alias{gini-class}
\alias{plot,gini,missing-method}
\alias{show,gini-method}

\alias{VarCorr}
\alias{VarCorr-methods}

\title{Classes and Methods for a Compound Poisson Linear Model Object}

\description{
Documented here are the \code{"cplm"} class and its derived classes \code{"cpglm"}, \code{"cpglmm"}, and \code{"bcplm"}. Several primitive methods and statistical methods are created to facilitate the extraction of specific slots and further statistical analysis. \code{"gini"} is a class that stores the Gini indices and associated standard errors that could be used to perform model comparison involving the compound Poisson distribution. \code{"NullNum"},  \code{"NullList"}, \code{"NullFunc"} and \code{"ListFrame"} are virtual classes for \code{c("NULL", "numeric")}, \code{c("NULL","list")}, \code{c("NULL","function")} and \code{c("list","data.frame")}, respectively. 
}

\section{Objects from the Class}{

\describe{
  \item{\code{"cplm"}}{Objects can be created by calls of the form \code{new("cplm", ...)}.} 
  \item{\code{"cpglm"}}{Objects can be created by calls from \code{new("cpglm", ...)} or \code{cpglm}.} 
  \item{\code{"cpglmm"}}{Objects can be created by calls of the form \code{new("cpglmm", ...)}, or a call to \code{cpglmm}.}
  \item{\code{"summary.cpglmm"}}{Objects can be created by calls of the form \code{new("summary.cpglmm", ...)}, or a call to \code{summary} on a \code{cpglmm} object. }
  \item{\code{"bcplm"}}{Objects can be created by calls from \code{new("bcplm", ...)} or \code{bcplm}.}
  \item{\code{"gini"}}{Objects can be created by calls from \code{new("gini", ...)} or \code{gini}.}  
  \item{\code{"NullNum"}, \code{"NullList"}, \code{"NullFunc"}}{These are  virtual classes and  no objects may be created from them.}
 } 
}

\section{Slots}{
  The \code{"cplm"} class defines the slots common in all the model classes in the \code{cplm} package, and thus the utility methods defined on the \code{"cplm"} class such as \code{[}, \code{names} and so on are applicable to all of the derived classes. 
  \describe{
    \item{\code{call}:}{the matched call. }
    \item{\code{formula}:}{the formula supplied, class \code{"formula"}}
    \item{\code{contrasts}:}{the contrasts used, class \code{"NullList"} }
    \item{\code{link.power}:}{index of power link function, class \code{"numeric"}. See \code{\link[statmod]{tweedie}}.}
    \item{\code{model.frame}:}{the data frame used. class \code{"ListFrame"}. }
    \item{\code{inits}:}{initial values used, class \code{"NullList"}.}
  }
  

 The \code{"cpglm"} class extends \code{"cplm"} directly. Most of the slots have the same definition as those in \code{\link{glm}}.  The following slots are in addition to those in \code{"cplm"}:
\describe{
    \item{\code{coefficients}:}{estimated mean parameters, class \code{"numeric"}. }
    \item{\code{residuals}:}{the working residuals, that is the residuals in the final iteration of the IWLS fit, class \code{"numeric"}}
    \item{\code{fitted.values}:}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function, class \code{"numeric"} }
     \item{\code{linear.predictors}:}{the fitted linear predictors, class \code{"numeric"}}    
    \item{\code{weights}:}{working weights from the last iteration of the iterative least square, class \code{"numeric"}}
    \item{\code{df.residual}:}{residual degrees of freedom, class \code{"integer"}}
    \item{\code{deviance}:}{up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero. This is computed using \code{\link[tweedie]{tweedie.dev}}.}
    \item{\code{aic}:}{a version of Akaike's Information Criterion, minus twice the maximized log-likelihood plus twice the number of mean parameters. This is computed using the tweedie density approximation as in \code{\link[tweedie]{dtweedie}}. }
    \item{\code{offset}:}{the offset vector used, class \code{"NullNum"},} 
    \item{\code{prior.weights}:}{the weights initially supplied, a vector of \code{1}s if none were, class \code{"NullNum"}}        
    \item{\code{y}:}{the response vector used.}
    \item{\code{control}:}{the value of the control argument used, class \code{"list"} }
    \item{\code{p}:}{the maximum likelihood estimate of the index parameter.}    
    \item{\code{phi}:}{the maximum likelihood estimate of the dispersion parameter.}    
    \item{\code{vcov}:}{estimated variance-covariance matrix, class \code{"matrix"}}    
    \item{\code{iter}:}{the number of Fisher's scoring iterations in the final GLM.}
    \item{\code{converged}:}{indicating whether the algorithm has converged, class \code{"logical"}.}
    \item{\code{na.action}:}{method of handling \code{NA}'s, class \code{"NullFunc"}.}
}    
    
  The  \code{"cpglmm"} class extends \code{"cplm"} and the old version of \code{"mer"} class from \code{lme4} directly, and has the following additional slots:
  \describe{
    \item{\code{p}:}{estimated value of the index parameter, class \code{"numeric"}  }
    \item{\code{phi}:}{estimated value of the dispersion parameter, class \code{"numeric"} }
    \item{\code{bound.p}:}{the specified bounds of the index parameter, class \code{"numeric"}   }
    \item{\code{vcov}:}{estimated variance-covariance matrix, class \code{"matrix"}}
    \item{\code{smooths}:}{a list of smooth terms}
    }
The slots it used from the old \code{"mer"} class has the following slots (copied from \code{lme4_0.999999-2}):
  \describe{

    \item{\code{env}:}{An environment (class \code{"environment"})
      created for the evaluation of the nonlinear model function. }

    \item{\code{nlmodel}:}{The nonlinear model function as an object of
      class \code{"call"}. }

    \item{\code{frame}:}{The model frame (class \code{"data.frame"}).}

    \item{\code{call}:}{The matched call to the function that
      created the object. (class \code{"call"}).}

    \item{\code{flist}:}{The list of grouping factors for the random
      effects.}

    \item{\code{X}:}{Model matrix for the fixed effects. }

    \item{\code{Zt}:}{The transpose of model matrix for the random
      effects, stored as a compressed column-oriented sparse matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}).}

    \item{\code{pWt}:}{Numeric prior weights vector. This may be of length
      zero (0), indicating unit prior weights.}

    \item{\code{offset}:}{Numeric offset vector. This may be of length
      zero (0), indicating no offset.}

    \item{\code{y}:}{The response vector (class \code{"numeric"}).}

    \item{\code{Gp}:}{Integer vector of group pointers within the random
      effects vector.  The elements of \code{Gp} are 0-based indices of
      the first element from each random-effects term.  Thus the first
      element is always 0.  The last element is the total length of the
      random effects vector.}

    \item{\code{dims}:}{A named integer vector of dimensions.  Some of
      the dimensions are \eqn{n}, the number of observations, \eqn{p}, the
      number of fixed effects, \eqn{q}, the total number of random
      effects, \eqn{s}, the number of parameters in the nonlinear model
      function and \eqn{nt}, the number of random-effects terms in the
      model.}

    \item{\code{ST}:}{A list of S and T factors in the TSST' Cholesky
      factorization of the relative variance matrices of the random
      effects associated with each random-effects term.  The unit lower
      triangular matrix, \eqn{T}, and the diagonal matrix, \eqn{S}, for
      each term are stored as a single matrix with diagonal elements
      from \eqn{S} and off-diagonal elements from \eqn{T}.}

    \item{\code{V}:}{Numeric gradient matrix (class \code{"matrix"}) of
      the nonlinear model function.}

    \item{\code{A}:}{Scaled sparse model matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
      the the unit, orthogonal random effects, \eqn{U}.}

    \item{\code{Cm}:}{Reduced, weighted sparse model matrix (class
      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for the
      unit, orthogonal random effects, U.  .}

    \item{\code{Cx}:}{The \code{"x"} slot in the weighted sparse model
      matrix (class \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"})
      for the unit, orthogonal random effects, \eqn{U}, in generalized
      linear mixed models.  For these models the matrices \eqn{A} and
      \eqn{C} have the same sparsity pattern and only the \code{"x"}
      slot of \eqn{C} needs to be stored.}

    \item{\code{L}:}{The sparse lower Cholesky factor of \eqn{P(AA'+I)P'}
      (class \code{"\link[Matrix:CHMfactor-class]{dCHMfactor}"}) where \eqn{P}
      is the fill-reducing permutation calculated from the pattern of
      nonzeros in \eqn{A}.}

    \item{\code{deviance}:}{Named numeric vector containing the deviance
      corresponding to the maximum likelihood (the \code{"ML"} element)
      and \code{"REML"} criteria and various components. The
      \code{"ldL2"} element is twice the logarithm of the determinant of
      the Cholesky factor in the \code{L} slot. The \code{"usqr"}
      component is the value of the random-effects quadratic form.}

    \item{\code{fixef}:}{Numeric vector of fixed effects.}

    \item{\code{ranef}:}{Numeric vector of random effects on the
      original scale.}

    \item{\code{u}:}{Numeric vector of orthogonal, constant variance,
      random effects.}

    \item{\code{eta}:}{The linear predictor at the current values of the
      parameters and the random effects.}

    \item{\code{mu}:}{The means of the responses at the current parameter
      values.}

    \item{\code{muEta}:}{The diagonal of the Jacobian of \eqn{\mu}{mu}
      by \eqn{\eta}{eta}.  Has length zero (0) except for generalized
      mixed models.}

    \item{\code{var}:}{The diagonal of the conditional variance of
      \eqn{Y} given the random effects, up to prior weights.  In
      generalized mixed models this is the value of the variance
      function for the \code{\link{glm}} family.}

    \item{\code{resid}:}{The residuals, \eqn{y - \mu}{y-mu}, weighted by
      the \code{sqrtrWt} slot (when its length is \eqn{>0}).}

    \item{\code{sqrtXWt}:}{The square root of the weights applied to the
      model matrices \eqn{X} and \eqn{Z}. This may be of length zero
      (0), indicating unit weights.}

    \item{\code{sqrtrWt}:}{The square root of the weights applied to the
      residuals to obtain the weighted residual sum of squares. This may
      be of length zero (0), indicating unit weights.}

    \item{\code{RZX}:}{The dense solution (class \code{"matrix"}) to
      \eqn{L RZX = ST'Z'X = AX}.}

    \item{\code{RX}:}{The upper Cholesky factor (class \code{"matrix"})
      of the downdated \eqn{X'X}.}

  }
  
  The \code{"summary.cpglmm"} class \emph{contains} the \code{"cpglmm"}
  class and has the following additional slots:
  \describe{
    \item{\code{methTitle}:}{character string specifying a method title}
    \item{\code{logLik}:}{the same as \code{logLik(object)}.}
    \item{\code{ngrps}:}{the number of levels per grouping factor in the
      \code{flist} slot.}
    \item{\code{sigma}:}{the scale factor for the variance-covariance estimates}
    \item{\code{coefs}:}{the matrix of estimates, standard errors,
      etc. for the fixed-effects coefficients}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{AICtab}:}{a named vector of values of AIC, BIC, log-likelihood
      and deviance}
  }    
  
  The  \code{"bcplm"} class extends \code{"cplm"} directly, and has the following additional slots: 
  \describe{
    \item{\code{dims}:}{a named integer vector of dimensions. }
   \item{\code{sims.list}:}{an object of class \code{"mcmc.list"}. It is a list of \code{n.chains} \code{mcmc} objects, each \code{mcmc} object storing the simulation result from a Markov chain. See \code{\link[coda]{mcmc}} and \code{\link[coda]{mcmc.convert}}. Since this is an \code{"mcmc.list"} object, most methods defined in the  \code{coda} package can be directly applied to it. }
    \item{\code{Zt}:}{the transpose of model matrix for the random effects, stored as a compressed column-oriented sparse matrix (class \code{"dgCMatrix"}).}
    \item{\code{flist}:}{the list of grouping factors for the random effects.}  
    \item{\code{prop.var}:}{a named list of proposal variance-covariance matrix used in the Metropolis-Hasting update.}
    
}

   The \code{"gini"} class has the following slots:
  \describe{
    \item{\code{call}:}{the matched call. }
    \item{\code{gini}:}{a matrix of the Gini indices. The row names are corresponding to the base while the column names are corresponding to the scores.}
    \item{\code{sd}:}{a matrix of standard errors for each computed Gini index. }
    \item{\code{lorenz}:}{a list of matrices that determine the graph of the ordered Lorenz curve associated with each base and score combination. For each base, there is an associated matrix.   }
  }
  
}


\section{Extends}{

Class \code{"cpglm"} extends class \code{"\linkS4class{cplm}"}, directly.

Class \code{"cpglmm"} extends class \code{"\linkS4class{cplm}"}, directly;

Class \code{"summary.cpglmm"}  extends class \code{"\linkS4class{cpglmm}"}, directly; 
class \code{"\linkS4class{cplm}"}, by class \code{"\linkS4class{cpglmm}"}, distance 2.

Class \code{"bcplm"} extends class \code{"\linkS4class{cplm}"}, directly.
}


\section{Methods}{
  The following methods are defined for the class \code{"cplm"}, which are also applicable to all of the derived classes:
  \describe{
    \item{$}{\code{signature(x = "cplm")}:  extract a slot of \code{x} with a specified slot name, just as in list. }
    \item{[[}{\code{signature(x = "cplm", i = "numeric", j = "missing")}:  extract the i-th slot of a \code{"cpglm"} object, just as in list. }
    \item{[[}{\code{signature(x = "cplm", i = "character", j = "missing")}: extract the slots of a \code{"cpglm"} object with names in \code{i}, just as in list.}
    \item{[}{\code{signature(x = "cplm", i = "numeric", j = "missing", drop="missing")}:  extract the i-th slot of a \code{"cpglm"} object, just as in list. \code{i} could be a vector. }
    \item{[}{\code{signature(x = "cplm", i = "character", j = "missing", drop="missing")}:  extract the slots of a \code{"cpglm"} object with names in \code{i}, just as in list. \code{i} could be a vector. }    
    \item{names}{\code{signature(x = "cplm")}:  return the slot names. }
    \item{terms}{\code{signature(x = "cplm")}: extract the \code{terms} object from the model frame.  See \code{\link[stats]{terms}}.}
    \item{formula}{\code{signature(x = "cplm")}: extract the \code{formula} slot.  See \code{\link[stats]{formula}}.}    
    \item{model.matrix}{\code{signature(object = "cplm")}: extract the design matrix.  }
    \item{show}{\code{signature(object = "cplm")}: method for \code{show}. }
    \item{vcov}{\code{signature(object = "cplm")}: extract the variance-covariance matrix of a \code{"cplm"} object.}
 }
 
 The following methods are defined for the \code{"cpglm"} class:
  \describe{    
    \item{coef}{\code{signature(object = "cpglm")}:  extract the estimated coefficients.}
    \item{fitted}{\code{signature(object = "cpglm")}: return the fitted values. }
    \item{residuals}{\code{signature(object = "cpglm")}: extract residuals  from a \code{cpglm} object. You can also specify a \code{type} argument to indicate the type of residuals to be computed. See \code{\link[stats]{glm.summaries}}.}
    \item{resid}{\code{signature(object = "cpglm")}: same as \code{residuals}.}   
    \item{AIC}{\code{signature(object = "cpglm",k="missing")}: extract the AIC information from the \code{"cpglm"} object.  See \code{\link[stats]{AIC}}.}    
    \item{deviance}{\code{signature(object = "cpglm")}: extract the deviance from the \code{"cpglm"} object.  See \code{\link[stats]{deviance}}.}
    \item{summary}{\code{signature(object = "cpglm")}:  the same as \code{\link[stats]{glm.summaries}} except that both the dispersion and the index parameter are estimated using maximum likelihood estimation.  }    
    \item{predict}{\code{signature(object = "cpglm")}: generate predictions for new data sets}

   }
   
   
   The following are written for \code{"cpglmm"}: 
  \describe{
    \item{print}{\code{signature(x = "cpglmm")}:  print the object }
    \item{summary}{\code{signature(object = "cpglmm")}: summary results}
    \item{predict}{\code{signature(object = "cpglmm")}: generate predictions for new data sets}
    \item{VarCorr}{\code{signature(x = "cpglmm")}: estimation for the variance components }
    \item{vcov}{\code{signature(object = "cpglmm")}: variance-covariance matrix for fixed effects }
   }
   
  The following methods are available for the class \code{"bcplm"}: 
  \describe{
    \item{plot}{\code{signature(x = "bcplm", y = "missing")}: summarize the \code{"bcplm"} object with a trace of the sampled output and a density estimate for each variable in the chain. See \code{\link[coda]{plot.mcmc}}. }
    \item{summary}{\code{signature(object = "bcplm")}:  produce two sets of summary statistics. See \code{\link[coda]{summary.mcmc}}. }
    \item{VarCorr}{\code{signature(x = "bcplm")}: estimation for the variance components if the random effects are present }
    \item{fixef}{\code{signature(object = "bcplm")}: extract fixed effects. Additional arguments include: \code{sd = FALSE}: extract standard errors; \code{quantiles = NULL}: compute empirical quantiles. These additional statistics are stored as attributes in the returned results.}
   }

The following methods are defined for the \code{"gini"} class:
  \describe{
    \item{plot}{\code{signature(x = "gini", y = "missing")}: plot the ordered Lorenz curve from each model comparison. If \code{overlay = TRUE} (the default), different curves are plotted on the same graph for each base. }
    \item{show}{\code{signature(object = "gini")}: print the computed Gini indices and standard errors.}
   } 
}

\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
  See also \code{\link{cpglm}}, \code{\link{cpglmm}}, \code{\link{bcplm}}, \code{\link[stats]{glm}}.   
}

\keyword{classes}
