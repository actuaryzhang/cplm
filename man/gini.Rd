\name{gini}
\alias{gini}
\title{
The Gini index 
}
\description{
Compute Gini indices and their standard errors.
}
\usage{
gini(loss, score, base = NULL, data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loss}{a character that contains the name of the response variable. 
}
  \item{score}{
a character vector that contains the list of the scores, which are the predictions from the set of models to be compared.  
}
  \item{base}{a character that contains the name of a baseline statistic. If \code{NULL} (the default), each score will be successively used as the base. 
}
  \item{data}{a data frame containing the variables listed in the above arguments. 
}
  \item{\dots}{
not used. 
}
}
\details{
For model comparison involving the compound Poisson distribution, the usual mean squared loss function is not quite informative for capturing the differences between predictions and observations, due to the high proportions of zeros and the skewed heavy-tailed distribution of the positive losses. For this reason, Frees et al. (2011) develop an ordered version of the Lorenz curve and the associated Gini index as a statistical measure of the association between distributions, through which different predictive models can be compared. The idea is that a score (model) with a greater Gini index produces a greater separation among the observations. In the insurance context, a higher Gini index indicates greater ability to distinguish good risks from bad risks. Therefore, the model with the highest Gini index is preferred. 


This function computes the Gini indices and their asymptotic standard errors based on the ordered Lorenz curve. These metrics are mainly used for model comparison. Depending on the problem, there are generally two ways to do this. Take insurance predictive modeling as an example. First, when there is a baseline premium, we can compute the Gini index for each score (predictions from the model), and select the model with the highest Gini index.  Second, when there is no baseline premium (\code{base = NULL}), we successively specify the prediction from each model as the baseline premium and use the remaining models as the scores. This results in a matrix of Gini indices, and we select the model that is least vulnerable to alternative models using a "mini-max" argument - we select the score that provides the smallest of the maximal Gini indices, taken over competing scores.  
}

\value{
  \code{gini} returns an object of class \code{"gini"}. See \code{\link{gini-class}} for details of the return values as well as various methods available for this class.
}

\references{
\cite{Frees, E. W., Meyers, G. and Cummings, D. A. (2011). Summarizing Insurance Scores Using
a Gini Index. \emph{Journal of the American Statistical Association}, 495, 1085 - 1098.
}
}
\author{
Yanwei (Wayne) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{gini-class}}  for related information.
}
\examples{
\dontrun{

# Let's fit a series of models and compare them using the Gini index
da <- subset(AutoClaim, IN_YY == 1)
da <- transform(da, CLM_AMT = CLM_AMT / 1000)
                    
P1 <- cpglm(CLM_AMT ~ 1, data = da, offset = log(NPOLICY))


P2 <- cpglm(CLM_AMT ~ factor(CAR_USE) + factor(REVOLKED) + 
              factor(GENDER) + factor(AREA) + 
              factor(MARRIED) + factor(CAR_TYPE),
            data = da, offset = log(NPOLICY))

P3 <- cpglm(CLM_AMT ~ factor(CAR_USE) + factor(REVOLKED) + 
              factor(GENDER) + factor(AREA) + 
              factor(MARRIED) + factor(CAR_TYPE) +
              TRAVTIME + MVR_PTS + INCOME,
            data = da, offset = log(NPOLICY))

da <- transform(da, P1 = fitted(P1), P2 = fitted(P2), P3 = fitted(P3))
                
# compute the Gini indices
gg <- gini(loss = "CLM_AMT", score  = paste("P", 1:3, sep = ""), 
           data = da)
gg
           
# plot the Lorenz curves 
theme_set(theme_bw())
plot(gg)
plot(gg, overlay = FALSE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models } 
