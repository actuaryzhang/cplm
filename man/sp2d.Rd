\name{sp2d}
\alias{sp2d}

\title{
2-dimentional Radial Spline
}
\description{
2-dimentional radial spline generator used in \code{cpglmm}
}
\usage{
sp2d(x1, x2, k = max(20, min(length(x1)/4, 150)), by = NULL, 
allPen = FALSE, varying = NULL, diag = FALSE, 
knots1 = quantile(x1, probs = 1:k/(k + 1)), 
knots2 = quantile(x1, probs = 1:k/(k + 1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
the first covariate in the coordinate 
}
  \item{x2}{
the second covariate in the coordinate
}
  \item{k}{
number of knots 
}
  \item{by}{
 not used. This is just for compatibility with \code{amer}. 
}
  \item{allPen}{
not used. This is just for compatibility with \code{amer}.
}
  \item{varying}{
not used. This is just for compatibility with \code{amer}.
}
  \item{diag}{
not used. This is just for compatibility with \code{amer}.
}
  \item{knots1}{
vector of knot locations for the first covariate in the coordinate
}
  \item{knots2}{
vector of knot locations for the second covariate in the coordinate
}
}

\author{
Fabian Scheipl \email{fabian.scheipl@googlemail.com}
}


