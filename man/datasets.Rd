\name{datasets}
\alias{FineRoot}
\alias{ClaimTriangle}
\alias{AutoClaim}
\docType{data}

\title{Data sets in the cplm pakcage}
\description{The data sets included in package is described here.}
  
\usage{
data(FineRoot)
data(ClaimTriangle)
data(AutoClaim)
}

\format{
\code{FineRoot}: a data set used for the study of the fine root length density of plants. It is a data frame with 511 records and 5 variables: 
  \describe{
    \item{\code{Plant}:}{identifier of the apple tree, 1-8}
    \item{\code{Stock}:}{root stokcing, one of three different root stocks: Mark, MM106 and M26}
    \item{\code{Spacing}:}{between-row \eqn{\times} within-row spacings, one of the following two: \eqn{4 \times 2} meters and \eqn{5 \times 3} meters}
    \item{\code{Zone}:}{inner or outer}
    \item{\code{RLD}:}{root length density}
    }

\code{ClaimTriangle}: a data set from an insurance loss reserving triangle. It is a data frame with 55 records and 3 variables: 
  \describe{
    \item{\code{year}:}{the year when the accident occurs}
    \item{\code{lag}:}{development lag}
    \item{\code{increLoss}:}{incremental insurance loss in 1000s}
    }

 
\code{AutoClaim}: a motor insurance data set retrieved from
the SAS Enterprise Miner database. It is a data frame with 10296 records and 29 variables: 
  \describe{
    \item{\code{POLICYNO}:}{"character", the policy number}
    \item{\code{PLCYDATE}:}{"Date", policy effective date}
    \item{\code{CLM_FREQ5}:}{"integer", the number of claims in the past 5 years}
    \item{\code{CLM_AMT5}:}{"integer", the total claim amount in the past 5 years}
    \item{\code{CLM_AMT}:}{"integer", the claim amount in the current insured period}
    \item{\code{KIDSDRIV}:}{"integer", the number of driving children}
    \item{\code{TRAVTIME}:}{"integer", the distance to work}
    \item{\code{CAR_USE}:}{"factor", the primary use of the vehicle: "Commercial", "Private".}
    \item{\code{BLUEBOOK}:}{"integer", the value of the vehicle}
    \item{\code{RETAINED}:}{"integer", the number of years as a customer}
    \item{\code{NPOLICY}:}{"integer", the number of policies}
    \item{\code{CAR_TYPE}:}{"factor", the type of the car: "Panel Truck", "Pickup", "Sedan", "Sports Car", "SUV", "Van".}
    \item{\code{RED_CAR}:}{"factor", whether the color of the car is red: "no", "yes".}
    \item{\code{REVOLKED}:}{"factor", whether the dirver's license was invoked in the past 7 years: "No", "Yes",}
    \item{\code{MVR_PTS}:}{"integer", MVR violation records}
    \item{\code{CLM_FLAG}:}{"factor", whether a claim is reported: "No", "Yes".}
    \item{\code{AGE}:}{"integer", the age of the driver}
    \item{\code{HOMEKIDS}:}{"integer", the number of children}
    \item{\code{YOJ}:}{"integer", years at current job}
    \item{\code{INCOME}:}{"integer", annual income}
    \item{\code{GENDER}:}{"factor", the gender of the driver: "F", "M".}
    \item{\code{MARRIED}:}{"factor", married or not: "No", "Yes".}
    \item{\code{PARENT1}:}{"factor", single parent: "No", "Yes".}
    \item{\code{JOBCLASS}:}{"factor": "Unknown", "Blue Collar", "Clerical", "Doctor", "Home Maker", "Lawyer", "Manager", "Professional", "Student".}
    \item{\code{MAX_EDUC}:}{"factor", max education level:"<High School", "Bachelors", "High School", "Masters", "PhD".}
    \item{\code{HOME_VAL}:}{"integer", the value of the insured's home}
    \item{\code{SAMEHOME}:}{"integer", years in the current address}
    \item{\code{DENSITY}:}{"factor", home/work area: "Highly Rural", "Highly Urban", "Rural", "Urban".}
    \item{\code{IN_YY}:}{"logical", whether the record is used in the Yip and Yau (2005) paper.}
    }
}

\source{
	\cite{de Silva, H. N., Hall, A. J., Tustin, D. S. and Gandar, P. W. (1999). Analysis of distribution
of root length density of apple trees on different dwarfing rootstocks. \emph{Annals of 
Botany}, 83: 335-345.}

\cite{ Dunn, P.K. and Smyth, G.K. (2005). Series evaluation of Tweedie exponential dispersionmodels densities. \emph{Statistics and Computing}, 15, 267-280.}

  \cite{Peters G. W., Shevchenko P. V. and Wuthrich M. V. (2009). Model Uncertainty in Claims Reserving within Tweedie's Compound Poisson Models. \emph{Astin Bulletin}, 39(1),  1-33.}

\cite{Yip, K. C. H. and Yau, K. K. W. (2005). On Modeling Claim Frequency Data In General
Insurance With Extra Zeros. \emph{Insurance: Mathematics and Economics}, 36(2), 153-163.}
}


\keyword{datasets}

