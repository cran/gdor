\name{anova.gdor}
\alias{anova.gdor}
\alias{anova.gdorlist}
\title{Analysis of Deviance for Generalized Linear Model Fits}
\description{
  Compute an analysis of deviance table for one or more generalized
  linear model fits including the case when the MLE does not exist in 
  the coventional sense.
}
\usage{
\method{anova}{gdor}(object, \dots)

\method{anova}{gdorlist}(object, \dots)
}

\arguments{
  \item{object, \dots}{objects of class \code{gdor}, typically
    the result of a call to \code{\link{gdor}}.}
} 
\details{
  Specifying a single object gives a sequential analysis of deviance
  table for that fit.  That is, the reductions in the residual deviance
  as each term of the formula is added in turn are given in as
  the rows of a table, plus the residual deviances themselves.

  If more than one object is specified, the table has a row for the
  residual degrees of freedom and deviance for each model.
  For all but the first model, the change in degrees of freedom and
  deviance is also given. (This only makes statistical sense if the
  models are nested.)  It is conventional to list the models from
  smallest to largest, but this is up to the user.

  The table will contain test statistics (and P values) comparing 
  the reduction in deviance for the row to the residuals. The likelihood
  ratio test is used exclusively here in accordance with the theory when
  the MLE for the null model does not exist and the resulting object is
  of class \code{gdor}.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The comparison between two or more models by \code{anova} or
  \code{anova.gdorlist} will only be valid if they
  are fitted to the same dataset. This may be a problem if there are
  missing values and \R's default of \code{na.action = na.omit} is used,
  and \code{anova.gdorlist} will detect this with an error.
}
\seealso{
  \code{\link{gdor}}, \code{\link{anova}}.
}
\examples{
#  library(gdor)

#  data(catrec)

#  m1 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^2, family=poisson,
#		data=catrec)
#  m2 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^3, family=poisson,
#		data=catrec)
#  m3 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^4, family=poisson,
#		data=catrec)
#  anova(m1,m2)
#  anova(m2,m3)
}
\keyword{models}
\keyword{regression}