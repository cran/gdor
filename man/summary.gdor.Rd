\name{summary.gdor}
\alias{summary.gdor}
\alias{print.summary.gdor}
\title{Summarizing Poisson and Binomial Fits}
\usage{
\method{summary}{gdor}(object, \dots)

\method{print}{summary.gdor}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{object}{an object of class \code{"gdor"}, usually, a result of a
    call to \code{\link{gdor}}.}

  \item{x}{an object of class \code{"summary.gdor"}, usually, a result of a
    call to \code{summary.gdor}.}

  \item{digits}{the number of significant digits to use when printing.}

  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
   These functions are all \code{\link{methods}} for class \code{gdor} or
  \code{summary.gdor} objects.
}
\details{
  \code{print.summary.gdor} tries to be smart about formatting the
  coefficients, standard errors, etc. for the limiting conditional model of
  \code{gdor} objects as explained by the supporting theory. The 
  \code{coefficients} component of the result gives the estimated coefficients 
  and their estimated standard errors, together with their ratio. A fourth 
  column gives the two-tailed p-value corresponding to the z ratio based 
  on a Normal reference distribution.

  Aliased coefficients are omitted in the returned object but restored
  by the \code{print} method.

  Correlations are printed to two decimal places (or symbolically): to
  see the actual correlations print \code{summary(object)$correlation}
  directly.
}
\value{
  \code{summary.gdor} returns an object of class \code{"summary.gdor"}, a
  list with components

  \item{call}{the component from \code{object}.}
  \item{family}{the component from \code{object}.}
  \item{deviance}{the component from \code{object}.}
  \item{contrasts}{the component from \code{object}.}
  \item{df.residual}{the component from \code{object}.}
  \item{null.deviance}{the component from \code{object}.}
  \item{df.null}{the component from \code{object}.}
  \item{deviance.resid}{the deviance residuals:
    see \code{\link{residuals.glm}}.}
  \item{coefficients}{the matrix of coefficients, standard errors,
    z-values and p-values.  Aliased coefficients are omitted.}
  \item{df}{a 3-vector of the rank of the model and the number of
    residual degrees of freedom, plus number of non-aliased coefficients.}
}
\seealso{
  \code{\link{gdor}}, \code{\link{summary}}.
}
\examples{
#  library(gdor)
#  x <- 1:30
#  y <- c(rep(0, 12), rep(1, 11), rep(0, 7)) 
#  m1 <- gdor(y ~ x + I(x^2), family=binomial)
#  summary(m1)

#  data <- read.table(url("http://www.stat.umn.edu/geyer/gdor/catrec.txt"),
#		    header = TRUE)
#  m2 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^3, family=poisson,
#               data=data)
#  summary(m2)
}
\keyword{models}
\keyword{regression}