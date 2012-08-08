\name{gdor}
\title{Generic Directions of Recession}
\alias{gdor}
\description{
  \code{gdor} uses the methodology of generic directions of recession in the 
  event that the MLE does not exist in the usual sense, i.e. observed value of
  mean value parameter on the boundary, for binomial and poisson fits.
}
\usage{
gdor(formula, family, data, subset, na.action,
     start = NULL, etastart, mustart, offset,
     control = list(\dots), model = TRUE, method = "glm.fit",
     y = TRUE, contrasts = NULL, \dots)
}

\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

  \item{family}{a description of the error distribution and link
    function to be used in the model. This can be a character string
    naming a family function, a family function or the result of a call
    to a family function. This version of this package
    is limited to the binomial and poisson families.	
    (See \code{\link{family}} for details of family functions.)}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{glm} or \code{gdor} 
    is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{start}{starting values for the parameters in the linear predictor.}

  \item{etastart}{starting values for the linear predictor.}

  \item{mustart}{starting values for the vector of means.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{control}{a list of parameters for controlling the fitting
    process.  For \code{glm.fit} this is passed to
    \code{\link{glm.control}}.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{method}{the method to be used in fitting the model.  The default
    method \code{"glm.fit"} uses iteratively reweighted least squares
    (IWLS): the alternative \code{"model.frame"} returns the model frame
    and does no fitting.

    User-supplied fitting functions can be supplied either as a function
    or a character string naming a function, with a function which takes
    the same arguments as \code{glm.fit}.  If specified as a character
    string it is looked up from within the \pkg{stats} namespace.}

  \item{y}{For \code{glm} and \code{gdor}:
    logical values indicating whether the response vector used in the fitting 
    process should be returned as components of the returned value.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{\dots}{arguments to be used to form the default
    \code{control} argument if it is not supplied directly.}
}
\value{
  \code{gdor} returns an object of class inheriting from \code{"gdor"}
  which inherits from the class \code{"glm"} and \code{"lm"}. The object
  of class \code{"gdor"} is used for poisson and logistic regression
  applications when the MLE fails to exist in the usual convention.
  If a non-standard \code{method} is used, the object will also inherit
  from the class (if any) returned by that function.
  
  The function \code{\link{summary}} (i.e., \code{\link{summary.gdor}}) can
  be used to obtain or print a summary of the results and the function
  \code{\link{anova}} (i.e., \code{\link{anova.gdor}})
  to produce an analysis of variance table in the circumstance that the MLE
  does not exist.
}
\seealso{
  \code{\link{glm}}, \code{\link{anova.gdor}}, \code{\link{summary.gdor}}.
}
\references{
  Geyer, Charles J. (2009)
  \emph{Likelihood inference in exponential families and directions of
    recession.}
  Electronic Journal of Statistics 3, (2009), 259-289 (electronic)
}
\examples{
 # library(gdor)
 # x <- 1:30
 # y <- c(rep(0, 12), rep(1, 11), rep(0, 7)) 
 # m1 <- gdor(y ~ x + I(x^2), family=binomial)
 # summary(m1)

 # data <- read.table(url("http://www.stat.umn.edu/geyer/gdor/catrec.txt"),
 #		    header = TRUE)
 # m2 <- gdor(y ~ (v1 + v2 + v3 + v4 + v5 + v6 + v7)^3, family=poisson,
 #	       	data=data)
}
\keyword{models}
\keyword{regression}