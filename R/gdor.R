gdor <- function(formula, family, data, subset, na.action, start = NULL,
		 etastart, mustart, offset, control = list(...), model = TRUE,
		 method = "glm.fit", y = TRUE, contrasts = NULL, ...)
{
  call <- match.call()
  ## family
  if(is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
    
  if(identical(method, "model.frame")) return(mf)

  if (!is.character(method) && !is.function(method))
      stop("invalid 'method' argument")

  ## for back-compatibility in return result
  if (identical(method, "glm.fit"))
      control <- do.call("glm.control", control)

  mt <- attr(mf, "terms") # allow model.frame to have updated it
  minY <- Y <- model.response(mf, "any") # e.g. factors are allowed
  minY <- as.matrix(minY)
  dim <- dim(minY)[2]

  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if(!is.null(nm)) names(Y) <- nm
  }
  if(dim == 2){
    sum <- sum(Y[,1]) - sum(Y[,2])
    if(sum < 0) minY <- Y[,1]
    else minY <- Y[,2]
  }
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)

  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
      if(length(offset) != NROW(Y))
          stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
	  length(offset), NROW(Y)), domain = NA)
  }

  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")

  ## extract x, y, etc from the model formula and frame in preperation for a
  ## check of the MLE
  mf2 <- match.call(expand.dots = FALSE)
  m2 <- match(c("formula", "family", "data", "subset", 
		"na.action", "etastart", "mustart", "offset"), 
		names(mf2), 0L)
  mf2 <- mf2[c(1L, m2)] 
  mf2[[1L]] <- as.name("glm")
  x <- TRUE
  mf2$x <- TRUE
  flag <- flag2 <- FALSE
  gDOR <- out <- subset <- nsing <- lCM <- NULL
  
  if(family$family == "binomial"){
    out <- suppressWarnings(eval(mf2, parent.frame()))
    M <- out$x
    len <- length(M[2,])
    tanv <- M
    tanv[minY == 0, ] <- (-tanv[minY == 0, ])
    vrep <- cbind(0, 0, tanv)
    if(dim == 2){
      linear <- Y[,1] > 0 & Y[,2] > 0
      vrep[linear, 1] <- 1
      lout <- linearity(vrep, representation = "V")
      linear[lout] <- TRUE
      p <- ncol(tanv)
      hrep <- cbind(0, 0, -tanv, 0)
      hrep[!linear, ncol(hrep)] <- (-1)
      hrep[linear, 1] <- 1
      hrep <- rbind(hrep, c(0, 1, rep(0, p), -1))
      objv <- c(rep(0, p), 1)
      pout <- lpcdd(hrep, objv, minimize = FALSE)
      gDOR <- pout$primal.solution[1:p]
      ## check validity of gdor
      if(all(tanv %*% gDOR <= 0) == TRUE && (tanv %*% gDOR < 0) == TRUE){
	flag <- TRUE
	data.cond <- as.matrix(M[linear,])
	Y.cond <- Y[linear,]
	lCM <- eval(glm(Y.cond ~ data.cond + 0, family = binomial))
	nsing <- length(Y[,1]) - length(Y.cond[,1])
      }
    }
    else{
      lout <- linearity(vrep, representation = "V")
      hrep <- cbind(-vrep, -1)
      hrep <- rbind(hrep, c(0, 1, rep(0, len), -1))
      objv <- c(rep(0, len), 1)
      pout <- lpcdd(hrep, objv, minimize = FALSE)
      gDOR <- pout$primal.solution[-length(pout$primal.solution)]
      ## check validity of gdor
	if(all(tanv %*% gDOR < 0) == TRUE){
	flag <- TRUE 
	data.cond <- M[((M %*% gDOR) == 0),]
	mf2$data <- data.cond
	subset <- mf2$data
	if(length(M[((M %*% gDOR) == 0),]) == 0){
	  flag2 <- TRUE
	}
	else{
	  lCM <- eval(mf2, parent.frame())
	  nsing <- length(Y) - length(Y.cond)
	}
      }
    }
  }

  else if(family$family == "poisson"){
    out <- suppressWarnings(eval(mf2, parent.frame()))
    M <- out$x
    tanv <- out$x
    vrep <- cbind(0, 0, tanv)
    vrep[data$y > 0, 1] <- 1
    lout <- linearity(d2q(vrep), representation = "V")
    linear <- data$y > 0
    linear[lout] <- TRUE
    ## check validity of gdor
    if(!all(linear) == TRUE){
      flag <- TRUE
      hrep <- cbind(0, 0, -tanv, 0)
      hrep[!linear, ncol(hrep)] <- (-1)
      hrep[linear, 1] <- 1
      hrep <- rbind(hrep, c(0, 1, rep(0, ncol(out$x)),-1))
      objv <- c(rep(0, ncol(out$x)), 1)
      pout <- lpcdd(d2q(hrep), d2q(objv), minimize = FALSE)
      gDOR <- pout$primal.solution[-length(pout$primal.solution)]
      subset <- as.data.frame(data[linear, ])
      nsing <- (sum(data$y == 0) - sum(!linear))
      mf2$data <- subset
      if(length(subset) == 0){
	flag2 <- TRUE
      }
      else lCM <- eval(mf2, parent.frame())
    }
  }
  else stop(" 'family' not recognized")


    fit <- suppressWarnings(eval(call(if(is.function(method)) "method" else method,
		     x = X, y = Y, start = start,
		     etastart = etastart, mustart = mustart,
		     offset = offset, family = family, control = control,
		     intercept = attr(mt, "intercept") > 0L)))

    ## This calculated the null deviance from the intercept-only model
    ## if there is one, otherwise from the offset-only model.
    ## We need to recalculate by a proper fit if there is intercept and
    ## offset.
    ##
    ## The glm.fit calculation could be wrong if the link depends on the
    ## observations, so we allow the null deviance to be forced to be
    ## re-calculated by setting an offset (provided there is an intercept).
    ## Prior to 2.4.0 this was only done for non-zero offsets.
    if(length(offset) && attr(mt, "intercept") > 0L) {
        fit$null.deviance <-
            eval(call(if(is.function(method)) "method" else method,
                      x = X[, "(Intercept)", drop=FALSE], y = Y,
                      offset = offset, family = family,
                      control = control, intercept = TRUE))$deviance
    }
    if(model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if(x) fit$x <- X
    if(!y) fit$y <- NULL

    fit <- list(call = call, method = method, formula = formula, 
		family = family, flag = flag, flag2 = flag2, 
		lcm = lCM, model = out, x = x, y = y, terms = mt, 
		data = data, control = control,
		xlevels = .getXlevels(mt, mf), subset = subset,
		offset = offset, gDOR = gDOR, nsing = nsing,
		contrasts = attr(X, "contrasts"))
    class(fit) <- c(fit$class, "gdor")
    fit
}

anova.gdor <- function(object, ...)
{    
    ## check for multiple objects
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
	rep(FALSE, length(dotargs)) else (names(dotargs) != "")
    if(any(named))
	warning("the following arguments to 'anova.gdor' are invalid and dropped: ",
		paste(deparse(dotargs[named]), collapse=", "))
    dotargs <- dotargs[!named]
    is.gdor <- unlist(lapply(dotargs,function(x) inherits(x,"gdor")))
    dotargs <- dotargs[is.gdor]
    if (length(dotargs))
	return(anova.gdorlist(c(list(object), dotargs)))    
    object <- object$model

    ## extract variables from model
    varlist <- attr(object$terms, "variables")

    ## must avoid partial matching here.
    x <-
	if (n <- match("x", names(object), 0L))
	    object[[n]]
	else model.matrix(object)
    varseq <- attr(x, "assign")
    nvars <- max(0, varseq)
    resdev <- resdf <- NULL

    ## if there is more than one explanatory variable then
    ## recall glm.fit to fit variables sequentially

    ## for score tests, we need to do so in any case
    if(nvars > 1) {
	method <- object$method
    
        ## allow for 'y = FALSE' in the call (PR#13098)
        y <- object$y
        if(is.null(y)) { ## code from residuals.glm
            mu.eta <- object$family$mu.eta
            eta <- object$linear.predictors
            y <-   object$fitted.values + object$residuals * mu.eta(eta)
        }
	for(i in seq_len(nvars-1L)) {
	    ## explanatory variables up to i are kept in the model
	    ## use method from glm to find residual deviance
	    ## and df for each sequential fit
	    fit <- eval(call(if(is.function(method)) "method" else method,
                             x=x[, varseq <= i, drop = FALSE],
                             y=y,
                             weights=object$prior.weights,
                             start  =object$start,
                             offset =object$offset,
                             family =object$family,
                             control=object$control))
    }
  }

    ## add values from null and full model

    resdf <- c(object$df.null, resdf, object$df.residual)
    resdev <- c(object$null.deviance, resdev, object$deviance)

    ## construct table and title

    table <- data.frame(c(NA, -diff(resdf)),
			c(NA, pmax(0, -diff(resdev))), resdf, resdev)
    tl <- attr(object$terms, "term.labels")
    if (length(tl) == 0L) table <- table[1,,drop=FALSE] # kludge for null model
    dimnames(table) <- list(c("NULL", tl),
			    c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
    title <- paste("Analysis of Deviance Table", "\n\nModel: ",
		   object$family$family, ", link: ", object$family$link,
		   "\n\nResponse: ", as.character(varlist[-1L])[1L],
		   "\n\nTerms added sequentially (first to last)\n\n", sep="")

    ## calculate test statistics if needed

    df.dispersion <- Inf
    if(is.null(dispersion)) {
	dispersion <- summary(object, dispersion=dispersion)$dispersion
	df.dispersion <- if (dispersion == 1) Inf else object$df.residual
    }
    structure(table, heading = title, class= c("anova", "data.frame"))
}

anova.gdorlist <- function(object, ...)
{
  ## find responses for all models and remove
  ## any models with a different response

  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2L]])} ))
      sameresp <- responses==responses[1L]
      if(!all(sameresp)) {
	object <- object[sameresp]
	warning("models with response ", deparse(responses[!sameresp]),
		" removed because response differs from model 1")
      }

  ns <- sapply(object, function(x) length(x$residuals))
  if(any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")

  ## calculate the number of models

  nmodels <- length(object)
  if(nmodels==1)
    return(anova.gdor(object[[1L]]))

  ## find existence of gdor in null model

  flag <- as.logical(lapply(object, function(x) x$flag))
  flag <- as.logical(length(flag) == length(flag[flag == TRUE]))
  resdev <- resdf <- NULL
    
  ## extract statistics

  if(flag == FALSE){
    resdf  <- as.numeric(lapply(object, function(x) x$model$df.residual))
    resdev <- as.numeric(lapply(object, function(x) x$model$deviance))
  }
    
  else if(flag == TRUE){
    lCM <- lapply(object, function(x) x$lcm)
    M <- lapply(lCM, function(x) x$x)
    df <- NULL
      
    for(i in 1:length(M)){
      suppressWarnings(df[i] <- redundant(M[[i]]))
      resdf[i] <- as.numeric(length(df[[i]][,1]))
    }
    resdev <- as.numeric(lapply(lCM, function(x) x$deviance))
  }

  ## construct table and title

  table <- data.frame(resdf, resdev, c(NA, abs(diff(resdf))),
		  c(NA, -diff(resdev)) )
  variables <- lapply(object, function(x)
		  paste(deparse(formula(x)), collapse="\n") )
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", "Df",
				  "Deviance"))

    title <- "Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1L:nmodels),": ",
		  variables, sep="", collapse="\n")

  ## calculate test statistic

  bigmodel <- object[[order(resdf)[1L]]]
  dispersion <- summary(bigmodel, dispersion=dispersion)$dispersion
  df.dispersion <- if (dispersion == 1) Inf else min(resdf)

  table <- stat.anova(table = table, test = "LRT",
		      scale = dispersion, df.scale = df.dispersion,
		      n = length(bigmodel$residuals))
 
  structure(table, heading = c(title, topnote),
		 class = c("anova", "data.frame"))
}

summary.gdor <- function(object, ...)
{
  keep <- ans <- NULL
  dispersion <-
	  if(object$family$family %in% c("poisson", "binomial")) 1
	  else if(df.r > 0) {
	    est.disp <- TRUE
	    if(any(object$weights==0))
	      warning("observations with zero weight not used for calculating dispersion")
	      sum((object$weights*object$residuals^2)[object$weights > 0])/ df.r
	  } 
	  else {
	    est.disp <- TRUE
	    NaN
	  }
  flag <- object$flag
  flag2 <- object$flag2
  if(flag == FALSE){
    print(eval(summary(object$model)))
  }
  else{
    gDOR <- object$gDOR
    nsing <- object$nsing
    if(class(object$data) != "environment") data <- as.data.frame(object$data)
    if(flag2 == FALSE){
      lCM <- object$lcm
      coef.table <- as.data.frame(eval(summary(lCM)$coefficients))
      est.disp <- FALSE
      df.r <- lCM$df.residual

      ## these need not all exist, e.g. na.action.
      keep <- match(c("call","terms","family","deviance", "aic",
		    "contrasts", "df.residual","null.deviance","df.null",
		    "iter", "na.action"), names(lCM), 0L)
      ans <- c(lCM[keep],
	      list(deviance.resid = residuals(lCM, type = "deviance"),
		  coef = coef.table, flag2 = flag2,
		  #aliased = aliased,
		  gDOR = gDOR, dispersion = dispersion,
		  nsing = nsing))
    }
    else{
      keep <- match(c("call","terms","family","deviance", "aic",
		    "contrasts", "df.residual","null.deviance","df.null",
		    "iter", "na.action"), names(object$model), 0L)
      ans <- c(object$model[keep],
	      list(flag2 = flag2, gDOR = gDOR, nsing = nsing))
    #aliased <- is.na(coef(lCM))
    }
    class(ans) <- "summary.gdor"
    return(ans)
  }
}

print.summary.gdor <- function(x, digits = max(3, getOption("digits") - 3),
			       ...)
{
  x$call$data <- "data.cond"
  cat("\nCall:\n", paste(deparse(x$call),sep = "\n",collapse = "\n"), "\n\n", sep=" ")
  cat(" MLE failed to exist, solution is 'at infinity'. A limiting conditional model\n",
      "had to be implemented to obtain relevant output. MLE exists 'at infinity'\n", 
      "in direction: \n")
  print.default(x$gDOR, digits=digits, na.print=" ", pint.gap = 2)

  if(x$flag2 == TRUE)
    cat("\n In addition, limiting conditional model is completely degenerate.\n")

  else{
  cat("\nDeviance Residuals: \n")
  if(x$df.residual > 5) {
    x$deviance.resid <- quantile(x$deviance.resid,na.rm=TRUE)
    names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
  }
  xx <- zapsmall(x$deviance.resid, digits + 1)
  print.default(xx, digits=digits, na.print = "", print.gap = 2)

  cat("\nCoefficients: (", x$nsing,
      " not defined because of singularities)\n", sep = "")
  printCoefmat(x$coef, digits=digits, ...)

  cat("\n(Dispersion parameter for ", x$family$family,
      " family taken to be ", format(x$dispersion), ")\n\n",
      apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                        "deviance:"),
		  format(unlist(x[c("null.deviance","deviance")]),
			 digits= max(5, digits+1)), " on",
		  format(unlist(x[c("df.null","df.residual")])),
		  " degrees of freedom\n"),
	    1L, paste, collapse=" "), sep="")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    cat("AIC: ", format(x$aic, digits= max(4, digits+1)),"\n\n",
	"Number of Fisher Scoring iterations: ", x$iter,
	"\n", sep="")
  }
}