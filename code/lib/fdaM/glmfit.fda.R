glmfit.fda  = function(x, y, weights = rep(1, nobs), start = NULL,
                       etastart = NULL, mustart = NULL, dataind = NULL,
                       offset = rep(0, nobs), family = gaussian(), 
                       control = list(), intercept = TRUE) 
{
    control = do.call("glm.control", control)
    x       = as.matrix(x)
    xnames  = dimnames(x)[[2L]]
    ynames  = if (is.matrix(y)) rownames(y)
               else              names(y)
    conv    = FALSE
    nobs    = NROW(y)
    nvars   = ncol(x)
    EMPTY   = nvars == 0
    if (is.null(weights)) weights = rep.int(1, nobs)
    if (is.null(offset))  offset  = rep.int(0, nobs)
    variance = family$variance
    linkinv  = family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
    dev.resids  = family$dev.resids
    aic         = family$aic
    mu.eta      = family$mu.eta
    unless.null = function(x, if.null) 
                     if (is.null(x)) if.null
                     else            x
    valideta    = unless.null(family$valideta, function(eta) TRUE)
    validmu     = unless.null(family$validmu,  function(mu)  TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep = mustart
        eval(family$initialize)
        mustart = mukeep
    }
    if (EMPTY) {
        #  initialization if no covariates
        eta       = rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                call. = FALSE)
        mu        = linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", call. = FALSE)
        dev       = sum(dev.resids(y, mu, weights))
        w         = ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals = (y - mu)/mu.eta(eta)
        good      = rep(TRUE, length(residuals))
        boundary  = conv = TRUE
        coef      = numeric(0L)
        iter      = 0L
    }
    else {
        #  --------------------------------------------------------------
        #           normal initialization for covariates present
        #  --------------------------------------------------------------
        coefold = NULL
        #  initial value of link function eta
        eta = 
            if      (!is.null(etastart)) etastart
            else if (!is.null(start)) 
                if (length(start) != nvars) 
                    errstr =paste("length of 'start' should equal %d",
                                  "and correspond to initial coefs for %s")
                    stop(gettextf(errstr, nvars, 
                                  paste(deparse(xnames), collapse = ", ")), 
                         domain = NA)
                else {
                    coefold = start
                    offset + as.vector(
                       if (NCOL(x) == 1) x * start
                       else              x %*% start)
                }
            else family$linkfun(mustart)
        eta[!dataind] = y[!dataind]
        #  initial value of expectation mu
        mu = linkinv(eta)
        mu[!dataind] = y[!dataind]
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some", 
                 call. = FALSE)
        #  initial previous deviance value
        devold   = sum(dev.resids(y, mu, weights))
        boundary = conv = FALSE
        #  --------------------------------------------------------------
        #                      iteration loop
        #  --------------------------------------------------------------
        for (iter in 1L:control$maxit) {
            #  indices for positive weights
            good  = weights > 0
            #  variance values
            varmu = variance(mu)[good]
            if (any(is.na(varmu))) stop("NAs in V(mu)")
            if (any(varmu == 0))   stop("0s in V(mu)")
            #  values of derivative of mu wrt eta
            mu.eta.val = mu.eta(eta)
            mu.eta.val[!dataind] = 1
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            #  indices both good and with nonzero mu.eta.val values
            good = (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv = FALSE
                warning("no observations informative at iteration ", 
                  iter)
                break
            }
            #  vector of predicted values
            z = (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            #  root-weights
            w = sqrt((weights[good] * mu.eta.val[good]^2) /
                      variance(mu)[good])
            w[!dataind] = 1
            ngoodobs = as.integer(nobs - sum(!good))
            #  compute LS fit
            fit = .Fortran("dqrls", qr = x[good, ] * w, n = ngoodobs, 
                p = nvars, y = w * z, ny = 1L, tol = min(1e-07, 
                  control$epsilon/1000), coefficients = double(nvars), 
                residuals = double(ngoodobs), effects = double(ngoodobs), 
                rank = integer(1L), pivot = 1L:nvars, qraux = double(nvars), 
                work = double(2 * nvars), PACKAGE = "base")
            #  check coefficients for infinite values
            if (any(!is.finite(fit$coefficients))) {
                conv = FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                  iter), domain = NA)
                break
            }
            #  check for under-identified model
            if (nobs < fit$rank) 
                errstr = "X matrix has rank %d, but only %d observations"
                stop(gettextf(errstr,fit$rank, nobs), domain = NA)
            start[fit$pivot] = fit$coefficients
            #  update eta values
            eta = drop(x %*% start)
            #  update mu values
            mu  = linkinv(eta = eta + offset)
            #  update deviance values
            dev = sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, "\n")
            boundary = FALSE
            #  control for infinite deviance value
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  errmsg = 
                      paste("no valid set of coefficients has been found:",
                            please supply starting values")
                  stop(errmsg, call. = FALSE)
                warning("step size truncated due to divergence", 
                        call. = FALSE)
                ii = 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; cannot correct step size", 
                      call. = FALSE)
                  ii = ii + 1
                  start = (start + coefold)/2
                  eta = drop(x %*% start)
                  mu = linkinv(eta = eta + offset)
                  dev = sum(dev.resids(y, mu, weights))
                }
                boundary = TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            #  control for invalid values of eta and mu
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii = 1
                #  control for invalid eta or invalid mu values
                while (!(valideta(eta) && validmu(mu))) {
                  #  check for too many step size corrections
                  if (ii > control$maxit) 
                    stop("inner loop 2; cannot correct step size", 
                         call. = FALSE)
                  ii    = ii + 1
                  start = (start + coefold)/2
                  eta   = drop(x %*% start)
                  mu    = linkinv(eta = eta + offset)
                }
                boundary = TRUE
                dev      = sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            #  convergence test
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv = TRUE
                coef = start
                break
            }
            else {
                devold = dev
                coef   = coefold = start
            }
        }
        #  --------------------------------------------------------------
        #       iteration loop completed:  wrap-up calculations
        #  --------------------------------------------------------------
        if (!conv) 
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary) 
            warning("glm.fit: algorithm stopped at boundary value", 
                    call. = FALSE)
        eps = 10 * .Machine$double.eps
        #  warning for binomial case of 0 or 1 mu values
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) {
                warnmsg = paste("glm.fit: fitted probabilities",
                                "0 or 1 occurred")
                warning(warnmsg, call. = FALSE)
            }
        }
        #  warning of poisson case of near zero values
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("glm.fit: fitted rates numerically 0 occurred", 
                  call. = FALSE)
        }
        #  control for singular design
        if (fit$rank < nvars) 
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] = NA
        xxnames   = xnames[fit$pivot]
        residuals = (y - mu)/mu.eta(eta)
        fit$qr    = as.matrix(fit$qr)
        nr        = min(sum(good), nvars)
        if (nr < nvars) {
            Rmat = diag(nvars)
            Rmat[1L:nr, 1L:nvars] = fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat = fit$qr[1L:nvars, 1L:nvars]
        Rmat = as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] = 0
        names(coef)      = xnames
        colnames(fit$qr) = xxnames
        dimnames(Rmat)   = list(xxnames, xxnames)
    }
    #  --------------------------------------------------------------
    #               assemble list of returned values
    #  --------------------------------------------------------------
    wt = rep.int(0, nobs)
    wt[good] = w^2
    names(residuals) = ynames
    names(mu)        = ynames
    names(eta)       = ynames
    names(wt)        = ynames
    names(weights)   = ynames
    names(y)         = ynames
    if (!EMPTY) 
        names(fit$effects) = c(xxnames[seq_len(fit$rank)], rep.int("", 
            sum(good) - fit$rank))
    wtdmu     = if (intercept) sum(weights * y)/sum(weights)
                 else           linkinv(offset)
    nulldev   = sum(dev.resids(y, wtdmu, weights))
    n.ok      = nobs - sum(weights == 0)
    nulldf    = n.ok - as.integer(intercept)
    rank      = if (EMPTY) 0
                 else       fit$rank
    resdf     = n.ok - rank
    aic.model = aic(y, n, mu, weights, dev) + 2 * rank
    #  set up list object to be returned
    list(coefficients = coef, 
         residuals = residuals, 
         effects   = if (!EMPTY) fit$effects, 
         R         = if (!EMPTY) Rmat, 
         rank      = rank, 
         qr        = if (!EMPTY) structure(
                         fit[c("qr", "rank", "qraux", "pivot", "tol")], 
                         class = "qr"), 
         family    = family, 
         deviance  = dev, 
         aic       = aic.model, 
         iter      = iter, 
         weights   = wt, 
         df.null   = nulldf, 
         y         = y, 
         converged = conv, 
         boundary  = boundary
         fitted.values     = mu, 
         linear.predictors = eta, 
         null.deviance     = nulldev, 
         prior.weights     = weights, 
         df.residual       = resdf, 
        )
}

___________________________________________________________________________

> glm.control
function (epsilon = 1e-08, maxit = 25, trace = FALSE) 
{
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}

