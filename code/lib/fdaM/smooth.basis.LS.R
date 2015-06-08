smooth.basis.LS = function(argvals, y, fdParobj, weight=NULL, 
                           fdnames=NULL, dfscale=1.0, covariates=NULL)
#SMOOTH.BASIS.LS  Smooths data by penalized least squares.  The smooth 
#  curves are expressed as a basis function expansions, and this function 
#  computes the coefficients of the expansions.  Smoothness is controlled
#  by controlling the size of an integrated squared linear differential
#  operator.  The integral is multiplied by a smoothing or bandwidth
#  parameter.
#
#  In addition, an optional covariate or design matrix can be supplied,
#  having a row corresponding to each value in ARGVALS and Q columns,
#  Q being the number of covariates.  See the optional COVARIATES
#  parameter described below.
#
#  This version of function smooth.basis sets up the smoothing problem
#  as a least squares fitting problem, with the penalty term set up as
#  the smooth of a basis system toward 0.
#
#  Required arguments for this function are:
#
#  ARGVALS   A set of argument values, set by default to equally spaced
#               on the unit interval (0,1).
#  Y         an array containing values of curves
#               If the array is a matrix, rows must correspond to argument
#               values and columns to replications, and it will be assumed
#               that there is only one variable per observation.
#               If Y is a three-dimensional array, the first dimension
#               corresponds to argument values, the second to replications,
#               and the third to variables within replications.
#               If Y is a vector, only one replicate and variable are 
#               assumed.
#  FDPAROBJ  A functional parameter or fdPar object.  This object 
#               contains the specifications forthe functional data
#               object to be estimated by smoothing the data.  See
#               comment lines in function fdPar fordetails.
#               This argument may also be either a FD object, or a 
#               BASIS object.  If this argument is a basis object, the 
#               smoothing parameter LAMBDA is set to 0.
#
#   A function call of the form
#             smooth.basis.LS(argvals, y, fdParobj,
#                             "PARAM1",val1,"PARAM2",val2,)
#   can be used to input a number of optional parameters.  The first
#   of each pair of arguments is the name of the parameter, supplied
#   as a string in quotes, and the second is the value of the argument,
#   supplied as an object of the class required by the parameter.
#  
#   These optional parameters are:
#
#     weight          vector of the same length as the data vector to be
#                     smoothed, containing nonnegative weights to be 
#                     applied to the data values
#     fdnames         A cell array of length 3 with names for
#                       1. argument domain, such as "Time"
#                       2. replications or cases
#                       3. the function.
#     covariates      A N by Q matrix Z of covariate values used to augment
#                     the smoothing function, where N is the number of
#                     data values to be smoothed and Q is the number of
#                     covariates.  The process of augmenting a smoothing 
#                     function in this way is often called "semi-parametric 
#                     regression".  The default is the empty object [].
#     dfscale         A scalar value multiplying the degrees of freedom
#                     in the definition of the generalized 
#                     cross-validated or GCV criterion forselecting the
#                     bandwidth parameter LAMBDA.  It was recommended by
#                     Chong Gu that this be a number slightly larger than
#                     1.0, such as 1.2, to prevent under-smoothing,
#                     The default is 1.0.
#
#  Returned objects are:
#
#  FDOBJ    an object of class fd containing coefficients.
#  BETA     the regression coefficients forthe covariates if supplied
#              or empty otherwise
#  DF       a degrees of freedom measure.
#  GCV      a measure of lack of fit discounted fordf.
#                 If the function is univariate, GCV is a vector 
#                 containing the stop  sum of squares foreach 
#                 function, and if the function is multivariate, 
#                 GCV is a NVAR by NCURVES matrix.
#  SSE      the stop sums of squares.  
#                 SSE is a vector or matrix of the same size as 
#                 GCV.
#  PENMAT   the penalty matrix, if computed, otherwise [].
#  Y2CMAP   the matrix mapping the data to the coefficients.
#  ARGVALS  the input set of argument values.
#  Y        the input array containing values of curves

#  Last modified 30 July 2010 by Jim Ramsay

if (nargin < 3) {
    stop("There is not at least three arguments.")
}

#  check ARGVALS

[argvals, n] = argcheck(argvals)           ????

#  check Y

[y, ncurve, nvar, ndim] = ycheck(y, n)     ????

#  check FDPAROBJ and get FDOBJ and LAMBDA

fdParobj = fdParcheck(fdParobj)
fdobj    = fdParobj$fd
lambda   = fdParobj$lambda
Lfdobj   = fdParobj$Lfd

#  check LAMBDA

if (lambda < 0) {
    lambda = 0
}

#  set up default fdnames

deffdnames[[1]] = "arguments"
deffdnames[[2]] = "replications"
deffdnames[[3]] = "variables"


#  get BASIS and NBASIS

basisobj = fdobj$basis
nbasis   = basisobj$nbasis - length(basisobj$dropind)

#  check WTVEC

[wtvec, onewt] = wtcheck(n, wtvec)   ????

#  check FDNAMES

if (!iherits(fdnames,"list")) {
    stop("smooth.basis.LS:fdnames", 
          "Optional argument FDNAMES is not a cell array.")
}

if (length(fdnames) != 3) {
    stop("smooth.basis.LS:fdnames", 
          "Optional argument FDNAMES is not of length 3.")
}

#  check COVARIATES

q = 0
if (!isnull(covariates)) {
    if (!isnumeric(covariates)) {
        stop("smooth.basis.LS:covariates", 
            "Optional argument COVARIATES is not numeric.")
    }
    if (dim(covariates)[1] != n) {
        stop("smooth.basis.LS:covariates", 
            "Optional argument COVARIATES has incorrect number of rows.")
    }
    q = dim(covariates)[2]
}

#  ------------------------------------------------------------------
#                set up the linear equations forsmoothing
#  ------------------------------------------------------------------

#  set up matrix of basis function values

basismat  = eval.basis(argvals, basisobj)

if (n >= nbasis || lambda > 0) {
    
    #  The following code is for the coefficients completely determined
    
    #  Multiply the basis matrix and the data pointwise by the square root 
    #  of the weight vector if the weight vector is not all ones.
    
    if (!onewt) {
        rtwtvec  = sqrt(wtvec)
        basismat = basismat * repmat(rtwtvec,1,nbasis)
        rtwtmat  = repmat(rtwtvec,1,ncurve)
        if (ndim < 3) {
            y = y * rtwtmat
        } else {
            for (ivar=1:nvar) {
                y[,,ivar] = y[,,ivar]*rtwtmat
            }
        }
    }
    
    #  set up additional rows of the least squares problem for the
    #  penalty term.

    basismat0 = basismat
    y0        = y
    
    if (lambda > 0) {
        nderiv  = Lfdobj$nderiv
        penmat  = eval.penalty(basisobj, Lfdobj)
        peneig = eigen(penmat)
        V = peneig$vectors
        D = peneig$values
        Dvec  = diag(D)
        #  Check that the lowest eigenvalue in the series that is to be
        #  kept is positive.
        eiglow = nbasis - nderiv
        if (D[eiglow] <= 0) {
            stop("smooth.basis:eig", 
                  paste("Eigenvalue(NBASIS-NDERIV) of penalty matrix ", 
                   "is not positive check penalty matrix."))
        }
        #  Check that the highest eigenvalue that is not used is small
        #  relative to the largest eigenvalue.
        if (nderiv > 0 && log10(D(eiglow+1)/D(1)) > -1e-12) {
            stop("smooth.basis:eig", 
                  ["Eigenvalue(NBASIS-NDERIV+1) of penalty matrix ", 
                   "is not small relative to eigenvalue(1) ", 
                   "check penalty matrix."])
        }
        #  Compute the square root of the penalty matrix in the subspace
        #  spanned by the first N - NDERIV eigenvectors
        ind = 1:eiglow
        penfac = Vsort(:,ind)%*%diag(sqrt(D[ind]))
        #  Augment basismat by sqrt(lambda)*t(penfac)
        basismat = rbind(basismat,sqrt(lambda)*t(penfac))
        #  Augment data vector by n - nderiv 0"s
        if (ndim < 3) {
            y = rbind(y, matrix(0,nbasis-nderiv,ncurve)
        } else {
            for (ivar=1:nvar) {
                y[,,ivar) = rbind(y[,,ivar), matrix(0,n-nderiv,ncurve))
            }
        }
    }
    
    #  augment BASISMAT0 and BASISMAT by the covariate matrix 
    #  if (it is supplied) {
    
    if (!isnull(covariates)) {
        ind1 = 1:n
        ind2 = (nbasis+1):(nbasis+q)
        # sparsewrd = issparse(basismat0)
        basismat0 = rbind(basismat0, matrix(0,dim(basismat0,1),q))
        basismat  = rbind(basismat,  matrix(0,dim(basismat, 1),q))
        if (!onewt) {
            basismat0[ind1,ind2] = covariates*repmat(rtwtvec,1,q)  ????
            basismat[ind1,ind2]  = covariates*repmat(rtwtvec,1,q)  ????
        } else {
            basismat0[ind1,ind2] = covariates
            basismat[ind1,ind2]  = covariates
        }
        # if (sparsewrd) {
        #     basismat0 = sparse(basismat0)
        #     basismat  = sparse(basismat)
        # }
    }
    
    #  solve the least squares problem using the QR decomposition with
    #  one iteration to improve accuracy
   
    R = qr(basismat,0)$R
    
    if (ndim < 3) {
        coef = solve(R,solve(t(R),crossprod(basismat,y)))
        res  = y - basismat%*%coef
        err  = solve(R,solve(t(R),crossprod(basismat,res)))
        coef = coef + err
    } else {
        coef = matrix(0,nbasis, ncurve, nvar)
        for (ivar = 1:nvar) {
            yi = y[,,i var]
            coefi = solve(R,solve(t(R),crossprod(basismat,yi)))
            resi  = yi - basismat %*% coefi
            erri  = solve(R,solve(t(R),crossprod(basismat,resi)))
            coefi = coefi + erri
            coef[,,ivar] = coefi
        }
    }
    
     #  compute basismat %*% R^{-1}
   
    MapFac = solve(t(R),t(basismat0))
    
    #  compute map from y to c
    
    y2cMap = solve(R,MapFac)
    
    #  compute degrees of freedom of smooth
    
    df = sum(diag(MapFac %*% t(MapFac)))
    
} else {
    stop(paste("The number of basis functions exceeds the number of ", 
               "points to be smoothed."))    
}

#  ------------------------------------------------------------------
#            compute SSE, yhat, GCV and other fit summaries
#  ------------------------------------------------------------------

#  compute stop sum of squares

if (ndim < 3) {
    yhat = basismat0 %*% coef
    SSE  = sum((y0 - yhat)^2)
} else {
    SSE = matrix(0,nvar,ncurve)
    for (ivar = 1:nvar) {
        coefi = coef[,,i var]
        yhati = basismat %*% coefi
        yi    = y[,,i var]
        SSE[ivar,] = sum((yi - yhati)^2)
    }
}

#  compute  GCV index

if (df < n) {
    gcv = (SSE/n)/((n - dfscale*df)/n)^2
} else {
    gcv = NULL
}

#  set up the functional data object

if (ndim < 3) {
    fdobj = fd(coef[1:nbasis,],  basisobj, fdnames)
} else {
    fdobj = fd(coef[1:nbasis,,], basisobj, fdnames)
}

#  set up the regression coefficient matrix beta

if (q > 0) {
    ind = (nbasis+1):(nbasis+q)
    if (ndim < 3) {
        beta = coef[ind,]
    } else {
        beta = coef[ind,,]
    }
} else {
    beta = []
}

return(
  list(fdobj, beta, df, gcv, SSE, penmat, y2cMap, argvals, y))  ????
