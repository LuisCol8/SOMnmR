#' Linear combination fitting solve function
#'
#' Quadratic programming solution function for linear combination fitting (LCF)
#' @param LCF.stds Standards for LCF
#' @param LCF.samp NMR Sample(s) for LCF
#' @param NMRmeth  Regions to be integrated, methods available include: "4region", "Bonanomi", "Smernik" and Molecular mixing model ("MMM").
#' @param FixNC TRUE or FALSE, for fixing or not the NC ratio on the sample fitting.
#' @returns A dataframe containing the result of the quadratic programming exercise, constrained or not by the Nc ratio (FixNC)
#' @keywords normalization correction
#' @importFrom quadprog solve.QP
#' @export

MMM_solve_QP <- function (LCF.stds, LCF.samp, NMRmeth = NULL, FixNC) {

  if (is.null(NMRmeth)) {

    stop("Please provide a method for Molecular mixing model constrain (NMRmet = MMM, FixNC = FALSE: no NC constrain, NMRmet = MMM, FixNC = TRUE: NC constrain)")

  } else if (FixNC == FALSE) {

    ## extract the names of the standards
    LC.standard.names <- colnames(LCF.stds)

    ## transform the raw absorption to numeric values
    X <- sapply(LCF.stds, as.numeric)
    Y <- sapply(LCF.samp, as.numeric)

    ## solve the Choleski factorization to find the covariance matrix to be minimized
    Rinv <- solve(chol(t(X) %*% X))

    ## create the vector to be minimized
    d <- t(Y) %*% X

    ## matrix with constrains under which to minimize the quadratic function (values between 0 and 1)
    C <- cbind(rep(1,length(LCF.stds)), diag(length(LCF.stds)))

    ## vector holding the values of sum to one constrain
    b <- c(1, rep(0, length(LCF.stds)))

    ## actual fit / solving the quadratic problem
    LCF.solve <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)

    ## extract the solution
    raw.coeff <- as.data.frame(t(LCF.solve$solution))
    colnames(raw.coeff) <- LC.standard.names

    ## create the fitted spectrum
    fit.spec <- rowSums(data.frame(mapply(`*`,LCF.stds,raw.coeff)))

    ## create the R-factor as fitting statistics
    r.fac <- sum((Y - fit.spec)^2)/sum(Y^2)

    ## create the Sum of Squares as fitting statistics
    ssq <- sum((Y - fit.spec)^2)

    ## combinde the coefficient with the R-factor
    result <- cbind(raw.coeff, R.fac = r.fac, SSQ = ssq)

    } else if (FixNC == TRUE) {

      ## extract the names of the standards
      LC.standard.names <- colnames(LCF.stds)

      ## transform the raw absorption to numeric values
      X <- sapply(LCF.stds, as.numeric)
      Y <- sapply(LCF.samp, as.numeric)

      ## solve the Choleski factorization to find the covariance matrix to be minimized
      Rinv <- solve(chol(t(X) %*% X))

      ## create the vector to be minimized
      d <- t(Y) %*% X

      ## matrix with constrains under which to minimize the quadratic function (values between 0 and 1)
      C <- cbind(rep(1,length(LCF.stds)), diag(length(LCF.stds)))

      ## vector holding the values of sum to one contrain
      ntot  <- as.numeric(X[1,1])
      nsamp  <- as.numeric(Y[1])
      n <- c(nsamp/ntot)
      #n  <- (Y[1]/0.275)
      b <- c(1, n, rep(0, length(LCF.stds)-1))

      ## actual fit / solving the quadratic problem
      LCF.solve <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 2)

      ## extract the solution
      raw.coeff <- as.data.frame(t(LCF.solve$solution))
      colnames(raw.coeff) <- LC.standard.names

      ## create the fitted spectrum
      fit.spec <- rowSums(data.frame(mapply(`*`,LCF.stds,raw.coeff)))

      ## create the R-factor as fitting statistics
      r.fac <- sum((Y - fit.spec)^2)/sum(Y^2)

      ## create the Sum of Squares as fitting statistics
      ssq <- sum((Y - fit.spec)^2)

      ## combinde the coefficient with the R-factor
      result <- cbind(raw.coeff, R.fac = r.fac, SSQ = ssq)
      #result <- cbind(raw.coeff, R.fac = r.fac)

  }

  ## return result
  return(result)
}
