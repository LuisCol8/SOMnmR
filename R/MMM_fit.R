#' All combination fitting of NMR spectra.
#'
#' The function wraps the Linear combination fit of the integrated regions of the molecular mixing model.
#' @param sample Sample Integrals
#' @param standards List of all standards
#' @param ex.smaller Exclude portions smaller than a given value (decimal form), default to NULL
#' @param NMRmeth Regions to be integrated.
#' @param FixNC TRUE or FALSE, for fixing or not the NC ratio on the sample fitting.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @keywords normalization correction
#' @export
#' @importFrom utils head combn setTxtProgressBar txtProgressBar write.csv2


MMM_fit <- function (sample, standards, ex.smaller = NULL, NMRmeth, FixNC) {

  ## set exclude to zero or stop if not set correctly
  if (is.null(ex.smaller)) {
    ex.smaller <- 0
  } else {
    if (ex.smaller >= 1 | ex.smaller < 0) stop("You can only exclude portions between 0 and 1, e.g. 0.02 for 2 %")
  }

  ## extract standards and sample in given range
  LC.sample <- data.frame(sample)
  #LC.sample <- rbind(LC.sample[1:7,1], sum(as.numeric(LC.sample[8:9,1])))
  #print(LC.sample)
  LC.standards <- standards

  ## solve the fitting as Quadratic Programming problem
  result <- MMM_solve_QP(LCF.stds = LC.standards, LCF.samp = LC.sample, NMRmeth = NMRmeth, FixNC = FixNC)

  ## extract the standard names
  LC.standard.names.start <- colnames(LC.standards)
  LC.standard.names <- colnames(LC.standards)

  ## check which coefficients are below exclution limit
  fit.vals <- which(result[LC.standard.names.start] < ex.smaller)

  ## loop to process fitting until no standards are excluded any more
  while (length(fit.vals) > 0) {

    ## subset the remaining standards and their names
    LC.standards <- LC.standards[-fit.vals]
    LC.standard.names <- colnames(LC.standards)

    ## solve the fitting as Quadratic Programming problem
    result <- MMM_solve_QP(LCF.stds = LC.standards,  LCF.samp = LC.sample, NMRmeth = NMRmeth, FixNC = FixNC)

    ## check for standards below smaller % value
    fit.vals <- which(result[LC.standard.names] < ex.smaller)

    ## close while loop
  }

  ## find the names of the now fitted standards
  fit.stds <- match(LC.standard.names, LC.standard.names.start)

  ## create a dummy coefficient vector and fill it with the new fitting values
  dum.coef <- as.data.frame(t(rep(0, length(LC.standard.names.start))))
  colnames(dum.coef) <- LC.standard.names.start
  dum.coef[fit.stds] <- as.numeric(result[LC.standard.names])


  ## fill a result vector with the proportions and the statistics
  end.result <- cbind(dum.coef, R.fac = result$R.fac, SSQ = result$SSQ)
  #end.result <- cbind(dum.coef, SSQ = result$SSQ)

  ## return result
  return(end.result)

  ## close function
}
