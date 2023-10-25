#' Porting for linear combination fitting
#'
#' The function can be used to check which combinations of standards produce a good fit.
#' @param all.samples List of all samples
#' @param all.standards List of all standards
#' @param LC.vals The fitting range values for the linear combination fitting
#' @param amoSTD Use at most X standards
#' @param ex.smaller Exclude portions smaller than a given value (decimal form), default to NULL
#' @param file.output Possibility to have a file output, default to NULL
#' @param best.fits Possibility to output more than the best fit (e.g. the first 10 best fits), default to 1
#' @keywords normalization correction
#' @export
#' @importFrom utils head combn setTxtProgressBar txtProgressBar write.csv2

fit_LCF <- function (all.samples, all.standards, ecosys =NULL,
                     amoSTD, ex.smaller = NULL,
                     file.output = NULL, best.fits = NULL,  NMRmeth, FixNC) {
  
  full.result <- NULL
  
  if(is.null(best.fits)) {
    best.fits <- 1
  }
  
  ## check if more than two standards have been set
  if(amoSTD < 2) stop("Use a minimum of two combinations of standards")
  if(amoSTD > length(all.standards[[1]])) stop(paste("Use a maximum of", length(all.standards[[1]]), "combinations of standards", sep = " "))
  
  ## check if file output is set, default is TRUE
  if(is.null(file.output)) {
    file.output <- FALSE
  }
  
  ## set exclude to FALSE if not set
  if (is.null(ex.smaller)) {
    
    ex.smaller <- 0
    
  } else {
    if (ex.smaller >= 1 | ex.smaller <= 0) stop("You can only exclude portions between 0 and 1, e.g. 0.02 for 2 %")
    
  }
  
  for (i in 1:length(all.samples)) {
    
    ## create correct data frame of standards, depending of the energy of the sample
    standards <- data.frame(all.standards[[1]][1:8,1:6])
    
    ## extract all names of the used standards
    STD.names <- colnames(standards)
    
    ## create list of all standard combinations
    STD.combs <- NULL
    
    ## loop of all standards and their combinations
    ## amoSTD sets the "use at most x Standards"
    for (m in 2:amoSTD) {
      
      ## recent combination of m standards
      temp.comb <- as.data.frame(combn(STD.names, m))
      colnames(temp.comb) <- NULL
      
      ## add recent combination to list of standard combinations
      #STD.combs <- c(STD.combs, as.list(temp.comb))
      
      if (FixNC == TRUE) {
        ## add recent combination to list of standard combinations
        STD.combs <- c(STD.combs, as.list(temp.comb))
        STD.combs <- STD.combs[grepl("Protein",STD.combs)]
        
      } else {
        ## add recent combination to list of standard combinations
        STD.combs <- c(STD.combs, as.list(temp.comb))
        
      }
      
    }
    
    ## create dummy list of new results
    new.results <- NULL
    
    ## paste which sample is used now
    print(paste("Sample: ", all.samples[[i]]$name, ", start date: ", Sys.time(), sep = ""))
    
    ## paste how much combinations will be fitted
    print(paste("Fitting ", length(STD.combs), " combinations", sep = ""))
    
    ## create progress bar for the standards combinations
    pb <- txtProgressBar(min = 1, max = length(STD.combs), style = 3)
    
    ## loop over all standards combinations
    for (j in 1:length(STD.combs)) {
      
      ## extract the standard combination and reduce the standards to only those that are used for fitting
      fit.standards <- standards[as.vector(STD.combs[[j]])]
      
      ## check names not in fit.standards
      res.STD <- STD.names[grep("FALSE", STD.names %in% colnames(fit.standards))]
      
      ## LC fit the sample with the standards choosen before
      result <- MMM_fit(sample = all.samples[[i]]$data$Integral, standards = fit.standards,  NMRmeth = NMRmeth, FixNC = FixNC)
      
      exp.zero <- as.data.frame(t(rep(0, length(res.STD))))
      colnames(exp.zero) <- res.STD
      result <- cbind(result, exp.zero)
      
      result <- result[c(names(standards), "R.fac", "SSQ")]
      new.results <- rbind(new.results, result)
      
      ## set progress bar again
      setTxtProgressBar(pb, j)
      
    }
    
    ## close progress bar
    close(pb)
    
    ## print finished time
    print(paste("Finished: ", Sys.time(), sep = ""))
    
    
    if (!is.null(new.results)) {
      
      new.results.sorted <- new.results[order(new.results$R.fac),]
      
    } else {
      
      
      new.results.sorted <- as.data.frame(t(rep(0, length(c(names(standards), "R.fac", "SSQ")))))
      colnames(new.results.sorted) <- c(names(standards), "R.fac", "SSQ")
      
    }
    
    if (file.output == TRUE) {
      write.csv2(x = new.results.sorted, file = paste("LCF.SOMnmR", all.samples[[i]]$name, "csv", sep = "."), row.names = FALSE)
    }
    
    best.fit <- head(new.results.sorted, best.fits)
    
    row.names(best.fit) <- NULL
    row.names(best.fit)[1] <- all.samples[[i]]$name
    
    full.result <- rbind(full.result, best.fit)
    
    write.csv2(full.result, "temp.result.SOMnmR.csv", row.names = TRUE)
    
    ## close all samples loop
  }
  
  return(full.result)
  
  ## close function
}

