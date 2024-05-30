#' N/C data merge function
#'
#' This function allows you import a .csv file and create a dataframe with the C and N data of the samples been processed.
#' The created dataframe will be merged with the spectral data during the fitting.
#' @param NCdata Raw spectrum
#' @keywords normalization, correction, flattening
#' @export


nc_data <- function (NCdata) {

  NC.end <- list()
  NC <- read.csv(NCdata, colClasses = "character")

  for(i in 1:nrow(NC)) {

    file.name <- NC[i,1]

    calc <- c((as.numeric(NC[i,3])/14.0067)/(as.numeric(NC[i,2])/12.0107))

    ## create a list with name, element, edge and the data of the spectrum
    NC.end[[i]] <- list("name" = file.name, "NC" =calc)
    ## close file loop
  }

  ## return resulting data frame
  return(NC.end)

  ## close function
}
