#' Create a data frame of standard NMR areas
#'
#' The function creates a data frame with all standards of the selected ecosystem (Terrestrial or Aquatic).
#' @param ecosys Standards from the ecosystem to be fitted. "Terr" for terrestrial, "Aqua" for aquatic.
#' @keywords standards
#' @export
#' @examples


std_nmr <- function (ecosys=NULL) {

  ## create dummy vectors for standard spectra and names
  if (is.null(ecosys)) {

    stop("Please choose an ecosystem model composition by typic 'Terr' for Terrestrial or 'Aqua' for Aquatic")

  } else if (ecosys == "Terr") {

    ## Area composition
    MMM <- list(data.frame("Protein" = c(0.32,39.6,21.9,2.1,0.0,7.5,2.5,26.4), "Carbohydrates" = c(0.0,0.0,4.3,79.0,15.7,1.0,0.0,0.0),
                           "Lignin" = c(0.0,10.5,13.8,12.5,8.6,30.6,19.5,4.6), "Lipid" =c(0.0,75.6,4.5,9.0,0.0,3.6,0.7,6.6),
                           "Carbonyl" = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.00), "Char" = c(0.0,0.0,1.7,1.8,5.3,72.1,15.2,3.9)))


    #MMM <- list(data.frame("Protein" = c(0.275,35.4,22.6,3.5,0.0,8.9,1.3,28.3), "Carbohydrates" = c(0.0,0.0,0.0,83.3,16.7,0.0,0.0,0.0),
    #           "Lignin" = c(0.0,10.5,13.8,12.5,8.6,30.6,19.5,4.6), "Lipid" =c(0.0,75.6,4.5,9.0,0.0,3.6,0.7,6.6),
    #           "Carbonyl" = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.00), "Char" = c(0.0,0.0,1.7,1.8,5.3,72.1,15.2,3.9)))

  } else if (ecosys == "Aqua") {

    ## Area composition
    MMM <- list(data.frame("Protein" = c(0.266,36.6,24.7,2.9,0.0,4.5,1.0,30.4),"Carbohydrates" = c(0.0,0.0,0.0,83.3,16.7,0.0,0.0,0.0),
                "Lignin" = c(0.0,10.5,13.8,12.5,8.6,30.6,19.5,4.6), "Lipid" =c(0.0,94.4,0.0,0.0,0.0,0.0,0.0,5.6),
                "Carbonyl" = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.00), "Char" = c(0.0,0.0,1.7,1.8,5.3,72.1,15.2,3.9)))

  }
  ## transform results to data frame and change its column names
  ## return resulting data frame
  return(MMM)
  ## close function
}
