#' Create a data frame of standard NMR areas
#'
#' The function creates a data frame with all standards of the selected ecosystem (Terrestrial or Aquatic).
#' @param NMRmeth Regions to be integrated.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @keywords standards
#' @export
#' @examples
#' see_NMR_table <- NMR_table(NMRmeth="4region")


NMR_table <- function (NMRmeth=NULL) {

  ## create dummy vectors for standard spectra and names
  if (is.null(NMRmeth)) {

    stop("Please choose an preset region model composition by typing 'MMM' for Molecular mixing model, 'Bonanomi' or '4region'")

  } else if (NMRmeth == "MMM") {

    ## Area composition
    int_table <- data.frame("Component" = c("Alkyl", "Methoxy/N-Alkyl", "O-Alkyl", "Di-O-Alkyl", "Aromatic",
                                          "Phenolic", "Amide/Carboxyl", "Ketone"),
                          "From" = c(-10.0, 45.0, 60.0, 95.0, 110.0, 145.0, 165.0, 190),
                          "To" = c(45.0, 60.0, 95.0, 110.0, 145.0, 165.0, 190.0, 220.0),
                          "ID" = c("Alkyl", "N-Alkyl/Methoxyl", "O-Alkyl", "Di-O-Alkyl", "Aromatic",
                                          "Phenolic","Amide/Carboxyl&Ketone","Amide/Carboxyl&Ketone"))

  } else if (NMRmeth == "Bonanomi") {

    ## Area composition
    int_table <- data.frame("Component" = c("Alkyl", "Methoxy/N-Alkyl", "O-Alkyl", "Di-O-Alkyl", "Aromatic",
                                            "Phenolic", "Amide/Carboxyl"),
                            "From" = c(0.0, 45.0, 60.0, 90.0, 110.0, 140.0, 160.0),
                            "To" = c(45.0, 60.0, 90.0, 110.0, 140.0, 160.0, 190.0))

  } else if (NMRmeth == "4region") {

    ## Area composition
    int_table <- data.frame("Component" = c("Alkyl", "O/N-Alkyl", "Aryl", "Carboxyl"),
                            "From" = c(0.0, 45.0, 110.0, 160.0),
                            "To" = c(45.0, 110, 160.0, 220.0))

  } else if (NMRmeth == "Smernik") {

  ## Area composition
  int_table <- data.frame("Component" = c("Alkyl", "O/N-Alkyl", "Aryl", "Carboxyl"),
                          "From" = c(0.0, 45.0, 110.0, 165.0),
                          "To" = c(45.0, 110, 165.0, 185.0))

}
  ## transform results to data frame and change its column names
  ## return resulting data frame
  return(int_table)
  ## close function
}
