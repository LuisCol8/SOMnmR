#' Integration function
#'
#' This function allows you to integrate the 13C-NMR spectra using diferent integration regions.
#' The loaded Raw spectra can be integrated using the spinning side bands regions(default),
#' the Bonanomi("Bonanomi") regions or the Molecular Mixing Model regions("MMM").
#' The function returns the corrected, normalized and flattened spectrum
#' @param raw.spec Raw spectrum
#' @param NMRmeth Regions to be integrated.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM" or "MMM-SSB").
#' @param NMR_field Magnetic field of the NMR
#' @param NMR_rotation Rotation frequency of the sample probe in the NMR
#' @keywords normalization, integration
#' @export
#' @importFrom cmna simp
#' @examples

int_nmr <- function(raw.spec, NMRmeth=NULL, NMR_field=NULL, NMR_rotation=NULL) {
  
  raw.spec.end <- NULL
  
  if (is.null(NMRmeth)) {
    
    
    stop("Please choose an preset region model composition by typing 'MMM' for Molecular mixing model, 'Bonanomi' or '4region'")
    
  } else if (is.null(NMR_field)){
    
    stop("Please add the NMR Magnetic field")
    
  } else if (is.null(NMR_rotation)){
    
    stop("Please add the NMR roation frequency")
    
  } else if (!is.null(NMRmeth)) {
    
    raw.spec.end <- NULL
    
    int_table <- ssb_ofset(NMRmeth=NMRmeth, NMR_field=NMR_field, NMR_rotation=NMR_rotation)

    for (i in 1:length(raw.spec)) {
      Integral <- NULL
      name <- raw.spec[[i]]$name
      raw.spec.end[[i]] <- raw.spec[[i]]
      spectrum <- raw.spec[[i]]$data$raw.spec
      
      ## Extract ppm (x) and the intesity (y) for predifined intervals
      for (j in 1:nrow(int_table)){
        
        Int.min.j <- which(abs(spectrum[[c("ppm")]]-(int_table$From[j])) == min(abs(spectrum[[c("ppm")]]-(int_table$From[j]))))
        #print(Int.min.j)
        Int.max.j <- which(abs(spectrum[[c("ppm")]]-(int_table$To[j])) == min(abs(spectrum[[c("ppm")]]-(int_table$To[j]))))
        #print(Int.max.j)
        Int.x.j <- c(spectrum[[c("ppm")]][(Int.min.j:Int.max.j)])
        #print(Int.x.j)
        Int.y.j <- c(spectrum[[c("raw.intensity")]][(Int.min.j:Int.max.j)])
        Integral <- append(Integral,trapz(Int.x.j,Int.y.j))
        #print(length((Int.x.1)))
      }
      
      norm <- sum(Integral)
      normalized.Int <- (Integral/norm)*100
      f.integral <- data.frame(normalized.Int)
      Integral <-setNames(cbind(f.integral,int_table),c("Integral","From", "To", "Component", "Component_index", "Component_ssb", "sbb_index", "ssb_ofset"))
      raw.spec.end[[i]] <- list("name" = name, "data" = list("raw.spec" = spectrum,"Integral" = Integral))
    }
  }
  return(raw.spec.end)
}

