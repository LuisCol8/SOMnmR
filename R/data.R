#' Hall sub data set from Hall et al. (2020)
#' 
#' Contains 17 CP MAS 13C NMR spectra.
#' 
#' The spectra were taken in a NMR spectrometer with field strength of 300 MHz and MAS rate of 12 kHz
#'
#'@format A nested list with 17 sub-lists:
#' \describe{
#'   \item{1 to 17}{A list containing the soil NMR spectrum of one of the following sites.
#'     \describe{
#'       \item{name}{"Calhoun", "CPER", "DCFS", "elve", "GRSM", "HARV", "icac", "JERC", "KONZ", "LENO", "MOAB", "NIWO", "ONAQ", "samt", "SCBI", "UNDE", "WOOD"}
#'       \item{raw.spec}{A data frame with 2 columns:
#'         \describe{
#'           \item{ppm}{Numeric vector.}
#'           \item{raw.intensity}{Numeric vector.}
#'         }
#'       }
#'     }
#'   }
#' }
#' 
#' @source https://portal.edirepository.org/nis/mapbrowse?packageid=edi.575.1
#' 
#' @examples
#' data(Hall300)
"Hall300"

#' Smernik200  data set from Smernik et al. (2008)
#' 
#' Contains 15 CP MAS 13C NMR spectra.
#' 
#' The spectra were taken in a NMR spectrometer with field strength of 200 MHz and MAS rate of 5 kHz
#'
#'@format A nested list with 15 sub-lists:
#' \describe{
#'   \item{1 to 15}{A list containing the soil NMR spectrum of one of the following sites.
#'     \describe{
#'       \item{name}{"Control", "Burnt", "Burnt 1 year", "Control", "Control", "Control", "Control", "Burnt", "Burnt", "Burnt", "Burnt", "Burnt 1 year", "Burnt 1 year", "Burnt 1 year", "Burnt 1 year"}
#'       \item{raw.spec}{A data frame with 2 columns:
#'         \describe{
#'           \item{ppm}{Numeric vector.}
#'           \item{raw.intensity}{Numeric vector.}
#'         }
#'       }
#'     }
#'   }
#' }
#' 
#' @source Smernik et al., (2008) DOI: 10.1071/SR07128
#' 
#' @examples
#' data(Smernik200)
"Smernik200"

#' Smernik400  data set from Smernik et al. (2008)
#' 
#' Contains 15 CP MAS 13C NMR spectra.
#' 
#' The spectra were taken in a NMR spectrometer with field strength of 400 MHz and MAS rate of 7 kHz
#'
#'@format A nested list with 15 sub-lists:
#' \describe{
#'   \item{1 to 15}{A list containing the soil NMR spectrum of one of the following sites.
#'     \describe{
#'       \item{name}{"Control", "Burnt", "Burnt 1 year", "Control", "Control", "Control", "Control", "Burnt", "Burnt", "Burnt", "Burnt", "Burnt 1 year", "Burnt 1 year", "Burnt 1 year", "Burnt 1 year"}
#'       \item{raw.spec}{A data frame with 2 columns:
#'         \describe{
#'           \item{ppm}{Numeric vector.}
#'           \item{raw.intensity}{Numeric vector.}
#'         }
#'       }
#'     }
#'   }
#' }
#' 
#' @source Smernik et al., (2008) DOI: 10.1071/SR07128
#' 
#' @examples
#' data(Smernik400)
"Smernik400"

#' GarciaF200  sub data set from Garcia-Franco et al. (2021)
#' 
#' Contains 3 CP MAS 13C NMR spectra.
#' 
#' The spectra were taken in a NMR spectrometer with field strength of 200 MHz and MAS rate of 6.8 kHz
#'
#'@format A nested list with 3 sub-lists:
#' \describe{
#'   \item{1 to 3}{A list containing the vegetation NMR spectrum of one of the following sites.
#'     \describe{
#'       \item{name}{"EB_Vegetation", "Fendt_Vegetation", "Graswang_Vegetation"}
#'       \item{raw.spec}{A data frame with 2 columns:
#'         \describe{
#'           \item{ppm}{Numeric vector.}
#'           \item{raw.intensity}{Numeric vector.}
#'         }
#'       }
#'     }
#'   }
#' }
#' 
#' @source Garcia-Franco et al. (2021) DOI: 10.1007/s00374-020-01518-0
#' 
#' @examples
#' data(GarciaF200)
"GarciaF200"
