% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Smernik200}
\alias{Smernik200}
\title{Smernik200  data set from Smernik et al. (2008)}
\format{
A nested list with 15 sub-lists:
\describe{
  \item{1 to 15}{A list containing the soil NMR spectrum of one of the following sites.
    \describe{
      \item{name}{"Control", "Burnt", "Burnt 1 year", "Control", "Control", "Control", "Control", "Burnt", "Burnt", "Burnt", "Burnt", "Burnt 1 year", "Burnt 1 year", "Burnt 1 year", "Burnt 1 year"}
      \item{raw.spec}{A data frame with 2 columns:
        \describe{
          \item{ppm}{Numeric vector.}
          \item{raw.intensity}{Numeric vector.}
        }
      }
    }
  }
}
}
\source{
Smernik et al., (2008) DOI: 10.1071/SR07128
}
\usage{
Smernik200
}
\description{
Contains 15 CP MAS 13C NMR spectra.
}
\details{
The spectra were taken in a NMR spectrometer with field strength of 200 MHz and MAS rate of 5 kHz
}
\examples{
data(Smernik200)
}
\keyword{datasets}
