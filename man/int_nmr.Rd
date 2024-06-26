% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/int_nmr.R
\name{int_nmr}
\alias{int_nmr}
\title{Integration function}
\usage{
int_nmr(raw.spec, NMRmeth = NULL, NMR_field = NULL, NMR_rotation = NULL)
}
\arguments{
\item{raw.spec}{Raw spectrum}

\item{NMRmeth}{Regions to be integrated. Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM" or "MMM").}

\item{NMR_field}{Magnetic field of the NMR}

\item{NMR_rotation}{Rotation frequency of the sample probe in the NMR}
}
\value{
A nested list containing in the first level a string (name) and a list (data) which contains two data frames
one the raw spectrum and another the output of table of the integration with the spinning side bands.
}
\description{
This function allows you to integrate the 13C-NMR spectra using diferent integration regions.
The loaded Raw spectra can be integrated using the spinning side bands regions(default),
the Bonanomi("Bonanomi") regions or the Molecular Mixing Model regions("MMM").
The function returns the corrected, normalized and flattened spectrum
}
\examples{
data(GarciaF200)
Integralregions <- int_nmr(GarciaF200, NMRmeth = "4region", NMR_field = 200, NMR_rotation = 6800)
}
\keyword{integration}
\keyword{normalization}
