% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region_calc.R
\name{region_calc}
\alias{region_calc}
\title{Functional groups calculation}
\usage{
region_calc(
  batch_nmr = NULL,
  file = NULL,
  NMRmeth = NULL,
  ecosys = NULL,
  cn_constrain = FALSE,
  cndata = NULL
)
}
\arguments{
\item{batch_nmr}{Vector with file names, default}

\item{NMRmeth}{Regions to be integrated.
Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").}

\item{ecosys}{Standards to be used for the MMM, can be Terrestrial("Terr") or Aquatic ("Aqua")}

\item{raw.spec}{List of files already loaded with read_raw_spec()}

\item{nc_constrain}{Boolean for the fitting of the MMM, TRUE for N/C constrained fit, FALS for non-constrained. Default is FALSE.}

\item{ncdata}{Initial correction and normalization parameters}
}
\description{
This function loads, integates and calculates the functional group distribution from the raw spectra.
Output is a list with the raw data, integrals and corrected spectra.
}
