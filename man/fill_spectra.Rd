% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fill_spectra.R
\name{fill_spectra}
\alias{fill_spectra}
\title{Fill spectra (with 0s)}
\usage{
fill_spectra(chromdata, x = "mz")
}
\arguments{
\item{chromdata}{Tibble with columns \code{rt}, \code{mz}/\code{wl} and \code{intensity}, pre-grouped by peak (e.g. by ROI.)}

\item{x}{Name of spectrum x-coordinate (e.g. \code{mz} or \code{wl}.)}
}
\value{
The input tibble with 0-intensity observations added.
}
\description{
Take tibble of full-scan chromatographic data grouped by scan,
fill every missing x-coord (m/z or wl) in each scan with 0.
This is important prior to single-ion or -wavelength peak detection and integration.
}
\examples{

chromdata_filled <- chromdata \%>\%
 group_by(scan_rt) \%>\%
 fill_spectra
}
\keyword{fill}
\keyword{gaps}
\keyword{spectra}
\keyword{zero}
\keyword{zeroes}
