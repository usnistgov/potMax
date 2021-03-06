% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_threshold_fits.R
\name{gumbelMultiFit}
\alias{gumbelMultiFit}
\alias{gumbelMultiFit.declustered_series}
\alias{gumbelMultiFit.default}
\title{Maximum Likelihood Estimation of the Gumble Model for Many Thresholds}
\usage{
gumbelMultiFit(x, lt, n_min, n_max, weight_scale, progress_tf = TRUE)

\method{gumbelMultiFit}{declustered_series}(x, lt, n_min, n_max,
  weight_scale, progress_tf = TRUE)

\method{gumbelMultiFit}{default}(x, lt, n_min, n_max, weight_scale,
  progress_tf = TRUE)
}
\arguments{
\item{x}{An S3 object of class \code{declustered_series} or a numeric vector.
If the latter, the values to be thresholded and used in fitting.}

\item{lt}{(numeric scalar) The length of the time series in units
of time (seconds, minutes, hours, etc.).}

\item{n_min}{(numeric scalar) The minimum number of thresholded observations
to include}

\item{n_max}{(numeric scalar) The maximum number of thresholded observations
to include}

\item{weight_scale}{(numeric scalar) The value of \eqn{\tau}}

\item{progress_tf}{(logical scalar) Display a progress bar if TRUE, else not.}
}
\value{
An S3 object of class \code{gumbel_multi_fit} with elements

\describe{
  \item{\code{$all_fits}}{An object of type \code{gumbel_pot_fit} for each
  threshold}

  \item{\code{$thresholds}}{The thresholds for the fits}

  \item{\code{$weights}}{The weights associated with the fitted model for
  each threshold}

  \item{\code{$lt}}{The value of the \code{lt} argument}

  \item{\code{$n_min}}{The value of the \code{n_min} argument}

  \item{\code{$n_max}}{The value of the \code{n_max} argument}

  \item{\code{$weight_scale}}{The value of the \code{weight_scale} argument}
}
}
\description{
Fit the Gumbel like 2D extremal Poisson process for many
  thresholds
}
\details{
\code{gumbelMLE} and \code{gumbelWPlot} are called for a sequence of
  thresholds.  Weights associated with each fit are also calculated.  Suppose
  that for threshold \eqn{u_i} the maximum vertical distance from a point on
  the W plot to the \eqn{45^\circ} line is \eqn{\delta_i} such that the
  \eqn{\delta_i} are scaled to the unit interval.  The weight
  associated with threshold \eqn{u_i} is then

  \deqn{\frac{\exp\{-\tau\delta_i\}}{\sum\exp\{-\tau\delta_i\}}}
}
\section{Methods (by class)}{
\itemize{
\item \code{declustered_series}: 

\item \code{default}: 
}}

\examples{

\dontrun{

ddat <- decluster(-jp1tap813wind315$value)

multi_est <- gumbelMultiFit(x = ddat, lt = 100, n_min = 10, n_max = 50, weight_scale = 5)

}

}
