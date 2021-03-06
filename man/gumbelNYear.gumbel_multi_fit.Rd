% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n_year_vals_multi_thresh_fits.R
\name{gumbelNYear.gumbel_multi_fit}
\alias{gumbelNYear.gumbel_multi_fit}
\title{Return Values for the Gumbel Model Using Many Thresholds}
\usage{
\method{gumbelNYear}{gumbel_multi_fit}(x, N)
}
\arguments{
\item{x}{An S3 object of type \code{gumbel_multi_fit}}

\item{N}{(numeric scalar) The N in N-year return value.  This is a bit of a
misnomer since the unit of time does not have to be years.  The function
can calculate N-second, N-minute, N-hour, etc. return values as well.  In
fact, the unit of time is the same unit of time passed in for the \code{lt}
argument of \code{gumbelMultiFit}.  Naming the function for the unit of
time year is simply due to my past experince with calculating return values
on the time scale of years.}
}
\value{
An S3 object of class \code{gumbel_N_year_val_multi_thresh} with
  elements


\describe{

  \item{\code{$mu}}{A numeric vector containing the estimated location
  parameters, one for each threshold considered}

  \item{\code{$sigma}}{A numeric vector containing the estimated scale
  parameters, one for each threshold considered}

  \item{\code{$thresh}}{A numeric vector containing the thresholds being
  considered}

  \item{\code{$N}}{The value inputted for \code{N}}

  \item{\code{$N_year_val}}{(numeric scalar) The return value}
}
}
\description{
The function calculates return values using the Gumbel like 2D
  extremal Poisson process, but combining the results over many thresholds
}
\details{
Solves the equation

\deqn{\int_{y_N}^\infty\int_0^1\lambda(t, y)dtdy = \frac{1}{N}}

where \eqn{\lambda(t, y)} is given in the documentation for \code{gumbelMLE}
for many thresholds and combines the results by weighted averaging.  The
weights are described in the documentation for \code{gumbelMultiFit}.
}
\examples{

\dontrun{

complete_series <- -jp1tap1715wind270$value

declustered_obs <- decluster(complete_series)

gumbel_multi_fit <- gumbelMultiFit(x = declustered_obs,
                                   lt = 100,
                                   n_min = 10,
                                   n_max = 50,
                                   weight_scale = 5)

500_second_val <- gumbelNYear(x = gumbel_multi_fit, N = 500)
}

}
