\name{calc_irrad}
\alias{calc_irrad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the parameters associated with solar irradience}
\description{
    An R wrapper for Larson's C implemenation of Astronomical Almanac solar irradience
}
\usage{
calc_solar(date_time, lat, lon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{date_time}{
        A character vector of UTC date-time values.
    }
    \item{lat}{
        A numeric vector of latitude values.
    }
    \item{lon}{
        A numeric vector of longitude values.
    }
}
\value{
    A list containing all of the input variables (described above) in addition to
    `solar` (the amount of downwelling solar irradience), `cza` (cosine of the zenith angle), 
    `fdir` (the percentage of solar irradience that is direct), and `status`
    (with 0 indicating convergence and -1 indicating non-convergence of the Liljegren algorithm).
}
\author{
    Dan Putler
}

\examples{
calc_solar("2022-06-30 19:00:00", 37.7667, -122.4167)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{manip}
