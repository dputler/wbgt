\name{calc_wind}
\alias{calc_wind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the parameters associated with estimated wind speed}
\description{
    Estimate the two meter wind speed from the N meter speed adjusted for stability.
}
\usage{
calc_wind(speed, zspeed, solar, dT, daytime, urban)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{speed}{
        A numeric vector of wind speeds at a reference altitude in m/s.
    }
    \item{zspeed}{
        A numeric vector of reference altitudes in m.
    }
    \item{solar}{
        A numeric vector of solar irradience values in W/m^2.
    }
    \item{dT}{
        A numeric vector of the temperature difference between altitudes in K.
    }
    \item{daytime}{
        An integer vector where 0 indicates non-daytime and 1 (or any non-zero int) is daytime.
    }
    \item{urban}{
        An integer vector where 0 indicates non-urban and 1 (or any non-zero int) is urban.
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
