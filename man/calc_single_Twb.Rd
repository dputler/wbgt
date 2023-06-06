\name{calc_single_Twb}
\alias{calc_single_Twb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A single pass calculation of the natural wet bulb temperature}
\description{
    The thermodynamic model of heat transfer of a natural wet bulb temperature sensor wick.
}
\usage{
calc_single_Twb(Tair, rh, Pair, speed, solar, fdir, cza)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Tair}{
        A numeric vector of two meter air temperatures in K.
    }
    \item{rh}{
        A numeric vector of relative humidities.
    }
    \item{Pair}{
        A numeric vector of atmospheric pressure values in millibars.
    }
    \item{speed}{
        A numeric vector of two meter wind speeds in m/s
    }
    \item{solar}{
        A numeric vector of total downwelling solar irradiance in W/m^2.
    }
    \item{fdir}{
        The percentages of total solar irradiance that is direct beam.
    }
}
\value{
    The first pass through the Argonne model natural wet bulb temperature along
    with its subcomponents to find the python translation issue.
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
