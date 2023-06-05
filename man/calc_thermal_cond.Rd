\name{calc_thermal_cond}
\alias{calc_thermal_cond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the thermal condition of the air}
\description{
    The thermodynamic model of heat transfer of a natural wet bulb temperature sensor wick.
}
\usage{
calc_thermal_cond(Tair)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Tair}{
        A numeric vector of two meter air temperatures in K.
    }
}
\value{
    The thermal condition of the air.
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
