\name{calc_cyl_air}
\alias{calc_cyl_air}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the theat transfer of a cylindar (the wick) in the air}
\description{
    The thermodynamic model of heat transfer of a natural wet bulb temperature sensor wick.
}
\usage{
calc_cyl_air(T_a, T_w, P_air, speed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{T_a}{
        A numeric vector of two meter air temperatures in K.
    }
    \item{T_w}{
        A numeric vector of wick temperature values in K.
    }
    \item{P_air}{
        A numeric vector of surface air pressure values in millibars.
    }
    \item{speed}{
        A numeric vector of estimated wind speeds at 2 meters in m/s.
    }
}
\value{
    The heat transfer coefficient of the wick of a QUESTemp 34 WBGT monitor.
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
