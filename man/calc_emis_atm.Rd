\name{calc_emis_atm}
\alias{calc_emis_atm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the emisivity of the atmosphere}
\description{
    The thermodynamic model of heat transfer of a natural wet bulb temperature sensor wick.
}
\usage{
calc_emis_atm(Tair, rh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Tair}{
        A numeric vector of two meter air temperatures in K.
    }
    \item{rh}{
        A numeric vector of decimal relative humidity values (all between 0 and 1).
    }
}
\value{
    The emisitvity of the atmosphere.
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
