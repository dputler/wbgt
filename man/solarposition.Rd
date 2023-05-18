\name{solarposition}
\alias{solarposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the web bulb globe temperature (wbgt)}
\description{
    An R wrapper for Larson's C implemenation of Astronomical Almanac solar pos.
}
\usage{
solarposition(date_time, latitude, longitude)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{date_time}{
        A character vector of UTC date-time values.
    }
    \item{latitude}{
        A numeric vector of latitude values.
    }
    \item{lon}{
        A numeric vector of longitude values.
}
\value{
    A list containing all of the input variables (described above) in addition to
    `Twbg` (the wet bulb globe temperature in degrees C), `Tpsy` (psychrometric wet bulb temperature in degrees C), 
    `Tnwb` (natural wet bulb temperature in degrees C), `Tg` (the globe temperature in degrees C), and `status`
    (with 0 indicating convergence and -1 indicating non-convergence of the Liljegren algorithm).
}
\author{
    Dan Putler
}

\examples{
solarposition("2022-06-30 19:00:00", 37.7667, -122.4167)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{manip}
