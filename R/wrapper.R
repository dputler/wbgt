wbgt <- function(
    year,
    month,
    day,
    hour,
    minute,
    gmt,
    avg,
    lat,
    lon,
    solar,
    pres,
    Tair,
    relhum,
    speed,
    zspeed,
    dT,
    urban
) {
    num_obs <- length(year)
    Tg <- rep(0.0, num_obs)
    Tnwb <- rep(0.0, num_obs)
    Tpsy <- rep(0.0, num_obs)
    Twbg <- rep(0.0, num_obs)
    status <- rep(0, num_obs)

    out <- .C(
        "wbgt",
        num_obs = as.integer(num_obs),
        year = as.integer(year),
        month = as.integer(month),
        day = as.integer(day),
        hour = as.integer(hour),
        minute = as.integer(minute),
        gmt = as.integer(gmt),
        avg = as.integer(avg),
        lat = as.double(lat),
        lon = as.double(lon),
        solar = as.double(solar),
        pres = as.double(pres),
        Tair = as.double(Tair),
        relhum = as.double(relhum),
        speed = as.double(speed),
        zspeed = as.double(zspeed),
        dT = as.double(dT),
        urban = as.integer(urban),
        Tg = as.double(Tg),
        Tnwb = as.double(Tnwb),
        Tpsy = as.double(Tpsy),
        Twbg = as.double(Twbg),
        status = as.integer(status),
        PACKAGE="wbgt"
    )

    # any calculation based on pressure outside of normal range yields NA
    out[pres > 2000 | pres < 800] <- NA
    out
}

# Apply the wbgt function above to a data frame, assumed to contain columns
# with names equal to the names of the inputs to wbgt
# Returns copy of the original data frame with the wbgt appended as a new column
wbgt_df <- function(data) {
    if (!requireNamespace("dplyr", quietly=TRUE)) {
        stop("Please install dplyr in order to use this function.")
    }
    with(data, dplyr::mutate(data, wbgt = wbgt(year, month, day, hour, minute,
        gmt, avg, lat, lon, solar, pres, 
             Tair, relhum, speed, zspeed, dT, urban)$Twbg))
}

calc_solar <- function(date_time, lat, lon) {
    # Calculate the needed date-time variables
    time_1900 <- as.POSIXct("1900-01-01 00:00:00")
    date_time_posix <- as.POSIXct(date_time)
    year <- as.integer(format(date_time_posix, format = "%Y"))
    month <- as.integer(format(date_time_posix, format = "%m"))
    dt_day <- as.numeric(format(date_time_posix, format = "%d"))
    dt_hour <- as.numeric(format(date_time_posix, format = "%H"))
    dt_minute <- as.numeric(format(date_time_posix, format = "%M"))
    decimal_hour <- dt_hour + dt_minute / 60.0
    day <- dt_day + decimal_hour / 24.0
    days_1900 <- as.numeric(date_time_posix - time_1900) + decimal_hour / 24.0
    # Prep the output variables from the C call
    num_obs <- length(year)
    ap_ra <- rep(0.0, num_obs)
    ap_dec <- rep(0.0, num_obs)
    altitude <- rep(0.0, num_obs)
    refraction <- rep(0.0, num_obs)
    azimuth <- rep(0.0, num_obs)
    distance <- rep(0.0, num_obs)
    status <- rep(0, num_obs)
    # Call the C function
    out <- .C(
        "calc_solar",
        num_obs = as.integer(num_obs),
        year = as.integer(year),
        month = as.integer(month),
        day = as.double(day),
        days_1900 = as.double(days_1900),
        latitude = as.double(latitude),
        longitude = as.double(longitude),
        ap_ra = as.double(ap_ra),
        ap_dec = as.double(ap_dec),
        altitude = as.double(altitude),
        refraction = as.double(refraction),
        azimuth = as.double(azimuth),
        distance = as.integer(distance),
        status = as.integer(status),
        PACKAGE = "wbgt"
    )
    out
}