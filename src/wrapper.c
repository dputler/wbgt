/*
 * The MIT License (MIT)
 * 
 * Copyright (c) 2016 Max Lieblich
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *   
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include "wbgt.h"
#include <R.h>

/*
 * All parameters below are the same as those defined in wbgt.c, with the
 * exception of num_obs, which counts the number of observations and is used
 * for properly vectorizing the R code
 */
void wbgt(int *num_obs, int *year, int *month, int *day, int *hour, int *minute, int *gmt, int *avg, 
    double *lat, double *lon, double *solar, double *pres, double *Tair, double *relhum, double *speed, double *zspeed, 
    double *dT, int *urban, 
    double *Tg, double *Tnwb, double *Tpsy, double *Twbg,
    int *status)
{
  int n = *num_obs;
  double est_speed = 0.;
  for (int i = 0; i < n; ++i)
  {
    status[i] = calc_wbgt(year[i], month[i], day[i], hour[i], minute[i], gmt[i], avg[i],
        lat[i], lon[i], solar[i], pres[i], Tair[i], relhum[i], speed[i], zspeed[i],
        dT[i], urban[i],
        &est_speed, Tg + i, Tnwb + i, Tpsy + i, Twbg + i);
  }
}
void calc_solar(int *num_obs, int *year, int *month, double *day,
    double *days_1900, double *lat, double *lon,
    double *ap_ra, double *ap_dec, double *altitude, double *refraction,
    double *azimuth, double *distance, int *status)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    status[i] = solarposition(year[i], month[i], day[i], days_1900[i],
        lat[i], lon[i],
        ap_ra + i, ap_dec + i, altitude + i, refraction + i, azimuth + i, distance + i);
  }
}

void calc_wind(int *num_obs, double *speed, double *zspeed, double *solar,
    double *dT, int *daytime, int *urban, int *stb_cls, double *est_wind)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    stb_cls[i] = stab_srdt(daytime[i], speed[i], solar[i], dT[i]);
    est_wind[i] = est_wind_speed(speed[i], zspeed[i], stb_cls[i], urban[i]);
  }
}

void calc_cyl_air(int *num_obs, double diameter, double length, double *Tair,
    double *P_air, double *speed, double *hw)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    hw[i] = h_cylinder_in_air(diameter, length, Tair[i], P_air[i], speed[i]);
  }
}

void calc_viscosity(int *num_obs, double *Tair, double *visc)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    visc[i] = viscosity(Tair[i]);
  }
}

void calc_thermal_cond(int *num_obs, double *Tair, double *thrmcond)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    thrmcond[i] = thermal_cond(Tair[i]);
  }
}

void calc_h_evap(int *num_obs, double *Tair, double *hevap)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    hevap[i] = evap(Tair[i]);
  }
}

void calc_esat(int *num_obs, double *tk, double *esat_out)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    esat_out[i] = esat(tk[i], 0);
  }
}

void calc_dew_point(int *num_obs, double *e, double *dew_pt)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    dew_pt[i] = dew_point(tk[i], 0);
  }
}