/*
 * actpol/astro.h : libactpol header file
 *
 * Mike Nolta <mike@nolta.net>
 */

#pragma once

#include <math.h>
#include "quaternion.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

static inline double secs2days( double s ) { return s/86400.; }
static inline double deg2rad( double deg ) { return deg*M_PI/180.; }
static inline double rad2deg( double rad ) { return rad*180./M_PI; }
static inline double arcsec2rad( double sec ) { return deg2rad(sec/3600.); }

typedef struct
{
    double temperature_K;
    double pressure_mbar;
    double relative_humidity;   /* percentage */
    double tropospheric_lapse_rate_K_per_m;
}
ACTpolWeather;

void
ACTpolWeather_default(ACTpolWeather *weather);

double
actpol_refraction(ACTpolWeather *weather, double freq_GHz, double alt);

int
observed_altaz_to_mean_radec( const ACTpolWeather *weather, double freq_GHz,
        int n, const double ctime[], const double alt[], const double az[],
        double ra[], double dec[] );

void actpol_ang2vec(double alt, double az, double r[3]);
void actpol_vec2ang(double r[3], double *alt, double *az);

void
actpol_diurnal_aberration(double r[3], Quaternion q);

void
actpol_NWU_to_ITRS_quaternion(Quaternion q);

void
actpol_ITRS_to_GCRS_quaternion(double unixtime, Quaternion q);

#ifdef __cplusplus
}
#endif

