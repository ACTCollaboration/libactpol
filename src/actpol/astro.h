//
// actpol/astro.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#include "constants.h"
#include "math.h"
#include "quaternion.h"
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline void
actpol_ang2vec(double a, double d, double r[3])
{
    double cos_d = cos(d);
    r[0] = cos_d*cos(a);
    r[1] = cos_d*sin(a);
    r[2] = sin(d);
}

static inline void
actpol_vec2ang(const double r[3], double *a, double *d)
{
    *a = atan2(r[1], r[0]);
    *d = atan2(r[2], hypot(r[0],r[1]));
}

void
actpol_aberration(const double u[3], const double beta[3], Quaternion q);

void
actpol_diurnal_aberration(const double r[3], Quaternion q);

void
actpol_unixtime_to_jd_tt(double unixtime, double jd_tt[2]);

void
actpol_earth_orbital_beta(const double jd_tdb[2], double beta[3]);

void
actpol_annual_aberration(const double jd_tdb[2], const double r[3], Quaternion q);

void
actpol_rotate_NWU_to_ITRS(Quaternion q);

void
actpol_rotate_ITRS_to_GCRS(double unixtime, Quaternion q);

void
actpol_NWU_to_GCRS_rotation(double unixtime, Quaternion q);

int
actpol_altaz_to_radec(const ACTpolWeather *weather, double freq_GHz, double unixtime, double alt, double az, double *ra, double *dec);

int
actpol_radec_to_crude_altaz(double unixtime, double ra, double dec, double *alt, double *az);

#ifdef __cplusplus
}
#endif

