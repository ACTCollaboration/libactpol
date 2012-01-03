//
// actpol/astro.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

// get definition of M_PI
#include <math.h>
#include "constants.h"

#include "quaternion.h"
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline double secs2days( double s ) { return s/86400.; }
static inline double deg2rad( double deg ) { return deg*M_PI/180.; }
static inline double rad2deg( double rad ) { return rad*180./M_PI; }
static inline double arcsec2rad( double sec ) { return sec*M_PI/(180*3600); }
static inline double rad2arcsec( double rad ) { return rad*180*3600/M_PI; }

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
actpol_diurnal_aberration(const double r[3], Quaternion q);

void
actpol_rotate_NWU_to_ITRS(Quaternion q);

void
actpol_rotate_ITRS_to_GCRS(double unixtime, Quaternion q);

void
actpol_NWU_to_GCRS_rotation(double unixtime, Quaternion q);

#ifdef __cplusplus
}
#endif

