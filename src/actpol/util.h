//
// actpol/util.h : libactpol header file
//
// 2012 Mike Nolta <mike@nolta.net>
//

#pragma once

#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline double secs2days( double s ) { return s/86400.; }
static inline double deg2rad( double deg ) { return deg*M_PI/180.; }
static inline double rad2deg( double rad ) { return rad*180./M_PI; }
static inline double arcsec2rad( double sec ) { return sec*M_PI/(180*3600); }
static inline double rad2arcsec( double rad ) { return rad*180*3600/M_PI; }
static inline double jd2mjd( double jd ) { return jd - 2400000.5; }

#ifdef __cplusplus
}
#endif

