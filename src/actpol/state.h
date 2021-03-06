//
// actpol/state.h : libactpol header file
//
// 2011 Mike Nolta <mike@nolta.net>
//

#pragma once

#include "quaternion.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double mean_alt;
    double mean_az;
    double mean_throw;
}
ACTpolScan;

ACTpolScan *
ACTpolScan_alloc(void);

void
ACTpolScan_init(ACTpolScan *scan, double mean_alt, double mean_az, double mean_throw);

typedef struct
{
    double temperature_C;
    double pressure_mbar;
    double relative_humidity; // percentage
    double tropospheric_lapse_rate_K_per_m;
}
ACTpolWeather;

ACTpolWeather *
ACTpolWeather_alloc(void);

void
ACTpolWeather_default(ACTpolWeather *weather);

double
actpol_refraction(const ACTpolWeather *weather, double freq_GHz, double alt);
double
actpol_refraction_alma366(const ACTpolWeather *weather, double freq_GHz, double alt);
double
actpol_refraction_ulich(const ACTpolWeather *weather, double freq_GHz, double alt);

void
actpol_rotate_focalplane_to_NWU(double boresight_alt, double boresight_az, Quaternion q);

typedef struct
{
    double boresight_alt;
    double boresight_az;
    double unixtime;

    Quaternion NWU_to_GCRS_q;
    QuaternionSlerp slerp;
    double slerp_unixtime0;
    double slerp_length;

    // earth orbital beta (v/c)
    double eob_unixtime;
    double earth_orbital_beta[3];
}
ACTpolState;

ACTpolState *
ACTpolState_alloc(void);

void
ACTpolState_free(ACTpolState *state);

void
ACTpolState_init(ACTpolState *state);

void
ACTpolState_update(ACTpolState *state, double unixtime, double boresight_alt, double boresight_az);

void
ACTpolState_update_unixtime(ACTpolState *state, double unixtime);

void
ACTpolState_update_unixtime_fast(ACTpolState *state, double unixtime);

void
ACTpolState_update_boresight(ACTpolState *state, double boresight_alt, double boresight_az);

#ifdef __cplusplus
}
#endif

