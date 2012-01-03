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
    double temperature_C;
    double pressure_mbar;
    double relative_humidity; // percentage
    double tropospheric_lapse_rate_K_per_m;
}
ACTpolWeather;

void
ACTpolWeather_default(ACTpolWeather *weather);

double
actpol_refraction(const ACTpolWeather *weather, double freq_GHz, double alt);

typedef struct
{
    double boresight_alt;
    double boresight_az;
    double unixtime;
    ACTpolWeather weather;

    Quaternion q;
    QuaternionSlerp slerp;
    double slerp_unixtime0;
    double slerp_length;
}
ACTpolState;

ACTpolState *
ACTpolState_alloc(void);

void
ACTpolState_free(ACTpolState *state);

void
ACTpolState_update_unixtime(ACTpolState *state, double unixtime);

void
ACTpolState_update_unixtime_fast(ACTpolState *state, double unixtime);

void
ACTpolState_update_boresight(ACTpolState *state, double alt, double az);

#ifdef __cplusplus
}
#endif

