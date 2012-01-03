
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <slalib.h>

#include "actpol/constants.h"
#include "actpol/astro.h"
#include "actpol/state.h"

void
ACTpolWeather_default(ACTpolWeather *w)
{
    w->temperature_C = 0.;
    w->pressure_mbar = 550.;
    w->relative_humidity = 0.2;
    w->tropospheric_lapse_rate_K_per_m = 0.0065;
}

double
actpol_refraction(const ACTpolWeather *weather, double freq_GHz, double alt)
{
    double ref;
    double temperature_K = weather->temperature_C + 273.15;
    slaRefro(M_PI_2 - alt,
        ACTPOL_ELEVATION_METERS,
        temperature_K,
        weather->pressure_mbar,
        weather->relative_humidity,
        SPEED_OF_LIGHT_M_PER_S*1e-3/freq_GHz, // wavelength [microns]
        ACTPOL_LATITUDE,
        weather->tropospheric_lapse_rate_K_per_m,
        1e-8,
        &ref);
    return ref;
}

ACTpolState *
ACTpolState_alloc(void)
{
    ACTpolState *state = (ACTpolState *)malloc(sizeof(ACTpolState));
    assert(state);
    ACTpolWeather_default(&state->weather);
    state->slerp_unixtime0 = 0.;
    state->slerp_length = 660.; // 11 minutes
    return state;
}

void
ACTpolState_free(ACTpolState *state)
{
    free(state);
}


void
ACTpolState_update_boresight(ACTpolState *state, double alt, double az)
{
    state->boresight_alt = alt;
    state->boresight_az = az;
}

void
ACTpolState_update_unixtime(ACTpolState *state, double unixtime)
{
    state->unixtime = unixtime;
    actpol_NWU_to_GCRS_rotation(unixtime, state->q);
}

void
ACTpolState_update_unixtime_fast(ACTpolState *state, double unixtime)
{
    state->unixtime = unixtime;
    double t = (unixtime - state->slerp_unixtime0)/state->slerp_length;
    if (0. <= t && t <= 1.) {
        QuaternionSlerp_interpolate(&state->slerp, t, state->q);
    } else {
        Quaternion q1;
        state->slerp_unixtime0 = unixtime;
        actpol_NWU_to_GCRS_rotation(unixtime, state->q);
        actpol_NWU_to_GCRS_rotation(unixtime + state->slerp_length, q1);
        QuaternionSlerp_init(&state->slerp, state->q, q1);
    }
}

