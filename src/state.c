
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
    w->temperature_K = 273.;
    w->pressure_mbar = 550.;
    w->relative_humidity = 0.2;
    w->tropospheric_lapse_rate_K_per_m = 0.0065;
}

double
actpol_refraction(const ACTpolWeather *weather, double freq_GHz, double alt)
{
    double ref;
    slaRefro(M_PI_2 - alt,
        ACTPOL_ELEVATION_METERS,
        weather->temperature_K,
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
    ACTpolWeather_default(&state->weather);
    return state;
}

void
ACTpolState_free(ACTpolState *state)
{
    free(state);
}

