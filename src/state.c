
#include <assert.h>
#include <stdlib.h>

#include <slalib.h>
#include <sofa.h>

#include "actpol/astro.h"
#include "actpol/constants.h"
#include "actpol/math.h"
#include "actpol/state.h"

void
ACTpolScan_init(ACTpolScan *scan, double mean_alt, double mean_az, double mean_throw)
{
    scan->mean_alt = mean_alt;
    scan->mean_az = mean_az;
    scan->mean_throw = mean_throw;
}

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
    return state;
}

void
ACTpolState_free(ACTpolState *state)
{
    free(state);
}

void
ACTpolState_init(ACTpolState *state)
{
    state->slerp_unixtime0 = 0.;
    state->slerp_length = 660.; // 11 minutes
    state->eob_unixtime = 0.;
}

void
actpol_rotate_focalplane_to_NWU(double boresight_alt, double boresight_az, Quaternion q)
{
    Quaternion_r3_mul(M_PI, q);
    Quaternion_r2_mul(M_PI_2 - boresight_alt, q);
    Quaternion_r3_mul(-boresight_az, q);
}

void
ACTpolState_update_boresight(ACTpolState *state, double boresight_alt, double boresight_az)
{
    state->boresight_alt = boresight_alt;
    state->boresight_az = boresight_az;
}

static void
ACTpolState_tt_jd(const ACTpolState *state, double jd_tt[2])
{
    double jd_utc[2], jd_tai[2];
    int stat;

    jd_utc[0] = UNIX_JD_EPOCH;
    jd_utc[1] = secs2days(state->unixtime);

    // utc -> tai
    stat = iauUtctai(jd_utc[0], jd_utc[1], jd_tai+0, jd_tai+1);
    assert(stat == 0);

    // tai -> tt
    stat = iauTaitt(jd_tai[0], jd_tai[1], jd_tt+0, jd_tt+1);
    assert(stat == 0);
}

static void
ACTpolState_update_earth_orbital_velocity(ACTpolState *state)
{
    double jd_tdb[2], pvb[2][3];

    ACTpolState_tt_jd(state, jd_tdb); // tt ~ tdb
    iauEpv00(jd_tdb[0], jd_tdb[1], pvb, pvb);

    for (int i = 0; i < 3; i++)
        state->earth_orbital_beta[i] = pvb[1][i]/SPEED_OF_LIGHT_AU_PER_D;

    state->eob_unixtime = state->unixtime;
}

void
ACTpolState_update_unixtime(ACTpolState *state, double unixtime)
{
    state->unixtime = unixtime;
    actpol_NWU_to_GCRS_rotation(unixtime, state->NWU_to_GCRS_q);
    ACTpolState_update_earth_orbital_velocity(state);
}

void
ACTpolState_update_unixtime_fast(ACTpolState *state, double unixtime)
{
    state->unixtime = unixtime;

    double t = (unixtime - state->slerp_unixtime0)/state->slerp_length;
    if (0. <= t && t <= 1.) {
        QuaternionSlerp_interpolate(&state->slerp, t, state->NWU_to_GCRS_q);
    } else {
        Quaternion q1;
        state->slerp_unixtime0 = unixtime;
        actpol_NWU_to_GCRS_rotation(unixtime, state->NWU_to_GCRS_q);
        actpol_NWU_to_GCRS_rotation(unixtime + state->slerp_length, q1);
        QuaternionSlerp_init(&state->slerp, state->NWU_to_GCRS_q, q1);
    }

    // recalc every 2.4 secs (time it takes to move ~0.1")
    if (unixtime > state->eob_unixtime + 2.4)
        ACTpolState_update_earth_orbital_velocity(state);
}

void
ACTpolState_update(ACTpolState *state, double unixtime,
    double boresight_alt, double boresight_az)
{
    ACTpolState_update_unixtime_fast(state, unixtime);
    ACTpolState_update_boresight(state, boresight_alt, boresight_az);
}

