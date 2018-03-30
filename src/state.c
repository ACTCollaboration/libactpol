
#include <assert.h>
#include <stdlib.h>

#include "../config.h"
#if HAVE_LIBSLALIB
#include <slalib.h>
#elif HAVE_LIBSLAREFRO
#include <slarefro.h>
#else
#error("Must have one of slalib or slarefro.")
#endif

#include "actpol/astro.h"
#include "actpol/constants.h"
#include "actpol/math.h"
#include "actpol/state.h"
#include "actpol/util.h"

ACTpolScan *
ACTpolScan_alloc(void)
{
    ACTpolScan *scan = (ACTpolScan *)malloc(sizeof(ACTpolScan));
    assert(scan);
    return scan;
}

void
ACTpolScan_init(ACTpolScan *scan, double mean_alt, double mean_az, double mean_throw)
{
    scan->mean_alt = mean_alt;
    scan->mean_az = mean_az;
    scan->mean_throw = mean_throw;
}

ACTpolWeather *
ACTpolWeather_alloc(void)
{
    ACTpolWeather *weather = (ACTpolWeather *)malloc(sizeof(ACTpolWeather));
    assert(weather);
    return weather;
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
#if HAVE_LIBSLALIB
    slaRefro(
#elif HAVE_LIBSLAREFRO
    slaf_refro(
#endif
        M_PI_2 - alt,
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

// ALMA technical memo 366
// doesn't agree w/ figure 3 yet
double
actpol_refraction_alma366(const ACTpolWeather *weather, double freq_GHz, double alt)
{
    const double cos_alt = cos(alt), sin_alt = sin(alt);

    /* surface temperature [K] */
    const double Ts = weather->temperature_C + 273.15;

    /* total surface barometric pressure [mb] */
    const double Ps = weather->pressure_mbar;

    /* surface saturated water vapor pressure [mb], eqn 16 */
    const double esat = (1.0007 + 3.46e-6*Ps) * 6.1121 * exp(17.502*(Ts - 273.15)/(Ts - 32.18));

    /* partial surface pressure of water vapor [mb], eqn 14 */
    const double Pw = weather->relative_humidity * esat;

    const double N0 = (77.6*Ps + (-5.6 + 3.75e5/Ts)*Pw)/Ts; /* eqn 12 */
    const double R0 = 1e-6*N0; /* [radians], eqn 11 */

    /* eqn 30 */
    const double A1 = 0.6306849
       + 0.6069e-4 * (Ps - 1013.25)
       - 0.2532e-4 * Pw
       - 0.9881e-6 * Pw * Pw
       - 0.5154e-3 * (Ts - 258.15)
       + 0.2880e-5 * (Ts - 258.15) * (Ts - 258.15);
    const double A2 = 1.302642;

    /* effective height of the atmosphere [km] */
    const double H = 8.31434 * Ts / (28.97 * 9.784);

    /* eqn 27 for I */
    const double I2_csc_alt = (0.5*6378./H)*sin_alt/(cos_alt*cos_alt);
    const double inv_mprime =
        (sin_alt + A1 /
          (I2_csc_alt + A2 /
            (sin_alt + 13.24969 /
              (I2_csc_alt + 173.4233))));

    return R0*cos_alt/inv_mprime;
}

double
actpol_refraction_ulich(const ACTpolWeather *weather, double freq_GHz, double alt)
{
    const double cos_alt = cos(alt), sin_alt = sin(alt);

    /* surface temperature [K] */
    const double Ts = weather->temperature_C + 273.15;

    /* total surface barometric pressure [mb] */
    const double Ps = weather->pressure_mbar;

    /* surface saturated water vapor pressure [mb], eqn 16 */
    const double esat = (1.0007 + 3.46e-6*Ps) * 6.1121 * exp(17.502*(Ts - 273.15)/(Ts - 32.18));

    /* partial surface pressure of water vapor [mb], eqn 14 */
    const double Pw = weather->relative_humidity * esat;

    const double R0 = 16.01*Ps/Ts - 1.15*Pw/Ts + 7.734937e4*Pw/(Ts*Ts);
    const double f = cos_alt/(sin_alt + 0.00175*tan(deg2rad(87.5) - alt));
    return R0*f;
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
ACTpolState_update_earth_orbital_velocity(ACTpolState *state)
{
    double jd_tdb[2];
    actpol_unixtime_to_jd_tt(state->unixtime, jd_tdb); // tt ~ tdb
    actpol_earth_orbital_beta(jd_tdb, state->earth_orbital_beta);
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

