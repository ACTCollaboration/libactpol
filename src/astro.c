
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <slalib.h>
#include <sofa.h>

#include "actpol/astro.h"
#include "actpol/constants.h"
#include "actpol/vec3.h"

#include "iers_bulletin_a.h"

// Julian date for unixtime = 0
#define UNIX_JD_EPOCH 2440587.5

static double inline mjd(double jd) { return jd - 2400000.5; }

void
actpol_ang2vec(double a, double d, double r[3])
{
    double cos_d = cos(d);
    r[0] = cos_d*cos(a);
    r[1] = cos_d*sin(a);
    r[2] = sin(d);
}

void
actpol_vec2ang(double r[3], double *a, double *d)
{
    *a = atan2(r[1], r[0]);
    *d = atan2(r[2], hypot(r[0],r[1]));
}

void
actpol_diurnal_aberration(double r[3], Quaternion q)
{
    const double v = 0.295043*M_PI/(180.*3600.);
    double n[3], east[3] = {0.,-1.,0.}; // NWU
    vec3_cross_product(r, east, n);
    double angle = v*vec3_norm(n)/vec3_norm(r);
    Quaternion_rot(q, -angle, n);
}

void
actpol_NWU_to_ITRS_quaternion(Quaternion q)
{
    Quaternion_r3_mul(M_PI, q);
    Quaternion_r2_mul(M_PI_2 - ACTPOL_LATITUDE, q);
    Quaternion_r3_mul(ACTPOL_LONGITUDE_EAST, q);
}

/*
void
actpol_UEN_to_ITRS_quaternion(Quaternion q)
{
    Quaternion_r2(-deg2rad(ACTPOL_LATITUDE_DEG), q);
    Quaternion_r3_mul(deg2rad(ACTPOL_LONGITUDE_EAST_DEG), q);
}
*/

void
actpol_ITRS_to_GCRS_quaternion(double unixtime, Quaternion q)
{
    const double sprime = 0.; // ~ 1e-4 "
    const double dX06 = 0.; // ~ 1e-4 "
    const double dY06 = 0.; // ~ 1e-4 "
    double X, Y, Z, s, E, d, theta;
    double ut1_minus_utc, xp, yp;
    double jd_utc[2], jd_ut1[2], jd_tai[2], jd_tt[2];
    int stat;

    jd_utc[0] = UNIX_JD_EPOCH;
    jd_utc[1] = secs2days(unixtime);

    double mjd_utc = mjd(jd_utc[0]) + jd_utc[1];
    stat = get_iers_bulletin_a(mjd_utc, &ut1_minus_utc, &xp, &yp);
    assert(stat == 0);

    // utc -> tai
    stat = iauUtctai(jd_utc[0], jd_utc[1], jd_tai+0, jd_tai+1);
    assert(stat == 0);

    // tai -> tt
    stat = iauTaitt(jd_tai[0], jd_tai[1], jd_tt+0, jd_tt+1);
    assert(stat == 0);

    // utc -> ut1
    stat = iauUtcut1(jd_utc[0], jd_utc[1], ut1_minus_utc, jd_ut1+0, jd_ut1+1);
    assert(stat == 0);

    // precession, nutation, frame bias
    iauXys06a(jd_tt[0], jd_tt[1], &X, &Y, &s);
    X += dX06;
    Y += dY06;
    Z = sqrt(1. - X*X - Y*Y);

    // X = sin d cos E; Y = sin d sin E; Z = cos d
    E = atan2(Y, X);
    d = acos(Z);

    // earth rotation
    theta = iauEra00(jd_ut1[0], jd_ut1[1]);

    Quaternion_r1_mul(-arcsec2rad(yp), q);
    Quaternion_r2_mul(-arcsec2rad(xp), q);
    Quaternion_r3_mul(-E-s+theta+sprime, q);
    Quaternion_r2_mul(d, q);
    Quaternion_r3_mul(E, q);
    /*
    Quaternion_r3_mul(-E, q);
    Quaternion_r2_mul(-d, q);
    Quaternion_r3_mul(E+s-theta-sprime, q);
    Quaternion_r2_mul(xp, q);
    Quaternion_r1_mul(yp, q);
    */
}

void
ACTpolWeather_default(ACTpolWeather *w)
{
    w->temperature_K = 273.;
    w->pressure_mbar = 550.;
    w->relative_humidity = 0.2;
    w->tropospheric_lapse_rate_K_per_m = 0.0065;
}

double
actpol_refraction(ACTpolWeather *weather, double freq_GHz, double alt)
{
    double ref;
    slaRefro(M_PI_2 - alt,
        ACTPOL_ELEVATION_METERS,
        weather->temperature_K,
        weather->pressure_mbar,
        weather->relative_humidity,
        299792.458/freq_GHz, // wavelength [micrometers]
        ACTPOL_LATITUDE,
        weather->tropospheric_lapse_rate_K_per_m,
        1e-8,
        &ref);
    return ref;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// old ACT pointing code

static double
convert_ctime_to_utc_mjd( double ctime )
{
    // Not strictly correct, because of leap seconds.
    return secs2days(ctime) + 40587.0;
}

static double
convert_utc_to_tt( double utc )
{
    return utc + secs2days(slaDtt(utc));
}

int
observed_altaz_to_mean_radec( const ACTpolWeather *weather, double freq_ghz,
        int n, const double ctime[],
        const double alt[], const double az[],
        double ra[], double dec[] )
{
    assert( n > 0 );
    assert( ctime != NULL );
    assert( alt != NULL );
    assert( az != NULL );
    assert( ra != NULL );
    assert( dec != NULL );

    int stat;
    double dut1, x, y;
    double amprms[21], aoprms[14];

    double utc = convert_ctime_to_utc_mjd( ctime[0] );

    stat = get_iers_bulletin_a( utc, &dut1, &x, &y );
    if ( stat != 0 )
        return stat;

    double wavelength_um = 299792.458/freq_ghz;

    slaAoppa( utc, dut1,
            ACTPOL_LONGITUDE_EAST,
            ACTPOL_LATITUDE,
            ACTPOL_ELEVATION_METERS,
            arcsec2rad(x),
            arcsec2rad(y),
            weather->temperature_K,
            weather->pressure_mbar,
            weather->relative_humidity,
            wavelength_um,
            weather->tropospheric_lapse_rate_K_per_m,
            aoprms );

    double tt = convert_utc_to_tt( utc );
    // using tt instead of tdb
    slaMappa( 2000.0, tt, amprms );
    //printf( "freq = %g\n", freq_ghz );

    for ( int i = 0; i < n; i++ )
    {
        double observed_az = az[i];
        double observed_zenith = M_PI/2 - alt[i];
        double apparent_ra, apparent_dec;

        utc = convert_ctime_to_utc_mjd( ctime[i] );
        slaAoppat( utc, aoprms );

        slaOapqk( "A", observed_az, observed_zenith, aoprms,
                &apparent_ra, &apparent_dec );

        //printf( "apparent ra,dec = %g, %g\n", apparent_ra, apparent_dec );

        slaAmpqk( apparent_ra, apparent_dec, amprms, ra+i, dec+i );
    }

    return 0;
}

