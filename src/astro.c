
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <slalib.h>
#include <sofa.h>

#include "actpol/astro.h"
#include "actpol/constants.h"
#include "actpol/math.h"
#include "actpol/util.h"
#include "actpol/vec3.h"

#include "iers_bulletin_a.h"

void
actpol_aberration(const double u[3], const double beta[3], Quaternion q)
{
    double n[3];
    vec3_cross_product(n, u, beta);
    double angle = vec3_norm(n);
    Quaternion_rot(q, -angle, n);
}

void
actpol_diurnal_aberration(const double r[3], Quaternion q)
{
    /*
    const double v = 0.295043*M_PI/(180.*3600.);
    double n[3], east[3] = {0.,-1.,0.}; // NWU
    vec3_cross_product(n, r, east);
    double angle = v*vec3_norm(n)/vec3_norm(r);
    Quaternion_rot(q, -angle, n);
    */
    const double beta[3] = {0., arcsec2rad(-0.295043), 0.}; // NWU
    actpol_aberration(r, beta, q);
}

void
actpol_unixtime_to_jd_tt(double unixtime, double jd_tt[2])
{
    double jd_utc[2], jd_tai[2];
    int stat;

    jd_utc[0] = UNIX_JD_EPOCH;
    jd_utc[1] = secs2days(unixtime);

    // utc -> tai
    stat = iauUtctai(jd_utc[0], jd_utc[1], jd_tai+0, jd_tai+1);
    assert(stat == 0);

    // tai -> tt
    stat = iauTaitt(jd_tai[0], jd_tai[1], jd_tt+0, jd_tt+1);
    assert(stat == 0);
}

void
actpol_earth_orbital_beta(const double jd_tdb[2], double beta[3])
{
    double pvb[2][3];
    iauEpv00(jd_tdb[0], jd_tdb[1], pvb, pvb);
    for (int i = 0; i < 3; i++)
        beta[i] = pvb[1][i]/SPEED_OF_LIGHT_AU_PER_D;
}

void
actpol_annual_aberration(const double jd_tdb[2], const double r[3], Quaternion q)
{
    double beta[3];
    actpol_earth_orbital_beta(jd_tdb, beta);
    actpol_aberration(r, beta, q);
}

void
actpol_rotate_NWU_to_ITRS(Quaternion q)
{
    Quaternion_r3_mul(M_PI, q);
    Quaternion_r2_mul(M_PI_2 - ACTPOL_LATITUDE, q);
    Quaternion_r3_mul(ACTPOL_LONGITUDE_EAST, q);
}

void
actpol_rotate_ITRS_to_GCRS(double unixtime, Quaternion q)
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

    double mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
    stat = actpol_get_iers_bulletin_a(mjd_utc, &ut1_minus_utc, &xp, &yp);
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
}

void
actpol_NWU_to_GCRS_rotation(double unixtime, Quaternion q)
{
    Quaternion_identity(q);
    actpol_rotate_NWU_to_ITRS(q);
    actpol_rotate_ITRS_to_GCRS(unixtime, q);
}

