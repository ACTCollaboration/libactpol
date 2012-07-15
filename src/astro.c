
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

int
actpol_altaz_to_radec(const ACTpolWeather *weather, double freq_GHz, double unixtime, double alt, double az, double *ra, double *dec)
{
    double r[3];
    double mat[3][3];
    double jd_tt[2];

    Quaternion focalplane_to_topo;
    Quaternion diurnal_aberration, focalplane_to_apparent;
    Quaternion focalplane_to_GCRS, NWU_to_GCRS;
    Quaternion GCRS_to_BCRS;
    Quaternion focalplane_to_BCRS;

    double ref = actpol_refraction(weather, freq_GHz, alt);

    // focalplane -> topo
    Quaternion_identity(focalplane_to_topo);
    actpol_rotate_focalplane_to_NWU(alt-ref, az, focalplane_to_topo);

    // diurnal aberration
    Quaternion_unit(focalplane_to_topo);
    Quaternion_to_matrix_col3(focalplane_to_topo, r);
    actpol_diurnal_aberration(r, diurnal_aberration);
    Quaternion_mul(focalplane_to_apparent, diurnal_aberration, focalplane_to_topo);

    // focalplane -> GCRS
    actpol_NWU_to_GCRS_rotation(unixtime, NWU_to_GCRS);
    Quaternion_mul(focalplane_to_GCRS, NWU_to_GCRS, focalplane_to_apparent);

    // center of array in GCRS
    Quaternion_unit(focalplane_to_GCRS);
    Quaternion_to_matrix_col3(focalplane_to_GCRS, r);

    // annual aberration correction
    actpol_unixtime_to_jd_tt(unixtime, jd_tt);
    actpol_annual_aberration(jd_tt, r, GCRS_to_BCRS);
    Quaternion_mul(focalplane_to_BCRS, GCRS_to_BCRS, focalplane_to_GCRS);

    Quaternion_conj(focalplane_to_BCRS); // transpose mat
    Quaternion_to_matrix(focalplane_to_BCRS, mat);
    double *rr = mat[2];
    *ra = atan2(rr[1], rr[0]);
    *dec = atan2(rr[2], hypot(rr[0],rr[1]));

    return 0;
}

int
actpol_radec_to_crude_altaz(double unixtime, double ra, double dec, double *alt, double *az)
{
    Quaternion NWU_to_GCRS;
    double mat[3][3];
    double r[3], p[3];

    p[0] = cos(dec)*cos(ra);
    p[1] = cos(dec)*sin(ra);
    p[2] = sin(dec);

    actpol_NWU_to_GCRS_rotation(unixtime, NWU_to_GCRS);
    Quaternion_conj(NWU_to_GCRS); // transpose mat
    Quaternion_to_matrix(NWU_to_GCRS, mat);

    matrix_times_vec3(r, mat, p);

    *az = -atan2(r[1], r[0]);
    *alt = atan2(r[2], hypot(r[0],r[1]));

    return 0;
}

