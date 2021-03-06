
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <slalib.h>
#include <sofa.h>
#include <wcslib/wcslib.h>

#include "actpol/actpol.h"
#include "actpol/iers_bulletin_a.h"
#include "actpol/math.h"
#include "actpol/oldact.h"
#include "actpol/vec3.h"
#include "actpol/util.h"

#define TEST_FILENAME "/tmp/libactpol_test_map.fits";

void
test_map(void)
{
    ACTpolMap *map;
    const char *filename = TEST_FILENAME;
    const char *ofilename = "!" TEST_FILENAME;
    const double test_ra = 1.2, test_dec = -2.3;
    const double test_val = 4.456;
    const double pixsize = 0.00825;
    long test_pix, pix, pix2;
    double ra, dec, ra2, dec2;
    int stat;

    wcserr_enable(1);

    map = ACTpolMap_new(-2, 4.3, -5, 5, pixsize);
    assert(map != NULL);
    assert(map->naxis1 > 0);
    assert(map->naxis2 > 0);
    assert(map->npix > 0);
    //printf("%ld %ld %ld\n", map->naxis1, map->naxis2, map->npix);

    test_pix = ACTpolMap_sky2pix(map, test_ra, test_dec);
    //printf("test_pix = %ld\n", test_pix);
    assert(test_pix >= 0);
    assert(test_pix < map->npix);
    map->data[test_pix] = test_val;

    stat = ACTpolMap_pix2sky(map, test_pix, &ra, &dec);
    //printf("ra, dec = %g, %g\n", ra, dec);
    assert(stat == 0);
    pix = ACTpolMap_sky2pix(map, ra, dec);
    //printf("pix = %ld\n", pix);
    assert(test_pix == pix);

    // detailed pixel test
    for (int i = 0; i != 20; ++i) {
        double ira = ra-i*pixsize*0.2, idec = dec+i*pixsize*0.2;
        pix = ACTpolMap_sky2pix(map, ira, idec);
        stat = ACTpolMap_pix2sky(map, pix, &ra2, &dec2);
        assert(stat == 0);
        pix2 = ACTpolMap_sky2pix(map, ra2, dec2);
        assert(pix == pix2);
        long pix3 = ACTpolMap_sky2pix_cea_fast(map, deg2rad(ira), sin(deg2rad(idec)));
        //printf("pix, pix3 = %ld, %ld\n", pix, pix3);
        assert(pix == pix3);
    }

    stat = ACTpolMap_write_to_fits(map, ofilename);
    assert(stat == 0);
    ACTpolMap_free(map);

    map = ACTpolMap_read_from_fits(filename);
    assert(map != NULL);
    assert(map->data[test_pix] == test_val);
    pix = ACTpolMap_sky2pix(map, test_ra, test_dec);
    assert(test_pix == pix);

    ACTpolMap_free(map);
}

void
print_mat(double mat[3][3])
{
    for (int i = 0; i < 3; i++)
        printf("%+.15f %+.15f %+.15f\n", mat[i][0], mat[i][1], mat[i][2]);
}

void
print_quat(Quaternion q)
{
    printf("q: %.15e %.15e %.15e %.15e\n", q[0], q[1], q[2], q[3]);
}

void
print_vec(const char *s, double q[3])
{
    printf("%s %+.15e %+.15e %+.15e\n", s, q[0], q[1], q[2]);
    printf("%s norm = %.15e\n", s, vec3_norm(q));
}

void
altaz2itrs(double alt, double az, double r[3])
{
    double uze[3], une[3], uwe[3];
    double coslat = ACTPOL_COS_LATITUDE;
    double sinlat = ACTPOL_SIN_LATITUDE;
    double coslon = ACTPOL_COS_LONGITUDE_EAST;
    double sinlon = ACTPOL_SIN_LONGITUDE_EAST;

    //Set up orthonormal basis vectors in local Earth-fixed system.
    //Define vector toward local zenith in Earth-fixed system (z axis).
    uze[0] = coslat * coslon;
    uze[1] = coslat * sinlon;
    uze[2] = sinlat;

    //Define vector toward local north in Earth-fixed system (x axis).
    une[0] = -sinlat * coslon;
    une[1] = -sinlat * sinlon;
    une[2] = coslat;

    //Define vector toward local west in Earth-fixed system (y axis).
    uwe[0] = sinlon;
    uwe[1] = -coslon;
    uwe[2] = 0.0;

    double n, e, u;
    u = sin(alt);
    n = cos(alt)*cos(az);
    e = cos(alt)*sin(az);

    for (int i = 0; i != 3; ++i)
        r[i] = u*uze[i] + n*une[i] - e*uwe[i];
}

void
check_focalplane_to_NWU(double alt, double az)
{
    double r[3], mat[3][3];
    Quaternion q;
    Quaternion_identity(q);
    actpol_rotate_focalplane_to_NWU(alt, az, q);

    Quaternion_to_matrix(q, mat);
    actpol_ang2vec(-az, alt, r);

    for (int i = 0; i != 3; ++i)
        assert(fabs(r[i] - mat[i][2]) < 2e-16);
}

void
check_horizon_to_itrs(double alt, double az)
{
    double t[3], h[3], r[3], mat[3][3];
    Quaternion q;

    altaz2itrs(alt, az, t);
    //print_vec("true", t);

    actpol_ang2vec(-az, alt, h);
    Quaternion_identity(q);
    actpol_rotate_NWU_to_ITRS(q);
    Quaternion_to_matrix(q, mat);
    matrix_times_vec3(r, mat, h);
    //print_vec("quat", r);

    for (int i = 0; i != 3; ++i)
        assert(fabs(r[i] - t[i]) < 3e-16);
}

void
check_itrs_to_gcrs(double unixtime)
{
    double djmjd0, date, time, utc, dat, tai, tt, tut, ut1;
    double X, Y, s, theta, xp, yp, dut1;

    double jd_utc[2];
    jd_utc[0] = 2440587.5;
    jd_utc[1] = secs2days(unixtime);
    double mjd_utc = jd2mjd(jd_utc[0]) + jd_utc[1];
    int stat = actpol_get_iers_bulletin_a(mjd_utc, &dut1, &xp, &yp);
    assert(stat == 0);
    xp = arcsec2rad(xp);
    yp = arcsec2rad(yp);

    djmjd0 = 2400000.5;
    time = modf(mjd_utc, &date);
    //iauCal2jd(2007, 4, 5, &djmjd0, &date);
    //time = 12./24.;
    utc = date + time;
    //iauDat(2007, 4, 5, time, &dat);
    dat = 34.;
    tai = utc + dat/86400.;
    tt = tai + 32.184/86400.;
    tut = time + dut1/86400.;
    ut1 = date + tut;

    iauXy06(djmjd0, tt, &X, &Y);
    s = iauS06(djmjd0, tt, X, Y);
    theta = iauEra00(djmjd0+date, tut);

    double rc2i[3][3], rc2ti[3][3], rpom[3][3], rc2it[3][3];
    iauC2ixys(X, Y, s, rc2i);
    iauCr(rc2i, rc2ti);
    iauRz(theta, rc2ti);
    iauPom00(xp, yp, 0., rpom);
    iauRxr(rpom, rc2ti, rc2it);
    //print_mat(rc2it);

    // quaternion code
    Quaternion q;
    Quaternion_identity(q);
    actpol_rotate_ITRS_to_GCRS(unixtime, q);
    double mat[3][3];
    //Quaternion_unit(q);
    Quaternion_to_matrix(q, mat);
    //print_mat(mat);

    for (int i = 0; i != 3; ++i)
        for (int j = 0; j != 3; ++j)
        {
            double x= fabs(rc2it[i][j]-mat[j][i]);
            //printf("%g\n", x);
            assert(x < 2e-11); // not sure why this is so large...
        }
}

void
ACTSite_init(ACTSite *s)
{
    s->latitude = ACTPOL_LATITUDE;
    s->east_longitude = ACTPOL_LONGITUDE_EAST;
    s->elevation_m = ACTPOL_ELEVATION_METERS;
    s->temperature_K = 273.;
    s->pressure_mb = 550.;
    s->relative_humidity = 0.2;
}

void
sofa_altaz_to_radec(ACTpolWeather *weather, double freq_GHz, double unixtime,
        double alt, double az, double *ra, double *dec)
{
    double zd = deg2rad(90) - alt;
    double utc1, utc2, dut1, elong, phi, hm, xp, yp, wl;

    utc1 = 2440587.5;
    utc2 = secs2days(unixtime);
    double mjd_utc = jd2mjd(utc1) + utc2;
    int stat = actpol_get_iers_bulletin_a(mjd_utc, &dut1, &xp, &yp);
    assert(stat == 0);
    xp = arcsec2rad(xp);
    yp = arcsec2rad(yp);

    elong = ACTPOL_LONGITUDE_EAST;
    phi = ACTPOL_LATITUDE;
    hm = ACTPOL_ELEVATION_METERS;
    wl = 299792.458/freq_GHz; // um

    stat = iauAtoc13("A", az, zd, utc1, utc2, dut1,
                       elong, phi, hm, xp, yp,
                       weather->pressure_mbar,
                       weather->temperature_C,
                       weather->relative_humidity,
                       wl, ra, dec);
    assert(stat == 0);
}

void
test_astro(void)
{
    Quaternion q, q1, q2, q3, q4;
    double alt=deg2rad(45.), az=deg2rad(123.), r[3], rp[3], alt2, az2;
    double unixtime = 1327019150.9773691;
    double mat[3][3], ra, dec, freq_GHz=150.;
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    check_focalplane_to_NWU(alt, az);
    check_horizon_to_itrs(alt, az);
    check_itrs_to_gcrs(unixtime);

    ACTpolArray *array = ACTpolArray_alloc(1);
    ACTpolArray_init(array, freq_GHz, 0., 0.);
    ACTpolFeedhorn_init(array->horn, 0., 0., 0.);
    ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
    ACTpolArrayCoords_init(coords, ACTPOL_COORDSYS_RA_SINDEC);
    ACTpolState *state = ACTpolState_alloc();
    ACTpolState_init(state);
    ACTpolState_update(state, unixtime, alt, az);
    ACTpolScan scan;
    ACTpolScan_init(&scan, alt, az, 5.);
    ACTpolArrayCoords_update_refraction(coords, &scan, &weather);
    ACTpolArrayCoords_update(coords, state);
    printf("noltaq ra, sin(dec) = %.8g, %.8g\n", rad2deg(coords->horn[0].a), coords->horn[0].b);

    ACTSite site;
    ACTSite_init(&site);
    observed_altaz_to_mean_radec( &site, freq_GHz, 1, &unixtime, &alt, &az, &ra, &dec );
    printf("slalib ra, sin(dec) = %.8g, %.8g\n", rad2deg(ra), sin(dec));

    double tol = arcsec2rad(0.05);
    assert(fabs(coords->horn[0].a - ra) < tol);
    assert(fabs(coords->horn[0].b - sin(dec)) < tol);

    tol = arcsec2rad(0.06);
    sofa_altaz_to_radec(&weather, freq_GHz, unixtime, alt, az, &ra, &dec);
    printf("  sofa ra, sin(dec) = %.8g, %.8g\n", rad2deg(ra), sin(dec));
    assert(fabs(coords->horn[0].a - ra) < tol);
    assert(fabs(coords->horn[0].b - sin(dec)) < tol);

    tol = 1e-15;
    actpol_altaz_to_radec(&weather, freq_GHz, unixtime, alt, az, &ra, &dec);
    printf("actpol ra, sin(dec) = %.8g, %.8g\n", rad2deg(ra), sin(dec));
    assert(fabs(coords->horn[0].a - ra) < tol);
    assert(fabs(coords->horn[0].b - sin(dec)) < tol);

    double gl, gb;
    actpol_altaz_to_galactic(&weather, freq_GHz, unixtime, alt, az, &gl, &gb);
    printf("               %.8g, %.8g\n", rad2deg(ra), rad2deg(dec));
    printf("galactic l,b = %.8g, %.8g\n", rad2deg(gl), rad2deg(gb));

    tol = deg2rad(0.01);
    for (int i = 0; i < 10; i++) {
        actpol_altaz_to_radec(&weather, freq_GHz, unixtime, alt, az, &ra, &dec);
        actpol_radec_to_crude_altaz(unixtime, ra, dec, &alt2, &az2);
        assert(fabs(alt - alt2) < tol);
        assert(fabs(az - az2) < tol);
        alt += deg2rad(2);
        az += deg2rad(2);
    }

    ACTpolState_free(state);
    ACTpolArrayCoords_free(coords);
    ACTpolArray_free(array);
}

void
test_astro2(void)
{
    Quaternion q, q1, q2, q3, q4;
    double alt=deg2rad(44.), az=deg2rad(128.), r[3], rp[3];
    double unixtime0 = 1337758400.;
    double mat[3][3], ra, dec, freq_GHz=150.;
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    check_focalplane_to_NWU(alt, az);
    check_horizon_to_itrs(alt, az);
    check_itrs_to_gcrs(unixtime0);

    ACTpolArray *array = ACTpolArray_alloc(1);
    ACTpolArray_init(array, freq_GHz, 0., 0.);
    ACTpolFeedhorn_init(array->horn, 0., 0., 0.);
    ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
    ACTpolArrayCoords_init(coords, ACTPOL_COORDSYS_RA_SINDEC);
    ACTpolState *state = ACTpolState_alloc();
    ACTpolState_init(state);
    ACTpolScan scan;
    ACTpolScan_init(&scan, alt, az, 5.);
    ACTpolArrayCoords_update_refraction(coords, &scan, &weather);

    ACTSite site;
    ACTSite_init(&site);
    double tol = arcsec2rad(0.06);

    for (int j = 0; j < 10; j++)
    {
        double unixtime = unixtime0 + 240.*j;
        ACTpolState_update(state, unixtime, alt, az);
        ACTpolArrayCoords_update(coords, state);
        //printf("noltaq ra, sin(dec) = %.8g, %.8g\n", rad2deg(coords->horn[0].a), coords->horn[0].b);

        observed_altaz_to_mean_radec( &site, freq_GHz, 1, &unixtime, &alt, &az, &ra, &dec );
        //printf("slalib ra, sin(dec) = %.8g, %.8g\n", rad2deg(ra), sin(dec));
        assert(fabs(sin(coords->horn[0].a) - sin(ra)) < tol);
        assert(fabs(cos(coords->horn[0].a) - cos(ra)) < tol);
        assert(fabs(coords->horn[0].b - sin(dec)) < tol);
    }

    ACTpolState_free(state);
    ACTpolArrayCoords_free(coords);
    ACTpolArray_free(array);
}

void
test_altaz(void)
{
    double freq_GHz = 150.;
    double unixtime0 = 1337758400.;
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    ACTpolArray *array = ACTpolArray_alloc(1);
    ACTpolArray_init(array, freq_GHz, 0., 0.);
    ACTpolFeedhorn_init(array->horn, 0., 0., 0.);
    ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
    ACTpolArrayCoords_init(coords, ACTPOL_COORDSYS_AZ_ALT);
    ACTpolState *state = ACTpolState_alloc();
    ACTpolState_init(state);
    ACTpolScan scan;
    ACTpolScan_init(&scan, 133, 45, 5.);
    ACTpolArrayCoords_update_refraction(coords, &scan, &weather);

    double boresight_alt = deg2rad(45);
    double boresight_az = deg2rad(-32);
    printf("%g %g\n", boresight_alt, boresight_az);
    ACTpolState_update(state, unixtime0, boresight_alt, boresight_az);
    ACTpolArrayCoords_update(coords, state);
    printf("%g %g\n", coords->horn[0].b, coords->horn[0].a);

    ACTpolState_free(state);
    ACTpolArrayCoords_free(coords);
    ACTpolArray_free(array);
}

int
main(int argc, char *argv[])
{
    test_altaz();
    test_astro2();
    test_astro();
    test_map();
    return 0;
}

