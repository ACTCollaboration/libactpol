
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <slalib.h>
#include <wcslib/wcslib.h>

#include "actpol/actpol.h"
#include "actpol/oldact.h"
#include "actpol/vec3.h"

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
        pix = ACTpolMap_sky2pix(map, ra-i*pixsize*0.2, dec+i*pixsize*0.2);
        stat = ACTpolMap_pix2sky(map, pix, &ra2, &dec2);
        assert(stat == 0);
        pix2 = ACTpolMap_sky2pix(map, ra2, dec2);
        //printf("pix, pix2 = %ld, %ld\n", pix, pix2);
        assert(pix == pix2);
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
test_astro(void)
{
    Quaternion q, q1, q2, q3, q4;
    double alt=deg2rad(45.), az=deg2rad(123.), r[3], rp[3];
    double unixtime = 1327019150.9773691;
    double mat[3][3], ra, dec, freq_GHz=150.;
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);

    check_focalplane_to_NWU(alt, az);
    check_horizon_to_itrs(alt, az);

    double ref = actpol_refraction(&weather, freq_GHz, alt);
    //altaz2itrs(alt-ref, az, r);
    actpol_ang2vec(-az, alt-ref, r);
    Quaternion_identity(q);
    //actpol_diurnal_aberration(r, q);
    actpol_rotate_NWU_to_ITRS(q);
    actpol_rotate_ITRS_to_GCRS(unixtime, q);
    Quaternion_to_matrix(q, mat);
    matrix_times_vec3(rp, mat, r);
    actpol_vec2ang(rp, &ra, &dec);
    printf("noltaq ra, dec = %g, %g\n", rad2deg(ra), rad2deg(dec));

    ACTSite site;
    ACTSite_init(&site);
    observed_altaz_to_mean_radec( &site, freq_GHz, 1, &unixtime, &alt, &az, &ra, &dec );
    printf("slalib ra, dec = %g, %g\n", rad2deg(ra), rad2deg(dec));
}

int
main(int argc, char *argv[])
{
    test_astro();
    test_map();
    return 0;
}

