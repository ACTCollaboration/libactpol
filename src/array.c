
#include <assert.h>
#include <malloc.h>
#include <stdlib.h>

#include "actpol/array.h"
#include "actpol/astro.h"
#include "actpol/math.h"
#include "actpol/state.h"
#include "actpol/vec3.h"
#include "actpol/util.h"

#include "debug.h"

void
actpol_detector_to_focalplane_rotation(double focalplane_x, double focalplane_y,
    double pol_angle, Quaternion q)
{
    Quaternion_r3(q, -pol_angle);
    Quaternion_r2_mul(focalplane_x, q);
    Quaternion_r1_mul(-focalplane_y, q);
    Quaternion_unit(q);
}

void
ACTpolFeedhorn_init(ACTpolFeedhorn *feedhorn, double focalplane_x,
    double focalplane_y, double pol_angle)
{
    actpol_detector_to_focalplane_rotation(focalplane_x,
        focalplane_y, pol_angle, feedhorn->focalplane_q);
}

ACTpolArray *
ACTpolArray_alloc(int nhorns)
{
    ACTpolArray *array;
    array = (ACTpolArray *)malloc(sizeof(ACTpolArray));
    assert(array);
    array->nhorns = nhorns;
    array->horn = (ACTpolFeedhorn *)memalign(32, nhorns*sizeof(ACTpolFeedhorn));
    assert(array->horn);
    return array;
}

void
ACTpolArray_free(ACTpolArray *array)
{
    free(array->horn);
    free(array);
}

void
ACTpolArray_init(ACTpolArray *array, double freq_GHz, double focalplane_x, double focalplane_y)
{
    array->freq_GHz = freq_GHz;
    actpol_detector_to_focalplane_rotation(focalplane_x,
        focalplane_y, 0., array->focalplane_q);
}

ACTpolFeedhorn *
ACTpolArray_get_feedhorn(ACTpolArray *array, int i)
{
    assert(i >= 0 && i < array->nhorns);
    return array->horn + i;
}

ACTpolArrayCoords *
ACTpolArrayCoords_alloc(const ACTpolArray *array)
{
    ACTpolArrayCoords *coords = (ACTpolArrayCoords *)malloc(sizeof(ACTpolArrayCoords));
    coords->array = array;
    coords->ref = (double *)malloc(sizeof(double) * array->nhorns);
    coords->horn = (ACTpolFeedhornCoords *)malloc(sizeof(ACTpolFeedhornCoords) * array->nhorns);
    return coords;
}

void
ACTpolArrayCoords_free(ACTpolArrayCoords *coords)
{
    free(coords->horn);
    free(coords->ref);
    free(coords);
}

void
ACTpolArrayCoords_init(ACTpolArrayCoords *coords)
{
    coords->mean_ref = arcsec2rad(30.);
}

ACTpolFeedhornCoords *
ACTpolArrayCoords_get_feedhorn_coords(ACTpolArrayCoords *coords, int i)
{
    assert(i >= 0 && i < coords->array->nhorns);
    return coords->horn + i;
}

void
ACTpolArrayCoords_update_refraction(ACTpolArrayCoords *coords,
    const ACTpolScan *scan, const ACTpolWeather *weather)
{
    const ACTpolArray *array = coords->array;

    Quaternion focalplane_to_NWU_q;
    Quaternion_identity(focalplane_to_NWU_q);
    actpol_rotate_focalplane_to_NWU(scan->mean_alt, scan->mean_az, focalplane_to_NWU_q);

    coords->mean_ref = 0.;
    for (int i = 0; i != coords->array->nhorns; ++i)
    {
        Quaternion q;
        Quaternion_mul(q, focalplane_to_NWU_q, array->horn[i].focalplane_q);

        double mat[3][3];
        Quaternion_to_matrix(q, mat);
        double alt = asin(mat[2][2]);

        coords->ref[i] = actpol_refraction(weather, array->freq_GHz, alt);
        assert(coords->ref[i] > 0. && coords->ref[i] < 1e-3);
        coords->mean_ref += coords->ref[i];
    }
    coords->mean_ref /= coords->array->nhorns;
    //DEBUG("mean_ref = %g\"\n", rad2arcsec(coords->mean_ref));
}

static void
compute_mean_focalplane_to_BCRS(const ACTpolArrayCoords *coords,
    const ACTpolState *state, Quaternion focalplane_to_BCRS)
{
    double r[3];
    Quaternion q;

    // focalplane -> topo
    Quaternion focalplane_to_topo;
    Quaternion_identity(focalplane_to_topo);
    actpol_rotate_focalplane_to_NWU(
        state->boresight_alt - coords->mean_ref,
        state->boresight_az,
        focalplane_to_topo);

    // diurnal aberration
    Quaternion diurnal_aberration, focalplane_to_apparent;
    Quaternion_mul(q, focalplane_to_topo, coords->array->focalplane_q);
    Quaternion_unit(q);
    Quaternion_to_matrix_col3(q, r);
    actpol_diurnal_aberration(r, diurnal_aberration);
    Quaternion_mul(focalplane_to_apparent, diurnal_aberration, focalplane_to_topo);

    // focalplane -> GCRS
    Quaternion focalplane_to_GCRS;
    Quaternion_mul(focalplane_to_GCRS, state->NWU_to_GCRS_q, focalplane_to_apparent);

    // center of array in GCRS
    Quaternion_mul(q, focalplane_to_GCRS, coords->array->focalplane_q);
    Quaternion_unit(q);
    Quaternion_to_matrix_col3(q, r);

    // annual aberration correction
    Quaternion GCRS_to_BCRS;
    actpol_aberration(r, state->earth_orbital_beta, GCRS_to_BCRS);

    Quaternion_mul(focalplane_to_BCRS, GCRS_to_BCRS, focalplane_to_GCRS);
    Quaternion_unit(focalplane_to_BCRS);
}

int
ACTpolArrayCoords_update(ACTpolArrayCoords *coords, const ACTpolState *state)
{
    const ACTpolArray *array = coords->array;

    Quaternion focalplane_to_BCRS;
    compute_mean_focalplane_to_BCRS(coords, state, focalplane_to_BCRS);

    #pragma omp parallel for
    for (int i = 0; i != array->nhorns; ++i)
    {
        Quaternion q;
        Quaternion_mul(q, focalplane_to_BCRS, array->horn[i].focalplane_q);

        double mat[3][3];
        Quaternion_conj(q); // transpose mat
        Quaternion_to_matrix(q, mat);
        double *p1 = mat[0];
        double *p2 = mat[1];
        double *r = mat[2];

        ACTpolFeedhornCoords *horn = coords->horn+i;
        horn->ra = atan2(r[1], r[0]);
        horn->dec = atan2(r[2], hypot(r[0],r[1]));
        horn->sindec = r[2];

        // w = r x z
        double w[3], z[3] = {0, 0, 1};
        vec3_cross_product(w, r, z);
        vec3_unit(w);

        // n = w x r
        double n[3];
        vec3_cross_product(n, w, r);

        double sin_g = vec3_dot_product(p1, w);
        double cos_g = vec3_dot_product(p1, n);
        horn->sin2gamma = 2*sin_g*cos_g;
        horn->cos2gamma = 2*cos_g*cos_g - 1;
    }

    return 0;
}

int
ACTpolArrayCoords_update_fast(ACTpolArrayCoords *coords, const ACTpolState *state)
{
    const ACTpolArray *array = coords->array;

    Quaternion focalplane_to_BCRS;
    compute_mean_focalplane_to_BCRS(coords, state, focalplane_to_BCRS);

    Quaternion q0;
    Quaternion_mul(q0, focalplane_to_BCRS, array->focalplane_q);
    double r0[3];
    Quaternion_to_matrix_col3(q0, r0);

    double x0 = r0[0];
    double y0 = r0[1];
    double ra0 = atan2(y0, x0);

    #pragma omp parallel for
    for (int i = 0; i != array->nhorns; ++i)
    {
        Quaternion q;
        Quaternion_mul(q, focalplane_to_BCRS, array->horn[i].focalplane_q);

        double p1[3], r[3];
        Quaternion_to_matrix_col1(q, p1);
        Quaternion_to_matrix_col3(q, r);

        ACTpolFeedhornCoords *horn = coords->horn+i;

        /*
         * Define z = y/x, z0 = y0/x0.
         *
         * atan(z) - atan(z0) = atan(z) + atan(-z0)
         *                    = atan((z - z0)/(1 + z*z0))
         *
         * Let p = (z - z0)/(1 + z*z0)
         *       = (x0*y - x*y0)/(x*x0 + y*y0).
         *
         * For |p| << 1, atan(p) ~ p/(1 + 0.28*p^2)
         * (Abramowitz & Stegun, 4.4.48)
         *
         * Approximation appears to be good to ~0.01". By contrast,
         * 2nd order Taylor expansion is only good to ~10".
         */

        double p = (x0*r[1] - r[0]*y0)/(x0*r[0] + y0*r[1]);
        double ra1 = p/(1. + 0.28*p*p);
        horn->ra = ra0 + ra1;
        horn->sindec = r[2];

        // w = r x z
        const double z[3] = {0., 0., 1.};
        double w[3];
        vec3_cross_product(w, r, z);

        // n = w x r
        double n[3];
        vec3_cross_product(n, w, r);

        double sin_g = vec3_dot_product(p1, w);
        double cos_g = vec3_dot_product(p1, n);
        double norm2 = sin_g*sin_g + cos_g*cos_g;
        /*
        horn->sin2gamma = 2.*sin_g*cos_g/norm2;
        horn->cos2gamma = 2.*cos_g*cos_g/norm2 - 1.;
        */
        double two_cos_g_div_norm2 = 2.*cos_g/norm2;
        horn->sin2gamma = sin_g*two_cos_g_div_norm2;
        horn->cos2gamma = cos_g*two_cos_g_div_norm2 - 1.;
    }

    return 0;
}

